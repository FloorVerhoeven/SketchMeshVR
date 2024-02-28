#include "MeshExtrusion.h"
#include "LaplacianRemesh.h"
#include "Plane.h"
#include <iostream>
#include <algorithm>
#include <igl/unproject_ray.h>
#include <igl/triangle/triangulate.h>
#include <igl/edge_topology.h>
#include <igl/slice.h>

using namespace std;
using namespace igl;

int MeshExtrusion::prev_vertex_count = -1;
int MeshExtrusion::ID = -1;

bool MeshExtrusion::extrude_prepare(Stroke& base, SurfacePath& surface_path) {
	base.counter_clockwise();
	
	bool success = surface_path.create_from_stroke_extrude(base);
	if (!success) {
		return false;
	}
	post_extrude_prepare_update_points(base, surface_path);
	return true;
}

/*Here stroke is the silhouette stroke, and base is the extrusion base stroke. */
bool MeshExtrusion::extrude_main(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers, Eigen::VectorXi &part_of_original_stroke, Eigen::VectorXi &new_mapped_indices, Eigen::VectorXi &sharp_edge, SurfacePath& surface_path, Stroke& stroke, Stroke& base, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport) {
	if(V.rows() != prev_vertex_count) {
		ID++;
		prev_vertex_count = V.rows();
	}

	Mesh m(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, ID);
	stroke.counter_clockwise();
	bool remesh_success = true;
	Eigen::VectorXi boundary_vertices = LaplacianRemesh::remesh_extrusion_remove_inside(m, surface_path, model, view, proj, viewport, remesh_success);

	if (!remesh_success) {
		//NOTE: remesh will have cleared m.new_mapped_indices but since this is cleared every time it is needed, there's no need to revert this
		return false;
	}


	vector<int> sharp_edge_indices;
	for(int i = 0; i < m.sharp_edge.rows(); i++) {
		if(m.sharp_edge[i]) {
			sharp_edge_indices.push_back(i);
		}
	}

	Eigen::MatrixXi EV, FE, EF;
	igl::edge_topology(m.V, m.F, EV, FE, EF);

    //Store the vertex indices for sharp edges in the original mesh
	Eigen::MatrixXi sharpEV;
	Eigen::VectorXi sharpEV_row_idx, sharpEV_col_idx(2);
	sharpEV_row_idx = Eigen::VectorXi::Map(sharp_edge_indices.data(), sharp_edge_indices.size());
	sharpEV_col_idx.col(0) << 0, 1;
	igl::slice(EV, sharpEV_row_idx, sharpEV_col_idx, sharpEV);

    //Compute the real 3D positions of the silhouette stroke
	Eigen::RowVector3d center(0, 0, 0), pr;
	Eigen::Vector3d normal(0, 0, 0), source, dir;
	Eigen::Vector2d tmp0;
	get_normal_and_center(center, normal, m, boundary_vertices);

	Eigen::Matrix4f modelview = stroke.viewer.core.view * stroke.viewer.core.model;
	igl::project(center, modelview, stroke.viewer.core.proj, stroke.viewer.core.viewport, pr);
	tmp0 = pr.leftCols(2);
	unproject_ray(tmp0, modelview, stroke.viewer.core.proj, stroke.viewer.core.viewport, source, dir);	

	Eigen::Vector3d camera_to_center = source + dir;
	Eigen::Vector3d normal2 = normal.cross(camera_to_center);
	normal2.normalize();

	double max, min;
	int most_left_vertex_idx = 0, most_right_vertex_idx = 0;
	find_left_and_right(most_left_vertex_idx, most_right_vertex_idx, min, max, m, center, boundary_vertices, normal2);

	Eigen::Vector3d v_tmp = (m.V.row(boundary_vertices[most_left_vertex_idx]) - m.V.row(boundary_vertices[most_right_vertex_idx])).transpose();
	Eigen::Vector3d pop_surface_normal = normal.cross(v_tmp);
	pop_surface_normal.normalize();


	Eigen::MatrixXd silhouette_vertices(0, 3);
	Eigen::RowVector3d pop_surface_point = m.V.row(boundary_vertices[most_left_vertex_idx]);
	compute_silhouette_positions(silhouette_vertices, stroke, modelview, center, normal2, min, max, pop_surface_point, pop_surface_normal);

	if(silhouette_vertices.rows() == 0) {
		cout << "No valid silhouette vertices. Try again." << endl;
		return false;
	}

	//Possibly reverse silhouette_vertices, so they start closest to the most_left_vertex
	if((m.V.row(boundary_vertices[most_left_vertex_idx]) - silhouette_vertices.row(0)).norm() > (m.V.row(boundary_vertices[most_right_vertex_idx]) - silhouette_vertices.row(0)).norm()) {
		silhouette_vertices = silhouette_vertices.colwise().reverse().eval();
	}

	
	vector<int> sil_original_indices = add_silhouette_vertices(m, stroke.get_ID(), silhouette_vertices);
	vector<int> sil_original_indices_looped = sil_original_indices;
	sil_original_indices_looped.push_back(sil_original_indices[0]);
	stroke.set_closest_vert_bindings(sil_original_indices_looped);

	//Create one of the two loops
	Eigen::MatrixXd front_loop3D = silhouette_vertices;
	vector<int> front_loop_base_original_indices;
	create_loop(m, front_loop3D, boundary_vertices, front_loop_base_original_indices, most_right_vertex_idx, most_left_vertex_idx);


	//Create the second loop
	Eigen::MatrixXd back_loop3D = silhouette_vertices.colwise().reverse();
	vector<int> back_loop_base_original_indices;
	create_loop(m, back_loop3D, boundary_vertices, back_loop_base_original_indices, most_left_vertex_idx, most_right_vertex_idx);

	Eigen::Vector3d x_vec = m.V.row(boundary_vertices[most_right_vertex_idx]) - m.V.row(boundary_vertices[most_left_vertex_idx]);
	x_vec.normalize();
	Eigen::Vector3d y_vec = normal.cross(x_vec);
	Eigen::Vector3d offset = x_vec.cross(y_vec);
	offset *= 0.05;

    //Generate the two half-domes
	generate_mesh(m, front_loop3D, center, x_vec*-1, y_vec, offset, silhouette_vertices.rows(), front_loop_base_original_indices, sil_original_indices);
	reverse(sil_original_indices.begin(), sil_original_indices.end());
	generate_mesh(m, back_loop3D, center, x_vec*-1, y_vec, offset, silhouette_vertices.rows(), back_loop_base_original_indices, sil_original_indices);

    //Update tracking variables with new vertex indices
	post_extrude_main_update_points(stroke, silhouette_vertices);
	post_extrude_main_update_bindings(base, surface_path);
	update_sharp_edges(m, sharpEV);

	return true;
}

/** Generates a mesh inside the loop3D vertices. **/
void MeshExtrusion::generate_mesh(Mesh& m, Eigen::MatrixXd loop3D, Eigen::Vector3d center, Eigen::Vector3d x_vec, Eigen::Vector3d y_vec, Eigen::Vector3d offset, int nr_silhouette_vert, vector<int> loop_base_original_indices, vector<int> sil_original_indices) {
	int size_before_gen = m.V.rows();
	Eigen::MatrixXd loop2D(loop3D.rows(), 2);
	Eigen::Vector3d vec;
	Eigen::MatrixXi loop_stroke_edges(loop3D.rows(), 2);

	for(int i = 0; i < loop3D.rows(); i++) {
		vec = (loop3D.row(i).transpose() - center);
		loop2D.row(i) << vec.dot(x_vec), vec.dot(y_vec);
		loop_stroke_edges.row(i) << i, ((i + 1) % loop3D.rows());
	}

	Eigen::MatrixXd V2;
	Eigen::MatrixXi F2;
	Eigen::MatrixXi vertex_markers, edge_markers;
	igl::triangle::triangulate(loop2D, loop_stroke_edges, Eigen::MatrixXd(0, 0), Eigen::MatrixXi::Constant(loop2D.rows(), 1, 1), Eigen::MatrixXi::Constant(loop_stroke_edges.rows(), 1, 1), "YqV", V2, F2, vertex_markers, edge_markers); //Capital Q silences triangle's output in cmd line. Also retrieves markers to indicate whether or not an edge/vertex is on the mesh boundary
																																																															  
	Eigen::RowVector3d vert;
	for(int i = 0; i < V2.rows(); i++) {
		if(i >= loop2D.rows()) { //Interior vertices
			vert = center;
			vert += x_vec.transpose()*V2(i, 0);
			vert += y_vec.transpose()*V2(i, 1);
			vert += offset.transpose();
			m.V.conservativeResize(m.V.rows() + 1, Eigen::NoChange);
			m.V.row(m.V.rows() - 1) = vert;
			m.part_of_original_stroke.conservativeResize(m.part_of_original_stroke.rows() + 1);
			m.vertex_boundary_markers.conservativeResize(m.vertex_boundary_markers.rows() + 1);
			m.part_of_original_stroke(m.part_of_original_stroke.rows() - 1) = 0;
			m.vertex_boundary_markers(m.vertex_boundary_markers.rows() - 1) = 0; //Interior mesh vertex, so non-boundary
		}
	}
	update_face_indices(m, F2, sil_original_indices, loop_base_original_indices, nr_silhouette_vert, size_before_gen, loop2D.rows());
}

void MeshExtrusion::get_normal_and_center(Eigen::RowVector3d& center, Eigen::Vector3d& normal, Mesh& m, Eigen::VectorXi &boundary_vertices) {
	for(int i = 0; i < boundary_vertices.rows(); i++) {
		center += m.V.row(boundary_vertices[i]);
	}
	center /= boundary_vertices.rows();

	Eigen::Vector3d vec0, vec1;
	for(int i = 0; i < boundary_vertices.rows(); i++) {
		vec0 = m.V.row(boundary_vertices[i]) - center;
		vec1 = m.V.row(boundary_vertices[(i + 1) % boundary_vertices.rows()]) - center;
		normal += vec1.cross(vec0);
	}
	normal.normalize();
}

/** Computes the left- and rightmost vertex index of the base stroke as seen from how the silhouette stroke was drawn. **/
void MeshExtrusion::find_left_and_right(int& most_left_vertex_idx, int& most_right_vertex_idx, double& min, double& max, Mesh& m, Eigen::RowVector3d& center, Eigen::VectorXi &boundary_vertices, Eigen::Vector3d& normal2) {
	Eigen::Vector3d center_to_vertex = m.V.row(boundary_vertices[0]) - center;
	double dot_prod = normal2.dot(center_to_vertex);
	max = dot_prod;
	min = dot_prod;
	for(int i = 1; i < boundary_vertices.rows(); i++) {
		center_to_vertex = m.V.row(boundary_vertices[i]) - center;
		dot_prod = normal2.dot(center_to_vertex);
		if(dot_prod > max) {
			max = dot_prod;
			most_left_vertex_idx = i;
		} else if(dot_prod < min) {
			min = dot_prod;
			most_right_vertex_idx = i;
		}
	}
}

/** Takes the silhouette's 2D positions as it was drawn and projects it onto the correct plane (more or less parallel to the base stroke). **/
void MeshExtrusion::compute_silhouette_positions(Eigen::MatrixXd& silhouette_vertices, Stroke& stroke, Eigen::Matrix4f& modelview, Eigen::RowVector3d& center, Eigen::Vector3d& normal2, int min, int max, Eigen::RowVector3d& pop_surface_point, Eigen::Vector3d& pop_surface_normal) {
	Eigen::RowVector3d v; 
	Eigen::Vector3d source, dir, center_to_vertex;
	Plane pop_surface(pop_surface_point, pop_surface_normal);

	for(int i = 0; i < stroke.get_stroke2DPoints().rows(); i++) {
		Eigen::Vector2d tmp = stroke.get_stroke2DPoints().row(i);
		unproject_ray(tmp, modelview, stroke.viewer.core.proj, stroke.viewer.core.viewport, source, dir);
		v = pop_surface.intersect_point(source, dir);
		center_to_vertex = v - center;
		if(normal2.dot(center_to_vertex) > max || normal2.dot(center_to_vertex) < min) { //Silhouette point is out of range
			continue;
		}
		silhouette_vertices.conservativeResize(silhouette_vertices.rows() + 1, Eigen::NoChange);
		silhouette_vertices.row(silhouette_vertices.rows() - 1) = v;
	}
}

/** Adds the silhouette vertices to the mesh and also inserts their tracking variables. **/
vector<int> MeshExtrusion::add_silhouette_vertices(Mesh& m, int stroke_ID, Eigen::MatrixXd& silhouette_vertices) {
	vector<int> sil_original_indices;
	int size_before_silhouette = m.V.rows();
	m.V.conservativeResize(m.V.rows() + silhouette_vertices.rows(), Eigen::NoChange);
	m.part_of_original_stroke.conservativeResize(m.part_of_original_stroke.rows() + silhouette_vertices.rows());
	m.vertex_boundary_markers.conservativeResize(m.vertex_boundary_markers.rows() + silhouette_vertices.rows());
	for(int i = 0; i < silhouette_vertices.rows(); i++) {
		m.V.row(size_before_silhouette + i) = silhouette_vertices.row(i);
		m.part_of_original_stroke[size_before_silhouette + i] = 0;
		m.vertex_boundary_markers[size_before_silhouette + i] = stroke_ID;
		sil_original_indices.push_back(size_before_silhouette + i);
	}
	return sil_original_indices;
}

/** Creates a loop consisting of the silhouette vertices and half of the base stroke vertices. **/
void MeshExtrusion::create_loop(Mesh& m, Eigen::MatrixXd& loop3D, Eigen::VectorXi& boundary_vertices, vector<int> &loop_base_original_indices, int start_idx, int end_idx) {
	int idx = start_idx;
	while(true) {
		loop3D.conservativeResize(loop3D.rows() + 1, Eigen::NoChange);
		loop3D.row(loop3D.rows() - 1) = m.V.row(boundary_vertices[idx]);
		loop_base_original_indices.push_back(boundary_vertices[idx]);
		if(idx == end_idx) {
			break;
		}
		idx++;
		if(idx == boundary_vertices.rows()) {
			idx = 0;
		}
	}

}

/** Updates the sharp_edge tracker with new edge indices after the mesh topology changed. **/
void MeshExtrusion::update_sharp_edges(Mesh& m, Eigen::MatrixXi sharpEV) {
	Eigen::MatrixXi EV, FE, EF;
	igl::edge_topology(m.V, m.F, EV, FE, EF);
	m.sharp_edge.resize(EV.rows());
	m.sharp_edge.setZero();

	int start, end, equal_pos;
	Eigen::VectorXi col1Equals, col2Equals;
	for(int i = 0; i < sharpEV.rows(); i++) {
		start = sharpEV(i, 0); //No vertices were removed in the meantime, so use indices as-is
		end = sharpEV(i, 1);
		if(start == -1 || end == -1) { //Sharp edge no longer exists
			continue;
		}

		col1Equals = EV.col(0).cwiseEqual(std::min(start, end)).cast<int>();
		col2Equals = EV.col(1).cwiseEqual(std::max(start, end)).cast<int>();
		(col1Equals + col2Equals).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge

		m.sharp_edge[equal_pos] = 1; //Set this edge to be sharp
	}
}

/** Updates the old face indices with the new indexing and inserts new faces that are made in triangulation. **/
void MeshExtrusion::update_face_indices(Mesh& m, Eigen::MatrixXi& F2, vector<int> sil_original_indices, vector<int> loop_base_original_indices, int nr_silhouette_vert, int size_before_gen, int loop2D_size) {
	int original_F_size = m.F.rows();
	m.F.conservativeResize(m.F.rows() + F2.rows(), Eigen::NoChange);
	int vert_idx_in_mesh;
	for(int i = 0; i < F2.rows(); i++) {
		for(int j = 0; j < 3; j++) {
			if(F2(i, j) < nr_silhouette_vert) { //Silhouette vertex
				vert_idx_in_mesh = sil_original_indices[F2(i, j)];
			} else if(F2(i, j) < loop2D_size) { //Extrusion base vertex
				vert_idx_in_mesh = loop_base_original_indices[F2(i, j) - nr_silhouette_vert];
			} else {
				vert_idx_in_mesh = size_before_gen + F2(i, j) - loop2D_size; //Interior point
			}
			m.F(original_F_size + i, j) = vert_idx_in_mesh;
		}
	}
}

/** Updates the base stroke's 3DPoints with the new vertices. **/
void MeshExtrusion::post_extrude_prepare_update_points(Stroke& stroke, SurfacePath& surface_path) {
	vector<PathElement> path = surface_path.get_path();
	Eigen::MatrixX3d new_3DPoints(path.size(), 3);

	for(int i = 0; i < path.size(); i++) {
		new_3DPoints.row(i) = path[i].get_vertex().transpose();
	}

	stroke.set3DPoints(new_3DPoints);
}

/** Updates the silhouette stroke's 3DPoints with the new vertices. **/
void MeshExtrusion::post_extrude_main_update_points(Stroke &stroke, Eigen::MatrixXd new_positions) {
	Eigen::MatrixX3d new_3DPoints(new_positions.rows() + 1, 3);
	for(int i = 0; i < new_positions.rows(); i++) {
		new_3DPoints.row(i) = new_positions.row(i);
	}
	new_3DPoints.row(new_3DPoints.rows() - 1) = new_positions.row(0); //Make looped

	stroke.set3DPoints(new_3DPoints); 
}

/** Update the vertex bindings for the extrusion base stroke. Assumes indexing in m.V before any vertices are removed, so it will become the correct indexing after stroke.update_vert_bindings() is called. **/
void MeshExtrusion::post_extrude_main_update_bindings(Stroke& base, SurfacePath& surface_path) {
	vector<PathElement> path = surface_path.get_path();
	vector<int> new_closest_vertex_indices(path.size());
	for(int i = 0; i < path.size() - 1; i++) {
		new_closest_vertex_indices[i] = path[i].get_v_idx();
	}
	new_closest_vertex_indices[path.size() - 1] = new_closest_vertex_indices[0]; //In laplacianRemesh, the last path vertex is skipped because it is a duplicate of the first one. So need to set it manually

	base.set_closest_vert_bindings(new_closest_vertex_indices);
}
