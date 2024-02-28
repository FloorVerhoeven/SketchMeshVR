#include "MeshExtrusion.h"
#include "LaplacianRemesh.h"
#include "CleanStroke3D.h"
#include <iostream>
#include <algorithm>
#include <igl/triangle/triangulate.h>
#include <igl/edge_topology.h>
#include <igl/slice.h>
#include <igl/cat.h>

using namespace std;
using namespace igl;

bool MeshExtrusion::extrude_prepare(Stroke& base, SurfacePath& surface_path) {
	base.counter_clockwise();
	
	bool success = surface_path.create_from_stroke_extrude(base);
	if (!success) {
		return false;
	}
	post_extrude_prepare_update_points(base, surface_path);
	return true;
}

bool MeshExtrusion::extrude_main(Mesh& m, SurfacePath& surface_path, Stroke& stroke, Stroke& base, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, Eigen::MatrixXi& replacing_vertex_bindings) {
	Eigen::MatrixXi start_F = m.F;
	Eigen::MatrixXd start_V = m.V;
	Eigen::VectorXi start_edge_boundary_markers = m.edge_boundary_markers;
	Eigen::VectorXi start_vertex_is_fixed = m.vertex_is_fixed;
	Eigen::VectorXi start_sharp_edge = m.sharp_edge;
	Eigen::VectorXi start_new_mapped_indices = m.new_mapped_indices;

	stroke.counter_clockwise();
	int remesh_success = 1;
	Eigen::VectorXi boundary_vertices = LaplacianRemesh::remesh_extrusion_remove_inside(m, surface_path, model, view, proj, viewport, remesh_success, replacing_vertex_bindings);

	if (!remesh_success) {
		//NOTE: remesh will have cleared m.new_mapped_indices but since this is cleared every time it is needed, there's no need to revert this
		return false;
	}

	Eigen::MatrixXi EV, FE, EF;
	igl::edge_topology(m.V, m.F, EV, FE, EF);

	Eigen::MatrixXi original_sharp_or_boundary_edges(0, 4);
	for (int i = 0; i < m.sharp_edge.rows(); i++) {
		if (m.sharp_edge[i] || m.edge_boundary_markers[i]) {
			original_sharp_or_boundary_edges.conservativeResize(original_sharp_or_boundary_edges.rows() + 1, Eigen::NoChange);
			original_sharp_or_boundary_edges.bottomRows(1) << EV(i, 0), EV(i, 1), m.edge_boundary_markers[i], m.sharp_edge[i];
		}
	}

    //Compute the real 3D positions of the silhouette stroke
	Eigen::RowVector3d center(0, 0, 0);
	Eigen::Vector3d normal(0, 0, 0);
	get_normal_and_center(center, normal, m, boundary_vertices);

	stroke.project_with_PCA_given_target(normal);
	Eigen::MatrixX3d sil_3D_points = stroke.get3DPoints();
	Eigen::MatrixXd silhouette_vertices(sil_3D_points.rows() - 1, 3);
	if (stroke.has_been_reversed) {
		silhouette_vertices = sil_3D_points.bottomRows(sil_3D_points.rows() - 1);
	} else {
		silhouette_vertices = sil_3D_points.topRows(sil_3D_points.rows()-1);
	}


	int most_left_vertex_idx = 0, most_right_vertex_idx = 0;
	find_left_and_right(most_left_vertex_idx, most_right_vertex_idx, m, boundary_vertices, silhouette_vertices.row(0), silhouette_vertices.row(silhouette_vertices.rows()-1), center);
	
	Eigen::RowVector3d dir = (silhouette_vertices.row(silhouette_vertices.rows()-1) - silhouette_vertices.row(0)).normalized();
	remove_out_of_bounds_silhouette(silhouette_vertices, center, m.V.row(boundary_vertices[most_left_vertex_idx]), m.V.row(boundary_vertices[most_right_vertex_idx]), dir);

	if(silhouette_vertices.rows() == 0) {
		m.F = start_F;
		m.V = start_V;
		m.new_mapped_indices = start_new_mapped_indices;
		m.edge_boundary_markers = start_edge_boundary_markers;
		m.vertex_is_fixed = start_vertex_is_fixed;
		m.sharp_edge = start_sharp_edge;
		std::cerr << "No valid silhouette vertices. Try again." << endl;
		return false;
	}


	//This makes sure that we are inserting the correct base vertices at the correct spot before/after the silhouette vertices
	if ((m.V.row(boundary_vertices[most_right_vertex_idx]) - silhouette_vertices.row(0)).norm() > (m.V.row(boundary_vertices[most_right_vertex_idx]) - silhouette_vertices.row(silhouette_vertices.rows()-1)).norm()) {
		Eigen::MatrixXd tmp_sil_vert = silhouette_vertices.colwise().reverse();
		silhouette_vertices = tmp_sil_vert;
	}

	//double sil_unit_length = get_stroke_unit_length(silhouette_vertices); //Distance used for resampling later. Calculate this before adding the left- and right-most vertices (because they might be far away from the silhouette endpoints)

	silhouette_vertices.conservativeResize(silhouette_vertices.rows() + 1, Eigen::NoChange);
	Eigen::MatrixXd tmp_sil = silhouette_vertices.topRows(silhouette_vertices.rows() - 1);
	silhouette_vertices.block(1, 0, silhouette_vertices.rows() - 1, 3) = tmp_sil; 
	silhouette_vertices.row(0) = m.V.row(boundary_vertices[most_right_vertex_idx]); //Add the right-most vertex
	silhouette_vertices.conservativeResize(silhouette_vertices.rows() + 1, Eigen::NoChange);
	silhouette_vertices.row(silhouette_vertices.rows() - 1) = m.V.row(boundary_vertices[most_left_vertex_idx]); //Add the left-most vertex

	double sil_unit_length = get_base_unit_length(m.V, boundary_vertices); //Use the intersample distance for the extrusion's base stroke, since this gives a more uniform sampling when we make the loops
	Eigen::MatrixXd silhouette_vertices_smoothed = smooth_stroke(silhouette_vertices);
	silhouette_vertices_smoothed = CleanStroke3D::resample_by_length_sub(silhouette_vertices_smoothed, 0, silhouette_vertices_smoothed.rows()-1, sil_unit_length); //Resample, including the space between the outer silhouette vertices and the left- and right-most base vertices
	tmp_sil = silhouette_vertices_smoothed.middleRows(1, silhouette_vertices_smoothed.rows() - 2); //Remove the left- and right-most vertices again
	silhouette_vertices = tmp_sil;

	Eigen::MatrixXi added_edges(0, 4);
	int size_before_adding_sil = m.V.rows();
	vector<int> sil_original_indices = add_silhouette_vertices(m, stroke.get_ID(), silhouette_vertices, added_edges);
	vector<int> sil_original_indices_looped = sil_original_indices;
	sil_original_indices_looped.insert(sil_original_indices_looped.begin(), boundary_vertices[most_right_vertex_idx]);
	sil_original_indices_looped.push_back(boundary_vertices[most_left_vertex_idx]);
	sil_original_indices_looped.push_back(sil_original_indices[0]);
	stroke.set_closest_vert_bindings(sil_original_indices_looped);

	//Create the first loop
	Eigen::MatrixXd front_loop3D = silhouette_vertices;
	vector<int> front_loop_base_original_indices;
	create_loop(m, front_loop3D, boundary_vertices, front_loop_base_original_indices, most_left_vertex_idx, most_right_vertex_idx);
	//Create the second loop
	Eigen::MatrixXd back_loop3D = silhouette_vertices.colwise().reverse();
	vector<int> back_loop_base_original_indices;
	create_loop(m, back_loop3D, boundary_vertices, back_loop_base_original_indices, most_right_vertex_idx, most_left_vertex_idx);

	added_edges.conservativeResize(added_edges.rows() + 1, Eigen::NoChange); 
	added_edges.bottomRows(1) <<  boundary_vertices[most_right_vertex_idx], size_before_adding_sil, stroke.get_ID(), 0; //Add connection between first silhouette vertex and closest base vertex

	added_edges.conservativeResize(added_edges.rows() + 1, Eigen::NoChange);
	added_edges.bottomRows(1) << boundary_vertices[most_left_vertex_idx], size_before_adding_sil + silhouette_vertices.rows() - 1, stroke.get_ID(), 0; //Add connection between last silhouette vertex and closest base vertex
	
	//Get center and normal for first half-dome and generate mesh
	Eigen::RowVector3d new_center_front(0,0,0);
	for (int i = 0; i < front_loop3D.rows(); i++) {
		new_center_front += front_loop3D.row(i);
	}
	for (int i = 0; i < silhouette_vertices.rows(); i++) {
		new_center_front -= silhouette_vertices.row(i);
	}
	new_center_front /= (front_loop3D.rows() - silhouette_vertices.rows());

	Eigen::Vector3d normal_front(0,0,0);
	Eigen::Vector3d vec0, vec1;
	for (int i = 0; i < front_loop3D.rows(); i++) {
		vec0 = front_loop3D.row(i) - new_center_front;
		vec1 = front_loop3D.row((i + 1) % front_loop3D.rows()) - new_center_front;
		normal_front += vec0.cross(vec1);
	}
	normal_front.normalize();

	Eigen::Vector3d x_vec = (m.V.row(boundary_vertices[most_right_vertex_idx]) - m.V.row(boundary_vertices[most_left_vertex_idx])).normalized();
	Eigen::Vector3d y_vec = normal_front.cross(x_vec);
	Eigen::Vector3d offset = x_vec.cross(y_vec);
	offset *= 0.05;

	generate_mesh(m, front_loop3D, new_center_front, x_vec, y_vec, offset, silhouette_vertices.rows(), front_loop_base_original_indices, sil_original_indices);

	//Get center and normal for second half-dome and generate mesh
	Eigen::RowVector3d new_center_back(0, 0, 0);
	for (int i = 0; i < back_loop3D.rows(); i++) {
		new_center_back += back_loop3D.row(i);
	}
	for (int i = 0; i < silhouette_vertices.rows(); i++) {
		new_center_back -= silhouette_vertices.row(i);
	}
	new_center_back /= (back_loop3D.rows() - silhouette_vertices.rows());


	Eigen::Vector3d normal_back(0, 0, 0);
	for (int i = 0; i < back_loop3D.rows(); i++) {
		vec0 = back_loop3D.row(i) - new_center_back;
		vec1 = back_loop3D.row((i + 1) % back_loop3D.rows()) - new_center_back;
		normal_back += vec0.cross(vec1);
	}
	normal_back.normalize();

	y_vec = normal_back.cross(x_vec);
	offset = x_vec.cross(y_vec);
	offset *= 0.05;

	reverse(sil_original_indices.begin(), sil_original_indices.end());
	generate_mesh(m, back_loop3D, new_center_back, x_vec, y_vec, offset, silhouette_vertices.rows(), back_loop_base_original_indices, sil_original_indices);


    //Update tracking variables with new vertex indices
	post_extrude_main_update_bindings(base, surface_path);

	Eigen::MatrixXi new_edge_indicators = original_sharp_or_boundary_edges;
	new_edge_indicators = igl::cat(1, new_edge_indicators, added_edges); //Both are already correctly indexed into the new V
	try {
		update_edge_indicators(m, new_edge_indicators);
	}
	catch(int ex){
		if (ex == -1) {
			std::cerr << "Extrusion results in a non edge-manifold mesh, which cannot be processed. Please try again. " << std::endl;
			m.F = start_F;
			m.V = start_V;
			m.new_mapped_indices = start_new_mapped_indices;
			m.edge_boundary_markers = start_edge_boundary_markers;
			m.vertex_is_fixed = start_vertex_is_fixed;
			m.sharp_edge = start_sharp_edge;
			return false;
		}
	}
	return true;
}

Eigen::MatrixXd MeshExtrusion::smooth_stroke(Eigen::MatrixXd& input_points) {
	Eigen::MatrixXd new_points = Eigen::MatrixXd::Zero(input_points.rows(), 3);
	int nr_iterations = max(5.0, input_points.rows() / 4.0);
	for (int i = 0; i < nr_iterations; i++) {
		move_to_middle(input_points, new_points);
		move_to_middle(new_points, input_points);
	}
	return new_points;
}

void MeshExtrusion::move_to_middle(Eigen::MatrixXd &positions, Eigen::MatrixXd &new_positions) {
	int n = positions.rows();
	Eigen::Vector3d prev, cur, next;
	new_positions.row(0) = positions.row(0);
	new_positions.bottomRows(1) = positions.bottomRows(1);

	for (int i = 1; i < n - 1; i++) { //Skip the first and last vertex, because we don't have a loop here and want to keep them fixed
		prev = positions.row(((i - 1) + n) % n);
		cur = positions.row(i % n);
		next = positions.row(((i + 1) + n) % n);

		new_positions(i, 0) = (cur[0] * 2 + prev[0] + next[0]) / 4;
		new_positions(i, 1) = (cur[1] * 2 + prev[1] + next[1]) / 4;
		new_positions(i, 2) = (cur[2] * 2 + prev[2] + next[2]) / 4;
	}
}

double MeshExtrusion::get_stroke_unit_length(Eigen::MatrixXd& vertices) {
	//Assumes that the matrix of vertices can be used as is (i.e. if we need the unit length of a closed curve, vertices needs to be looped already)
	double total = 0.0;
	for (int i = 0; i < vertices.rows() - 1; i++) {
		total += (vertices.row(i) - vertices.row(i + 1)).norm();
	}
	return total / (vertices.rows() - 1);
}

//Boundary_vertices here won't be looped, but we also don't necessarily need a loop to determine the average sampling distance
double MeshExtrusion::get_base_unit_length(Eigen::MatrixXd& V, Eigen::VectorXi& boundary_vertices) {
	double total = 0.0;
	for (int i = 0; i < boundary_vertices.rows() - 1; i++) {
		total += (V.row(boundary_vertices[i]) - V.row(boundary_vertices[i + 1])).norm();
	}
	return total / (boundary_vertices.rows() - 1);
}

void MeshExtrusion::remove_out_of_bounds_silhouette(Eigen::MatrixXd& silhouette_vertices, Eigen::RowVector3d& center, const Eigen::RowVector3d& left_most, const Eigen::RowVector3d& right_most, Eigen::RowVector3d& dir) {
	vector<int> keep_indices;
	Eigen::Vector3d center_to_vertex;
	double most_left_dot_prod = dir.dot(left_most - center);
	double most_right_dot_prod = dir.dot(right_most - center);
	double dot_prod;

	for (int i = 0; i < silhouette_vertices.rows(); i++) {
		center_to_vertex = silhouette_vertices.row(i) - center;
		dot_prod = dir.dot(center_to_vertex);
		if (dot_prod > most_left_dot_prod || dot_prod < most_right_dot_prod) { //Skip silhouette point because it's out of range
			continue;
		}
		else {
			keep_indices.push_back(i);
		}
	}

	Eigen::MatrixXd result;
	Eigen::VectorXi keep_row_idx, keep_col_idx(3);
	keep_row_idx = Eigen::VectorXi::Map(keep_indices.data(), keep_indices.size());
	keep_col_idx.col(0) << 0, 1, 2;
	igl::slice(silhouette_vertices, keep_row_idx, keep_col_idx, result);
	silhouette_vertices = result;
}

/** Generates a mesh inside the loop3D vertices. **/
void MeshExtrusion::generate_mesh(Mesh& m, Eigen::MatrixXd loop3D, Eigen::Vector3d center, Eigen::Vector3d x_vec, Eigen::Vector3d y_vec, Eigen::Vector3d offset, int nr_silhouette_vert, vector<int> loop_base_original_indices, vector<int> sil_original_indices) {
	int size_before_gen = m.V.rows();
	Eigen::MatrixXd loop2D(loop3D.rows(), 2);
	Eigen::Vector3d vec;
	Eigen::MatrixXi loop_stroke_edges(loop3D.rows(), 2);

	for(int i = 0; i < loop3D.rows(); i++) {
		vec = 1000*(loop3D.row(i).transpose() - center);
		loop2D.row(i) << vec.dot(x_vec), vec.dot(y_vec);
		loop_stroke_edges.row(i) << i, ((i + 1) % loop3D.rows());
	}
	double mean_squared_sample_dist = 0.0;
	for (int i = 0; i < loop2D.rows(); i++) {
		mean_squared_sample_dist += (loop2D.row(i) - loop2D.row((i + 1) % loop2D.rows())).squaredNorm();
	}
	mean_squared_sample_dist /= loop2D.rows();

	Eigen::MatrixXd V2;
	Eigen::MatrixXi F2;
	Eigen::MatrixXi vertex_markers, edge_markers;
	igl::triangle::triangulate(loop2D, loop_stroke_edges, Eigen::MatrixXd(0, 0), Eigen::MatrixXi::Constant(loop2D.rows(), 1, 1), Eigen::MatrixXi::Constant(loop_stroke_edges.rows(), 1, 1), "Yq25Qa" + to_string(0.5*mean_squared_sample_dist), V2, F2, vertex_markers, edge_markers); //Capital Q silences triangle's output in cmd line. Also retrieves markers to indicate whether or not an edge/vertex is on the mesh boundary
	V2 /= 1000.0;

	Eigen::RowVector3d vert;
	for(int i = 0; i < V2.rows(); i++) {
		if(i >= loop2D.rows()) { //Interior vertices
			vert = center;
			vert += x_vec.transpose()*V2(i, 0);
			vert += y_vec.transpose()*V2(i, 1);
			vert += offset.transpose();
			m.V.conservativeResize(m.V.rows() + 1, Eigen::NoChange);
			m.V.row(m.V.rows() - 1) = vert;
			m.vertex_is_fixed.conservativeResize(m.vertex_is_fixed.rows() + 1);
			m.vertex_is_fixed(m.vertex_is_fixed.rows() - 1) = 0;
		}
	}
	update_face_indices(m, F2, sil_original_indices, loop_base_original_indices, nr_silhouette_vert, size_before_gen, loop2D.rows());
}

/** Gets the center point in the extrusion base and the normal of its average plane. **/
void MeshExtrusion::get_normal_and_center(Eigen::RowVector3d& center, Eigen::Vector3d& normal, Mesh& m, Eigen::VectorXi &boundary_vertices) {
	for(int i = 0; i < boundary_vertices.rows(); i++) {
		center += m.V.row(boundary_vertices[i]);
	}
	center /= boundary_vertices.rows();

	Eigen::Vector3d vec0, vec1;
	for(int i = 0; i < boundary_vertices.rows(); i++) {
		vec0 = m.V.row(boundary_vertices[i]) - center;
		vec1 = m.V.row(boundary_vertices[(i + 1) % boundary_vertices.rows()]) - center;
		normal += vec0.cross(vec1);
	}
	normal.normalize();
}

/** Finds the leftmost and rightmost vertex of the extrusion base (according to some local set of axes). **/
void MeshExtrusion::find_left_and_right(int& most_left_vertex_idx, int& most_right_vertex_idx, Mesh& m, Eigen::VectorXi &boundary_vertices, const Eigen::RowVector3d& sil_start, const Eigen::RowVector3d& sil_end, Eigen::RowVector3d& center) {
	Eigen::Vector3d dir = sil_end - sil_start;
	dir.normalize();
	Eigen::Vector3d center_to_vertex;
	double dot_prod = dir.dot((m.V.row(boundary_vertices[0]) - center).normalized());
	double min = dot_prod, max = dot_prod;
	most_left_vertex_idx = 0;
	most_right_vertex_idx = 0;

	for (int i = 1; i < boundary_vertices.rows(); i++) {
		center_to_vertex = (m.V.row(boundary_vertices[i]) - center).normalized();
		dot_prod = dir.dot(center_to_vertex);
		if (dot_prod > max) {
			max = dot_prod;
			most_left_vertex_idx = i;
		}
		else if (dot_prod < min) {
			min = dot_prod;
			most_right_vertex_idx = i;
		}
	}
}

/** Adds the silhouette vertices to the mesh and also inserts their tracking variables. **/
vector<int> MeshExtrusion::add_silhouette_vertices(Mesh& m, int stroke_ID, Eigen::MatrixXd& silhouette_vertices, Eigen::MatrixXi& added_edges) {
	vector<int> sil_original_indices;
	int size_before_silhouette = m.V.rows();
	m.V.conservativeResize(m.V.rows() + silhouette_vertices.rows(), Eigen::NoChange);
	m.vertex_is_fixed.conservativeResize(m.vertex_is_fixed.rows() + silhouette_vertices.rows());

	for (int i = 0; i < silhouette_vertices.rows(); i++) {
		m.V.row(size_before_silhouette + i) = silhouette_vertices.row(i);
		m.vertex_is_fixed[size_before_silhouette + i] = 1;
		sil_original_indices.push_back(size_before_silhouette + i);
	}

	for (int i = 0; i < silhouette_vertices.rows() - 1; i++) {
		added_edges.conservativeResize(added_edges.rows() + 1, Eigen::NoChange);
		added_edges.bottomRows(1) << size_before_silhouette + i, size_before_silhouette + i + 1, stroke_ID, 0;
	}
	return sil_original_indices;
}

/** Creates a loop consisting of the silhouette vertices and half of the base stroke vertices. **/
void MeshExtrusion::create_loop(Mesh& m, Eigen::MatrixXd& loop3D, Eigen::VectorXi& boundary_vertices, vector<int> &loop_base_original_indices, int start_idx, int end_idx) {
	int idx = start_idx;
	while (true) {
		loop3D.conservativeResize(loop3D.rows() + 1, Eigen::NoChange);
		loop3D.row(loop3D.rows() - 1) = m.V.row(boundary_vertices[idx]);
		loop_base_original_indices.push_back(boundary_vertices[idx]);
		if (idx == end_idx) {
			break;
		}
		idx++;
		if (idx == boundary_vertices.rows()) {
			idx = 0;
		}
	}
}

/** Updates the old indicators for sharp edges and edge_boundary_markers after the mesh topology changed and inserts newly generated sharp edges & boundary edges. **/
void MeshExtrusion::update_edge_indicators(Mesh& m, Eigen::MatrixXi& edges_to_update) {
	if (!igl::is_edge_manifold(m.F)) {
		throw - 1;
		return;
	}
	Eigen::MatrixXi EV, FE, EF;
	igl::edge_topology(m.V, m.F, EV, FE, EF);
	m.sharp_edge.resize(EV.rows());
	m.sharp_edge.setZero();
	m.edge_boundary_markers.resize(EV.rows());
	m.edge_boundary_markers.setZero();

	int start, end, equal_pos;
	Eigen::VectorXi col1Equals, col2Equals;
	for (int i = 0; i < edges_to_update.rows(); i++) {
		start = edges_to_update(i, 0);
		end = edges_to_update(i, 1); 

		col1Equals = EV.col(0).cwiseEqual(min(start, end)).cast<int>();
		col2Equals = EV.col(1).cwiseEqual(max(start, end)).cast<int>();
		(col1Equals + col2Equals).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge

		m.edge_boundary_markers[equal_pos] = edges_to_update(i, 2);
		m.sharp_edge[equal_pos] = edges_to_update(i, 3);
	}
}

/** Updates the old face indices with the new indexing and inserts new faces that are made in triangulation. **/
void MeshExtrusion::update_face_indices(Mesh& m, Eigen::MatrixXi& F2, vector<int> sil_original_indices, vector<int> loop_base_original_indices, int nr_silhouette_vert, int size_before_gen, int loop2D_size) {
	int original_F_size = m.F.rows();
	m.F.conservativeResize(m.F.rows() + F2.rows(), Eigen::NoChange);
	int vert_idx_in_mesh;
	for (int i = 0; i < F2.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			if (F2(i, j) < nr_silhouette_vert) { //Silhouette vertex
				vert_idx_in_mesh = sil_original_indices[F2(i, j)];
			}
			else if (F2(i, j) < loop2D_size) { //Extrusion base vertex
				vert_idx_in_mesh = loop_base_original_indices[F2(i, j) - nr_silhouette_vert];
			}
			else {
				vert_idx_in_mesh = size_before_gen + F2(i, j) - loop2D_size; //Interior point
			}
			m.F(original_F_size + i, j) = vert_idx_in_mesh;
		}
	}
}

/** Updates the base stroke's 3DPoints with the new vertices. PathElements of the type FACE stay in the same position, but PathElements of the type EDGE get added. **/
void MeshExtrusion::post_extrude_prepare_update_points(Stroke& stroke, SurfacePath& surface_path) {
	vector<PathElement> path = surface_path.get_path();
	Eigen::MatrixX3d new_3DPoints(path.size(), 3);

	for (int i = 0; i < path.size(); i++) {
		new_3DPoints.row(i) = path[i].get_vertex().transpose();
	}

	stroke.set3DPoints(new_3DPoints);
}

/** Updates the silhouette stroke's 3DPoints with the new vertices. **/
void MeshExtrusion::post_extrude_main_update_points(Stroke &stroke, Eigen::MatrixXd new_positions) {
	Eigen::MatrixX3d new_3DPoints(new_positions.rows() + 1, 3);
	for (int i = 0; i < new_positions.rows(); i++) {
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
