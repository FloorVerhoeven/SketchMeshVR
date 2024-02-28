#include "LaplacianRemesh.h"
#include <igl/edge_topology.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/slice.h>
#include <igl/adjacency_list.h>
#include <igl/project.h>
#include <igl/cat.h>


using namespace std;
using namespace igl;

bool LaplacianRemesh::is_front_loop = true; //Determines extrusion or cut
bool LaplacianRemesh::remove_inside_faces = true;

Eigen::MatrixXi LaplacianRemesh::EV, LaplacianRemesh::FE, LaplacianRemesh::EF;
vector<vector<int>> LaplacianRemesh::VV;

Eigen::VectorXi LaplacianRemesh::remesh_cut_remove_inside(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, bool& remesh_success) {
	is_front_loop = false;
	remove_inside_faces = true;
	adjacency_list(m.F, VV);
	return remesh(m, surface_path, model, view, proj, viewport, remesh_success);
}

Eigen::VectorXi LaplacianRemesh::remesh_extrusion_remove_inside(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, bool& remesh_success) {
	is_front_loop = true;
	remove_inside_faces = true;
	adjacency_list(m.F, VV);
	return remesh(m, surface_path, model, view, proj, viewport, remesh_success);
}

Eigen::VectorXi LaplacianRemesh::remesh(Mesh& m, SurfacePath& surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, bool& remesh_success) {
	vector<bool> dirty_face(m.F.rows());
	vector<int> dirty_vertices(m.V.rows());
	m.new_mapped_indices.resize(m.V.rows()); //Resize the map from old to new (clean) vertex indices to allow it to contain the number of vertices that are in the mesh at the start

	vector<int> sharp_edge_indices;
	for(int i = 0; i < m.sharp_edge.rows(); i++) {
		if(m.sharp_edge[i]) {
			sharp_edge_indices.push_back(i);
		}
	}

	vector<vector<int>> VF, VI;
	igl::edge_topology(m.V, m.F, EV, FE, EF);
	igl::vertex_triangle_adjacency(m.V.rows(), m.F, VF, VI);

	Eigen::MatrixXi sharpEV;
	Eigen::VectorXi sharpEV_row_idx, sharpEV_col_idx(2);
	sharpEV_row_idx = Eigen::VectorXi::Map(sharp_edge_indices.data(), sharp_edge_indices.size());
	sharpEV_col_idx.col(0) << 0, 1;
	igl::slice(EV, sharpEV_row_idx, sharpEV_col_idx, sharpEV); //Keep only the sharp edges in the original mesh

	int face = -1, edge = -1;
	vector<PathElement> path = surface_path.get_path();
	if(!is_front_loop) { //Append first element at the end in case of cutting, to make it a loop and compatible with extrusion's loop
		path.push_back(path[0]);
	}

	//Find outside and inside vertices. +1 means that the vertex will stay, -1 means it will be destroyed
	for(int i = 0; i < path.size(); i++) {
		if(path[i].get_type() == PathElement::FACE) {
			face = path[i].get_ID();
		} else {
			edge = path[i].get_ID();
			if(EF(edge, 0) == face) {
				dirty_vertices[EV(edge, 0)] += 1;
				dirty_vertices[EV(edge, 1)] -= 1;
			} else {
				dirty_vertices[EV(edge, 0)] -= 1;
				dirty_vertices[EV(edge, 1)] += 1;
			}
			face = (EF(edge, 0) == face) ? EF(edge, 1) : EF(edge, 0); //get the polygon on the other side of the edge
		}
	}

	//Collect vertices on the boundary
	vector<int> outer_boundary_vertices;
	vector<int> inner_boundary_vertices;
	for(int i = 0; i < m.V.rows(); i++) {
		if(dirty_vertices[i] < 0) { //vertices that are removed
			inner_boundary_vertices.push_back(i);
		} else if(dirty_vertices[i] > 0) { //vertices that stay
			outer_boundary_vertices.push_back(i);
		}
	}

	if (outer_boundary_vertices.size() == 0) { //There are no vertices left behind (cut that removes everthing or extrude that doesn't include at least 1 vertex)
		remesh_success = false;
		return Eigen::VectorXi::Zero(1);
	}

	//Collect faces along the path
	for(int i = 0; i < path.size(); i++) {
		if(path[i].get_type() == PathElement::EDGE) {
			dirty_face[EF(path[i].get_ID(), 0)] = true;
			dirty_face[EF(path[i].get_ID(), 1)] = true;
		}
	}

	if(remove_inside_faces) {
		for(int i = 0; i < m.V.rows(); i++) {
			if(dirty_vertices[i] < 0) {
				for(int j = 0; j < VF[i].size(); j++) {
					if(!dirty_face[VF[i][j]]) {
						propagate_dirty_faces(VF[i][j], dirty_face);
					}
				}
			}
		}
	}

	//Remove dirty faces
	vector<int> clean_faces;
	for(int i = 0; i < m.F.rows(); i++) {
		if(!dirty_face[i]) {
			clean_faces.push_back(i);
		}
	}

	if (clean_faces.size() == 0) { //There are no faces left behind (cut that removes everthing or extrude that doesn't include at least 1 vertex)
		remesh_success = false;
		return Eigen::VectorXi::Zero(1);
	}

	//Determine which vertices are clean
	Eigen::VectorXi vertex_is_clean = Eigen::VectorXi::Zero(m.V.rows());
	vector<int> clean_vertices;
	for(int i = 0; i < m.V.rows(); i++) {
		if(dirty_vertices[i] == 0) { //Vertices with a - value are already known to be removed and vertices with a + value are already known to stay
			for(int j = 0; j < VF[i].size(); j++) {
				if(!dirty_face[VF[i][j]]) { //If the vertex is adjacent to a clean face, then the vertex has to be clean
					vertex_is_clean(i) = 1;
					clean_vertices.push_back(i);
					break; //Can skip checking all other adjacent faces
				}
			}
		} else if(dirty_vertices[i] > 0) {
			vertex_is_clean(i) = 1;
			clean_vertices.push_back(i);
		}
	}


	Eigen::MatrixXi tmp_F;
	Eigen::VectorXi row_idx, col_idx(3);
	row_idx = Eigen::VectorXi::Map(clean_faces.data(), clean_faces.size());
	col_idx.col(0) << 0, 1, 2;
	igl::slice(m.F, row_idx, col_idx, tmp_F); //Keep only the clean faces in the mesh
	m.F = tmp_F;

	//NOTE: Output from sort_boundary_vertices is not necessarily counter-clockwise
	outer_boundary_vertices = sort_boundary_vertices(path[0].get_vertex(), outer_boundary_vertices, m);
	if(!remove_inside_faces) {
		inner_boundary_vertices = sort_boundary_vertices(path[0].get_vertex(), inner_boundary_vertices, m);
	}

	Eigen::RowVector3d mean_viewpoint = compute_mean_viewpoint(m, inner_boundary_vertices);
	
	Eigen::MatrixXd tmp_V;
	Eigen::VectorXi tmp_part_of, tmp_markers;
	int size_before_removing = m.V.rows();
	row_idx = Eigen::VectorXi::Map(clean_vertices.data(), clean_vertices.size());
	col_idx.col(0) << 0, 1, 2;
	igl::slice(m.V, row_idx, col_idx, tmp_V);
	m.V = tmp_V;

	col_idx.resize(1);
	col_idx.col(0) << 0;
	igl::slice(m.part_of_original_stroke, row_idx, col_idx, tmp_part_of);
	igl::slice(m.vertex_boundary_markers, row_idx, col_idx, tmp_markers);
	m.part_of_original_stroke = tmp_part_of; //Compact part_of_original_stroke by removing the values for vertices that are now removed
	m.vertex_boundary_markers = tmp_markers; //Compact vertex_boundary_markers by removing the values for vertices that are now removed


	//Map indices in m.V at the start of remeshing to indices in m.V after remeshing
	int count = 0;
	m.new_mapped_indices.setConstant(m.new_mapped_indices.rows(), -1); //Initialize all values to -1
	for(int i = 0; i < vertex_is_clean.rows(); i++) {
		if(vertex_is_clean[i]) {
			m.new_mapped_indices[i] = count; 
			count++;
		}
	}
	
	update_face_indices(m);


	vector<int> path_vertices;
	int original_V_size = m.V.rows();
	for(int i = 0; i < path.size() - 1; i++) { //Do not use the last path vertex, as it is a copy of the first and creates unwanted behaviour
		path_vertices.push_back(original_V_size + i); //Store the vertex indices (in m.V) of the path vertices
	}


	update_mesh_values(m, path, surface_path.get_origin_stroke_ID(), vertex_is_clean.rows());

	Eigen::MatrixXi added_sharpEV(path.size() - 1, 2);
	for(int i = 0; i < path.size() - 1; i++) {
		surface_path.get_path_element(i).set_v_idx(size_before_removing + i);
		added_sharpEV.row(i) << vertex_is_clean.rows() + i, vertex_is_clean.rows() + ((i + 1) + path.size() - 1) % (path.size() - 1); //Use mod (path.size() - 1) to compensate for start point that is in path double
	}

	//Update the outer_boundary_vertices (which are indices into V before it was sliced to keep only the clean vertices)
	for(int i = 0; i < outer_boundary_vertices.size(); i++) {
		outer_boundary_vertices[i] = m.new_mapped_indices[outer_boundary_vertices[i]];
	}

	//Make the outer_boundary_vertices CCW
	Eigen::Matrix4f modelview = view * model;
	Eigen::VectorXi row_idx2, col_idx2(3);
	row_idx2 = Eigen::VectorXi::Map(outer_boundary_vertices.data(), outer_boundary_vertices.size());
	col_idx2.col(0) << 0, 1, 2;
	igl::slice(m.V, row_idx2, col_idx2, tmp_V); //Keep only the clean "boundary" vertices in the mesh
	if(!is_counter_clockwise_boundaries(tmp_V, modelview, proj, viewport, mean_viewpoint, !is_front_loop)) {
		reverse(outer_boundary_vertices.begin(), outer_boundary_vertices.end());
	}


	stitch(path_vertices, outer_boundary_vertices, m);
	if(!remove_inside_faces) {
		reverse_path(path_vertices);
		stitch(path_vertices, inner_boundary_vertices, m);
	}


	Eigen::MatrixXi sharpEVcat = igl::cat(1, sharpEV, added_sharpEV);
	update_sharp_edges(m, sharpEVcat);

	return Eigen::VectorXi::Map(path_vertices.data(), path_vertices.size());
}

/** Compute the mean of the boundary vertices on the side of the mesh that is removed.
	In the case of CUT this is used as a viewpoint to determine the orientation (CW or CCW) of the boundary vertices that remain.
	Do this before slicing m.V and removing the vertices that we need. **/
Eigen::RowVector3d LaplacianRemesh::compute_mean_viewpoint(Mesh &m, vector<int> inner_boundary_vertices) {
	Eigen::MatrixXd removed_V;
	Eigen::VectorXi row_idx_removed = Eigen::VectorXi::Map(inner_boundary_vertices.data(), inner_boundary_vertices.size());
	Eigen::VectorXi col_idx_removed(3);
	col_idx_removed.col(0) << 0, 1, 2;
	igl::slice(m.V, row_idx_removed, col_idx_removed, removed_V);

	return removed_V.colwise().mean();
} 

/** Updates the faces with new vertex indexing. **/
void LaplacianRemesh::update_face_indices(Mesh &m) {
	for(int i = 0; i < m.F.rows(); i++) {
		for(int j = 0; j < 3; j++) {
			m.F(i, j) = m.new_mapped_indices[m.F(i, j)];
		}
	}
}

/** Adds values for a bunch of tracking variables for the new mesh vertices. **/
void LaplacianRemesh::update_mesh_values(Mesh& m, vector<PathElement> path, int stroke_ID, int new_mapped_start) {
	int size_before = m.V.rows();
	m.V.conservativeResize(m.V.rows() + path.size() - 1, Eigen::NoChange);
	m.part_of_original_stroke.conservativeResize(m.part_of_original_stroke.rows() + path.size() - 1);
	m.vertex_boundary_markers.conservativeResize(m.vertex_boundary_markers.rows() + path.size() - 1);
	m.new_mapped_indices.conservativeResize(m.new_mapped_indices.rows() + path.size() - 1, Eigen::NoChange);
	for(int i = 0; i < path.size() - 1; i++) {
		m.V.row(size_before + i) << path[i].get_vertex().transpose();
		m.part_of_original_stroke[size_before + i] = 0;
		m.vertex_boundary_markers[size_before + i] = stroke_ID;
		m.new_mapped_indices(new_mapped_start + i) = size_before + i;
	}
}

/** Updates the old sharp edges after the mesh topology changed and inserts newly generated sharp edges. **/
void LaplacianRemesh::update_sharp_edges(Mesh& m, Eigen::MatrixXi& all_sharpEV) {
	igl::edge_topology(m.V, m.F, EV, FE, EF);
	m.sharp_edge.resize(EV.rows());
	m.sharp_edge.setZero();

	int start, end, equal_pos;
	Eigen::VectorXi col1Equals, col2Equals;
	for(int i = 0; i < all_sharpEV.rows(); i++) {
		start = m.new_mapped_indices(all_sharpEV(i, 0));
		end = m.new_mapped_indices(all_sharpEV(i, 1));
		if(start == -1 || end == -1) { //Sharp edge no longer exists
			continue;
		}

		col1Equals = EV.col(0).cwiseEqual(min(start, end)).cast<int>();
		col2Equals = EV.col(1).cwiseEqual(max(start, end)).cast<int>();
		(col1Equals + col2Equals).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge

		m.sharp_edge[equal_pos] = 1; //Set this edge to be sharp
	}
}

/** Floodfills all dirty faces. **/
void LaplacianRemesh::propagate_dirty_faces(int face, vector<bool>& dirty_face) {
	dirty_face[face] = true;
	for(int i = 0; i < 3; i++) {
		int next_face = (EF(FE(face, i), 0) == face) ? EF(FE(face, i), 1) : EF(FE(face, i), 0); //get the adjacent polygon
		if(!dirty_face[next_face]) {
			propagate_dirty_faces(next_face, dirty_face);
		}
	}
}

/** Sorts the boundary vertices such that they make up a continuous path that starts with the vertex that's closest to start_vertex. 
	Takes care of backtracking over original stroke boundary vertices that is sometimes necessary. **/
vector<int> LaplacianRemesh::sort_boundary_vertices(Eigen::Vector3d start_vertex, vector<int> boundary_vertices, Mesh& m) {
	igl::edge_topology(m.V, m.F, EV, FE, EF);

	vector<int> sorted_boundary_vertices;
	int start_v = find_closest(boundary_vertices, start_vertex, m);
	int v = start_v;
	int prev = start_v;

	Eigen::VectorXi tmp1, tmp2;
	int min_idx, max_idx, equal_pos, v_pos, max_val, only_option_idx;
	bool have_one_option = false;
	while(true) {
		sorted_boundary_vertices.push_back(v);
		for(int i = 0; i < VV[v].size(); i++) {
			if(std::find(boundary_vertices.begin(), boundary_vertices.end(), VV[v][i]) != boundary_vertices.end()) { //Only continue if the adjacent vertex is also a boundary vertex
				if(v < VV[v][i]) {
					min_idx = v;
					v_pos = 0;
					max_idx = VV[v][i];
				} else {
					min_idx = VV[v][i];
					v_pos = 1;
					max_idx = v;
				}
				tmp1 = EV.col(0).cwiseEqual(min_idx).cast<int>();
				tmp2 = EV.col(1).cwiseEqual(max_idx).cast<int>();
				max_val = (tmp1 + tmp2).maxCoeff(&equal_pos); //Tmp1 and tmp2 will contain a 1 on the positions where they equal the min and max_idx. When adding them, we'll get 2 on the row that contains both of them
				if((max_val != 2 && VV[v][i] != prev) || (max_val == 2 && EF(equal_pos, v_pos) == -1)) { //The first condition means that both vertices were only connected to faces that have been removed, so it is a valid move. The second condition means that the edge has a NULL face on the correct side of the edge, so proceed to adjacent vertex
					prev = v;
					v = VV[v][i];
					have_one_option = false; //Remove our backtracking option
					break;
				} else if(max_val != 2 && VV[v][i] == prev) { //Check if we have any other options before we start backtracking over a boundary
					only_option_idx = VV[v][i]; //Remember the vertex that we could backtrack to
					have_one_option = true;
					continue;
				}
			}
		}

		if(have_one_option) { //If we still have the backtracking vertex buffered, use it.
			prev = v;
			v = only_option_idx;
			have_one_option = false;
		}

		if(v == start_v) {
			break;
		}
	}
	return sorted_boundary_vertices;
}

/** Finds the closest vertex to point base. **/
int LaplacianRemesh::find_closest(vector<int> vertices, Eigen::Vector3d base, Mesh& m) {
	int closest = vertices[0];
	double min = (m.V.row(vertices[0]) - base.transpose()).norm();
	double d;

	for(int i = 1; i < vertices.size(); i++) {
		d = (m.V.row(vertices[i]) - base.transpose()).norm();
		if(d < min) {
			min = d;
			closest = vertices[i];
		}
	}
	return closest;
}

/** Stitches the boundary vertices (that are old, already existing mesh vertices) to the path vertices of the newly drawn stroke. **/
void LaplacianRemesh::stitch(std::vector<int> path_vertices, std::vector<int> boundary_vertices, Mesh& m) {
	int path_idx = 0, outer_idx = 0;

	boundary_vertices = reorder(boundary_vertices, m.V.row(path_vertices[0]), m);

	Eigen::RowVector3d start_path_v = m.V.row(path_vertices[0]), start_outer_v = m.V.row(boundary_vertices[0]);
	Eigen::RowVector3d path_v = start_path_v, outer_v = start_outer_v, next_path_v, next_outer_v;
	int path_v_idx = path_vertices[0], outer_v_idx = boundary_vertices[0], next_path_v_idx, next_outer_v_idx;

	bool proceed_outer_v;

	while(true) {
		next_path_v = m.V.row(path_vertices[(path_idx + 1) % path_vertices.size()]);
		next_outer_v = m.V.row(boundary_vertices[(outer_idx + 1) % boundary_vertices.size()]);
		next_path_v_idx = path_vertices[(path_idx + 1) % path_vertices.size()];
		next_outer_v_idx = boundary_vertices[(outer_idx + 1) % boundary_vertices.size()];

		proceed_outer_v = true;
		if(outer_idx == boundary_vertices.size()) {
			proceed_outer_v = false;
		} else if(path_idx == path_vertices.size()) {
			proceed_outer_v = true;
		} else if((path_v - next_outer_v).norm() < (outer_v - next_path_v).norm()) {
			proceed_outer_v = true;
		} else {
			proceed_outer_v = false;
		}

		//Add faces
		if(proceed_outer_v) {
			m.F.conservativeResize(m.F.rows() + 1, Eigen::NoChange);
			m.F.row(m.F.rows() - 1) << path_v_idx, outer_v_idx, next_outer_v_idx;
			outer_v = next_outer_v;
			outer_v_idx = next_outer_v_idx;
			outer_idx++;
		} else {
			m.F.conservativeResize(m.F.rows() + 1, Eigen::NoChange);
			m.F.row(m.F.rows() - 1) << next_path_v_idx, path_v_idx, outer_v_idx;
			path_v = next_path_v;
			path_v_idx = next_path_v_idx;
			path_idx++;
		}

		if(path_idx == path_vertices.size() && outer_idx == boundary_vertices.size()) {
			break;
		}
	}
}

/** Reorders the boundary vertices so they start from the vertex closest to start_v. Keeps their original order **/
vector<int> LaplacianRemesh::reorder(vector<int> boundary_vertices, Eigen::Vector3d start_v, Mesh& m) {
	int index = 0;
	double min = INFINITY;
	for(int i = 0; i < boundary_vertices.size(); i++) {
		Eigen::Vector3d vert = m.V.row(boundary_vertices[i]);
		double d = (vert - start_v).norm();
		if(d < min) {
			min = d;
			index = i;
		}
	}

	vector<int> reordered;
	for(int i = 0; i < boundary_vertices.size(); i++) {
		reordered.push_back(boundary_vertices[(i + index) % boundary_vertices.size()]);
	}
	return reordered;
}

void LaplacianRemesh::reverse_path(vector<int> path_vertices) {
	reverse(path_vertices.begin(), path_vertices.end());
	int v = path_vertices[path_vertices.size() - 1];
	path_vertices.insert(path_vertices.begin() + 1, path_vertices.begin(), path_vertices.end() - 1);
	path_vertices[0] = v;
}

/** For cut we determine if the boundary is CCW by comparing the direction of the normal of the remaining boundary points against the vector from the center of the removed boundary points to the center of the remaining boundary points.
	For extrusion we simply look at the sign of the inner area of the projected points (since we're looking from the "correct" viewpoint)
**/
bool LaplacianRemesh::is_counter_clockwise_boundaries(Eigen::MatrixXd boundary_points, Eigen::Matrix4f modelview, Eigen::Matrix4f proj, Eigen::Vector4f viewport, Eigen::RowVector3d mean_viewpoint, bool cut) {
	if(cut) {
		Eigen::RowVector3d center = boundary_points.colwise().mean();
		Eigen::Vector3d normal(0, 0, 0);
		Eigen::Vector3d vec0, vec1;
		for(int i = 0; i < boundary_points.rows(); i++) {
			vec0 = boundary_points.row(i) - center;
			vec1 = boundary_points.row((i + 1) % boundary_points.rows()) - center;
			normal += vec1.cross(vec0);
		}
		normal.normalize();

		if(normal.dot(center - mean_viewpoint) < 0) {
			return false;
		} else {
			return true;
		}
	} else {
		double total_area = 0.0;
		Eigen::RowVector3d pt, vert;
		Eigen::RowVector2d prev, next;
		vert = boundary_points.row(boundary_points.rows() - 1);
		igl::project(vert, modelview, proj, viewport, pt);
		prev = pt.leftCols(2);

		for(int i = 0; i < boundary_points.rows(); i++) {
			vert = boundary_points.row(i);
			igl::project(vert, modelview, proj, viewport, pt);
			next = pt.leftCols(2);
			total_area += (prev[1] + next[1]) * (next[0] - prev[0]);
			prev = next;
		}

		if(total_area > 0) {
			return false;
		}

		return true;
	}
}