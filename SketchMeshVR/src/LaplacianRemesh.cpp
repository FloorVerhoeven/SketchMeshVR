#include "LaplacianRemesh.h"
#include "CleanStroke3D.h"
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

Eigen::VectorXi LaplacianRemesh::remesh_cut_remove_inside(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, int& remesh_success, int cut_clicked_face, Eigen::MatrixXi& replacing_vertex_bindings) {
	is_front_loop = false;
	remove_inside_faces = true;
	adjacency_list(m.F, VV);
	return remesh(m, surface_path, model, view, proj, viewport, remesh_success, cut_clicked_face, replacing_vertex_bindings);
}

Eigen::VectorXi LaplacianRemesh::remesh_extrusion_remove_inside(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, int& remesh_success, Eigen::MatrixXi& replacing_vertex_bindings) {
	is_front_loop = true;
	remove_inside_faces = true;
	adjacency_list(m.F, VV);
	return remesh(m, surface_path, model, view, proj, viewport, remesh_success, -1, replacing_vertex_bindings);
}

bool LaplacianRemesh::remesh_open_path(Mesh& m, Stroke& open_path_stroke, Eigen::MatrixXi& replacing_vertex_bindings, igl::opengl::glfw::Viewer &viewer) {
	replacing_vertex_bindings.resize(0, 4);
	remove_inside_faces = false;
	SurfacePath surface_path;
	surface_path.create_from_open_path(open_path_stroke);

	Eigen::MatrixXi start_F = m.F;
	Eigen::MatrixXd start_V = m.V;
	Eigen::VectorXi start_edge_boundary_markers = m.edge_boundary_markers;
	Eigen::VectorXi start_vertex_is_fixed = m.vertex_is_fixed;
	Eigen::VectorXi start_sharp_edge = m.sharp_edge;
	Eigen::VectorXi start_new_mapped_indices = m.new_mapped_indices;

	vector<bool> dirty_face(m.F.rows());
	Eigen::VectorXi dirty_vertices(m.V.rows());
	dirty_vertices.setZero();

	int face = -1, edge = -1;
	vector<PathElement> path = surface_path.get_path();

	m.new_mapped_indices.resize(m.V.rows()); //Resize the map from old to new (clean) vertex indices to allow it to contain the number of vertices that are in the mesh at the start

	igl::edge_topology(m.V, m.F, EV, FE, EF);
	Eigen::MatrixXi startEV = EV;
	adjacency_list(m.F, VV);

	Eigen::MatrixXi original_sharp_or_boundary_edges(0, 4);
	for (int i = 0; i < m.sharp_edge.rows(); i++) {
		if (m.sharp_edge[i] || m.edge_boundary_markers[i]) {
			original_sharp_or_boundary_edges.conservativeResize(original_sharp_or_boundary_edges.rows() + 1, Eigen::NoChange);
			original_sharp_or_boundary_edges.bottomRows(1) << EV(i, 0), EV(i, 1), m.edge_boundary_markers[i], m.sharp_edge[i];
		}
	}

	for (int i = 0; i < path.size(); i++) {
		if (path[i].get_type() == PathElement::FACE) {
			face = path[i].get_ID();
			dirty_vertices[m.F(face, 0)] += 1; //These will ensure that all vertices of the end faces will also be in outer_boundary_vertices
			dirty_vertices[m.F(face, 1)] += 1;
			dirty_vertices[m.F(face, 2)] += 1;
		}
		else {
			edge = path[i].get_ID();
			dirty_vertices[EV(edge, 0)] += 1;
			dirty_vertices[EV(edge, 1)] += 1;
			dirty_face[EF(edge, 0)] = true;
			dirty_face[EF(edge, 1)] = true;
			face = (EF(edge, 0) == face) ? EF(edge, 1) : EF(edge, 0); //get the polygon on the other side of the edge
		}
	}

	//Collect vertices on the boundary
	vector<int> outer_boundary_vertices;
	for (int i = 0; i < m.V.rows(); i++) {
		if (dirty_vertices[i] > 0) {
			outer_boundary_vertices.push_back(i);
		}
		m.new_mapped_indices[i] = i;
	}

	//Remove dirty faces
	vector<int> clean_faces;
	for (int i = 0; i < m.F.rows(); i++) {
		if (!dirty_face[i]) {
			clean_faces.push_back(i);
		}
	}

	Eigen::MatrixXi tmp_F;
	Eigen::MatrixXi prev_F = m.F;
	Eigen::VectorXi row_idx, col_idx(3);
	row_idx = Eigen::VectorXi::Map(clean_faces.data(), clean_faces.size());
	col_idx.col(0) << 0, 1, 2;
	igl::slice(m.F, row_idx, col_idx, tmp_F); //Keep only the clean faces in the mesh
	m.F = tmp_F;

	//NOTE: Output from sort_boundary_vertices is not necessarily counter-clockwise
	try {
		outer_boundary_vertices = sort_boundary_vertices(path[0].get_vertex(), outer_boundary_vertices, m);
	}
	catch (int ex) {
		if (ex == -1) {
			std::cerr << "Could not process adding stroke. Please try again. " << std::endl;
			m.F = prev_F;
			return false;
		}
		else {
			std::cout << "Something weird happened. Check what's going on." << std::endl;
			m.F = prev_F;
			return false;
		}
	}

	Eigen::MatrixXi replacing_edges(0, 4); //Contains edge start, edge end, edge_marker and edge_sharp
	Eigen::MatrixXd edge_split_positions(0, 3);

	for (int i = 0; i < path.size(); i++) {
		if (path[i].get_type() == PathElement::EDGE && (m.sharp_edge[path[i].get_ID()] || m.edge_boundary_markers[path[i].get_ID()])) {
			path[i].fixed = true;
			replacing_edges.conservativeResize(replacing_edges.rows() + 2, Eigen::NoChange);
			replacing_edges(replacing_edges.rows() - 2, 1) = startEV(path[i].get_ID(), 0); //Pair of original edge start index and middle vertex index (middle vertex index gets set after resampling)
			replacing_edges(replacing_edges.rows() - 2, 2) = m.edge_boundary_markers[path[i].get_ID()];
			replacing_edges(replacing_edges.rows() - 2, 3) = m.sharp_edge[path[i].get_ID()];
			replacing_edges(replacing_edges.rows() - 1, 1) = startEV(path[i].get_ID(), 1); //Pair of original edge end index and middle vertex index (middle vertex index gets set after resampling)
			replacing_edges(replacing_edges.rows() - 1, 2) = m.edge_boundary_markers[path[i].get_ID()];
			replacing_edges(replacing_edges.rows() - 1, 3) = m.sharp_edge[path[i].get_ID()];
			edge_split_positions.conservativeResize(edge_split_positions.rows() + 1, Eigen::NoChange);
			edge_split_positions.bottomRows(1) << path[i].get_vertex().transpose();
			replacing_vertex_bindings.conservativeResize(replacing_vertex_bindings.rows() + 1, Eigen::NoChange);
			replacing_vertex_bindings.bottomRows(1) << m.edge_boundary_markers[path[i].get_ID()], startEV(path[i].get_ID(), 0), startEV(path[i].get_ID(), 1), -1; //Middle vertex index gets set after resampling
		}
		else {
			path[i].fixed = false;
		}
	}



	Eigen::RowVector3d mean_viewpoint = compute_mean_viewpoint(m, outer_boundary_vertices);
	Eigen::Matrix4f modelview = viewer.oculusVR.get_start_action_view() * viewer.core.get_model();

	double unit_length = compute_average_length_of_crossing_edges(path, start_V, startEV);
	Eigen::MatrixXd resampled_path = CleanStroke3D::resample_by_length_with_fixes(path, unit_length);
	if (!is_counter_clockwise_boundaries(resampled_path, modelview, viewer.core.get_proj(), viewer.core.viewport, mean_viewpoint, true)) {
		resampled_path = resampled_path.colwise().reverse().eval();
	}

	vector<int> path_vertices;
	int original_V_size = m.V.rows();
	for (int i = 0; i < resampled_path.rows(); i++) {
		path_vertices.push_back(original_V_size + i); //Store the vertex indices (in m.V) of the path vertices
	}


	//We need the old vertex count for this, so do it before update_mesh_values (which will increase the vertex count)
	vector<PathElement> new_surface_path;
	for (int i = 0; i < resampled_path.rows(); i++) {
		new_surface_path.push_back(PathElement(resampled_path.row(i), m.V.rows() + i));

		for (int j = 0; j < edge_split_positions.rows(); j++) {
			if (resampled_path.row(i) == edge_split_positions.row(j)) {
				replacing_edges(j * 2 + 0, 0) = m.V.rows() + i; //Middle vertex is connected to 2 other vertices
				replacing_edges(j * 2 + 1, 0) = m.V.rows() + i;
				replacing_vertex_bindings(j, 3) = m.V.rows() + i;
				break;
			}
		}
	}
	surface_path.set_path(new_surface_path);


	Eigen::MatrixXi added_edges(0, 4);
	update_mesh_values(m, resampled_path, surface_path.get_origin_stroke_ID(), m.V.rows(), false, added_edges);

	//Add interior path vertices in reversed direction to the path, so we can complete the loop around it
	for (int i = path_vertices.size() - 2; i > 0; i--) {
		path_vertices.push_back(path_vertices[i]);
	}

	//Make the outer_boundary_vertices CCW
	Eigen::MatrixXd tmp_V;
	Eigen::VectorXi row_idx2, col_idx2(3);
	row_idx2 = Eigen::VectorXi::Map(outer_boundary_vertices.data(), outer_boundary_vertices.size());
	col_idx2.col(0) << 0, 1, 2;
	igl::slice(m.V, row_idx2, col_idx2, tmp_V); //Keep only the clean "boundary" vertices in the mesh

	if (!is_counter_clockwise_boundaries(tmp_V, modelview, viewer.core.get_proj(), viewer.core.viewport, mean_viewpoint, !is_front_loop)) {
		reverse(outer_boundary_vertices.begin(), outer_boundary_vertices.end());
	}

	try {
		stitch(path_vertices, outer_boundary_vertices, m, false, replacing_edges.leftCols(2));
	}
	catch (int ex) {
		if (ex == -1) {
				std::cerr << "Stitching the new stroke failed. Please try again" << std::endl;
				m.F = start_F;
				m.V = start_V;
				m.new_mapped_indices = start_new_mapped_indices;
				m.edge_boundary_markers = start_edge_boundary_markers;
				m.vertex_is_fixed = start_vertex_is_fixed;
				m.sharp_edge = start_sharp_edge;
				return false;
		}
	}
	Eigen::MatrixXi new_edge_indicators = igl::cat(1, original_sharp_or_boundary_edges, replacing_edges);
	new_edge_indicators = igl::cat(1, new_edge_indicators, added_edges);
	try {
		update_edge_indicators(m, new_edge_indicators);
	}
	catch (int ex) {
		if (ex == -1) {
			std::cerr << "Remeshing results in a non edge-manifold mesh, which cannot be processed. Please try again. " << std::endl;
			m.F = start_F;
			m.V = start_V;
			m.new_mapped_indices = start_new_mapped_indices;
			m.edge_boundary_markers = start_edge_boundary_markers;
			m.vertex_is_fixed = start_vertex_is_fixed;
			m.sharp_edge = start_sharp_edge;
			return false;
		}
	}

	update_add_path_points_and_bindings(open_path_stroke, surface_path);
	return true;
}

int LaplacianRemesh::remesh_cutting_path(Mesh& m, Stroke& cut_path_stroke, Eigen::MatrixXi& replacing_vertex_bindings, igl::opengl::glfw::Viewer &viewer) {
	replacing_vertex_bindings.resize(0, 4);
	int remesh_success = 1;
	is_front_loop = false;
	remove_inside_faces = false;
	adjacency_list(m.F, VV);

	SurfacePath surface_path;
	bool success = surface_path.create_from_stroke_cut(cut_path_stroke);
	if (!success) {
		return 0;
	}
	remesh(m, surface_path, viewer.core.get_model(), viewer.oculusVR.get_start_action_view(), viewer.core.get_proj(), viewer.core.viewport, remesh_success, -1, replacing_vertex_bindings);
	update_add_path_points_and_bindings(cut_path_stroke, surface_path);
	return remesh_success;
}

Eigen::VectorXi LaplacianRemesh::remesh(Mesh& m, SurfacePath& surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, int& remesh_success, int cut_clicked_face, Eigen::MatrixXi& replacing_vertex_bindings) {
	replacing_vertex_bindings.resize(0, 4);
	Eigen::MatrixXi start_F = m.F;
	Eigen::MatrixXd start_V = m.V;
	Eigen::VectorXi start_edge_boundary_markers = m.edge_boundary_markers;
	Eigen::VectorXi start_vertex_is_fixed = m.vertex_is_fixed;
	Eigen::VectorXi start_sharp_edge = m.sharp_edge;
	Eigen::VectorXi start_new_mapped_indices = m.new_mapped_indices;

	vector<vector<int>> VF, VI;
	igl::edge_topology(m.V, m.F, EV, FE, EF);
	igl::vertex_triangle_adjacency(m.V.rows(), m.F, VF, VI);
	Eigen::MatrixXi startEV = EV;
	Eigen::MatrixXd startV = m.V;


	vector<bool> dirty_face(m.F.rows());
	Eigen::VectorXi dirty_vertices(m.V.rows());
	dirty_vertices.setZero();
	m.new_mapped_indices.resize(m.V.rows()); //Resize the map from old to new (clean) vertex indices to allow it to contain the number of vertices that are in the mesh at the start

	Eigen::MatrixXi original_sharp_or_boundary_edges(0, 4);
	for (int i = 0; i < m.sharp_edge.rows(); i++) {
		if (m.sharp_edge[i] || m.edge_boundary_markers[i]) {
			original_sharp_or_boundary_edges.conservativeResize(original_sharp_or_boundary_edges.rows() + 1, Eigen::NoChange);
			original_sharp_or_boundary_edges.bottomRows(1) << EV(i, 0), EV(i, 1), m.edge_boundary_markers[i], m.sharp_edge[i];
		}
	}

	int face = -1, edge = -1;
	vector<PathElement> path = surface_path.get_path();

	if (!is_front_loop) { //Append first element at the end in case of cutting, to make it a loop and compatible with extrusion's loop
		path.push_back(path[0]);
	}

	//Find outside and inside vertices. +1 means that the vertex will stay, -1 means it will be destroyed
	for (int i = 0; i < path.size(); i++) {
		if (path[i].get_type() == PathElement::FACE) {
			face = path[i].get_ID();
		}
		else {
			edge = path[i].get_ID();
			if (EF(edge, 0) == face) {
				dirty_vertices[EV(edge, 0)] += 1;
				dirty_vertices[EV(edge, 1)] -= 1;
			}
			else {
				dirty_vertices[EV(edge, 0)] -= 1;
				dirty_vertices[EV(edge, 1)] += 1;
			}
			face = (EF(edge, 0) == face) ? EF(edge, 1) : EF(edge, 0); //get the polygon on the other side of the edge
		}
	}

	//Collect vertices on the boundary
	vector<int> outer_boundary_vertices;
	vector<int> inner_boundary_vertices;
	for (int i = 0; i < m.V.rows(); i++) {
		if (dirty_vertices[i] < 0) { //vertices that are removed
			inner_boundary_vertices.push_back(i);
		}
		else if (dirty_vertices[i] > 0) { //vertices that stay
			outer_boundary_vertices.push_back(i);
		}
	}

	if (outer_boundary_vertices.size() == 0 || inner_boundary_vertices.size() == 0) { //There are no vertices left behind (cut that removes everthing or extrude that doesn't include at least 1 vertex) or none are removed
		std::cerr << "This action either removes all or no vertices, which is not possible. Please try again. " << std::endl;
		remesh_success = false;
		return Eigen::VectorXi::Zero(1);
	}

	//Collect faces along the path
	for (int i = 0; i < path.size(); i++) {
		if (path[i].get_type() == PathElement::EDGE) {
			dirty_face[EF(path[i].get_ID(), 0)] = true;
			dirty_face[EF(path[i].get_ID(), 1)] = true;
		}
	}

	//Catches case where a surface path enters and exits a triangle through the same edge. This triangle will be marked a dirty face but one of its vertices won't be added to the outer boundary vertices, and thus won't be stitched -> leads to holes.
	for (int i = 0; i < dirty_face.size(); i++) {
		if (dirty_face[i]) {
			for (int j = 0; j < 3; j++) {
				if (dirty_vertices[m.F(i, j)] == 0) {
					dirty_vertices[m.F(i, j)] = 1;
					outer_boundary_vertices.push_back(m.F(i, j)); 
				}
			}					
		}
	}

	//User clicked on one of the faces that the cutline goes through. Ignore & let them click again
	if (cut_clicked_face >= 0 && dirty_face[cut_clicked_face] && remove_inside_faces) {
		std::cerr << "You clicked too close to the cut line. Please try again. " << std::endl;
		remesh_success = -1;
		return Eigen::VectorXi::Zero(1);
	}

	if (remove_inside_faces) {
		for (int i = 0; i < m.V.rows(); i++) {
			if (dirty_vertices[i] < 0) {
				for (int j = 0; j < VF[i].size(); j++) {
					if (!dirty_face[VF[i][j]]) {
						propagate_dirty_faces(VF[i][j], dirty_face);
					}
				}
			}
		}
	}

	//Revert dirty/clean faces and vertices in case the user clicked on the "clean" part when cutting
	if (cut_clicked_face >= 0 && !dirty_face[cut_clicked_face] && remove_inside_faces) {
		for (int i = 0; i < dirty_face.size(); i++) {
			dirty_face[i] = !dirty_face[i];
		}

		dirty_vertices *= -1;

		inner_boundary_vertices.swap(outer_boundary_vertices);

		//Collect faces along the path
		for (int i = 0; i < path.size(); i++) {
			if (path[i].get_type() == PathElement::EDGE) {
				dirty_face[EF(path[i].get_ID(), 0)] = true;
				dirty_face[EF(path[i].get_ID(), 1)] = true;
			}
		}
	}

	//Remove dirty faces
	vector<int> clean_faces;
	for (int i = 0; i < m.F.rows(); i++) {
		if (!dirty_face[i]) {
			clean_faces.push_back(i);
		}
	}


	if (clean_faces.size() == 0) { //There are no faces left behind (cut that removes everything or extrude that doesn't include at least 1 vertex)
		if (cut_clicked_face == -1) { //Extrusion
			std::cerr << "The base does not include at least 1 vertex.  Please try again." << std::endl;
		}
		else { //Cut
			std::cerr << "This removes all vertices. Please try again." << std::endl;
		}
		remesh_success = false;
		return Eigen::VectorXi::Zero(1);
	}

	//Determine which vertices are clean
	Eigen::VectorXi vertex_is_clean = Eigen::VectorXi::Zero(m.V.rows());
	vector<int> clean_vertices;
	if (!remove_inside_faces) { //In the case of a path-adding "cut" stroke, we don't remove any vertices, so they're all clean
		dirty_vertices.setConstant(1);
	}

	for (int i = 0; i < m.V.rows(); i++) {
		if (dirty_vertices[i] == 0) { //Vertices with a - value are already known to be removed and vertices with a + value are already known to stay
			for (int j = 0; j < VF[i].size(); j++) {
				if (!dirty_face[VF[i][j]]) { //If the vertex is part of a clean face, then the vertex has to be clean
					vertex_is_clean(i) = 1;
					clean_vertices.push_back(i);
					break; //Can skip checking all other adjacent faces for this vertex
				}
			}
		}
		else if (dirty_vertices[i] > 0) {
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
	try {
		outer_boundary_vertices = sort_boundary_vertices(path[0].get_vertex(), outer_boundary_vertices, m);
		if (!remove_inside_faces) {
			inner_boundary_vertices = sort_boundary_vertices(path[0].get_vertex(), inner_boundary_vertices, m);
		}
	}
	catch (int ex) {
		if (ex == -1) {
			std::cerr << "Could not process cut stroke or extrusion base. Please try again. " << std::endl;
			m.F = start_F;
			remesh_success = false;
			return Eigen::VectorXi::Zero(1);
		}
		else {
			std::cout << "Something weird happened. Check what's going on." << std::endl;
			m.F = start_F;
			remesh_success = false;
			return Eigen::VectorXi::Zero(1);
		}
	}

	Eigen::RowVector3d mean_viewpoint_inner = compute_mean_viewpoint(m, inner_boundary_vertices);
	Eigen::RowVector3d mean_viewpoint_outer = compute_mean_viewpoint(m, outer_boundary_vertices);

	Eigen::MatrixXd tmp_V;
	Eigen::VectorXi tmp_markers, tmp_fixed;
	int size_before_removing = m.V.rows();
	row_idx = Eigen::VectorXi::Map(clean_vertices.data(), clean_vertices.size());
	col_idx.col(0) << 0, 1, 2;
	igl::slice(m.V, row_idx, col_idx, tmp_V);
	m.V = tmp_V; //Takes care of removing vertices that no longer exist due to all their adjacent faces being removed

	col_idx.resize(1);
	col_idx.col(0) << 0;
	igl::slice(m.vertex_is_fixed, row_idx, col_idx, tmp_fixed);
	m.vertex_is_fixed = tmp_fixed; //Compact vertex_is_fixed by removing the values for vertices that are now removed

	//Map indices in m.V at the start of remeshing to indices in m.V after remeshing
	int count = 0;
	m.new_mapped_indices.setConstant(m.new_mapped_indices.rows(), -1);
	for (int i = 0; i < vertex_is_clean.rows(); i++) {
		if (vertex_is_clean[i]) {
			m.new_mapped_indices[i] = count;
			count++;
		}
	}
	update_face_indices(m);


	Eigen::MatrixXi replacing_edges(0, 4);
	Eigen::MatrixXd edge_split_positions(0, 3);
	for (int i = 0; i < path.size(); i++) {
		if (path[i].get_type() == PathElement::EDGE && (m.sharp_edge[path[i].get_ID()] || m.edge_boundary_markers[path[i].get_ID()])) {
			path[i].fixed = true;
			replacing_edges.conservativeResize(replacing_edges.rows() + 2, Eigen::NoChange);
			replacing_edges(replacing_edges.rows() - 2, 1) = startEV(path[i].get_ID(), 0); //Pair of original edge start index and middle vertex index (middle vertex index gets set after resampling)
			replacing_edges(replacing_edges.rows() - 2, 2) = m.edge_boundary_markers[path[i].get_ID()];
			replacing_edges(replacing_edges.rows() - 2, 3) = m.sharp_edge[path[i].get_ID()];
			replacing_edges(replacing_edges.rows() - 1, 1) = startEV(path[i].get_ID(), 1); //Pair of original edge end index and middle vertex index (middle vertex index gets set after resampling)
			replacing_edges(replacing_edges.rows() - 1, 2) = m.edge_boundary_markers[path[i].get_ID()];
			replacing_edges(replacing_edges.rows() - 1, 3) = m.sharp_edge[path[i].get_ID()];
			edge_split_positions.conservativeResize(edge_split_positions.rows() + 1, Eigen::NoChange);
			edge_split_positions.bottomRows(1) << path[i].get_vertex().transpose();
			replacing_vertex_bindings.conservativeResize(replacing_vertex_bindings.rows() + 1, Eigen::NoChange);
			replacing_vertex_bindings.bottomRows(1) << m.edge_boundary_markers[path[i].get_ID()], startEV(path[i].get_ID(), 0), startEV(path[i].get_ID(), 1), -1; //Middle vertex index gets set after resampling
		}
		else {
			path[i].fixed = false;
		}
	}

	double unit_length = (is_front_loop) ? compute_average_distance_between_onPolygon_vertices(path, true) : compute_average_length_of_crossing_edges(path, startV, startEV);

	if (CleanStroke3D::get_stroke_length(path, 0, path.size() - 1) / unit_length < 12) {
		unit_length = CleanStroke3D::get_stroke_length(path, 0, path.size() - 1) / 12; //Adapt length such that we get 12 samples
	}


	Eigen::MatrixXd resampled_path = CleanStroke3D::resample_by_length_with_fixes(path, unit_length);

	Eigen::Matrix4f modelview = view * model;

	//For non-extrusion base strokes, check if they are counter clockwise (extrusion bases are made CCW when created)
	if (!is_front_loop && !is_counter_clockwise_boundaries(resampled_path, modelview, proj, viewport, mean_viewpoint_inner, !is_front_loop)) {
		resampled_path = resampled_path.colwise().reverse().eval();
	}

	vector<int> path_vertices;
	int original_V_size = m.V.rows();
	for (int i = 0; i < resampled_path.rows() - 1; i++) { //Do not use the last path vertex, as it is a copy of the first and creates unwanted behaviour
		path_vertices.push_back(original_V_size + i); //Store the vertex indices (in m.V) of the path vertices
	}


	Eigen::MatrixXi added_edges(0, 4);
	update_mesh_values(m, resampled_path, surface_path.get_origin_stroke_ID(), vertex_is_clean.rows(), true, added_edges);

	vector<PathElement> new_surface_path;
	for (int i = 0; i < resampled_path.rows(); i++) { //Make the new surfacePath looped (unlike the path_vertices indices)
		new_surface_path.push_back(PathElement(resampled_path.row(i), size_before_removing + i));

		for (int j = 0; j < edge_split_positions.rows(); j++) {
			if (resampled_path.row(i) == edge_split_positions.row(j) && i<resampled_path.rows()-1) { //Don't check the last vertex for looped paths that have the last vertex as a copy of the 0-th vertex (will give out-of-bounds vertex index)
				replacing_edges(j * 2 + 0, 0) = vertex_is_clean.rows() + i; //Middle vertex is connected to 2 other vertices
				replacing_edges(j * 2 + 1, 0) = vertex_is_clean.rows() + i;
				replacing_vertex_bindings(j, 3) = vertex_is_clean.rows() + i;
				break;
			}
		}
	}

	surface_path.set_path(new_surface_path);

	//Update the outer_boundary_vertices (which are indices into V before it was sliced to keep only the clean vertices)
	for (int i = 0; i < outer_boundary_vertices.size(); i++) {
		outer_boundary_vertices[i] = m.new_mapped_indices[outer_boundary_vertices[i]];
	}


	//Make the outer_boundary_vertices CCW
	Eigen::VectorXi row_idx2, col_idx2(3);
	row_idx2 = Eigen::VectorXi::Map(outer_boundary_vertices.data(), outer_boundary_vertices.size());
	col_idx2.col(0) << 0, 1, 2;
	igl::slice(m.V, row_idx2, col_idx2, tmp_V); //Keep only the clean "boundary" vertices in the mesh

	//Only use the regular check for non-extrusion generated boundaries (check will go wrong if user rotated the mesh between drawing the base and silhouette
	if (!is_front_loop && !is_counter_clockwise_boundaries(tmp_V, modelview, proj, viewport, mean_viewpoint_inner, !is_front_loop)) {
		reverse(outer_boundary_vertices.begin(), outer_boundary_vertices.end());
	}
	else if (is_front_loop) { //Extrusion
		Eigen::Vector3d path_normal = get_normal_from_curve(resampled_path);
		Eigen::Vector3d boundary_normal = get_normal_from_curve(tmp_V);
		if (path_normal.dot(boundary_normal) < 0) { //Normal based on the surface path and boundary vertices point in different directions so we need to reverse the order of the boundary vertices
			reverse(outer_boundary_vertices.begin(), outer_boundary_vertices.end());
		}
	}

	try {
		stitch(path_vertices, outer_boundary_vertices, m, false, replacing_edges.leftCols(2));
	}
	catch (int ex) {
		if (ex == -1) {
			std::cerr << "Stitching the new stroke failed. Please try again" << std::endl;
			m.F = start_F;
			m.V = start_V;
			m.new_mapped_indices = start_new_mapped_indices;
			m.edge_boundary_markers = start_edge_boundary_markers;
			m.vertex_is_fixed = start_vertex_is_fixed;
			m.sharp_edge = start_sharp_edge;
			remesh_success = false;
			return Eigen::VectorXi::Zero(1);
		}
	}


	if (!remove_inside_faces) {
		for (int i = 0; i < inner_boundary_vertices.size(); i++) {
			inner_boundary_vertices[i] = m.new_mapped_indices[inner_boundary_vertices[i]];
		}
		row_idx2 = Eigen::VectorXi::Map(inner_boundary_vertices.data(), inner_boundary_vertices.size());
		igl::slice(m.V, row_idx2, col_idx2, tmp_V); //Keep only the clean "boundary" vertices in the mesh
		if (is_counter_clockwise_boundaries(tmp_V, modelview, proj, viewport, mean_viewpoint_outer, !is_front_loop)) { //We want the inner_boundary_vertices to be clockwise when seen from the other side
			reverse(inner_boundary_vertices.begin(), inner_boundary_vertices.end());
		}
		try {
			stitch(path_vertices, inner_boundary_vertices, m, true, replacing_edges.leftCols(2));
		}
		catch (int ex) {
			if (ex == -1) {
				std::cerr << "Stitching the new stroke failed. Please try again" << std::endl;
				m.F = start_F;
				m.V = start_V;
				m.new_mapped_indices = start_new_mapped_indices;
				m.edge_boundary_markers = start_edge_boundary_markers;
				m.vertex_is_fixed = start_vertex_is_fixed;
				m.sharp_edge = start_sharp_edge;
				remesh_success = false;
				return Eigen::VectorXi::Zero(1);
			}
		}
	}

	Eigen::MatrixXi new_edge_indicators;
	new_edge_indicators = igl::cat(1, original_sharp_or_boundary_edges, added_edges); //Add sharp and boundary edge indicators that are created due to the newly added curve
	new_edge_indicators = igl::cat(1, new_edge_indicators, replacing_edges); //Add sharp edge & boundary indicators that are created due to old sharp or boundary edges being broken

	try {
		update_edge_indicators(m, new_edge_indicators);
	}
	catch (int ex) {
		if (ex == -1) {
			std::cerr << "Remeshing results in a non edge-manifold mesh, which cannot be processed. Please try again. " << std::endl;
			m.F = start_F;
			m.V = start_V;
			m.new_mapped_indices = start_new_mapped_indices;
			m.edge_boundary_markers = start_edge_boundary_markers;
			m.vertex_is_fixed = start_vertex_is_fixed;
			m.sharp_edge = start_sharp_edge;
			remesh_success = false;
			return Eigen::VectorXi::Zero(1);
		}
	}
	return Eigen::VectorXi::Map(path_vertices.data(), path_vertices.size());
}

Eigen::Vector3d LaplacianRemesh::get_normal_from_curve(Eigen::MatrixXd curve_points) {
	Eigen::RowVector3d center = curve_points.colwise().mean();
	Eigen::Vector3d normal(0, 0, 0);
	Eigen::Vector3d vec0, vec1;

	for (int i = 0; i < curve_points.rows(); i++) {
		vec0 = curve_points.row(i) - center;
		vec1 = curve_points.row((i + 1) % curve_points.rows()) - center;
		normal += vec1.cross(vec0);
	}
	normal.normalize();
	return normal;
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
	for (int i = 0; i < m.F.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			m.F(i, j) = m.new_mapped_indices[m.F(i, j)];
		}
	}
}

/** Adds values for a bunch of tracking variables for the new mesh vertices. **/
void LaplacianRemesh::update_mesh_values(Mesh& m, Eigen::MatrixXd path, int stroke_ID, int new_mapped_start, bool hold_back_due_to_loop, Eigen::MatrixXi& added_edges) {
	int size_before = m.V.rows();
	m.V.conservativeResize(m.V.rows() + path.rows() - hold_back_due_to_loop, Eigen::NoChange);
	m.vertex_is_fixed.conservativeResize(m.vertex_is_fixed.rows() + path.rows() - hold_back_due_to_loop);
	m.new_mapped_indices.conservativeResize(m.new_mapped_indices.rows() + path.rows() - hold_back_due_to_loop, Eigen::NoChange);
	
	for (int i = 0; i < path.rows() - hold_back_due_to_loop; i++) {
		m.V.row(size_before + i) << path.row(i);
		m.vertex_is_fixed[size_before + i] = 1;
		m.new_mapped_indices(new_mapped_start + i) = size_before + i;
	}

	for (int i = 0; i < path.rows() - 1; i++) { //Added cutting paths have path.row(0)==path.row(N), added open paths have path.row(0)!=path.row(N). We want to loop the added edges for cutting paths but not for open paths, so always just go to path.rows()-1
		added_edges.conservativeResize(added_edges.rows() + 1, Eigen::NoChange);
		added_edges.bottomRows(1) << new_mapped_start + i, new_mapped_start + ((i + 1) + path.rows() - hold_back_due_to_loop) % (path.rows() - hold_back_due_to_loop), stroke_ID, remove_inside_faces; //Remove_inside_faces is true for extrusion and cut (which generate sharp edges) and false for adding curves (which generate smooth faces)
	}
}

/** Updates the old indicators for sharp edges and edge_boundary_markers after the mesh topology changed and inserts newly generated sharp edges & boundary edges. **/
void LaplacianRemesh::update_edge_indicators(Mesh& m, Eigen::MatrixXi& edges_to_update) {
	if (!igl::is_edge_manifold(m.F)) {
		throw - 1;
		return;
	}
	igl::edge_topology(m.V, m.F, EV, FE, EF);
	m.sharp_edge.resize(EV.rows());
	m.sharp_edge.setZero();
	m.edge_boundary_markers.resize(EV.rows());
	m.edge_boundary_markers.setZero();

	int start, end, equal_pos, val;
	Eigen::VectorXi col1Equals, col2Equals;
	for (int i = 0; i < edges_to_update.rows(); i++) {
		start = m.new_mapped_indices(edges_to_update(i, 0));
		end = m.new_mapped_indices(edges_to_update(i, 1));
		
		if (start == -1 || end == -1) { //Edge no longer exists
			continue;
		}

		col1Equals = EV.col(0).cwiseEqual(min(start, end)).cast<int>();
		col2Equals = EV.col(1).cwiseEqual(max(start, end)).cast<int>();
		val = (col1Equals + col2Equals).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge
		if (val == 2) { //When adding an open path, we might have the case where both vertices still exist, but the edge between them disappeared
			m.edge_boundary_markers[equal_pos] = edges_to_update(i, 2);
			m.sharp_edge[equal_pos] = edges_to_update(i, 3);
		}
	}
}

/** Floodfills all dirty faces. **/
void LaplacianRemesh::propagate_dirty_faces(int face, vector<bool>& dirty_face) {
	dirty_face[face] = true;
	for (int i = 0; i < 3; i++) {
		int next_face = (EF(FE(face, i), 0) == face) ? EF(FE(face, i), 1) : EF(FE(face, i), 0); //Get the adjacent polygon
		if (!dirty_face[next_face]) {
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
	int iter = 0;

	Eigen::VectorXi tmp1, tmp2;
	int min_idx, max_idx, equal_pos, v_pos, max_val, backtracking_option_idx, disconnected_option_idx;
	bool backtracking_option_available = false, disconnected_option_available = false;
	while (true) {
		iter++;
		if (iter == 10 * boundary_vertices.size()) { //We're not finding a solution
			throw - 1;
			return sorted_boundary_vertices; //Will be discarded anyway
		}

		sorted_boundary_vertices.push_back(v);
		for (int i = 0; i < VV[v].size(); i++) {
			if (std::find(boundary_vertices.begin(), boundary_vertices.end(), VV[v][i]) != boundary_vertices.end()) { //Only continue if the adjacent vertex is also a boundary vertex
				if (v < VV[v][i]) {
					min_idx = v;
					v_pos = 0;
					max_idx = VV[v][i];
				}
				else {
					min_idx = VV[v][i];
					v_pos = 1;
					max_idx = v;
				}
				tmp1 = EV.col(0).cwiseEqual(min_idx).cast<int>();
				tmp2 = EV.col(1).cwiseEqual(max_idx).cast<int>();
				max_val = (tmp1 + tmp2).maxCoeff(&equal_pos); //Tmp1 and tmp2 will contain a 1 on the positions where they equal the min and max_idx. When adding them, we'll get 2 on the row that contains both of them
				if (max_val == 2 && EF(equal_pos, v_pos) == -1 && EF(equal_pos, !v_pos) != -1) { //An edge between these vertices still exists and has a NULL face on the correct side of the edge (but not on both sides), so proceed to adjacent vertex
					prev = v;
					v = VV[v][i];
					disconnected_option_available = false; //Remove the disconnected edge option
					backtracking_option_available = false; //Remove the backtracking option
					break;
				}
				else if (max_val != 2 && VV[v][i] != prev) { //At least one of the vertices was only connected to faces that have been removed (e.g. the edge was along the mesh boundary and is cut off)
					//First check if we can find an adjacent vertex for which a connecting edge still exists, otherwise choose the current vertex
					disconnected_option_idx = VV[v][i]; //Remember the vertex along a disconnected edge that we could fall back to
					disconnected_option_available = true;
					continue;
				}
				else if (max_val != 2 && VV[v][i] == prev && !disconnected_option_available) { //Backtracking over a disconnected edge (only set this option when we're out of non-backtracking disconnected options
					backtracking_option_idx = VV[v][i]; //Remember the vertex that we could backtrack to
					backtracking_option_available = true;
					continue;
				}
			}
		}

		if (disconnected_option_available) { //If we still have an adjacent vertex along a disconnected edge buffered, use it.
			prev = v;
			v = disconnected_option_idx;
			disconnected_option_available = false;
		}
		else if (backtracking_option_available) { //If we still have the backtracking vertex buffered, use it.
			prev = v;
			v = backtracking_option_idx;
			backtracking_option_available = false;
		}

		if (v == start_v) {
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

	for (int i = 1; i < vertices.size(); i++) {
		d = (m.V.row(vertices[i]) - base.transpose()).norm();
		if (d < min) {
			min = d;
			closest = vertices[i];
		}
	}
	return closest;
}

/** Stitches the boundary vertices (that are old, already existing mesh vertices) to the path vertices of the newly drawn stroke. **/
void LaplacianRemesh::stitch(std::vector<int> path_vertices, std::vector<int> boundary_vertices, Mesh& m, bool reverse, Eigen::MatrixXi replacing_edges) {
	int path_idx = 0, outer_idx = 0;
	boundary_vertices = reorder(boundary_vertices, m.V.row(path_vertices[0]), m);

	Eigen::RowVector3d start_path_v = m.V.row(path_vertices[0]), start_outer_v = m.V.row(boundary_vertices[0]);
	Eigen::RowVector3d path_v = start_path_v, outer_v = start_outer_v, next_path_v, next_outer_v;
	int path_v_idx = path_vertices[0], outer_v_idx = boundary_vertices[0], next_path_v_idx, next_outer_v_idx;

	bool proceed_outer_v, prev_proceed_path = false;
	int iter = 0;
	while (true) {
		if (iter == 4000) {
			throw - 1;
			return;
		}
		next_path_v = m.V.row(path_vertices[(path_idx + 1) % path_vertices.size()]);
		next_outer_v = m.V.row(boundary_vertices[(outer_idx + 1) % boundary_vertices.size()]);
		next_path_v_idx = path_vertices[(path_idx + 1) % path_vertices.size()];
		next_outer_v_idx = boundary_vertices[(outer_idx + 1) % boundary_vertices.size()];

		proceed_outer_v = true;
		if (is_replacing_endpoint(outer_v_idx, replacing_edges) && is_replacing_midpoint(next_path_v_idx, replacing_edges)) { //To ensure that split edges will be connected again
			proceed_outer_v = false;
		}
		else if (is_replacing_endpoint(next_outer_v_idx, replacing_edges) && is_replacing_midpoint(path_v_idx, replacing_edges)) { //To ensure that split edges will be connected again
			proceed_outer_v = true;
		}
		else if (outer_idx == boundary_vertices.size()) {
			proceed_outer_v = false;
		}
		else if (path_idx == path_vertices.size()) {
			proceed_outer_v = true;
		}
		else if (path_idx == path_vertices.size() / 2 && prev_proceed_path == true) { //The first thing after we created a triangle using the "last" path vertex ALWAYS has to be moving the outer 
			proceed_outer_v = true;
		}
		else if ((path_v - next_path_v).dot(outer_v - next_outer_v) < 0 && path_vertices.size() / 2 == path_idx) {
			proceed_outer_v = true;
		}
		else if ((path_v - next_outer_v).norm() < (outer_v - next_path_v).norm()) {
			proceed_outer_v = true;
		}
		else {
			proceed_outer_v = false;
		}

		//Add faces
		if (proceed_outer_v) {
			m.F.conservativeResize(m.F.rows() + 1, Eigen::NoChange);
			m.F.row(m.F.rows() - 1) << path_v_idx, outer_v_idx, next_outer_v_idx;
			if (reverse) {
				m.F.row(m.F.rows() - 1) = m.F.row(m.F.rows() - 1).reverse().eval();
			}
			outer_v = next_outer_v;
			outer_v_idx = next_outer_v_idx;
			outer_idx++;
			prev_proceed_path = false;
		}
		else {
			m.F.conservativeResize(m.F.rows() + 1, Eigen::NoChange);
			m.F.row(m.F.rows() - 1) << next_path_v_idx, path_v_idx, outer_v_idx;
			if (reverse) {
				m.F.row(m.F.rows() - 1) = m.F.row(m.F.rows() - 1).reverse().eval();
			}
			path_v = next_path_v;
			path_v_idx = next_path_v_idx;
			path_idx++;
			prev_proceed_path = true;
		}

		if (path_idx == path_vertices.size() && outer_idx == boundary_vertices.size()) {
			break;
		}
		iter++;
	}
}

/** Reorders the boundary vertices so they start from the vertex closest to start_v. Keeps their original order **/
vector<int> LaplacianRemesh::reorder(vector<int> boundary_vertices, Eigen::Vector3d start_v, Mesh& m) {
	int index = 0;
	double min = INFINITY;
	for (int i = 0; i < boundary_vertices.size(); i++) {
		Eigen::Vector3d vert = m.V.row(boundary_vertices[i]);
		double d = (vert - start_v).norm();
		if (d < min) {
			min = d;
			index = i;
		}
	}

	vector<int> reordered;
	for (int i = 0; i < boundary_vertices.size(); i++) {
		reordered.push_back(boundary_vertices[(i + index) % boundary_vertices.size()]);
	}
	return reordered;
}

/** For cut we determine if the boundary is CCW by comparing the direction of the normal of the remaining boundary points against the vector from the center of the removed boundary points to the center of the remaining boundary points.
	For extrusion we simply look at the sign of the inner area of the projected points (since we're looking from the "correct" viewpoint) **/
bool LaplacianRemesh::is_counter_clockwise_boundaries(Eigen::MatrixXd boundary_points, Eigen::Matrix4f modelview, Eigen::Matrix4f proj, Eigen::Vector4f viewport, Eigen::RowVector3d mean_viewpoint, bool cut) {
	if (cut) {
		Eigen::RowVector3d center = boundary_points.colwise().mean();
		Eigen::Vector3d normal(0, 0, 0);
		Eigen::Vector3d vec0, vec1;
		for (int i = 0; i < boundary_points.rows(); i++) {
			vec0 = boundary_points.row(i) - center;
			vec1 = boundary_points.row((i + 1) % boundary_points.rows()) - center;
			normal += vec1.cross(vec0);
		}
		normal.normalize();

		if (normal.dot(center - mean_viewpoint) < 0) {
			return false;
		}
		else {
			return true;
		}
	}
	else {
		double total_area = 0.0;
		Eigen::RowVector3d pt, vert;
		Eigen::RowVector2d prev, next;
		vert = boundary_points.row(boundary_points.rows() - 1);
		igl::project(vert, modelview, proj, viewport, pt);
		prev = pt.leftCols(2);

		for (int i = 0; i < boundary_points.rows(); i++) {
			vert = boundary_points.row(i);
			igl::project(vert, modelview, proj, viewport, pt);
			next = pt.leftCols(2);
			total_area += (prev[1] + next[1]) * (next[0] - prev[0]);
			prev = next;
		}

		if (total_area > 0) {
			return false;
		}

		return true;
	}
}

/** Computes the average distance between PathElements that are on faces (and not on edges). **/
double LaplacianRemesh::compute_average_distance_between_onPolygon_vertices(std::vector<PathElement> path, bool ignore_first) {
	Eigen::MatrixX3d onPoly_vertices = Eigen::MatrixX3d::Zero(0, 3);
	int start = ignore_first; //For extrusion bases we want to rotate the vector so that we start at the first point that was drawn (because path actually starts with the last one for extrusion bases). For open paths we don't want to do this
	for (int i = start; i < path.size(); i++) {
		if (path[(i + path.size()) % path.size()].get_type() == PathElement::ElementType::FACE) {
			onPoly_vertices.conservativeResize(onPoly_vertices.rows() + 1, Eigen::NoChange);
			onPoly_vertices.row(onPoly_vertices.rows() - 1) = path[(i + path.size()) % path.size()].get_vertex().transpose();
		}
	}

	double total = 0.0;
	for (int i = 0; i < onPoly_vertices.rows() - 1; i++) {
		total += (onPoly_vertices.row(i) - onPoly_vertices.row(i + 1)).norm();
	}

	return total / (onPoly_vertices.rows() - 1 - start);
}

/** Computes the average length of the edges that are crossed by the SurfacePath**/
double LaplacianRemesh::compute_average_length_of_crossing_edges(std::vector<PathElement> path, Eigen::MatrixXd startV, Eigen::MatrixXi startEV) {
	double total = 0.0;
	int count = 0;
	for (int i = 0; i < path.size(); i++) {
		if (path[i].get_type() == PathElement::ElementType::EDGE) {
			total += (startV.row(startEV(path[i].get_ID(), 0)) - startV.row(startEV(path[i].get_ID(), 1))).norm();
			count++;
		}
	}
	return total / count;
}

/** Will update the stroke points that formed the origin of the SurfacePath. Will for example add path vertices on the crossed edges. **/
void LaplacianRemesh::update_add_path_points_and_bindings(Stroke& stroke, SurfacePath& surface_path) {
	vector<PathElement> path = surface_path.get_path();
	int nr_to_remove = 0;
	if (path[0].get_vertex() == path[path.size() - 1].get_vertex()) {
		nr_to_remove = 1;
	}

	Eigen::MatrixX3d new_3DPoints(path.size() + 1 - nr_to_remove, 3); //Increase size by 1 (if it's not a loop) because we want it to become a loop again
	vector<int> new_closest_vertex_indices(path.size() + 1 - nr_to_remove);
	for (int i = 0; i < path.size() - nr_to_remove; i++) {
		new_3DPoints.row(i) = path[i].get_vertex().transpose();
		new_closest_vertex_indices[i] = path[i].get_v_idx();
	}
	new_3DPoints.row(new_3DPoints.rows() - 1) = new_3DPoints.row(0);
	new_closest_vertex_indices[new_closest_vertex_indices.size() - 1] = new_closest_vertex_indices[0];

	stroke.set3DPoints(new_3DPoints);
	stroke.set_closest_vert_bindings(new_closest_vertex_indices);
}

bool LaplacianRemesh::is_replacing_endpoint(int v_idx, Eigen::MatrixXi& replacing_edges) {
	return (replacing_edges.col(1).cwiseEqual(v_idx)).any();
}

bool LaplacianRemesh::is_replacing_midpoint(int v_idx, Eigen::MatrixXi& replacing_edges) {
	return (replacing_edges.col(0).cwiseEqual(v_idx)).any();
}