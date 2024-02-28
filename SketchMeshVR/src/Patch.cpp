#include "Patch.h"
#include <igl/edge_topology.h>
#include <igl/internal_angles.h>
#include <igl/upsample.h>
#include <igl/slice.h>
#include <igl/adjacency_list.h>
#include <igl/cat.h>
#include <unordered_set>

static Eigen::MatrixXi EV, FE, EF, patchEV, patchFE, patchEF;
static int next_mesh_ID;

std::vector<Patch*> Patch::init_patches(Mesh& h) {
	next_mesh_ID = h.ID + 1;
	igl::edge_topology(h.V, h.F, EV, FE, EF);

	std::vector<Patch*> patches;
	h.face_patch_map.clear();
	h.face_patch_map.resize(h.F.rows(), nullptr);

	for (int i = 0; i < h.F.rows(); i++) {
		Patch* face_patch = h.face_patch_map[i];
		if (face_patch == nullptr) {
			Patch* new_patch = new Patch();
			Eigen::VectorXi faces;
			propagate_patch(new_patch, i, faces, h.face_patch_map, h.sharp_edge);
			(*new_patch).create_mesh_structure(h, faces);
			patches.push_back(new_patch);
		}
	}
	return patches;
}

void Patch::propagate_patch(Patch* patch, int face, Eigen::VectorXi& faces, std::vector<Patch*> &face_patch_map, Eigen::VectorXi& sharp_edge) {
	face_patch_map[face] = patch;
	faces.conservativeResize(faces.rows() + 1);
	faces.tail(1) << face;
	for (int i = 0; i < 3; i++) {
		int edge = FE(face, i);
		int face2 = (EF(edge, 0) == face) ? EF(edge, 1) : EF(edge, 0); //get the adjacent polygon
		if (face2 != -1 && (face_patch_map[face2] == nullptr) && (!sharp_edge[edge])) {
			propagate_patch(patch, face2, faces, face_patch_map, sharp_edge);
		}
	}
}

void Patch::create_mesh_structure(Mesh& m, Eigen::VectorXi& faces) {
	patch_vertices.resize(0, 3);
	patch_vertex_is_init = Eigen::VectorXi::Zero(m.V.rows());
	patch_faces.resize(0, 3);
	patch_edge_is_init = Eigen::VectorXi::Zero(EF.rows());
	mesh_to_patch_indices = Eigen::VectorXi::Constant(m.V.rows(), -1);
	Eigen::MatrixXi tmp_sharp_edge(0, 2), tmp_edge_boundary_markers(0, 3);

	for (int i = 0; i < faces.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			int edge = FE(faces[i], j);
			get_patch_edge(edge, patch_edge_is_init, patch_vertex_is_init, faces[i], patch_vertices, m.sharp_edge, tmp_sharp_edge, m.edge_boundary_markers, tmp_edge_boundary_markers, m.V, m.vertex_is_fixed, m.new_mapped_indices, mesh_to_patch_indices);
		}
	}

	for (int i = 0; i < faces.rows(); i++) {
		patch_faces.conservativeResize(patch_faces.rows() + 1, Eigen::NoChange);
		patch_faces.bottomRows(1) << mesh_to_patch_indices[m.F(faces[i], 0)], mesh_to_patch_indices[m.F(faces[i], 1)], mesh_to_patch_indices[m.F(faces[i], 2)];
	}

	igl::edge_topology(patch_vertices, patch_faces, patchEV, patchFE, patchEF);

	Eigen::VectorXi col1Equals, col2Equals;
	int start, end, val, equal_pos;
	new_sharp_edge.resize(patchEV.rows());
	new_sharp_edge.setZero();
	for (int i = 0; i < tmp_sharp_edge.rows(); i++) { //Go over all sharp edges in the patch and set them
		start = mesh_to_patch_indices[tmp_sharp_edge(i, 0)];
		end = mesh_to_patch_indices[tmp_sharp_edge(i, 1)];
		col1Equals = patchEV.col(0).cwiseEqual(std::min(start, end)).cast<int>();
		col2Equals = patchEV.col(1).cwiseEqual(std::max(start, end)).cast<int>();
		val = (col1Equals + col2Equals).maxCoeff(&equal_pos);
		if (val == 2) {
			new_sharp_edge[equal_pos] = 1;
		}
	}

	new_edge_boundary_markers.resize(patchEV.rows());
	for (int i = 0; i < tmp_edge_boundary_markers.rows(); i++) {
		start = mesh_to_patch_indices[tmp_edge_boundary_markers(i, 0)];
		end = mesh_to_patch_indices[tmp_edge_boundary_markers(i, 1)];
		col1Equals = patchEV.col(0).cwiseEqual(std::min(start, end)).cast<int>();
		col2Equals = patchEV.col(1).cwiseEqual(std::max(start, end)).cast<int>();
		val = (col1Equals + col2Equals).maxCoeff(&equal_pos);
		if (val == 2) {
			new_edge_boundary_markers[equal_pos] = tmp_edge_boundary_markers(i, 2);
		}
	}

	mesh.V = patch_vertices;
	mesh.F = patch_faces;
	mesh.edge_boundary_markers = new_edge_boundary_markers;
	mesh.vertex_is_fixed = new_vertex_is_fixed;
	mesh.new_mapped_indices = new_mapped_indices;
	mesh.sharp_edge = new_sharp_edge;
	mesh.ID = next_mesh_ID;
	mesh.mesh_to_patch_indices = mesh_to_patch_indices; //Maps a vertex index from the base mesh to a vertex index in the mesh of a underlying patch
	next_mesh_ID++;
}

void Patch::get_patch_edge(int edge, Eigen::VectorXi& patch_edge_is_init, Eigen::VectorXi& patch_vertex_is_init, int face, Eigen::MatrixXd& patch_vertices, Eigen::VectorXi& sharp_edge, Eigen::MatrixXi& tmp_sharp_edge, Eigen::VectorXi &edge_boundary_markers, Eigen::MatrixXi &tmp_edge_boundary_markers, Eigen::MatrixXd& V_orig, Eigen::VectorXi & vertex_is_fixed_orig, Eigen::VectorXi& new_mapped_indices_orig, Eigen::VectorXi& mesh_to_patch_indices) {
	if (patch_edge_is_init[edge]) {
		return;
	}

	int start = EV(edge, 0);
	int end = EV(edge, 1);
	get_patch_vertex(start, face, patch_vertices, patch_vertex_is_init, V_orig, vertex_is_fixed_orig, new_mapped_indices_orig, mesh_to_patch_indices);
	get_patch_vertex(end, face, patch_vertices, patch_vertex_is_init, V_orig, vertex_is_fixed_orig, new_mapped_indices_orig, mesh_to_patch_indices);
	if (sharp_edge[edge]) {
		tmp_sharp_edge.conservativeResize(tmp_sharp_edge.rows() + 1, Eigen::NoChange);
		tmp_sharp_edge.bottomRows(1) << start, end;
	}
	tmp_edge_boundary_markers.conservativeResize(tmp_edge_boundary_markers.rows() + 1, Eigen::NoChange);
	tmp_edge_boundary_markers.bottomRows(1) << start, end, edge_boundary_markers[edge];
	patch_edge_is_init[edge] = 1;
}

void Patch::get_patch_vertex(int v_idx, int face, Eigen::MatrixXd& patch_vertices, Eigen::VectorXi& patch_vertex_is_init, Eigen::MatrixXd& V_orig, Eigen::VectorXi& vertex_is_fixed_orig, Eigen::VectorXi& new_mapped_indices_orig, Eigen::VectorXi& mesh_to_patch_indices) {
	if (patch_vertex_is_init[v_idx]) {
		return;
	}

	patch_vertices.conservativeResize(patch_vertices.rows() + 1, Eigen::NoChange);
	patch_vertices.bottomRows(1) = V_orig.row(v_idx);
	parent_vertices.conservativeResize(parent_vertices.rows() + 1);
	parent_vertices.tail(1) << v_idx;
	mesh_to_patch_indices[v_idx] = parent_vertices.rows() - 1;
	new_vertex_is_fixed.conservativeResize(new_vertex_is_fixed.rows() + 1);
	new_vertex_is_fixed.tail(1) << vertex_is_fixed_orig.row(v_idx);
	new_mapped_indices.conservativeResize(new_mapped_indices.rows() + 1);
	new_mapped_indices.tail(1) << -1;
	patch_vertex_is_init[v_idx] = 1;
}

void Patch::update_parent_vertex_positions(Eigen::MatrixXd& base_V) {
	for (int i = 0; i < mesh.V.rows(); i++) {
		base_V.row(parent_vertices[i]) = mesh.V.row(i);
	}
}

void Patch::update_patch_vertex_positions(Eigen::MatrixXd& base_V) {
	for (int i = 0; i < base_V.rows(); i++) {
		if (mesh_to_patch_indices[i] != -1) {
			mesh.V.row(mesh_to_patch_indices[i]) = base_V.row(i);
		}
	}
}

//This will upsample *this* Patch, and it will also split the triangles that are connected to this patch via a sharp edge (because the upsampling of this patch will insert points on the patch boundary). Every triangle in the patch will be split into 4 triangles (by inserting new vertices at the edge midpoints).
void Patch::upsample_patch(Eigen::MatrixXd& base_V, Eigen::MatrixXi& base_F, std::vector<Patch*> base_face_patch_map, Eigen::VectorXi& base_edge_boundary_markers, Eigen::VectorXi& base_sharp_edge, Eigen::VectorXi& base_vertex_is_fixed, Eigen::MatrixXi& replacing_vertex_bindings) {
	replacing_vertex_bindings.resize(0, 4);
	Eigen::MatrixXi old_baseF = base_F;
	int start_size_baseV = base_V.rows();

	//Store EV of base_mesh before upsampling
	Eigen::MatrixXi startEV, startFE, startEF;
	igl::edge_topology(base_V, base_F, startEV, startFE, startEF);

	Eigen::MatrixXd new_patchV;
	Eigen::MatrixXi new_patchF;
	igl::upsample(mesh.V, mesh.F, new_patchV, new_patchF);
	//Assuming that each face is only part of 1 Patch, we can take the base_mesh's face_patch_map, and remove the faces belonging to this patch with a slice (keeping only the faces that do not map to the current patch)
	std::vector<int> clean_faces;
	for (int i = 0; i < base_face_patch_map.size(); i++) {
		if (base_face_patch_map[i] == this) {
			//Face should be removed from base_F		
		}
		else {
			clean_faces.push_back(i);
		}
	}

	//Store all clean faces from the base mesh (that are not part of this Patch) in tmp_F
	Eigen::MatrixXi tmp_F;
	Eigen::VectorXi row_idx = Eigen::VectorXi::Map(clean_faces.data(), clean_faces.size());
	Eigen::VectorXi col_idx(3);
	col_idx.col(0) << 0, 1, 2;
	igl::slice(base_F, row_idx, col_idx, tmp_F);

	//Append newly created vertices to the vertex list of the base (don't care about adding them to the patch, because this will be reinitialized afterwards anyway)
	int nr_added_points = new_patchV.rows() - mesh.V.rows();
	Eigen::MatrixXd newly_addedV = new_patchV.bottomRows(nr_added_points);
	base_V.conservativeResize(base_V.rows() + nr_added_points, Eigen::NoChange);
	base_V.bottomRows(nr_added_points) = newly_addedV;

	int start_size = mesh.V.rows();
	parent_vertices.conservativeResize(parent_vertices.rows() + nr_added_points); //Add mapping from patch face indices to base_mesh face indices for the newly added vertices
	for (int i = 0; i < nr_added_points; i++) {
		parent_vertices[start_size + i] = start_size_baseV + i;
	}

	Eigen::MatrixXi new_patchF_reindexed(new_patchF.rows(), 3);
	for (int i = 0; i < new_patchF.rows(); i++) {
		new_patchF_reindexed.row(i) << parent_vertices[new_patchF(i, 0)], parent_vertices[new_patchF(i, 1)], parent_vertices[new_patchF(i, 2)];
	}

	//Append newly created faces to the cleaned-up face list of the base
	base_F.resize(tmp_F.rows() + new_patchF_reindexed.rows(), 3);
	base_F.topRows(tmp_F.rows()) = tmp_F;
	base_F.bottomRows(new_patchF_reindexed.rows()) = new_patchF_reindexed;


	//Construct EV for the new base_V and base_F
	Eigen::MatrixXi newEV, newFE, newEF;
	igl::edge_topology(base_V, base_F, newEV, newFE, newEF);


	//Loop over all edges in the new mesh (containing the upsampled faces for *this* patch and previous faces for all other patches.
	//If there was a connection between the 2 edge vertices before upsampling, copy the edge_boundary_marker for the specific edge. 
	//Else it is a subdivided edge. If exactly one of the edge vertices already existed in the start_V, then we have an edge that is half of a broken up edge (if both vertices are new, it is an internal edge)
	//For the subdivided edge, find the other subdivided edge that attaches to the middle vertex. Check in the original vertex_boundary_markers if the edge between these 2 endpoints was marked, and copy that marker onto both edges
	Eigen::MatrixXi faces_to_add(0, 3);
	Eigen::VectorXi new_base_edge_boundary_markers(newEV.rows()), new_base_sharp_edge(newEV.rows()), new_base_vertex_is_fixed(base_vertex_is_fixed.rows() + nr_added_points), dirty_faces(base_F.rows());
	new_base_edge_boundary_markers.setConstant(-1);
	new_base_sharp_edge.setConstant(-1);
	new_base_vertex_is_fixed.topRows(base_vertex_is_fixed.rows()) = base_vertex_is_fixed;
	dirty_faces.setZero();
	int dirty_count = 0;
	int start, end, old_edge, equal_pos, maxval, old_edge_start, old_edge_end, original_edge, other_edge, marker_val, sharp_val, tmp1, tmp2;
	Eigen::VectorXi right_end_point, connects_old_point, not_current_edge;
	for (int i = 0; i < newEV.rows(); i++) {
		if (new_base_sharp_edge[i] != -1) { //Edge has already been covered
			continue;
		}

		start = newEV(i, 0);
		end = newEV(i, 1);
		if (start >= start_size_baseV && end >= start_size_baseV) { //Both vertices are newly added, so they form an interior edge (which is not sharp and not part of a boundary).
			new_base_edge_boundary_markers[i] = 0;
			new_base_sharp_edge[i] = 0;
		}
		else if (start < start_size_baseV && end < start_size_baseV) { //Both vertices were existing already, copy their existing markers
			old_edge = find_edge(start, end, startEV);
			new_base_edge_boundary_markers[i] = base_edge_boundary_markers[old_edge];
			new_base_sharp_edge[i] = base_sharp_edge[old_edge];
		}
		else { //One vertex is old (smallest index), one is new (larger index)
			//Find all edges in newEV that connect to the vertex with the larger index (new vertex). Only one of these (except for the current edge i) will be with a vertex that has an index < start_size_baseV
			right_end_point = (newEV.col(1).cwiseEqual(std::max(start, end))).cast<int>();
			connects_old_point = (newEV.col(0).array() < start_size_baseV).cast<int>();
			not_current_edge = (newEV.col(0).array() != std::min(start, end)).cast<int>();
			maxval = (right_end_point + connects_old_point + not_current_edge).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge
			if (maxval == 3) {
				old_edge_start = std::min(start, end);
				old_edge_end = std::min(newEV(equal_pos, 0), newEV(equal_pos, 1));
				original_edge = find_edge(old_edge_start, old_edge_end, startEV);
				other_edge = find_edge(old_edge_end, std::max(start, end), newEV);

				marker_val = base_edge_boundary_markers[original_edge];
				sharp_val = base_sharp_edge[original_edge];

				//Copy the markers of the edge that is split
				new_base_edge_boundary_markers[i] = marker_val;
				new_base_sharp_edge[i] = sharp_val;
				new_base_edge_boundary_markers[other_edge] = marker_val;
				new_base_sharp_edge[other_edge] = sharp_val;

				tmp1 = std::min(old_edge_start, old_edge_end);
				tmp2 = std::max(old_edge_start, old_edge_end);

				replacing_vertex_bindings.conservativeResize(replacing_vertex_bindings.rows() + 1, Eigen::NoChange);
				replacing_vertex_bindings.bottomRows(1) << marker_val, tmp1, tmp2, std::max(start, end);

				//Only set the newly added midpoint vertices to fixed if both endpoints of the original edge are fixed and the edge has a boundary marker (to prevent fixed edges from forming when 2 endpoints from the original edge each are fixed but are not connected by a boundary marked edge)
				if (base_vertex_is_fixed[old_edge_start] && base_vertex_is_fixed[old_edge_end] && base_edge_boundary_markers[original_edge]) {
					new_base_vertex_is_fixed[std::max(start, end)] = 1;
				}
				else {
					new_base_vertex_is_fixed[std::max(start, end)] = 0;
				}

				//Find face that is adjacent to the patch boundary edge that we just split, remove that face and replace it with 2 half faces (outside the patch's domain)
				if (sharp_val) {
					int edge = find_edge(old_edge_start, old_edge_end, startEV);
					int cur_edge = find_edge(old_edge_start, old_edge_end, newEV); //Find the edge between the original triangle "base" in the current EV
					int face1 = startEF(edge, 0);
					int face2 = startEF(edge, 1);

					int mid_vertex = std::max(start, end);
					int vert_to_connect;
					if (base_face_patch_map[face2] != this) { //Only add 2 new faces if the face is not on the same patch (we don't want to add split faces across sharp edges that are within 1 patch)
						for (int j = 0; j < 3; j++) {
							if ((old_baseF(face2, j) != old_edge_start) && (old_baseF(face2, j) != old_edge_end)) {
								vert_to_connect = old_baseF(face2, j);
								faces_to_add.conservativeResize(faces_to_add.rows() + 1, Eigen::NoChange);
								faces_to_add.bottomRows(1) << old_edge_start, vert_to_connect, mid_vertex;
								faces_to_add.conservativeResize(faces_to_add.rows() + 1, Eigen::NoChange);
								faces_to_add.bottomRows(1) << mid_vertex, vert_to_connect, old_edge_end;
								int face2_cur = (newEF(cur_edge, 0) == -1) ? newEF(cur_edge, 1) : newEF(cur_edge, 0);
								dirty_faces[face2_cur] = 1;
							}
						}
					}
					else if (base_face_patch_map[face1] != this) {
						for (int j = 0; j < 3; j++) {
							if ((old_baseF(face1, j) != old_edge_start) && (old_baseF(face1, j) != old_edge_end)) {
								vert_to_connect = old_baseF(face1, j);
								faces_to_add.conservativeResize(faces_to_add.rows() + 1, Eigen::NoChange);
								faces_to_add.bottomRows(1) << old_edge_start, mid_vertex, vert_to_connect;
								faces_to_add.conservativeResize(faces_to_add.rows() + 1, Eigen::NoChange);
								faces_to_add.bottomRows(1) << mid_vertex, old_edge_end, vert_to_connect;
								int face1_cur = (newEF(cur_edge, 0) == -1) ? newEF(cur_edge, 1) : newEF(cur_edge, 0);
								dirty_faces[face1_cur] = 1;
							}
						}
					}
				}
			}
		}
	}

	//Remove the faces that are adjacent to *this* patch via a sharp edge (because they have been replaced by 2 half faces)
	clean_faces.clear();
	for (int i = 0; i < base_F.rows(); i++) {
		if (!dirty_faces[i]) {
			clean_faces.push_back(i);
		}
	}
	row_idx = Eigen::VectorXi::Map(clean_faces.data(), clean_faces.size());
	igl::slice(base_F, row_idx, col_idx, tmp_F);

	base_F = tmp_F;
	base_F = igl::cat(1, base_F, faces_to_add); //Add the replacing split faces

	update_edge_indicators(base_V, base_F, new_base_sharp_edge, new_base_edge_boundary_markers, newEV);

	base_sharp_edge = new_base_sharp_edge;
	base_edge_boundary_markers = new_base_edge_boundary_markers;
	base_vertex_is_fixed = new_base_vertex_is_fixed;
}

//Specifically adapted version for use in Patch. Searches for an edge between vertices start and end in the given list of edges EV
int Patch::find_edge(int start, int end, Eigen::MatrixXi& EV) {
	Eigen::VectorXi col1Equals, col2Equals;
	int equal_pos;
	col1Equals = EV.col(0).cwiseEqual(std::min(start, end)).cast<int>();
	col2Equals = EV.col(1).cwiseEqual(std::max(start, end)).cast<int>();
	int val = (col1Equals + col2Equals).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge

	if (val != 2) {
		return -1; //Edge between start and end doesn't exist in EV
	}
	return equal_pos;
}

//Used to update the sharp_edge and edge_boundary_marker indicators when we go from a mesh that corresponds to the edges in oldEV to one with newly added edges in EV (but those newly added edges have edge_boundary_markers 0 and aren't sharp)
void Patch::update_edge_indicators(Eigen::MatrixXd& meshV, Eigen::MatrixXi& meshF, Eigen::VectorXi& sharp_edge, Eigen::VectorXi& edge_boundary_markers, Eigen::MatrixXi& oldEV) {
	Eigen::MatrixXi EV, FE, EF;
	igl::edge_topology(meshV, meshF, EV, FE, EF);
	Eigen::VectorXi new_sharp_edge(EV.rows());
	Eigen::VectorXi new_edge_boundary_markers(EV.rows());
	new_sharp_edge.setZero(); //Can init new_sharp_edge and new_edge_boundary_markers to zeroes, because edges that have not previously been initialized, won't be part of a sharp edge or boundary
	new_edge_boundary_markers.setZero();

	int new_edge;
	for (int i = 0; i < oldEV.rows(); i++) {
		if (sharp_edge[i] || edge_boundary_markers[i]) { //Only need to find the new edge when there is a non-zero indicator
			new_edge = find_edge(oldEV(i, 0), oldEV(i, 1), EV);
			if (new_edge == -1) { //Edge not found anymore. Likely a edge on the (sharp) patch boundary, that was part of a face that is now split (and thus removed)
				continue;
			}
			new_sharp_edge[new_edge] = sharp_edge[i];
			new_edge_boundary_markers[new_edge] = edge_boundary_markers[i];
		}
	}

	sharp_edge = new_sharp_edge;
	edge_boundary_markers = new_edge_boundary_markers;
}