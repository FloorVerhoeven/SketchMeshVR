#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/per_vertex_normals.h>
#include <igl/adjacency_list.h>
#include <igl/edge_topology.h>
#include "SurfaceSmoothing.h"
#include "Patch.h"

using namespace std;
using namespace igl;

typedef Eigen::Triplet<double> T;

int SurfaceSmoothing::prev_vertex_count = -1;

static int prev_mesh_ID = -1;
static int prev_patch_count = -1;
Eigen::MatrixX3d SurfaceSmoothing::vertex_normals(0, 3);
std::vector<Eigen::SparseLU<Eigen::SparseMatrix<double>>*> solver1_array;
std::vector<Eigen::SparseLU<Eigen::SparseMatrix<double>>*> solver2_array;

Eigen::VectorXd initial_curvature, average_edge_lengths;
Eigen::VectorXd SurfaceSmoothing::curvatures(ptrdiff_t(0));
int ID = 0, iteration = 0;
int no_boundary_vertices, no_boundary_adjacent_vertices;
double SurfaceSmoothing::vertex_weight = 1000.0;
double SurfaceSmoothing::edge_weight = 0.001;
double SurfaceSmoothing::factor = 0.1;
double SurfaceSmoothing::existing_curvature_weight = 0.01; //Controls how fast the surface inflates. A lower value will lead to quicker inflation, but thin regions will stay flatter
vector<vector<int>> neighbors;


//Map to the mesh IDs
std::vector<Eigen::MatrixXd> precomputed_L;
std::vector<Eigen::VectorXd> precomputed_laplacian_weights;
std::unordered_map<int, Eigen::SparseMatrix<double>> SurfaceSmoothing::precompute_matrix_for_LM_and_edges;
std::unordered_map<int, Eigen::SparseMatrix<double>> SurfaceSmoothing::AT_for_LM_and_edges;
std::unordered_map<int, Eigen::SparseMatrix<double>> SurfaceSmoothing::precompute_matrix_for_positions;
std::unordered_map<int, Eigen::SparseMatrix<double>> SurfaceSmoothing::AT_for_positions;

void SurfaceSmoothing::smooth(Mesh& base_mesh, bool& BOUNDARY_IS_DIRTY, bool force_update){
	if (base_mesh.patches.empty()) {
		base_mesh.patches = Patch::init_patches(base_mesh);
	}

	if(force_update || (prev_vertex_count != base_mesh.V.rows()) || (prev_patch_count != base_mesh.patches.size())) { //The mesh topology has changed, so we need to reset the precomputed matrices. If only boundary constraints got added, we can reuse part of the matrices so don't clear it all out but instead overwrite certain parts
		clear_precomputed_matrices(base_mesh.patches.size());
		prev_vertex_count = base_mesh.V.rows();
		prev_patch_count = base_mesh.patches.size();
		iteration = 0;
		for (int i = 0; i < base_mesh.patches.size(); i++) {
			solver1_array.push_back(new Eigen::SparseLU<Eigen::SparseMatrix<double>>);
			solver2_array.push_back(new Eigen::SparseLU<Eigen::SparseMatrix<double>>);
		}
	}
	if(BOUNDARY_IS_DIRTY) {
		iteration = 0;
	}

	for (int i = 0; i < base_mesh.patches.size(); i++) {
		Patch* patch = (base_mesh.patches[i]);
		std::cout << "Patch " << i << " has " << (*patch).mesh.F.rows() << " faces" << std::endl;
		smooth_main((*patch).mesh, BOUNDARY_IS_DIRTY);
		(*patch).update_parent_vertex_positions(base_mesh.V);
		prev_mesh_ID = (*patch).mesh.ID;
	}

	BOUNDARY_IS_DIRTY = false;
	iteration++;
}

void SurfaceSmoothing::clear_precomputed_matrices(int no_patches) {
	precomputed_L.clear();
	precomputed_L.resize(no_patches + 1);
	precomputed_laplacian_weights.clear();
	precomputed_laplacian_weights.resize(no_patches + 1);

	for(auto it = precompute_matrix_for_LM_and_edges.begin(); it != precompute_matrix_for_LM_and_edges.end(); ++it) {
		it->second.resize(0, 0); //Resize the referenced Eigen::MatrixXd to size 0
	}
	precompute_matrix_for_LM_and_edges.clear();

	for(auto it = precompute_matrix_for_positions.begin(); it != precompute_matrix_for_positions.end(); ++it) {
		it->second.resize(0, 0); //Resize the referenced Eigen::MatrixXd to size 0
	}
	precompute_matrix_for_positions.clear(); 
}

void SurfaceSmoothing::get_average_edge_lengths(Eigen::MatrixXd& V) {
	double total;
	average_edge_lengths.resize(V.rows());
	for (int i = 0; i < V.rows(); i++) {
		total = 0.0;
		for (int j = 0; j < neighbors[i].size(); j++) {
			total += (V.row(i) - V.row(neighbors[i][j])).norm();
		}
		average_edge_lengths[i] = (total / neighbors[i].size());
	}
}

/*void SurfaceSmoothing::performance_test(Mesh& m) {
	no_boundary_vertices = (m.vertex_is_fixed.array() > 0).count();
	no_boundary_adjacent_vertices = 0;
	for (int i = 0; i < m.V.rows(); i++) {
		if (m.vertex_is_fixed[i] > 0) {
			for (int j = 0; j < neighbors[i].size(); j++) {
				if (m.vertex_is_fixed[neighbors[i][j]] == 0) {
					no_boundary_adjacent_vertices++;
				}
			}
		}
	}
}*/

void SurfaceSmoothing::smooth_main(Mesh &m, bool BOUNDARY_IS_DIRTY) {
	Eigen::MatrixXd& L = precomputed_L[m.ID - 1]; //Patch mesh IDs start at 1 (0 is reserved for the base mesh). So take -1 to index into the vector
	
	if(L.rows()==0 && L.cols()==0) { //Stuff that gets computed once, when we get a new mesh topology
		L = compute_laplacian_matrix(m);
		precomputed_L[m.ID - 1] = L;
		adjacency_list(m.F, neighbors);
	}
	else if (m.ID != prev_mesh_ID){
		adjacency_list(m.F, neighbors);
	}
	

	if(BOUNDARY_IS_DIRTY || m.ID != prev_mesh_ID) {
		no_boundary_vertices = (m.vertex_is_fixed.array() > 0).count();
		no_boundary_adjacent_vertices = 0;
		for(int i = 0; i < m.V.rows(); i++) {
			if (m.vertex_is_fixed[i] > 0) {
				for(int j = 0; j < neighbors[i].size(); j++) {
					if (m.vertex_is_fixed[neighbors[i][j]] == 0) {
						no_boundary_adjacent_vertices++;
					}
				}
			}
		}
	}

	get_average_edge_lengths(m.V);
	igl::per_vertex_normals(m.V, m.F, PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, vertex_normals);
	Eigen::VectorXd& target_LMs = compute_target_LMs(m, L, BOUNDARY_IS_DIRTY);
	Eigen::VectorXd& target_edge_lengths = compute_target_edge_lengths(m, L);
	compute_target_vertices(m, L, target_LMs, target_edge_lengths, BOUNDARY_IS_DIRTY);
}

Eigen::MatrixXd SurfaceSmoothing::compute_laplacian_matrix(Mesh &m) {
	Eigen::SparseMatrix<double> A;
	adjacency_matrix(m.F, A);
	// sum each row 
	Eigen::SparseVector<double> Asum;
	sum(A, 1, Asum);
	// Convert row sums into diagonal of sparse matrix
	Eigen::SparseMatrix<double> Adiag;
	diag(Asum, Adiag);
	// Build uniform laplacian
	Eigen::SparseMatrix<double> U;
	U = (Adiag.cwiseInverse())*(Adiag - A);

	return (Eigen::MatrixXd) U;
}

Eigen::VectorXd SurfaceSmoothing::compute_initial_curvature(Mesh& m) {
	Eigen::VectorXd initial_curvatures(m.V.rows());
	initial_curvatures.setZero();
		for (int i = 0; i < m.V.rows(); i++) {
			if (m.vertex_is_fixed[i] > 0) {
				Eigen::RowVector3d center(0, 0, 0);
				int count = 0;
				double l = 0.0;
				for (int j = 0; j < neighbors[i].size(); j++) {
					if (m.vertex_is_fixed[neighbors[i][j]] == 0) {
						continue;
					}
					center += m.V.row(neighbors[i][j]);
					count++;
					l += (m.V.row(i) - m.V.row(neighbors[i][j])).norm();
				}
				if (count == 1) {
					continue; //Leave curvature at 0
				}
				double d = vertex_normals.row(i).dot(m.V.row(i) - (1.0 / count)*center);
				l /= count;

				double w = 1;
				w *= l*l;

				initial_curvatures[i] = d/w;
				if (d < 0) {
					initial_curvatures[i] = 0;
				}
			}
		}
	return initial_curvatures;
}

Eigen::VectorXd SurfaceSmoothing::compute_target_LMs(Mesh &m, Eigen::MatrixXd &L, bool BOUNDARY_IS_DIRTY) {
	Eigen::SparseMatrix<double>& A = get_precompute_matrix_for_LM_and_edges(m.ID);
    Eigen::SparseMatrix<double> AT;
	if(A.rows() == 0 && A.cols() == 0 || (BOUNDARY_IS_DIRTY && iteration == 0)) {
		std::vector<T> tripletList;
		tripletList.reserve(m.V.rows()*(m.V.rows() + 1));
		for(int i = 0; i < m.V.rows(); i++) {
			for (int j = 0; j < m.V.rows(); j++) {
				if (abs(L(i, j)) > 0) {
					tripletList.push_back(T(i, j, L(i, j)));
				}
			}
		
			if (m.vertex_is_fixed[i] > 0) { //Constrain only the boundary in the first iteration
				tripletList.push_back(T(m.V.rows() + i, i, 1));
			}
		}
		A = Eigen::SparseMatrix<double>(m.V.rows() * 2, m.V.rows());
		A.setFromTriplets(tripletList.begin(), tripletList.end());
		A.prune(0.0);
		AT = A.transpose();
		(*solver1_array[m.ID - 1]).compute(AT*A); //Take ID-1 because ID 0 is reserved for the basemesh and patches start from 1

		set_precompute_matrix_for_LM_and_edges(m.ID, A);
        set_AT_for_LM_and_edges(m.ID, AT);
	}
	
	if(iteration == 1) { //Constrain all vertices after the first iteration
		std::vector<T> tripletList;
		tripletList.reserve(m.V.rows()*(m.V.rows() + 1));
		for (int i = 0; i < m.V.rows(); i++) {
			for (int j = 0; j < m.V.rows(); j++) {
				if (abs(L(i, j)) > 0) {
					tripletList.push_back(T(i, j, L(i, j)));
				}
			}
			tripletList.push_back(T(m.V.rows() + i, i, existing_curvature_weight));
		}

		A.setFromTriplets(tripletList.begin(), tripletList.end());
		A.prune(0.0);
		AT = A.transpose();
		(*solver1_array[m.ID - 1]).compute(AT*A); //Take ID-1 because ID 0 is reserved for the basemesh and patches start from 1
		set_precompute_matrix_for_LM_and_edges(m.ID, A);
        set_AT_for_LM_and_edges(m.ID, AT);
	}
	
    Eigen::VectorXd current_curvatures;
	if(iteration == 0) {
		initial_curvature = compute_initial_curvature(m);
    }else{
        current_curvatures = get_curvatures(m);
		for (int i = 0; i < current_curvatures.rows(); i++) {
			current_curvatures[i] /= (average_edge_lengths[i] * average_edge_lengths[i]);
		}
    }

	Eigen::VectorXd b = Eigen::VectorXd::Zero(m.V.rows()*2);
    for(int i = 0; i < m.V.rows(); i++) {
		if(iteration == 0) {
			if (m.vertex_is_fixed[i] > 0) {
				b[m.V.rows() + i] = initial_curvature[i];
			}
		} else {
			b[m.V.rows() + i] = existing_curvature_weight*current_curvatures[i];
		}
	}

	AT = get_AT_for_LM_and_edges(m.ID);
	Eigen::VectorXd target_LM = (*solver1_array[m.ID - 1]).solve(AT*b);
	for (int i = 0; i < target_LM.rows(); i++) {
		if (target_LM(i) < 0) {
			target_LM(i) = 0;
		}
	}
	return target_LM;
}

Eigen::VectorXd SurfaceSmoothing::compute_target_edge_lengths(Mesh &m, Eigen::MatrixXd &L) {
	Eigen::SparseMatrix<double> A = get_precompute_matrix_for_LM_and_edges(m.ID);
    Eigen::SparseMatrix<double> AT = get_AT_for_LM_and_edges(m.ID);
	//A will be set up when computing the target LMs

	Eigen::VectorXd b = Eigen::VectorXd::Zero(m.V.rows() * 2);

	if(iteration == 0){
		for(int i = 0; i < m.V.rows(); i++) {
			if (m.vertex_is_fixed[i] > 0) {
				b[m.V.rows() + i] = average_edge_lengths[i]*0.1;
			}
		}
	} else {
		for(int i = 0; i < m.V.rows(); i++) {
			b[m.V.rows() + i] = average_edge_lengths[i] * 0.1;
		}
	}
	

	Eigen::VectorXd target_edge_lengths = (*solver1_array[m.ID - 1]).solve(AT*b);
	target_edge_lengths.cwiseMax(0.0); //Clamp edge lengths to non-negative

	return target_edge_lengths;
}

void SurfaceSmoothing::compute_target_vertices(Mesh &m, Eigen::MatrixXd &L, Eigen::VectorXd &target_LMs, Eigen::VectorXd &target_edge_lengths, bool BOUNDARY_IS_DIRTY) {
	Eigen::MatrixXd target_laplacians(m.V.rows(), 3);
	for (int i = 0; i < target_LMs.rows(); i++) {
		target_laplacians.row(i) = vertex_normals.row(i) * target_LMs[i] * (average_edge_lengths[i] * average_edge_lengths[i]);
	}
	
	Eigen::SparseMatrix<double> A = get_precompute_matrix_for_positions(m.ID);
    Eigen::SparseMatrix<double> AT = get_AT_for_positions(m.ID);
	Eigen::VectorXd& laplacian_weights = precomputed_laplacian_weights[m.ID - 1];
	if((A.rows() == 0 && A.cols() == 0) || BOUNDARY_IS_DIRTY) {
		laplacian_weights = Eigen::VectorXd::Ones(m.V.rows());
		Eigen::VectorXi& on_border = points_on_border(m);
		for (int i = 0; i < on_border.rows(); i++) {
			if (on_border[i]) {
				laplacian_weights[i] *= factor;
			}
		}
		precomputed_laplacian_weights[m.ID - 1] = laplacian_weights;

		std::vector<T> tripletList;
		tripletList.reserve(m.V.rows()*(m.V.rows() + 13));
		int count = 0, count2 = 0;

		for(int i = 0; i < m.V.rows(); i++) {
			for(int j = 0; j < m.V.rows(); j++) {
				if (abs(L(i, j)) > 0) {
					tripletList.push_back(T(i, j, L(i, j)*laplacian_weights[i]));
				}
			}
			if (m.vertex_is_fixed[i] > 0) { //Constrain only the boundary LMs and edges
				tripletList.push_back(T(m.V.rows() + count, i, vertex_weight));
				count++;

				for(int j = 0; j < neighbors[i].size(); j++) {
					if (m.vertex_is_fixed[neighbors[i][j]] == 0) {
						tripletList.push_back(T(m.V.rows() + no_boundary_vertices + count2, i, edge_weight));
						tripletList.push_back(T(m.V.rows() + no_boundary_vertices + count2, neighbors[i][j], -edge_weight));
						count2++;
					}
				}
			}
		}

		A = Eigen::SparseMatrix<double>(m.V.rows() + no_boundary_vertices + no_boundary_adjacent_vertices, m.V.rows());
		A.setFromTriplets(tripletList.begin(), tripletList.end());
		A.prune(0.0);
		AT = A.transpose();
		(*solver2_array[m.ID - 1]).compute(AT*A);

		set_precompute_matrix_for_positions(m.ID, A);
        set_AT_for_positions(m.ID, AT);
	}
    
    Eigen::VectorXd bx = Eigen::VectorXd::Zero(m.V.rows() + no_boundary_vertices + no_boundary_adjacent_vertices);
	Eigen::VectorXd by = Eigen::VectorXd::Zero(m.V.rows() + no_boundary_vertices + no_boundary_adjacent_vertices);
	Eigen::VectorXd bz = Eigen::VectorXd::Zero(m.V.rows() + no_boundary_vertices + no_boundary_adjacent_vertices);

	int count = 0, count2 = 0;
	Eigen::Vector3d delta;
	for(int i = 0; i < m.V.rows(); i++) {
		delta = target_laplacians.row(i);
		bx[i] = delta(0) * laplacian_weights[i];
		by[i] = delta(1) * laplacian_weights[i];
		bz[i] = delta(2) * laplacian_weights[i];

		if (m.vertex_is_fixed[i] > 0) {

			bx[m.V.rows() + count2] = m.V(i, 0) * vertex_weight;
			by[m.V.rows() + count2] = m.V(i, 1) * vertex_weight;
			bz[m.V.rows() + count2] = m.V(i, 2) * vertex_weight;
			count2++;

			for(int j = 0; j < neighbors[i].size(); j++) {
				if (m.vertex_is_fixed[neighbors[i][j]] == 0) {

					Eigen::Vector3d vec = (m.V.row(i) - m.V.row(neighbors[i][j]));
					vec.normalize();
					Eigen::Vector3d edge_vec = (target_edge_lengths[i] + target_edge_lengths[neighbors[i][j]]) / 2.0 * vec;

					bx[m.V.rows() + no_boundary_vertices + count] = edge_vec(0) * edge_weight;
					by[m.V.rows() + no_boundary_vertices + count] = edge_vec(1) * edge_weight;
					bz[m.V.rows() + no_boundary_vertices + count] = edge_vec(2) * edge_weight;
					count++;
				}
			}
		}
	}

	Eigen::VectorXd Vnewx = (*solver2_array[m.ID - 1]).solve(AT*bx);
	Eigen::VectorXd Vnewy = (*solver2_array[m.ID - 1]).solve(AT*by);
	Eigen::VectorXd Vnewz = (*solver2_array[m.ID - 1]).solve(AT*bz);

	for(int i = 0; i < m.V.rows(); i++) {
		if (m.vertex_is_fixed[i] == 0) { //Only update non-fixed points (curve points are fixed)
			m.V.row(i) << Vnewx[i], Vnewy[i], Vnewz[i];
		}
	}
}

Eigen::MatrixX3d SurfaceSmoothing::compute_vertex_laplacians(Mesh &m) {
	Eigen::MatrixX3d laplacians(m.V.rows(),3);
    Eigen::Vector3d vec;
    int nr_neighbors;
    for(int i = 0; i < m.V.rows(); i++) {
		vec = Eigen::Vector3d::Zero();
        nr_neighbors = neighbors[i].size();
		for(int j = 0; j < nr_neighbors; j++) {
			vec += m.V.row(neighbors[i][j]);
		}
		laplacians.row(i) = m.V.row(i) - (vec.transpose()/nr_neighbors);
	}
	return laplacians;
}

Eigen::VectorXd SurfaceSmoothing::get_curvatures(Mesh &m) {
	Eigen::VectorXd curvatures(m.V.rows());
	Eigen::MatrixX3d& laplacians = compute_vertex_laplacians(m);
    int sign;
	for (int i = 0; i < m.V.rows(); i++) {
		curvatures[i] = laplacians.row(i).dot(vertex_normals.row(i).transpose());
	}
	return curvatures;
}

Eigen::VectorXi SurfaceSmoothing::points_on_border(Mesh& m) {
	Eigen::MatrixXi EV, FE, EF;
	igl::edge_topology(m.V, m.F, EV, FE, EF);

	Eigen::VectorXi point_on_border(m.V.rows());
	point_on_border.setZero();

	std::vector<int> tmp(m.mesh_to_patch_indices.data(), m.mesh_to_patch_indices.data() + m.mesh_to_patch_indices.size());
	for (int i = 0; i < m.sharp_edge.rows(); i++) {
		if (m.sharp_edge[i]) {	
			point_on_border[EV(i, 0)] = 1;
			point_on_border[EV(i, 1)] = 1;
		}
	}
	return point_on_border;
}

Eigen::SparseMatrix<double> SurfaceSmoothing::get_precompute_matrix_for_LM_and_edges(int ID) {
	return precompute_matrix_for_LM_and_edges[ID];
}

void SurfaceSmoothing::set_precompute_matrix_for_LM_and_edges(int ID, Eigen::SparseMatrix<double> LM_edge) {
	precompute_matrix_for_LM_and_edges[ID] = LM_edge;
}

Eigen::SparseMatrix<double> SurfaceSmoothing::get_AT_for_LM_and_edges(int ID) {
	return AT_for_LM_and_edges[ID];
}

void SurfaceSmoothing::set_AT_for_LM_and_edges(int ID, Eigen::SparseMatrix<double> AT) {
	AT_for_LM_and_edges[ID] = AT;
}

Eigen::SparseMatrix<double> SurfaceSmoothing::get_precompute_matrix_for_positions(int ID) {
	return precompute_matrix_for_positions[ID];
}

void SurfaceSmoothing::set_precompute_matrix_for_positions(int ID, Eigen::SparseMatrix<double> pos) {
	precompute_matrix_for_positions[ID] = pos;
}

Eigen::SparseMatrix<double> SurfaceSmoothing::get_AT_for_positions(int ID) {
	return AT_for_positions[ID];
}

void SurfaceSmoothing::set_AT_for_positions(int ID, Eigen::SparseMatrix<double> AT) {
	AT_for_positions[ID] = AT;
}