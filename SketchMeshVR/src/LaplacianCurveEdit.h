#ifndef _Laplacian_Curve_Edit_H_
#define _Laplacian_Curve_Edit_H_
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <Eigen/OrderingMethods>

#include <vector>
#include <unordered_map>

class LaplacianCurveEdit {
public:
	void setup_for_update_curve(std::vector<int> vertices, std::vector<int> fixed_vertices, Eigen::MatrixXi&  edges, std::vector<int> fixed_edges, Eigen::MatrixXi vertex_triplets, Eigen::MatrixXi edge_triplets, Eigen::MatrixXd& V, Eigen::MatrixXi& EV_);
	void update_curve(Eigen::MatrixXd & V);
	LaplacianCurveEdit();
	LaplacianCurveEdit(const LaplacianCurveEdit &);
	LaplacianCurveEdit& operator=(LaplacianCurveEdit other);

private:
	int CONSTRAINT_WEIGHT = 10000;
	Eigen::SparseMatrix<double> A, A_L1_T;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solverL1; //Solver for final vertex positions (with L1)
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solverPosRot;
	Eigen::VectorXd B, PosRot;
	Eigen::MatrixXd original_L0, original_L1;
	Eigen::MatrixXi EV;
	Eigen::VectorXi is_fixed;
	std::vector<Eigen::Matrix3d> Rot;
	std::vector<std::vector<int>> vertex_triplet_to_rot_idx;

	std::vector<int> fixed_vertices_local;
	std::unordered_map<int, int> vertex_global_to_local;
	std::unordered_map<int, int> edge_global_to_local;

	std::vector<int> vertices, fixed_vertices, fixed_edges;
	Eigen::MatrixXi edges, edge_triplets, vertex_triplets;
	Eigen::MatrixXd V;

	void setup_for_L1_position_step();
	void solve_for_pos_and_rot(Eigen::MatrixXd & V);
	void update_rot();
	Eigen::Matrix3d compute_orthonormal(Eigen::Matrix3d & rot);
	void final_L1_pos(Eigen::MatrixXd & V);
	Eigen::Matrix3d average_rot(Eigen::Matrix3d & r0, Eigen::Matrix3d & r1);
	int find_edge(int start, int end);

	bool curve_structure_was_updated;
};


#endif