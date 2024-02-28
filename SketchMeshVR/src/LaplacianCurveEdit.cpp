#include "LaplacianCurveEdit.h"

#include <iostream>

typedef Eigen::Triplet<double> T;

LaplacianCurveEdit::LaplacianCurveEdit() {}
LaplacianCurveEdit::LaplacianCurveEdit(const LaplacianCurveEdit& other) {}
LaplacianCurveEdit& LaplacianCurveEdit::operator=(LaplacianCurveEdit other) { return *this; }

void LaplacianCurveEdit::setup_for_update_curve(std::vector<int> vertices_, std::vector<int> fixed_vertices_, Eigen::MatrixXi& edges_, std::vector<int> fixed_edges_, Eigen::MatrixXi vertex_triplets_, Eigen::MatrixXi edge_triplets_, Eigen::MatrixXd& V_, Eigen::MatrixXi& EV_) {
	V = V_;
	EV = EV_;
	vertices = vertices_;
	fixed_vertices = fixed_vertices_;
	edges = edges_;
	fixed_edges = fixed_edges_;
	vertex_triplets = vertex_triplets_;
	edge_triplets = edge_triplets_;

	B.resize(edges.rows() * 3 + edge_triplets.rows() * 9 + fixed_vertices.size() * 3 + fixed_edges.size() * 3);
	Rot.resize(edges.rows());
	original_L0.resize(edges.rows(), 3);
	vertex_global_to_local.clear();
	curve_structure_was_updated = true;

	for (int i = 0; i < edges.rows(); i++) {
		original_L0.row(i) = V.row(edges(i, 1)) - V.row(edges(i, 0));
	}

	for (int i = 0; i < vertices.size(); i++) {
		vertex_global_to_local.insert({vertices[i], i });
	}

	int idx;
	fixed_vertices_local.clear();
	fixed_vertices_local.resize(fixed_vertices.size());
	is_fixed.resize(vertices.size());
	is_fixed.setZero();
	for (int i = 0; i < fixed_vertices.size(); i++) { //Create a copy with local indexing for fixed_vertices, and also keep the global indexing
		auto loc = std::find(vertices.begin(), vertices.end(), fixed_vertices[i]);
		idx = distance(vertices.begin(), loc);
		fixed_vertices_local[i] = idx;
		is_fixed[idx] = 1;
	}
	//Fixed_edges already has local indexing into edges

	setup_for_L1_position_step();
}

void LaplacianCurveEdit::setup_for_L1_position_step() {
	std::vector<T> tripletList;
	tripletList.reserve(vertex_triplets.rows() * 3 + fixed_vertices.size());
	original_L1.resize(vertex_triplets.rows(), 3);
	vertex_triplet_to_rot_idx.resize(vertex_triplets.rows(), std::vector<int>(2));

	for (int i = 0; i < vertex_triplets.rows(); i++) {
		int v0 = vertex_triplets(i, 0);
		int v1 = vertex_triplets(i, 1);
		int v2 = vertex_triplets(i, 2);

		original_L1.row(i) = V.row(v1) - (0.5*V.row(v0) + 0.5*V.row(v2));
		vertex_triplet_to_rot_idx[i][0] = find_edge(v0, v1);
		vertex_triplet_to_rot_idx[i][1] = find_edge(v1, v2);

		tripletList.push_back(T(i, vertex_global_to_local.at(v0), -0.5));
		tripletList.push_back(T(i, vertex_global_to_local.at(v1), 1.0));
		tripletList.push_back(T(i, vertex_global_to_local.at(v2), -0.5));
	}

	for (int i = 0; i < fixed_vertices_local.size(); i++) {
		tripletList.push_back(T(vertex_triplets.rows() + i, fixed_vertices_local[i], CONSTRAINT_WEIGHT));
	}
	Eigen::SparseMatrix<double> A_L1(vertex_triplets.rows() + fixed_vertices_local.size(), vertices.size());
	A_L1.setFromTriplets(tripletList.begin(), tripletList.end());

	A_L1_T = A_L1.transpose();
	solverL1.compute(A_L1_T*A_L1);
}

void LaplacianCurveEdit::update_curve(Eigen::MatrixXd& V) {
	for (int i = 0; i < edges.rows(); i++) {
		Rot[i] = Eigen::Matrix3d::Identity();
	}

	for (int i = 0; i < 2; i++) {
		try {
			solve_for_pos_and_rot(V);
		}
		catch (int ex) {
			throw -1;
		}
		update_rot();
	}
	final_L1_pos(V);
}

void LaplacianCurveEdit::solve_for_pos_and_rot(Eigen::MatrixXd& V) {
	A.setZero();
	B.setZero();
	std::vector<T> tripletList;
	tripletList.reserve(edges.rows() * 15 + edge_triplets.rows() * 3 * 27 + fixed_vertices_local.size() * 3 + fixed_edges.size() * 3);

	int v0, v1;
	for (int i = 0; i < edges.rows(); i++) { //Laplacian of the position: v1' - v0' = dR * R * L(v1 - v0)
		v0 = vertex_global_to_local.at(edges(i, 0));
		v1 = vertex_global_to_local.at(edges(i, 1));

		tripletList.push_back(T(i * 3 + 0, v0 * 3, -1)); //L0
		tripletList.push_back(T(i * 3 + 0, v1 * 3, 1)); //L0
		tripletList.push_back(T(i * 3 + 0, vertices.size() * 3 + i * 3, 0));
		tripletList.push_back(T(i * 3 + 0, vertices.size() * 3 + i * 3 + 1, -(Rot[i].row(2).dot(original_L0.row(i)))));
		tripletList.push_back(T(i * 3 + 0, vertices.size() * 3 + i * 3 + 2, Rot[i].row(1).dot(original_L0.row(i))));
		B[i * 3] = Rot[i].row(0).dot(original_L0.row(i)); //Constant

		tripletList.push_back(T(i * 3 + 1, v0 * 3 + 1, -1)); //L0
		tripletList.push_back(T(i * 3 + 1, v1 * 3 + 1, 1)); //L0
		tripletList.push_back(T(i * 3 + 1, vertices.size() * 3 + i * 3, Rot[i].row(2).dot(original_L0.row(i))));
		tripletList.push_back(T(i * 3 + 1, vertices.size() * 3 + i * 3 + 1, 0));
		tripletList.push_back(T(i * 3 + 1, vertices.size() * 3 + i * 3 + 2, -Rot[i].row(0).dot(original_L0.row(i))));
		B[i * 3 + 1] = Rot[i].row(1).dot(original_L0.row(i)); //Constant

		tripletList.push_back(T(i * 3 + 2, v0 * 3 + 2, -1)); //L0
		tripletList.push_back(T(i * 3 + 2, v1 * 3 + 2, 1)); //L0
		tripletList.push_back(T(i * 3 + 2, vertices.size() * 3 + i * 3, -Rot[i].row(1).dot(original_L0.row(i)))); //Note: this is different from teddy, but I think teddy is wrong here
		tripletList.push_back(T(i * 3 + 2, vertices.size() * 3 + i * 3 + 1, Rot[i].row(0).dot(original_L0.row(i))));
		tripletList.push_back(T(i * 3 + 2, vertices.size() * 3 + i * 3 + 2, 0));
		B[i * 3 + 2] = Rot[i].row(2).dot(original_L0.row(i)); //Constant
	}

	int e0, e1, e2;
	for (int i = 0; i < edge_triplets.rows(); i++) { //Laplacian of the rotation matrix: L dR*R = 0 (r_i*R_i - r_j*R_j = 0 for i,j in Edges)
		e0 = edge_triplets(i, 0);
		e1 = edge_triplets(i, 1);
		e2 = edge_triplets(i, 2);

		for (int j = 0; j < 3; j++) {
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + j, vertices.size() * 3 + e0 * 3 + 0, 0));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + j, vertices.size() * 3 + e0 * 3 + 1, -0.5 * Rot[e0](2, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + j, vertices.size() * 3 + e0 * 3 + 2, -0.5 * -Rot[e0](1, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + j, vertices.size() * 3 + e1 * 3 + 0, 0));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + j, vertices.size() * 3 + e1 * 3 + 1, Rot[e1](2, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + j, vertices.size() * 3 + e1 * 3 + 2, -Rot[e1](1, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + j, vertices.size() * 3 + e2 * 3 + 0, 0));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + j, vertices.size() * 3 + e2 * 3 + 1, -0.5 * Rot[e2](2, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + j, vertices.size() * 3 + e2 * 3 + 2, -0.5 * -Rot[e2](1, j)));
			B[edges.rows() * 3 + i * 9 + j] = 0.5 * Rot[e0](0, j) - 1 * Rot[e1](0, j) + 0.5 * Rot[e2](0, j);

			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 3 + j, vertices.size() * 3 + e0 * 3 + 0, -0.5 * -Rot[e0](2, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 3 + j, vertices.size() * 3 + e0 * 3 + 1, 0));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 3 + j, vertices.size() * 3 + e0 * 3 + 2, -0.5 * Rot[e0](0, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 3 + j, vertices.size() * 3 + e1 * 3 + 0, -Rot[e1](2, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 3 + j, vertices.size() * 3 + e1 * 3 + 1, 0));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 3 + j, vertices.size() * 3 + e1 * 3 + 2, Rot[e1](0, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 3 + j, vertices.size() * 3 + e2 * 3 + 0, -0.5 * -Rot[e2](2, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 3 + j, vertices.size() * 3 + e2 * 3 + 1, 0));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 3 + j, vertices.size() * 3 + e2 * 3 + 2, -0.5 * Rot[e2](0, j)));
			B[edges.rows() * 3 + i * 9 + 3 + j] = 0.5 * Rot[e0](1, j) - 1 * Rot[e1](1, j) + 0.5 * Rot[e2](1, j);

			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 6 + j, vertices.size() * 3 + e0 * 3 + 0, -0.5*Rot[e0](1, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 6 + j, vertices.size() * 3 + e0 * 3 + 1, -0.5 * -Rot[e0](0, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 6 + j, vertices.size() * 3 + e0 * 3 + 2, 0));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 6 + j, vertices.size() * 3 + e1 * 3 + 0, Rot[e1](1, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 6 + j, vertices.size() * 3 + e1 * 3 + 1, -Rot[e1](0, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 6 + j, vertices.size() * 3 + e1 * 3 + 2, 0));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 6 + j, vertices.size() * 3 + e2 * 3 + 0, -0.5*Rot[e2](1, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 6 + j, vertices.size() * 3 + e2 * 3 + 1, -0.5 * -Rot[e2](0, j)));
			tripletList.push_back(T(edges.rows() * 3 + i * 9 + 6 + j, vertices.size() * 3 + e2 * 3 + 2, 0));
			B[edges.rows() * 3 + i * 9 + 6 + j] = 0.5 * Rot[e0](2, j) - 1 * Rot[e1](2, j) + 0.5 * Rot[e2](2, j);
		}
	}

	for (int i = 0; i < fixed_vertices_local.size(); i++) { //Position of fixed vertices: v_i - v_i' = 0
		tripletList.push_back(T(edges.rows() * 3 + edge_triplets.rows() * 9 + i * 3, fixed_vertices_local[i] * 3, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(edges.rows() * 3 + edge_triplets.rows() * 9 + i * 3 + 1, fixed_vertices_local[i] * 3 + 1, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(edges.rows() * 3 + edge_triplets.rows() * 9 + i * 3 + 2, fixed_vertices_local[i] * 3 + 2, CONSTRAINT_WEIGHT));

		B[edges.rows() * 3 + edge_triplets.rows() * 9 + i * 3] = CONSTRAINT_WEIGHT * V(fixed_vertices[i], 0);
		B[edges.rows() * 3 + edge_triplets.rows() * 9 + i * 3 + 1] = CONSTRAINT_WEIGHT * V(fixed_vertices[i], 1);
		B[edges.rows() * 3 + edge_triplets.rows() * 9 + i * 3 + 2] = CONSTRAINT_WEIGHT * V(fixed_vertices[i], 2);
	}

	for (int i = 0; i < fixed_edges.size(); i++) {
		tripletList.push_back(T(edges.rows() * 3 + edge_triplets.rows() * 9 + fixed_vertices_local.size() * 3 + i * 3 + 0, vertices.size() * 3 + fixed_edges[i] * 3, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(edges.rows() * 3 + edge_triplets.rows() * 9 + fixed_vertices_local.size() * 3 + i * 3 + 1, vertices.size() * 3 + fixed_edges[i] * 3 + 1, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(edges.rows() * 3 + edge_triplets.rows() * 9 + fixed_vertices_local.size() * 3 + i * 3 + 2, vertices.size() * 3 + fixed_edges[i] * 3 + 2, CONSTRAINT_WEIGHT));

		B[edges.rows() * 3 + edge_triplets.rows() * 9 + fixed_vertices_local.size() * 3 + i * 3 + 0] = 0;
		B[edges.rows() * 3 + edge_triplets.rows() * 9 + fixed_vertices_local.size() * 3 + i * 3 + 1] = 0;
		B[edges.rows() * 3 + edge_triplets.rows() * 9 + fixed_vertices_local.size() * 3 + i * 3 + 2] = 0;
	}

	A = Eigen::SparseMatrix<double>(edges.rows() * 3 + edge_triplets.rows() * 9 + fixed_vertices_local.size() * 3 + fixed_edges.size() * 3, vertices.size() *3 + edges.rows() * 3);
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	Eigen::SparseMatrix<double> A_before = A;
	A.prune(0.0);

	Eigen::SparseMatrix<double> AT = A.transpose();
	if (curve_structure_was_updated) {
		solverPosRot.analyzePattern(AT*A);
		curve_structure_was_updated = false;
	} 
	solverPosRot.factorize(AT*A);
	if (solverPosRot.info() != Eigen::Success) {
		throw -1;
	}
	PosRot = solverPosRot.solve(AT*B);
}

void LaplacianCurveEdit::update_rot() {
	Eigen::Matrix3d newRot;
	double rx, ry, rz;
	for (int i = 0; i < Rot.size(); i++) {
		rx = PosRot[vertices.size() * 3 + i * 3];
		ry = PosRot[vertices.size() * 3 + i * 3 + 1];
		rz = PosRot[vertices.size() * 3 + i * 3 + 2];

		newRot.row(0) << 1, -rz, ry;
		newRot.row(1) << rz, 1, -rx;
		newRot.row(2) << -ry, rx, 1;

		Rot[i] = newRot * Rot[i]; //Safe for matrix-matrix multiplication with Eigen
		Rot[i] = compute_orthonormal(Rot[i]);
	}

}

Eigen::Matrix3d LaplacianCurveEdit::compute_orthonormal(Eigen::Matrix3d& rot) {
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(rot, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d VT = svd.matrixV().transpose();
	return svd.matrixU()*VT;
}

void LaplacianCurveEdit::final_L1_pos(Eigen::MatrixXd &V) {
	Eigen::VectorXd Bx(vertex_triplets.rows() + fixed_vertices_local.size()), By(vertex_triplets.rows() + fixed_vertices_local.size()), Bz(vertex_triplets.rows() + fixed_vertices_local.size());
	Eigen::Vector3d Rl;
	Eigen::Matrix3d r0, r1;
	for (int i = 0; i < vertex_triplets.rows(); i++) {
		r0 = Rot[vertex_triplet_to_rot_idx[i][0]];
		r1 = Rot[vertex_triplet_to_rot_idx[i][1]];

		Rl = average_rot(r0, r1)* original_L1.row(i).transpose();
		Bx[i] = Rl(0);
		By[i] = Rl(1);
		Bz[i] = Rl(2);
	}

	for (int i = 0; i < fixed_vertices.size(); i++) {
		Bx[vertex_triplets.rows() + i] = V(fixed_vertices[i], 0) *CONSTRAINT_WEIGHT;
		By[vertex_triplets.rows() + i] = V(fixed_vertices[i], 1) *CONSTRAINT_WEIGHT;
		Bz[vertex_triplets.rows() + i] = V(fixed_vertices[i], 2) *CONSTRAINT_WEIGHT;
	}

	Eigen::VectorXd xpos = solverL1.solve(A_L1_T*Bx);
	Eigen::VectorXd ypos = solverL1.solve(A_L1_T*By);
	Eigen::VectorXd zpos = solverL1.solve(A_L1_T*Bz);
	for (int i = 0; i < vertices.size(); i++) { //update the position of non-fixed vertices of the stroke that is being pulled on
		if (!is_fixed[i]) {
			V.row(vertices[i]) << xpos[i], ypos[i], zpos[i];
		}
	}
}

Eigen::Matrix3d LaplacianCurveEdit::average_rot(Eigen::Matrix3d& r0, Eigen::Matrix3d& r1) {
	Eigen::Matrix3d rnew;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			rnew(i, j) = (r0(i, j) + r1(i, j)) / 2;
		}
		double length = rnew.row(i).norm();
		for (int j = 0; j < 3; j++) {
			rnew(i, j) /= length;
		}
	}
	return rnew;
}

int LaplacianCurveEdit::find_edge(int start, int end) {
	for (int i = 0; i < edges.rows(); i++) {
		if (edges(i, 0) == start) {
			if (edges(i, 1) == end) {
				return i;
			}
		}
		else if (edges(i, 0) == end) {
			if (edges(i, 1) == start) {
				return i;
			}
		}
	}
	return -1;
}