#include "CurveDeformation.h"

using namespace std;
using namespace igl;

int moving_vertex_ID, handle_ID, CONSTRAINT_WEIGHT = 10000, stroke_ID, no_vertices, no_ROI_vert = -1;
double current_ROI_size, curve_diag_length;
bool CurveDeformation::smooth_deform_mode, stroke_is_loop, prev_loop_type;

vector<Eigen::Matrix3d> Rot;
vector<int> vert_bindings;

Eigen::RowVector3d start_pos;
Eigen::VectorXi fixed_indices, fixed_indices_local, is_fixed; //Indicates for every vertex on the pulled curve whether or not it is fixed (aka it is 0 when it is in the ROI and 1 when it is not);
Eigen::VectorXd B, PosRot;
Eigen::MatrixXd original_L0, original_L1;

Eigen::SparseMatrix<double> A, A_L1_T;
Eigen::SparseLU<Eigen::SparseMatrix<double>> solverL1; //Solver for final vertex positions (with L1)
Eigen::SparseLU<Eigen::SparseMatrix<double>> solverPosRot;


void CurveDeformation::startPullCurve(Stroke& _stroke, int _handle_ID) {
	start_pos = (_stroke.get3DPoints()).row(_handle_ID);
	moving_vertex_ID = _stroke.get_vertex_idx_for_point(_handle_ID); //The global vertex index (in V) for the moving vertex
	current_ROI_size = 0.0;
	curve_diag_length = compute_curve_diag_length(_stroke);
	handle_ID = _handle_ID;
	no_vertices = _stroke.get3DPoints().rows() - 1; //We should ignore the last one, which is a copy of the first one
	no_ROI_vert = -1;
	Rot.resize(no_vertices);
	vert_bindings = _stroke.get_closest_vert_bindings();
	stroke_is_loop = _stroke.is_loop;
	stroke_ID = _stroke.get_ID();
}

//pos is the unprojection from the position to where the user dragged the vertex
void CurveDeformation::pullCurve(Eigen::RowVector3d& pos, Eigen::MatrixXd& V, Eigen::VectorXi& part_of_original_stroke) {
	double drag_size = (pos - start_pos).norm();
	drag_size /= curve_diag_length;
	bool ROI_is_updated = false;
	if(!current_ROI_size || (fabs(drag_size - current_ROI_size) > 0.01)) { //Take the deformation and roi size as percentages
		ROI_is_updated = update_ROI(drag_size);
	}
	if(no_ROI_vert == 0 || current_ROI_size == 0) { //If we have no other "free" vertices other than the handle vertex, we simply move it to the target position
		V.row(moving_vertex_ID) = pos;
	} else {
		if(ROI_is_updated) {
			setup_for_update_curve(V);
		}
		V.row(moving_vertex_ID) = pos;
		update_curve(V, part_of_original_stroke);
	}
}

double CurveDeformation::compute_curve_diag_length(Stroke& _stroke) {
	Eigen::Vector3d maxBB = _stroke.get3DPoints().colwise().maxCoeff();
	Eigen::Vector3d minBB = _stroke.get3DPoints().colwise().minCoeff();
	return (maxBB - minBB).norm();
}

bool CurveDeformation::update_ROI(double drag_size) {
	int no_ROI_vert_tmp;
	if(smooth_deform_mode && no_vertices > 1) {
		no_ROI_vert_tmp = max(round(0.17*no_vertices) + 1, min(round(drag_size * no_vertices) + 1, ceil(((no_vertices - 1) / 2) - 1))); //Determine how many vertices to the left and to the right to have free (at most half-1 of all vertices on each side, always at least 1/6th of the number of vertices + 1 vertex fixed, to take care of thin meshes with many vertices)
	} else {
		no_ROI_vert_tmp = max(0.0, min(round(drag_size / 4.0 * no_vertices), ceil(((no_vertices - 1) / 2) - 1))); //Determine how many vertices to the left and to the right to have free (at most half-1 of all vertices on each side)
	}

	if(((no_ROI_vert == no_ROI_vert_tmp) && prev_loop_type == stroke_is_loop) || no_ROI_vert_tmp == 0) { //number of vertices in ROI didn't change
		return false;
	}

	current_ROI_size = drag_size;
	no_ROI_vert = no_ROI_vert_tmp;
	prev_loop_type = stroke_is_loop;

	int ROI_1, ROI_2;
	compute_ROI_boundaries(ROI_1, ROI_2);
	setup_fixed_indices(ROI_1, ROI_2);
	return true;
}

void CurveDeformation::compute_ROI_boundaries(int& ROI_1, int& ROI_2) {
	if(stroke_is_loop) {
		ROI_1 = (((handle_ID - no_ROI_vert) + no_vertices) % no_vertices);
		ROI_2 = (((handle_ID + no_ROI_vert) + no_vertices) % no_vertices);
	} else {
		ROI_1 = max(0, handle_ID - no_ROI_vert); //For non-loop strokes (e.g. added strokes), don't wrap around but instead cap at the first and last vertex of the stroke
		ROI_2 = min(no_vertices - 1, handle_ID + no_ROI_vert);
	}
}

void CurveDeformation::setup_fixed_indices(int ROI_1, int ROI_2) {
	vector<int> fixed, fixed_local;
	is_fixed = Eigen::VectorXi::Zero(no_vertices);

	if(ROI_1 < ROI_2) {
		for(int i = 0; i < ROI_1; i++) {
			is_fixed[i] = 1;
			fixed.push_back(vert_bindings[i]);
			fixed_local.push_back(i);
		}
		for(int i = ROI_2 + 1; i < no_vertices; i++) {
			is_fixed[i] = 1;
			fixed.push_back(vert_bindings[i]);
			fixed_local.push_back(i);
		}
	} else {
		for(int i = ROI_2 + 1; i < ROI_1; i++) {
			is_fixed[i] = 1;
			fixed.push_back(vert_bindings[i]);
			fixed_local.push_back(i);
		}
	}

	fixed.push_back(moving_vertex_ID);
	fixed_local.push_back(handle_ID);
	is_fixed[handle_ID] = 1;

	fixed_indices = Eigen::VectorXi::Map(fixed.data(), fixed.size());
	fixed_indices_local = Eigen::VectorXi::Map(fixed_local.data(), fixed_local.size());

}

void CurveDeformation::setup_for_update_curve(Eigen::MatrixXd& V) {
	A.resize(no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + fixed_indices.size() * 3, no_vertices * 3 + no_vertices * 3); //Solve for x,y,z simultaneously since we cannot reuse A anyway
	B.resize(no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + fixed_indices.size() * 3);
	original_L0.resize(no_vertices, 3);

	if(stroke_is_loop) {
		for(int i = 0; i < no_vertices; i++) {
			original_L0.row(i) = V.row(vert_bindings[i]) - V.row(((vert_bindings[((i - 1) + no_vertices) % no_vertices]))); //This assumes that the stroke is looped, which might not always be true for added control strokes.
		}
	} else {
		original_L0.row(0) = V.row(vert_bindings[0]) - V.row(vert_bindings[1]);
		for(int i = 1; i < no_vertices; i++) {
			original_L0.row(i) = V.row(vert_bindings[i]) - V.row(vert_bindings[i - 1]);
		}

	}
	setup_for_L1_position_step(V);
}

void CurveDeformation::setup_for_L1_position_step(Eigen::MatrixXd& V) {
	Eigen::SparseMatrix<double> A_L1(no_vertices + fixed_indices.size(), no_vertices);
	original_L1.resize(no_vertices, 3);
	int cur, next, prev;
	if(stroke_is_loop) {
		for(int i = 0; i < no_vertices; i++) {
			cur = vert_bindings[i];
			next = vert_bindings[(i + 1) % no_vertices];
			prev = vert_bindings[(((i - 1) + no_vertices) % no_vertices)];

			A_L1.insert(i, (((i - 1) + no_vertices) % no_vertices)) = -0.5;
			A_L1.insert(i, i) = 1;
			A_L1.insert(i, (i + 1) % no_vertices) = -0.5;

			original_L1.row(i) = V.row(cur) - (0.5*V.row(next) + 0.5*V.row(prev));
		}
	} else {
		A_L1.insert(0, 0) = 1;
		A_L1.insert(0, 1) = -1;
		original_L1.row(0) = V.row(vert_bindings[0]) - V.row(vert_bindings[1]);

		A_L1.insert(no_vertices - 1, no_vertices - 1) = 1;
		A_L1.insert(no_vertices - 1, no_vertices - 2) = -1;
		original_L1.row(no_vertices - 1) = V.row(vert_bindings[no_vertices - 1]) - V.row(vert_bindings[no_vertices - 2]);

		for(int i = 1; i < no_vertices - 1; i++) {
			cur = vert_bindings[i];
			next = vert_bindings[(i + 1) % no_vertices];
			prev = vert_bindings[(((i - 1) + no_vertices) % no_vertices)];

			A_L1.insert(i, i - 1) = -0.5;
			A_L1.insert(i, i) = 1;
			A_L1.insert(i, i + 1) = -0.5;

			original_L1.row(i) = V.row(cur) - (0.5*V.row(next) + 0.5*V.row(prev));
		}
	}

	for(int i = 0; i < fixed_indices.size(); i++) {
		A_L1.insert(no_vertices + i, fixed_indices_local[i]) = CONSTRAINT_WEIGHT;
	}

	A_L1_T = A_L1.transpose();
	solverL1.compute(A_L1_T*A_L1);
}

void CurveDeformation::update_curve(Eigen::MatrixXd& V, Eigen::VectorXi & part_of_original_stroke) {
	for(int i = 0; i < no_vertices; i++) {
		Rot[i] = Eigen::Matrix3d::Identity();
	}

	for(int i = 0; i < 2; i++) {
		solve_for_pos_and_rot(V);
		update_rot();
	}
	final_L1_pos(V, part_of_original_stroke);
}

void CurveDeformation::solve_for_pos_and_rot(Eigen::MatrixXd& V) {
	A.setZero();
	B.setZero();

	int prev, cur;
	int start = !stroke_is_loop; //Make sure to skip wrapping around for index 0
	for(int i = 0; i < no_vertices; i++) {
		prev = (((i - 1) + no_vertices) % no_vertices);
		cur = i;
		if(i > 0 || stroke_is_loop) {
			A.insert(i * 3 + 0, prev * 3) = -1; //L0
			A.insert(i * 3 + 0, cur * 3) = 1; //L0
		}
		A.insert(i * 3, no_vertices * 3 + cur * 3) = 0;
		A.insert(i * 3, no_vertices * 3 + cur * 3 + 1) = -(Rot[i].row(2).dot(original_L0.row(i)));
		A.insert(i * 3, no_vertices * 3 + cur * 3 + 2) = Rot[i].row(1).dot(original_L0.row(i));
		B[i * 3] = Rot[i].row(0).dot(original_L0.row(i)); //constant

		if(i > 0 || stroke_is_loop) {
			A.insert(i * 3 + 1, prev * 3 + 1) = -1; //L0
			A.insert(i * 3 + 1, cur * 3 + 1) = 1; //L0
		}
		A.insert(i * 3 + 1, no_vertices * 3 + cur * 3) = Rot[i].row(2).dot(original_L0.row(i));
		A.insert(i * 3 + 1, no_vertices * 3 + cur * 3 + 1) = 0;
		A.insert(i * 3 + 1, no_vertices * 3 + cur * 3 + 2) = -Rot[i].row(0).dot(original_L0.row(i));
		B[i * 3 + 1] = Rot[i].row(1).dot(original_L0.row(i));

		if(i > 0 || stroke_is_loop) {
			A.insert(i * 3 + 2, prev * 3 + 2) = -1;//L0
			A.insert(i * 3 + 2, cur * 3 + 2) = 1;//L0
		}
		A.insert(i * 3 + 2, no_vertices * 3 + cur * 3) = -Rot[i].row(1).dot(original_L0.row(i));
		A.insert(i * 3 + 2, no_vertices * 3 + cur * 3 + 1) = Rot[i].row(0).dot(original_L0.row(i));
		A.insert(i * 3 + 2, no_vertices * 3 + cur * 3 + 2) = 0;
		B[i * 3 + 2] = Rot[i].row(2).dot(original_L0.row(i));
	}

	//Setup for r_i*R_i - r_j*R_j = 0 for i,j in Edges
	int prev_rot;
	for(int i = start; i < no_vertices; i++) { //For non-closed strokes, start at index 1 (we can't fill in a partial entry since they depend on eachother, unlike the previous block)
		prev = (((i - 1) + no_vertices) % no_vertices);
		prev_rot = (((i - 1) + no_vertices) % no_vertices);
		for(int j = 0; j < 3; j++) {
			A.insert(no_vertices * 3 + i * 9 + j, no_vertices * 3 + prev * 3) = 0;
			A.insert(no_vertices * 3 + i * 9 + j, no_vertices * 3 + prev * 3 + 1) = -1 * Rot[prev_rot](2, j);
			A.insert(no_vertices * 3 + i * 9 + j, no_vertices * 3 + prev * 3 + 2) = -1 * -Rot[prev_rot](1, j);
			A.insert(no_vertices * 3 + i * 9 + j, no_vertices * 3 + i * 3) = 0;
			A.insert(no_vertices * 3 + i * 9 + j, no_vertices * 3 + i * 3 + 1) = 1 * Rot[i](2, j);
			A.insert(no_vertices * 3 + i * 9 + j, no_vertices * 3 + i * 3 + 2) = 1 * -Rot[i](1, j);
			B[no_vertices * 3 + i * 9 + j] = -1 * Rot[i](0, j) + 1 * Rot[prev_rot](0, j);

			A.insert(no_vertices * 3 + i * 9 + 3 + j, no_vertices * 3 + prev * 3) = -1 * -Rot[prev_rot](2, j);
			A.insert(no_vertices * 3 + i * 9 + 3 + j, no_vertices * 3 + prev * 3 + 1) = 0;
			A.insert(no_vertices * 3 + i * 9 + 3 + j, no_vertices * 3 + prev * 3 + 2) = -1 * Rot[prev_rot](0, j);
			A.insert(no_vertices * 3 + i * 9 + 3 + j, no_vertices * 3 + i * 3) = 1 * -Rot[i](2, j);
			A.insert(no_vertices * 3 + i * 9 + 3 + j, no_vertices * 3 + i * 3 + 1) = 0;
			A.insert(no_vertices * 3 + i * 9 + 3 + j, no_vertices * 3 + i * 3 + 2) = 1 * Rot[i](0, j);
			B[no_vertices * 3 + i * 9 + 3 + j] = -1 * Rot[i](1, j) + 1 * Rot[prev_rot](1, j);

			A.insert(no_vertices * 3 + i * 9 + 6 + j, no_vertices * 3 + prev * 3) = -1 * Rot[prev_rot](1, j);
			A.insert(no_vertices * 3 + i * 9 + 6 + j, no_vertices * 3 + prev * 3 + 1) = -1 * -Rot[prev_rot](0, j);
			A.insert(no_vertices * 3 + i * 9 + 6 + j, no_vertices * 3 + prev * 3 + 2) = 0;
			A.insert(no_vertices * 3 + i * 9 + 6 + j, no_vertices * 3 + i * 3) = 1 * Rot[i](1, j);
			A.insert(no_vertices * 3 + i * 9 + 6 + j, no_vertices * 3 + i * 3 + 1) = 1 * -Rot[i](0, j);
			A.insert(no_vertices * 3 + i * 9 + 6 + j, no_vertices * 3 + i * 3 + 2) = 0;
			B[no_vertices * 3 + i * 9 + 6 + j] = -1 * Rot[i](2, j) + 1 * Rot[prev_rot](2, j);
		}
	}

	//Setup for v_i - v_i' = 0
	for(int i = 0; i < fixed_indices.size(); i++) {
		A.insert(no_vertices * 3 + no_vertices * 9 + i * 3, fixed_indices_local[i] * 3) = CONSTRAINT_WEIGHT;
		A.insert(no_vertices * 3 + no_vertices * 9 + i * 3 + 1, fixed_indices_local[i] * 3 + 1) = CONSTRAINT_WEIGHT;
		A.insert(no_vertices * 3 + no_vertices * 9 + i * 3 + 2, fixed_indices_local[i] * 3 + 2) = CONSTRAINT_WEIGHT;

		B[no_vertices * 3 + no_vertices * 9 + i * 3] = CONSTRAINT_WEIGHT * V(fixed_indices[i], 0);
		B[no_vertices * 3 + no_vertices * 9 + i * 3 + 1] = CONSTRAINT_WEIGHT * V(fixed_indices[i], 1);
		B[no_vertices * 3 + no_vertices * 9 + i * 3 + 2] = CONSTRAINT_WEIGHT * V(fixed_indices[i], 2);
	}

	//Setup for r_i*R_i - R_i' = 0 (note that we want all off-diagonal elements in r_i to be equal to 0 in order to multiply R_i with the identity matrix)
	for(int i = 0; i < fixed_indices.size() - 1; i++) { //Skip the last fixed index (the handle index)
		A.insert(no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + i * 3, no_vertices * 3 + fixed_indices_local[i] * 3) = CONSTRAINT_WEIGHT;
		A.insert(no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + i * 3 + 1, no_vertices * 3 + fixed_indices_local[i] * 3 + 1) = CONSTRAINT_WEIGHT;
		A.insert(no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + i * 3 + 2, no_vertices * 3 + fixed_indices_local[i] * 3 + 2) = CONSTRAINT_WEIGHT;
		B[no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + i * 3] = 0;
		B[no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + i * 3 + 1] = 0;
		B[no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + i * 3 + 2] = 0;
	}
	Eigen::SparseMatrix<double> AT = A.transpose();
	solverPosRot.compute(AT*A);
	PosRot = solverPosRot.solve(AT*B);
}

void CurveDeformation::update_rot() {
	Eigen::Matrix3d newRot;
	double rx, ry, rz;
	for(int i = 0; i < Rot.size(); i++) {
		rx = PosRot[no_vertices * 3 + i * 3];
		ry = PosRot[no_vertices * 3 + i * 3 + 1];
		rz = PosRot[no_vertices * 3 + i * 3 + 2];

		newRot.row(0) << 1, -rz, ry;
		newRot.row(1) << rz, 1, -rx;
		newRot.row(2) << -ry, rx, 1;
		Rot[i] = newRot*Rot[i]; //Safe for matrix-matrix multiplication with Eigen
		Rot[i] = compute_orthonormal(Rot[i]);
	}
}

Eigen::Matrix3d CurveDeformation::compute_orthonormal(Eigen::Matrix3d& rot) {
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(rot, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d VT = svd.matrixV().transpose();
	return svd.matrixU()*VT;
}

void CurveDeformation::final_L1_pos(Eigen::MatrixXd &V, Eigen::VectorXi & part_of_original_stroke) {
	Eigen::VectorXd Bx(no_vertices + fixed_indices.size()), By(no_vertices + fixed_indices.size()), Bz(no_vertices + fixed_indices.size());
	Eigen::Vector3d Rl;
	for(int i = 0; i < no_vertices; i++) {
		Rl = Rot[i] * original_L1.row(i).transpose();
		Bx[i] = Rl(0);
		By[i] = Rl(1);
		Bz[i] = Rl(2);
	}

	for(int i = 0; i < fixed_indices.size(); i++) {
		Bx[no_vertices + i] = V(fixed_indices[i], 0) * CONSTRAINT_WEIGHT;
		By[no_vertices + i] = V(fixed_indices[i], 1) * CONSTRAINT_WEIGHT;
		Bz[no_vertices + i] = V(fixed_indices[i], 2) * CONSTRAINT_WEIGHT;
	}

	Eigen::VectorXd xpos = solverL1.solve(A_L1_T*Bx);
	Eigen::VectorXd ypos = solverL1.solve(A_L1_T*By);
	Eigen::VectorXd zpos = solverL1.solve(A_L1_T*Bz);
	if(stroke_ID == 0) { //We're pulling on the original stroke
		for(int i = 0; i < no_vertices; i++) { //update the position of non-fixed vertices of the stroke that is being pulled on
			if(!is_fixed[i]) {
				V.row(vert_bindings[i]) << xpos[i], ypos[i], zpos[i];
			}
		}
	} else { //We're pulling on an added stroke, and we want to avoid deforming the original curve, so don't update positions of the original stroke
		for(int i = 0; i < no_vertices; i++) { //update the position of non-fixed stroke vertices
			if(!is_fixed[i] && !part_of_original_stroke[vert_bindings[i]]) {
				V.row(vert_bindings[i]) << xpos[i], ypos[i], zpos[i];
			}
		}
	}

}
