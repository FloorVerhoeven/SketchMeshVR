#ifndef _Curve_Deformation_H_
#define _Curve_Deformation_H_
#include "Stroke.h"

class CurveDeformation {
public:
	static void startPullCurve(Stroke & _stroke, int handle_ID);
	static void pullCurve(Eigen::RowVector3d& pos, Eigen::MatrixXd& V, Eigen::VectorXi & part_of_original_stroke);

	static bool smooth_deform_mode; //Used in the viewer menu to decide on smooth/sharp deformation

private:

	static double compute_curve_diag_length(Stroke & _stroke);
	static bool update_ROI(double drag_size);
	static void compute_ROI_boundaries(int & ROI_1, int & ROI_2);
	static void setup_fixed_indices(int ROI_1, int ROI_2);
	static void setup_for_update_curve(Eigen::MatrixXd& V);
	static void setup_for_L1_position_step(Eigen::MatrixXd& V);
	static void update_curve(Eigen::MatrixXd & V, Eigen::VectorXi & part_of_original_stroke);
	static void solve_for_pos_and_rot(Eigen::MatrixXd& V);
	static void update_rot();
	static Eigen::Matrix3d compute_orthonormal(Eigen::Matrix3d & rot);
	static void final_L1_pos(Eigen::MatrixXd & V, Eigen::VectorXi & part_of_original_stroke);
};
#endif

