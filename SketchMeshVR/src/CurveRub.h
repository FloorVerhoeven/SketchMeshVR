#ifndef _Curve_Rub_H_
#define _Curve_Rub_H_

#include <Eigen/Core>
#include <vector>

class CurveRub {

public:
	static void set_prev_target(int target);
	static void start_rubbing(int handleID, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& edge_boundary_markers);
	static void rubbing(Eigen::Vector3d& pos, int seam, Eigen::MatrixXd& V);
private:
	static std::vector<int> select_vertices_along_path(Eigen::Vector3d& pos, int seam, Eigen::MatrixXd& V);
	static void smooth(int v, int edge, Eigen::MatrixXd& V);
	static Eigen::RowVector3d get_target_position(int v0, int v1, int v2, int v3, int v4, Eigen::MatrixXd& V);
	static int get_next_vertex(int v0, int v1);
	static int find_next_v(int vert, Eigen::Vector3d& pos, int seam, Eigen::MatrixXd& V);
	static int find_edge(int start, int end);
	static bool on_seam_intersection(int v, int edge0);

};

#endif
