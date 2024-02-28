#ifndef _Curve_Deformation_H_
#define _Curve_Deformation_H_
#include "LaplacianCurveEdit.h"

class CurveDeformation {

public:
	static void startPullCurve(int _handle_ID, Eigen::MatrixXd & V, Eigen::MatrixXi & F, Eigen::VectorXi& edge_boundary_markers);
	static void pullCurve(const Eigen::RowVector3d & pos, Eigen::MatrixXd & V, Eigen::VectorXi& edge_boundary_markers);

	static Eigen::MatrixXi EF, EV, FE;
	static std::vector<std::vector<int>> neighbors;

	struct PulledCurve {
		std::vector<int> vertices, fixed_vertices, fixed_edges;
		Eigen::MatrixXi edges, edge_triplets, vertex_triplets;

		LaplacianCurveEdit laplacian_curve_edit;

		PulledCurve() {
			edges.resize(0, 2);
			edge_triplets.resize(0, 3);
			vertex_triplets.resize(0, 3);
		}
	};
private:

	static bool update_ROI_test(double drag_size, Eigen::MatrixXd & V, Eigen::VectorXi& edge_boundary_markers);
	static std::vector<int> collect_vertices_within_drag_length(double drag_size, Eigen::MatrixXd & V, Eigen::VectorXi& edge_boundary_markers);
	static void propagate(int prev_edge, int vert, double remaining_distance, std::vector<int>& vertices_in_range, Eigen::MatrixXd & V, Eigen::VectorXi& edge_boundary_markers);
	static void propagate(PulledCurve & curve, int prev_edge, int vert, int edge_marker, Eigen::VectorXi & edge_consumed, Eigen::VectorXi & edge_boundary_markers);
	static int find_next_edge(int vert, int prev_edge, int edge_marker, Eigen::VectorXi & edge_boundary_markers);
	static double find_furthest_edge_vertex(Eigen::MatrixXd & V, Eigen::VectorXi & edge_boundary_markers);
	static void collect_vertex_distances(double drag_size, Eigen::MatrixXd & V, Eigen::VectorXi & edge_boundary_markers);
	static void propagate_distances(int prev_edge, int vert, double remaining_distance, std::vector<int>& vertices_in_range, Eigen::MatrixXd & V, Eigen::VectorXi & edge_boundary_markers);
	static int find_edge(int start, int end);
	static void sort_by_distance(std::vector<int>& vertices_in_range);
	static bool distance_sorter(int a, int b);
	static bool create_pulled_curve_by_propagation(int vert, Eigen::VectorXi & edge_consumed, PulledCurve & curve, Eigen::VectorXi& edge_boundary_markers);
	static std::vector<int> get_adjacent_seam_vertices(int vert, PulledCurve & curve);
	static int get_adjacent_seam_edge(int vert, int edge, Eigen::VectorXi & edge_boundary_markers);
	
};
#endif

