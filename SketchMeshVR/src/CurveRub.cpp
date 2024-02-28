#include <igl/edge_topology.h>
#include <igl/adjacency_list.h>
#include "CurveRub.h"

 int prev_target;
 Eigen::MatrixXi rubEV, rubEF, rubFE;
 std::vector<std::vector<int>> rub_neighbors;
 Eigen::VectorXi rub_edge_boundary_markers;

void CurveRub::set_prev_target(int target) {
	prev_target = target;
}

void CurveRub::start_rubbing(int handleID, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& edge_boundary_markers) {
	set_prev_target(handleID);
	igl::edge_topology(V, F, rubEV, rubFE, rubEF);
	igl::adjacency_list(F, rub_neighbors);
	rub_edge_boundary_markers = edge_boundary_markers;
}

void CurveRub::rubbing(Eigen::Vector3d& pos, int seam, Eigen::MatrixXd& V) {
	std::vector<int> rubbed_vertices = select_vertices_along_path(pos, seam, V);
	if (rubbed_vertices.size() == 1) {
		return;
	}

	int v, v_next, edge;
	for (int i = 0; i < rubbed_vertices.size() - 1; i++) {
		v = rubbed_vertices[i];
		v_next = rubbed_vertices[i + 1];
		edge = find_edge(v, v_next);
		if (on_seam_intersection(v, edge)) {
			continue;
		}
		smooth(v, edge, V);
	}
	set_prev_target(rubbed_vertices.back());
}

std::vector<int> CurveRub::select_vertices_along_path(Eigen::Vector3d& pos, int seam, Eigen::MatrixXd& V) {
	std::vector<int> rubbed_vertices;
	int vert = prev_target;
	while (true) {
		rubbed_vertices.push_back(vert);
		vert = find_next_v(vert, pos, seam, V);
		if (vert == -1) {
			break;
		}
	}
	return rubbed_vertices;
}

void CurveRub::smooth(int v, int edge, Eigen::MatrixXd& V) {
	int v2 = v;
	int v3 = rubEV(edge, 0) == v ? rubEV(edge, 1) : rubEV(edge, 0);
	int v4 = get_next_vertex(v2, v3);
	int v1 = get_next_vertex(v3, v2);
	int v0 = get_next_vertex(v2, v1);

	if (v0 == -1 || v1 == -1 || v2 == -1 || v3 == -1 || v4 == -1) {
		return;
	}

	V.row(v2) = get_target_position(v0, v1, v2, v3, v4, V);
}

Eigen::RowVector3d CurveRub::get_target_position(int v0, int v1, int v2, int v3, int v4, Eigen::MatrixXd& V) {
	Eigen::RowVector3d vec0 = V.row(v1) - (0.5*V.row(v0) + 0.5*V.row(v2));
	Eigen::RowVector3d vec2 = V.row(v3) - (0.5*V.row(v2) + 0.5*V.row(v4));
	Eigen::RowVector3d vec = (vec0 + vec2)*0.5;
	vec.normalize();
	vec *= (vec0.norm() + vec2.norm())*0.5;
	Eigen::RowVector3d target = (0.5*V.row(v1) + 0.5*V.row(v3)) + vec;
	return (0.5*V.row(v2) + 0.5*target);
}

int CurveRub::get_next_vertex(int v0, int v1) {
	if (v0 == -1 || v1 == -1)
		return -1;

	int edge0 = find_edge(v0, v1);
	if (edge0 == -1)
		return -1;


	int edge1 = -1;
	int edge;
	for (int i = 0; i< rub_neighbors[v1].size(); i++) {
		edge = find_edge(v1, rub_neighbors[v1][i]);
		if (edge != edge0 && rub_edge_boundary_markers[edge] == rub_edge_boundary_markers[edge0]) {
			edge1 = edge;
			break;
		}
	}
	if (edge1 == -1)
		return -1;

	int v2 = rubEV(edge1, 0) == v1 ? rubEV(edge1, 1) : rubEV(edge1, 0);
	return v2;
}

int CurveRub::find_next_v(int vert, Eigen::Vector3d& pos, int seam, Eigen::MatrixXd& V) {
	int next_v = -1;
	double min_dist = (pos - V.row(vert).transpose()).norm();
	for (int i = 0; i < rub_neighbors[vert].size(); i++) {
		int end = rub_neighbors[vert][i];
		int edge = find_edge(vert, end);
		if (!rub_edge_boundary_markers[edge]) {
			continue;
		}
		if (rub_edge_boundary_markers[edge] != seam) {
			continue;
		}
		double dist = (pos - V.row(end).transpose()).norm();
		if (dist < min_dist) {
			next_v = end;
		}
	}
	return next_v;
}

int CurveRub::find_edge(int start, int end) {
	Eigen::VectorXi col1Equals, col2Equals;
	int equal_pos;
	col1Equals = rubEV.col(0).cwiseEqual(std::min(start, end)).cast<int>();
	col2Equals = rubEV.col(1).cwiseEqual(std::max(start, end)).cast<int>();
	(col1Equals + col2Equals).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge

	return equal_pos;
}

bool CurveRub::on_seam_intersection(int	v, int edge0) {
	int edge;
	for (int i = 0; i < rub_neighbors[v].size(); i++) {
		edge = find_edge(v, rub_neighbors[v][i]);
		if (edge != edge0 && rub_edge_boundary_markers[edge] && rub_edge_boundary_markers[edge] != rub_edge_boundary_markers[edge0]) {
			return true;
		}
	}
	return false;
}