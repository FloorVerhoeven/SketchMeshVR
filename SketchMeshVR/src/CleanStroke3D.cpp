#include <CleanStroke3D.h>
#include <vector>

using namespace std;
using namespace igl;

Eigen::MatrixXd CleanStroke3D::resample_by_length_with_fixes(std::vector<PathElement> path_vertices, double unit_length) {
	if (path_vertices.size() <= 1) {
		Eigen::MatrixXd result(path_vertices.size(), 3);
		for (int i = 0; i < path_vertices.size(); i++) {
			result.row(i) = path_vertices[i].get_vertex();
		}
		return result;
	}

	if (path_vertices[0].get_vertex() == path_vertices[path_vertices.size() - 1].get_vertex()) {
		for (int i = 0; i < path_vertices.size(); i++) {
			if (path_vertices[i].fixed) {
				path_vertices.pop_back();
				std::vector<PathElement> reordered;
				for (int j = 0; j < path_vertices.size(); j++) {
					reordered.push_back(path_vertices[(i + j) % path_vertices.size()]);
				}
				reordered.push_back(reordered[0]);
				path_vertices = reordered;
				break;
			}
		}
	}

	Eigen::MatrixXd resampled(1, 3), resample_sub(0, 3);
	resampled.row(0) = path_vertices[0].get_vertex();

	int idx0 = 0;
	while (true) {
		int idx1 = find_next_fixed_vertex(path_vertices, idx0);
		resample_sub = resample_by_length_sub(path_vertices, idx0, idx1, unit_length);
		resampled.conservativeResize(resampled.rows() + resample_sub.rows(), Eigen::NoChange);
		resampled.bottomRows(resample_sub.rows()) = resample_sub;
		if (idx1 == path_vertices.size() - 1) {
			break;
		}
		idx0 = idx1;
	}
	return resampled;
}

int CleanStroke3D::find_next_fixed_vertex(std::vector<PathElement> path_vertices, int idx) {
	for (int i = idx + 1; i < path_vertices.size(); i++) {
		if (path_vertices[i].fixed) {
			return i;
		}
	}
	return path_vertices.size() - 1;
}

Eigen::MatrixXd CleanStroke3D::resample_by_length_sub(std::vector<PathElement> path_vertices, int start_index, int end_index, double unit_length) {
	double length = get_stroke_length(path_vertices, start_index, end_index);
	int n = 1 + (int)(length / unit_length);

	Eigen::RowVector3d v0 = path_vertices[start_index].get_vertex();
	Eigen::RowVector3d v1 = path_vertices[end_index].get_vertex();

	Eigen::MatrixXd resampled_points(0, 3);

	if (length < unit_length) { //Actual stroke is shorter than requested inter-sample length
		resampled_points.conservativeResize(resampled_points.rows() + 1, Eigen::NoChange);
		resampled_points.row(resampled_points.rows() - 1) = v1;
		return resampled_points;
	}

	double total = 0.0, prev_total = 0.0, next_spot = unit_length;
	Eigen::RowVector3d prev = v0, next = path_vertices[start_index + 1].get_vertex();
	int index = start_index + 1, count = 0;

	while (true) {
		next = path_vertices[index].get_vertex();
		total += (prev - next).norm();
		while (total >= next_spot) { //The along-path distance_to_vert to the next path_vertex is bigger than where we would want the next sample, so we create an interpolated sample
			double t = (next_spot - prev_total) / (total - prev_total);
			resampled_points.conservativeResize(resampled_points.rows() + 1, Eigen::NoChange);
			resampled_points.row(resampled_points.rows() - 1) = prev*(1 - t) + next*t;
			next_spot += unit_length;
			count++;
			if (count == n - 1) {
				break;
			}

		}
		if (count == n - 1) {
			break;
		}
		prev = next;
		prev_total = total;
		index++;

		if (index == end_index) {
			break;
		}
	}

	resampled_points.conservativeResize(resampled_points.rows() + 1, Eigen::NoChange);
	resampled_points.row(resampled_points.rows() - 1) = v1;
	return resampled_points;
}

Eigen::MatrixXd CleanStroke3D::resample_by_length_sub(Eigen::MatrixXd path_vertices, int start_index, int end_index, double unit_length) {
	double length = get_stroke_length(path_vertices, start_index, end_index);
	int n = 1 + (int)(length / unit_length);

	Eigen::RowVector3d v0 = path_vertices.row(start_index);
	Eigen::RowVector3d v1 = path_vertices.row(end_index);

	Eigen::MatrixXd resampled_points(0, 3);

	if (length < unit_length) { //Actual stroke is shorter than requested inter-sample length
		resampled_points.conservativeResize(resampled_points.rows() + 1, Eigen::NoChange);
		resampled_points.row(resampled_points.rows() - 1) = v1;
		return resampled_points;
	}

	double total = 0.0, prev_total = 0.0, next_spot = unit_length;
	Eigen::RowVector3d prev = v0, next = path_vertices.row(start_index + 1);
	int index = start_index + 1, count = 0;

	while (true) {
		next = path_vertices.row(index);
		total += (prev - next).norm();
		while (total >= next_spot) { //The along-path distance_to_vert to the next path_vertex is bigger than where we would want the next sample, so we create an interpolated sample
			double t = (next_spot - prev_total) / (total - prev_total);
			resampled_points.conservativeResize(resampled_points.rows() + 1, Eigen::NoChange);
			resampled_points.row(resampled_points.rows() - 1) = prev*(1 - t) + next*t;
			next_spot += unit_length;
			count++;
			if (count == n - 1) {
				break;
			}

		}
		if (count == n - 1) {
			break;
		}
		prev = next;
		prev_total = total;
		index++;

		if (index == end_index + 1) {
			break;
		}
	}

	resampled_points.conservativeResize(resampled_points.rows() + 1, Eigen::NoChange);
	resampled_points.row(resampled_points.rows() - 1) = v1;
	return resampled_points;
}

double CleanStroke3D::get_stroke_length(std::vector<PathElement> path_vertices, int start_index, int end_index) {
	double length = 0.0;
	for (int i = start_index; i < end_index; i++) {
		length += vertex_distance(path_vertices[i], path_vertices[i + 1]);
	}
	return length;
}

double CleanStroke3D::get_stroke_length(Eigen::MatrixXd stroke) {
	double length = 0.0;
	for (int i = 0; i < stroke.rows() - 1; i++){
		length += (stroke.row(i) - stroke.row(i+1)).norm();
	}
	return length;
}

double CleanStroke3D::get_stroke_length(Eigen::MatrixXd stroke, int start, int end) {
	double length = 0.0;
	for (int i = start; i < end; i++) {
		length += (stroke.row(i) - stroke.row(i + 1)).norm();
	}
	return length;
}

double CleanStroke3D::vertex_distance(PathElement prev, PathElement next) {
	return (prev.get_vertex() - next.get_vertex()).norm();
}

Eigen::MatrixXd CleanStroke3D::resample(Eigen::MatrixXd stroke) {
	int n = stroke.rows() - 1;
	return resample_by_number(stroke, n);
}

Eigen::MatrixXd CleanStroke3D::resample_by_number(Eigen::MatrixXd stroke, int n) {
	if (stroke.rows() <= 1) {
		return stroke;
	}

	double length = get_stroke_length(stroke);
	double unit = length / n;

	Eigen::RowVector3d v0 = stroke.row(0);
	Eigen::RowVector3d v1 = stroke.row(stroke.rows() - 1);

	Eigen::MatrixXd resampled_points(0, 3);
	resampled_points.conservativeResize(resampled_points.rows() + 1, Eigen::NoChange);
	resampled_points.row(resampled_points.rows() - 1) = v0;

	double total = 0.0, prev_total = 0.0, next_spot = unit, t;
	Eigen::RowVector3d prev = v0, next;
	int index = 1, count = 0;

	while (true) {
		if (index == stroke.rows()) {
			break;
		}
		next = stroke.row(index);

		total += (prev - next).norm();
		while (total >= next_spot) {
			t = (next_spot - prev_total) / (total - prev_total);
			resampled_points.conservativeResize(resampled_points.rows() + 1, Eigen::NoChange);
			resampled_points.row(resampled_points.rows() - 1) = prev*(1 - t) + next*t;
			next_spot += unit;
			count++;
			if (count == n - 1) {
				break;
			}
		}
		if (count == n - 1) {
			break;
		}
		prev = next;
		prev_total = total;
		index++;

	}

	resampled_points.conservativeResize(resampled_points.rows() + 1, Eigen::NoChange);
	resampled_points.row(resampled_points.rows() - 1) = v1;

	return resampled_points;
}