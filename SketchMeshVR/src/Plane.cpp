#include "Plane.h"
#include <Eigen/Geometry>

Plane::Plane(Eigen::Vector3d v0, Eigen::Vector3d v1, Eigen::Vector3d v2) {
	Eigen::Vector3d vec0 = v1 - v0;
	Eigen::Vector3d vec1 = v2 - v0;
	normal = vec0.cross(vec1);
	base = v0;
	normal.normalize();
}

Plane::Plane(Eigen::Vector3d point, Eigen::Vector3d normal_) :
	base(point),
	normal(normal_) {
}

/** Computes and returns the 3D point where the line segment from edge_start to edge_end crosses the plane. **/
Eigen::RowVector3d Plane::cross_point(Eigen::RowVector3d edge_start, Eigen::RowVector3d edge_end) {
	double start_to_surface = signed_distance(edge_start);
	double end_to_surface = signed_distance(edge_end);
	double start_to_end = start_to_surface - end_to_surface;
	Eigen::RowVector3d result = edge_start + (edge_end - edge_start)*(start_to_surface / start_to_end);
	return result;
}


Eigen::RowVector3d Plane::cross_point(Eigen::RowVector3d edge_start, Eigen::RowVector3d edge_end, double& t_val) {
	double start_to_surface = signed_distance(edge_start);
	double end_to_surface = signed_distance(edge_end);
	double start_to_end = start_to_surface - end_to_surface;
	Eigen::RowVector3d result = edge_start + (edge_end - edge_start)*(start_to_surface / start_to_end);
	t_val = start_to_surface / start_to_end;
	return result;
}

/** Computes and returns the 3D point where the ray coming from ray_source in direction ray_dir crosses the plane. **/
Eigen::RowVector3d Plane::intersect_point(Eigen::Vector3d ray_source, Eigen::Vector3d ray_dir) {
	double t = ((base - ray_source).dot(normal)) / (ray_dir.dot(normal));
	return (ray_source + ray_dir * t).transpose();
}

double Plane::signed_distance(Eigen::RowVector3d v) {
	Eigen::Vector3d vec = v.transpose() - base;
	return normal.dot(vec);
}
