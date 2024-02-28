#ifndef _VR_Plane_H_
#define _VR_Plane_H_
#include <Eigen/Core>

class Plane {
public:
    Eigen::Vector3d base;
    Eigen::Vector3d normal;

	Plane(Eigen::Vector3d v0, Eigen::Vector3d v1, Eigen::Vector3d v2);
	Plane(Eigen::Vector3d point, Eigen::Vector3d normal);
	Eigen::RowVector3d cross_point(Eigen::RowVector3d edge_start, Eigen::RowVector3d edge_end);
	Eigen::RowVector3d cross_point(Eigen::RowVector3d edge_start, Eigen::RowVector3d edge_end, double & t_val);
	Eigen::RowVector3d intersect_point(Eigen::Vector3d ray_source, Eigen::Vector3d ray_dir);

private:
	double signed_distance(Eigen::RowVector3d v);
};
#endif
