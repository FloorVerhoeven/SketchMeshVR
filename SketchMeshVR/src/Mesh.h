#ifndef _MESH_H_
#define _MESH_H_
#include <Eigen/Core>
#include <vector>

class Patch;

class Mesh {

public:
	Mesh();
    Mesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi & edge_boundary_markers, Eigen::VectorXi &vertex_is_fixed, Eigen::VectorXi &new_mapped_indices, Eigen::VectorXi & sharp_edge, int ID);
	Mesh& operator=(Mesh other);
	~Mesh();
	Eigen::MatrixXd &V;
	Eigen::MatrixXi &F;
	Eigen::VectorXi &edge_boundary_markers;
	Eigen::VectorXi &vertex_is_fixed;
	Eigen::VectorXi &new_mapped_indices;
	Eigen::VectorXi &sharp_edge;
	std::vector<Patch*> patches;

	std::vector<Patch*> face_patch_map;
	Eigen::VectorXi mesh_to_patch_indices;
	int ID;
};

#endif
