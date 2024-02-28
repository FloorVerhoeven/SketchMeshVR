#ifndef _Laplacian_Remesh_H_
#define _Laplacian_Remesh_H_
#include "Mesh.h"
#include "SurfacePath.h"

class LaplacianRemesh {

public:
	static Eigen::VectorXi remesh_cut_remove_inside(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, bool& remesh_success);
	static Eigen::VectorXi remesh_extrusion_remove_inside(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, bool& remesh_success);
	static Eigen::VectorXi remesh(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, bool& remesh_success);

private:
	static bool is_front_loop;
	static bool remove_inside_faces;

	static Eigen::MatrixXi EV, FE, EF;
	static std::vector<std::vector<int>> VV;

	static int find_closest(std::vector<int> vertices, Eigen::Vector3d base, Mesh & m);
	static bool is_counter_clockwise_boundaries(Eigen::MatrixXd boundary_points, Eigen::Matrix4f modelview, Eigen::Matrix4f proj, Eigen::Vector4f viewport, Eigen::RowVector3d mean_viewpoint, bool cut);
	static std::vector<int> sort_boundary_vertices(Eigen::Vector3d start_vertex, std::vector<int> boundary_vertices, Mesh & m);
	static std::vector<int> reorder(std::vector<int> boundary_vertices, Eigen::Vector3d start_v, Mesh & m);
	static Eigen::RowVector3d compute_mean_viewpoint(Mesh & m, std::vector<int> inner_boundary_vertices);
	static void propagate_dirty_faces(int face, std::vector<bool>& dirty_face);
	static void stitch(std::vector<int> path_vertices, std::vector<int> boundary_vertices, Mesh & m);
	static void reverse_path(std::vector<int> path_vertices);
	static void update_sharp_edges(Mesh & m, Eigen::MatrixXi & all_sharpEV);
	static void update_mesh_values(Mesh & m, std::vector<PathElement> path, int stroke_ID, int new_mapped_start);
	static void update_face_indices(Mesh & m);


};

#endif