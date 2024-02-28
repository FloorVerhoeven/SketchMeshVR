#ifndef _Mesh_Extrusion_H_
#define _Mesh_Extrusion_H_
#include "Stroke.h"
#include "Mesh.h"
#include "SurfacePath.h"

class MeshExtrusion {

public:
	static bool extrude_prepare(Stroke & base, SurfacePath & surface_path);
	static bool extrude_main(Mesh & m, SurfacePath & surface_path, Stroke & stroke, Stroke & base, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, Eigen::MatrixXi & replacing_vertex_bindings);


private:

	static void generate_mesh(Mesh & m, Eigen::MatrixXd front_loop3D, Eigen::Vector3d center, Eigen::Vector3d x_vec, Eigen::Vector3d y_vec, Eigen::Vector3d offset, int nr_silhouette_vert, std::vector<int> loop_base_original_indices, std::vector<int> sil_original_indices);
	static void get_normal_and_center(Eigen::RowVector3d & center, Eigen::Vector3d & normal, Mesh & m, Eigen::VectorXi & boundary_vertices);
	static void find_left_and_right(int & most_left_vertex_idx, int & most_right_vertex_idx, Mesh & m, Eigen::VectorXi & boundary_vertices, const Eigen::RowVector3d & sil_start, const Eigen::RowVector3d & sil_end, Eigen::RowVector3d & center);
	static std::vector<int> add_silhouette_vertices(Mesh & m, int stroke_ID, Eigen::MatrixXd & silhouette_vertices, Eigen::MatrixXi & boundary_markers);
	static void create_loop(Mesh & m, Eigen::MatrixXd & loop3D, Eigen::VectorXi & boundary_vertices, std::vector<int> &loop_base_original_indices, int start_idx, int end_idx);
	static void update_edge_indicators(Mesh & m, Eigen::MatrixXi & edges_to_update);
	static void update_face_indices(Mesh & m, Eigen::MatrixXi & F2, std::vector<int> sil_original_indices, std::vector<int> loop_base_original_indices, int nr_silhouette_vert, int size_before_gen, int loop2D_size);
	static void post_extrude_prepare_update_points(Stroke& stroke, SurfacePath& surface_path);
	static void post_extrude_main_update_points(Stroke & stroke, Eigen::MatrixXd new_positions);
	static void post_extrude_main_update_bindings(Stroke & base, SurfacePath & surface_path);
	static Eigen::MatrixXd smooth_stroke(Eigen::MatrixXd & input_points);
	static void move_to_middle(Eigen::MatrixXd & positions, Eigen::MatrixXd & new_positions);
	static Eigen::MatrixXd TaubinFairing3D(Eigen::MatrixXd & original_3DPoints, int n);
	static void smooth_sub3D(Eigen::MatrixXd & points, double direction);
	static Eigen::RowVector3d to_sum_of_vectors3D(Eigen::RowVector3d vert, Eigen::RowVector3d prev, Eigen::RowVector3d next, double direction);
	static double get_stroke_unit_length(Eigen::MatrixXd & vertices);
	static double get_base_unit_length(Eigen::MatrixXd & V, Eigen::VectorXi & boundary_vertices);
	static void remove_out_of_bounds_silhouette(Eigen::MatrixXd & silhouette_vertices, Eigen::RowVector3d & center, const Eigen::RowVector3d & left_most, const Eigen::RowVector3d & right_most, Eigen::RowVector3d & dir);

};

#endif