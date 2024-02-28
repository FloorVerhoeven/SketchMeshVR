#ifndef _Mesh_Cut_H_
#define _Mesh_Cut_H_
#include "Stroke.h"
#include "Mesh.h"
#include "SurfacePath.h"

class MeshCut {

public:
	static bool cut(Eigen::MatrixXd & V, Eigen::MatrixXi & F, Eigen::VectorXi & vertex_boundary_markers, Eigen::VectorXi & part_of_original_stroke, Eigen::VectorXi & new_mapped_indices, Eigen::VectorXi & sharp_edge, Stroke & stroke);

private:
	static int prev_vertex_count;
	static int ID;

	static bool cut_main(Mesh& m, SurfacePath& surface_path, Stroke& stroke);
	static void mesh_open_hole(Eigen::VectorXi & boundary_vertices, Mesh & m);
	static void update_sharp_edges(Mesh & m, Eigen::MatrixXi & sharpEV);
	static void update_face_indices(Mesh & m, Eigen::MatrixXi & F2, Eigen::VectorXi & boundary_vertices, int original_v_size);
	static void project_points_to_2D(Eigen::VectorXi & boundary_vertices, Mesh & m, Eigen::MatrixXd & boundary_vertices_2D, Eigen::MatrixXi & stroke_edges, Eigen::RowVector3d & center, Eigen::Vector3d & x_vec, Eigen::Vector3d & y_vec);
	static void post_cut_update_points(Stroke& stroke, SurfacePath& surface_path);

};


#endif

