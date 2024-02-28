#include "Mesh.h"

Mesh::Mesh(Eigen::MatrixXd & V_, Eigen::MatrixXi & F_, Eigen::VectorXi & vertex_boundary_markers_, Eigen::VectorXi &part_of_original_stroke_, Eigen::VectorXi &new_mapped_indices_, Eigen::VectorXi &sharp_edge_, int ID_): V(V_), F(F_), vertex_boundary_markers(vertex_boundary_markers_), part_of_original_stroke(part_of_original_stroke_), new_mapped_indices(new_mapped_indices_), sharp_edge(sharp_edge_), ID(ID_)  {

}

Mesh::~Mesh() {
}