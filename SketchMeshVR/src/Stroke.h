#ifndef _VR_Stroke_H_
#define _VR_Stroke_H_
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>

#include <vector>
#include <unordered_map>
class Stroke {
public:
	
	Stroke();
	Stroke(Eigen::MatrixXd *V,  Eigen::MatrixXi *F, int stroke_ID_);
	Stroke(const Stroke& origin);
    Stroke& operator=(Stroke);
	void swap(Stroke & tmp);
	~Stroke();

	bool addSegment(Eigen::Vector3f & pos, igl::opengl::glfw::Viewer & viewer);
	void addSegmentAdd(Eigen::Vector3f & pos, igl::opengl::glfw::Viewer & viewer);
	void addSegmentCut(Eigen::Vector3f & pos, igl::opengl::glfw::Viewer & viewer);
	bool addSegmentExtrusionBase(Eigen::Vector3f & pos, igl::opengl::glfw::Viewer & viewer);
	void addSegmentExtrusionSilhouette(Eigen::Vector3f & pos, igl::opengl::glfw::Viewer & viewer);

	void prepend_first_point(igl::opengl::glfw::Viewer & viewer);
	void append_final_point(igl::opengl::glfw::Viewer & viewer);
    void counter_clockwise();
	bool toLoop();

	bool update_vert_bindings(Eigen::VectorXi & new_mapped_indices, Eigen::MatrixXi & replacing_vertex_bindings);
	void undo_stroke_add(Eigen::VectorXi & edge_boundary_markers, Eigen::VectorXi & sharp_edge, Eigen::VectorXi & vertex_is_fixed);
	void switch_stroke_edges_type(Eigen::VectorXi & sharp_edge);
	void update_Positions(Eigen::MatrixXd V, bool structure_changed);
	void generate3DMeshFromStroke(Eigen::VectorXi & edge_boundary_markers, Eigen::VectorXi & vertex_is_fixed, Eigen::MatrixXd & mesh_V, Eigen::MatrixXi & mesh_F, igl::opengl::glfw::Viewer & viewer);
	void project_with_PCA_given_target(Eigen::Vector3d target_vec);
	void resample_and_smooth_3DPoints(Eigen::Matrix4f& model, Eigen::Matrix4f& view, Eigen::Matrix4f& proj, Eigen::Vector4f& viewport);
	bool has_self_intersection(bool make_looped);
	int selectClosestVertex(Eigen::Vector3f pos, double & closest_distance);
    
    int get_vertex_idx_for_point(int i);
    int get_ID();
    Eigen::MatrixXd get_V() const;
    Eigen::MatrixXi get_F() const;
    Eigen::MatrixXd get_stroke2DPoints() const;
    Eigen::MatrixX3d get3DPoints();
	Eigen::MatrixX3d get3DPointsBack();
	Eigen::MatrixX3d get_hand_pos();
	Eigen::MatrixXi get_hit_faces();
	Eigen::MatrixXi get_stroke_edges();
	::std::vector<int> get_closest_vert_bindings();
    void set3DPoints(Eigen::MatrixX3d new_3DPoints);
    void set_closest_vert_bindings(::std::vector<int> new_vert_bindings);
	void setV(Eigen::MatrixXd* V_) { V = V_; };
	void setF(Eigen::MatrixXi* F_) { F = F_; };

	
	Eigen::MatrixXd* V;
	Eigen::MatrixXi* F;
	int stroke_ID; //Non-const for the sake of copy assignment operator 
	bool is_loop;
	bool has_points_on_mesh;
	bool has_been_outside_mesh;
	bool has_been_reversed;
	bool starts_on_mesh;
	bool ends_on_mesh;
	Eigen::RowVector3d stroke_color;

	Eigen::MatrixX3d stroke3DPoints;
	Eigen::MatrixX3d stroke3DPointsBack;
	Eigen::MatrixXd stroke2DPoints; //Used for early checking if point is new (in screen coordinates)
	Eigen::MatrixXi stroke_edges; //Indices into closest_vert_bindings, where every row contains the 2 end points of a stroke edge (looped)
	Eigen::MatrixXi faces_hit;
	Eigen::MatrixX3d hand_pos_at_draw; //Only used for extrusion base strokes.
	Eigen::Vector3d pos_before_cut;
	Eigen::Vector3d dir_before_cut;
	Eigen::Vector3d pos_after_cut;
	Eigen::Vector3d dir_after_cut;
	Eigen::VectorXd dep;
	bool prev_point_was_on_mesh; //Used for cutting strokes only. Indicates whether this is the first point outside of the mesh after we've been drawing on the mesh


	static const int MAX_NR_TRIANGLES = 10000; //For the entire mesh (front + backside)
	static constexpr double min_inter_point_distance = 0.0001125;


	std::vector<int> closest_vert_bindings;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:

	static Eigen::MatrixXd resample_stroke2D(Eigen::MatrixXd & original_stroke2DPoints);
	static void move_to_middle(Eigen::MatrixXd &positions, Eigen::MatrixXd &new_positions);
	static void TaubinFairing2D(Eigen::MatrixXd & original_stroke2DPoints, int n);
	static void smooth_sub2D(Eigen::MatrixXd & points, double direction);
	static Eigen::RowVector2d to_sum_of_vectors2D(Eigen::RowVector2d vert, Eigen::RowVector2d prev, Eigen::RowVector2d next, double direction);
	static void generate_backfaces(Eigen::MatrixXi &faces, Eigen::MatrixXi &back_faces);
	bool line_segments_intersect(Eigen::RowVector2d& p1, Eigen::RowVector2d& p2, Eigen::RowVector2d& p3, Eigen::RowVector2d& p4);
	double get_2D_length();
	Eigen::MatrixXd resample_stroke2D(Eigen::MatrixXd & original_2D, double unit_length, double length);
	Eigen::MatrixXd move_to_middle_smoothing(Eigen::MatrixXd & stroke2DPoints);
	double compute_stroke_diag();
	bool empty2D() const { return stroke2DPoints.isZero(); }

};

#include <igl/serialize.h>
namespace igl {
	namespace serialization {

		inline void serialization(bool s, Stroke& obj, std::vector<char>& buffer)
		{
			SERIALIZE_MEMBER(stroke_ID);
			SERIALIZE_MEMBER(stroke3DPoints);
			SERIALIZE_MEMBER(stroke3DPointsBack);
			SERIALIZE_MEMBER(stroke2DPoints);
			SERIALIZE_MEMBER(stroke_edges);
			SERIALIZE_MEMBER(faces_hit);
			SERIALIZE_MEMBER(hand_pos_at_draw);
			SERIALIZE_MEMBER(pos_before_cut);
			SERIALIZE_MEMBER(dir_before_cut);
			SERIALIZE_MEMBER(pos_after_cut);
			SERIALIZE_MEMBER(dir_after_cut);
			SERIALIZE_MEMBER(prev_point_was_on_mesh);
			SERIALIZE_MEMBER(dep);
			SERIALIZE_MEMBER(closest_vert_bindings);
			SERIALIZE_MEMBER(is_loop);
			SERIALIZE_MEMBER(has_points_on_mesh);
			SERIALIZE_MEMBER(has_been_outside_mesh);
			SERIALIZE_MEMBER(has_been_reversed);
			SERIALIZE_MEMBER(starts_on_mesh);
			SERIALIZE_MEMBER(ends_on_mesh);
			SERIALIZE_MEMBER(stroke_color);
		}

		template<>
		inline void serialize(const Stroke& obj, std::vector<char>& buffer)
		{
			serialization(true, const_cast<Stroke&>(obj), buffer);
		}

		template<>
		inline void deserialize(Stroke& obj, const std::vector<char>& buffer)
		{
			serialization(false, obj, const_cast<std::vector<char>&>(buffer));
		}
	}
}

#endif
