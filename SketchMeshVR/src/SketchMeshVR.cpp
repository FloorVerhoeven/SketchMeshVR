#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#else
#include <unistd.h>
#include <curses.h>
#endif

#define STB_IMAGE_IMPLEMENTATION
#include "igl/opengl/glfw/imgui/stb_image.h"

#include <igl/readOFF.h>
#include <igl/edge_topology.h>
#include <igl/cat.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/internal_angles.h>
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_corner_normals.h>

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/material_colors.h>
#include <igl/slice.h>
#include <igl/avg_edge_length.h>
#include <SketchMeshVR.h>
#include <Stroke.h>
#include "SurfaceSmoothing.h"
#include "MeshCut.h"
#include "CurveDeformation.h"
#include "MeshExtrusion.h"
#include <igl/opengl/oculusVR.h>
#include "Mesh.h"
#include "Patch.h"
#include "LaplacianRemesh.h"
#include "AppState.h"
#include "CurveRub.h"
#include <chrono>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;
using OculusVR = igl::opengl::OculusVR;

AppState appState;

OculusVR tmp;
Viewer viewer(tmp);
igl::opengl::glfw::imgui::ImGuiMenu menu;

Eigen::MatrixXd V_floor(4, 3), V_roof(4,3), V_left_wall(4,3), V_front_wall(4,3), V_right_wall(4,3), V_back_wall(4,3);
Eigen::MatrixXi F_floor(2, 3), F_roof(2,3), F_left_wall(2,3), F_front_wall(2,3), F_right_wall(2,3), F_back_wall(2,3);

// Vertex array, #V x3
Eigen::MatrixXd V(0, 3);
// Face array, #F x3
Eigen::MatrixXi F(0, 3);;
//Per vertex indicator of whether vertex position is fixed
Eigen::VectorXi vertex_is_fixed;
//Per edge indicator of what curve (or none if == 0) each edge belongs to
Eigen::VectorXi edge_boundary_markers;
//Per edge indicator of whether the edge is sharp (if == 1 then sharp, otherwise smooth)
Eigen::VectorXi sharp_edge;
//Takes care of index mapping from before a cut/extrusion action to after (since some vertices are removed)
Eigen::VectorXi new_mapped_indices;
//Tracker needed for updating the strokes' closest vertex bindings after a stroke edge has been split (e.g. due to stroke add, cut or extrude)
Eigen::MatrixXi replacing_vertex_bindings(0, 4);

Mesh* base_mesh;

int base_mesh_index, pointer_mesh_index;

std::function<void()> current_tool_HUD;
std::function<void(void)> menu_HUD;
//General
enum ToolMode { DRAW, ADD, CUT, EXTRUDE, PULL, REMOVE, CHANGE, SMOOTH, UPSAMPLE, NAVIGATE, TRANSLATE, NONE, DEFAULT, FAIL }; //NONE is used to indicate that a button was released, whereas DEFAULT indicates that one of the toggle buttons was pressed within its cooldown period, FAIL is used to indicate that something went wrong (e.g. user clicked too far away for PULL)

ToolMode tool_mode = DRAW, selected_tool_mode = DRAW, prev_tool_mode = NONE;
Stroke* added_stroke;
Stroke* extrusion_base;
vector<Stroke> stroke_collection;

//For smoothing
int initial_smooth_iter = 10;

//For selecting vertices
int handleID = -1;
int rub_seamID = -1;

//Variables for pulling a curve (and removing added control curves)
bool dirty_boundary = false;
int closest_stroke_ID, prev_closest_stroke_ID, prev_closest_stroke_idx = -1;

//Keeps track of the stroke IDs
int next_added_stroke_ID = 2; //Start at 2 because marker 1 belongs to the original boundary

//Variables for removing a control curve
bool stroke_was_removed = false;
int remove_stroke_clicked = 0;
Eigen::RowVector3d old_selected_stroke_color;

//Variables for cutting
bool cut_stroke_already_drawn = false;

//Variables for extrusion
bool extrusion_base_already_drawn = false;
SurfacePath base_surface_path, cut_surface_path;
Eigen::Matrix4f base_model, base_view, base_proj;
Eigen::Vector4f base_viewport;

bool has_recentered = false;
bool draw_should_block = false;
bool extrude_base_should_block = false;
bool upsample_activated = false;
bool upsample_init = false;


Eigen::RowVector3d red(1, 0, 0);
Eigen::RowVector3d black(0, 0, 0);
Eigen::RowVector3d blue(0, 0, 1);

Eigen::RowVector3d laser_end_point;
bool prev_laser_show;

Eigen::Vector3d prev_translation_point;
double prev_scale_size;

Eigen::VectorXd pull_path_starting_angles;

void sound_error_beep() {
#ifdef _WIN32
	Beep(500, 200);
#else
	beep();
#endif		
}

void draw_all_strokes() {
	viewer.data().set_edges(Eigen::MatrixXd(), Eigen::MatrixXi(), red); //Clear the non-original stroke edges
	Eigen::MatrixXd added_points;
	int points_to_hold_back;
	for (int i = 0; i < stroke_collection.size(); i++) {
		Eigen::VectorXi col_slice(3);
		col_slice.col(0) << 0, 1, 2;
		Eigen::MatrixXd edges_start, edges_end;
		Eigen::MatrixXi stroke_edges = stroke_collection[i].get_stroke_edges();
		igl::slice(V, stroke_edges.col(0).topRows(stroke_edges.rows() - !stroke_collection[i].is_loop), col_slice, edges_start);
		igl::slice(V, stroke_edges.col(1).topRows(stroke_edges.rows() - !stroke_collection[i].is_loop), col_slice, edges_end);

		viewer.data().add_edges(edges_start, edges_end, stroke_collection[i].stroke_color);
	}
}

void draw_extrusion_base() {
	//The extrusion base does not get resampled, so use simple old method
	Eigen::MatrixXd added_points = extrusion_base->get3DPoints();
	int points_to_hold_back = 1 + !extrusion_base->is_loop;
	viewer.data().add_edges(added_points.topRows(added_points.rows() - points_to_hold_back), added_points.middleRows(1, added_points.rows() - points_to_hold_back), red);
}

void rotate_unpushed_extrusion() {
	Eigen::MatrixXd new_points = extrusion_base->get3DPoints();
	viewer.data().rotate_points(new_points);
	extrusion_base->set3DPoints(new_points);
	draw_extrusion_base();

	new_points = base_surface_path.get3DPoints();
	viewer.data().rotate_points(new_points);
	base_surface_path.set_rotated_points(new_points);
}

void rotate_unpushed_cut() {
	Eigen::MatrixXd new_points = added_stroke->get3DPoints();
	viewer.data().rotate_points(new_points);
	added_stroke->set3DPoints(new_points);
	int nr_edges = new_points.rows() - 1;
	viewer.data().add_edges(new_points.topRows(nr_edges), new_points.middleRows(1, nr_edges), Eigen::RowVector3d(1, 0, 1));

	new_points = cut_surface_path.get3DPoints();
	viewer.data().rotate_points(new_points);
	cut_surface_path.set_rotated_points(new_points);
}

void handle_failed_cut() {
	sound_error_beep();
	prev_tool_mode = NONE;
	draw_all_strokes(); //Removes the drawn cut stroke
	Eigen::MatrixXd drawn_points = added_stroke->get3DPoints();
	int nr_edges = drawn_points.rows() - 1;
	viewer.data().add_edges(drawn_points.topRows(nr_edges), drawn_points.middleRows(1, nr_edges), black); //Display the stroke in black to show that it went wrong
	viewer.update_screen_while_computing = false;
}

void handle_failed_extrusion_base(bool hold_back) {
	sound_error_beep();
	prev_tool_mode = NONE;
	draw_all_strokes(); //Removes the drawn base stroke
	Eigen::MatrixXd drawn_points = extrusion_base->get3DPoints();
	int nr_edges = drawn_points.rows() - hold_back;
	viewer.data().add_edges(drawn_points.topRows(nr_edges), drawn_points.middleRows(1, nr_edges), black); //Display the stroke in black to show that it went wrong
	viewer.update_screen_while_computing = false;
}

void select_dragging_handle(Eigen::Vector3f& pos) {
	double closest_dist = INFINITY;
	double current_closest = closest_dist;
	closest_stroke_ID = -1;
	handleID = -1;
	int tmp_handleID;
	for (int i = 0; i < stroke_collection.size(); i++) { //Additional strokes that cross the original stroke will never be selected as the pulled curve when the user clicks a vertex that also belongs to the original boundary, since their vertex positions are the same and we check for SMALLER distances
		tmp_handleID = stroke_collection[i].selectClosestVertex(pos, closest_dist); //Returns the index into stroke3DPoints of the closest point
		if ((closest_dist < current_closest) && (tmp_handleID != -1)) {
			current_closest = closest_dist;
			handleID = tmp_handleID;
			closest_stroke_ID = i;
		}
	}
}

//TODO: just send right pos into here right away and remove use_right flag
void set_laser_points(Eigen::Vector3f& pos, bool use_right_direction) {
	if (!use_right_direction) {
		pos = viewer.oculusVR.get_left_hand_pos(); //Use the left hand position for the laser start point when requested
	}

	Eigen::MatrixX3d LP(2, 3);
	Eigen::MatrixXd laser_color(2, 3);
	vector<igl::Hit> hits;
	Eigen::Vector3f dir = use_right_direction ? viewer.oculusVR.get_right_touch_direction() : viewer.oculusVR.get_left_touch_direction();

	if (igl::ray_mesh_intersect(pos, dir, V, F, hits)) { //Intersect the ray from the Touch controller with the mesh to get the 3D point
		laser_end_point = (V.row(F(hits[0].id, 0))*(1.0 - hits[0].u - hits[0].v) + V.row(F(hits[0].id, 1))*hits[0].u + V.row(F(hits[0].id, 2))*hits[0].v);
	}
	else { //First check for intersections with the mesh, then with the floor and finally just set the end point at a far distance_to_vert
		viewer.selected_data_index = 0;
		bool hit_found = false;
		for (int i = 0; i < 6; i++) {
			if (igl::ray_mesh_intersect(pos, dir, viewer.data().V, viewer.data().F, hits)) {
				laser_end_point = (viewer.data().V.row(viewer.data().F(hits[0].id, 0))*(1.0 - hits[0].u - hits[0].v) + viewer.data().V.row(viewer.data().F(hits[0].id, 1))*hits[0].u + viewer.data().V.row(viewer.data().F(hits[0].id, 2))*hits[0].v);
				hit_found = true;
				break;
			}
			viewer.selected_data_index++;
		}

		if(!hit_found){
			laser_end_point = (pos + 10 * dir).cast<double>();
		}
		viewer.selected_data_index = base_mesh_index;
	}
	viewer.selected_data_index = pointer_mesh_index;
	LP.row(0) = pos.cast<double>();
	LP.row(1) = laser_end_point;
	laser_color.setZero();
	viewer.data().set_laser_points(LP, laser_color);
	viewer.selected_data_index = base_mesh_index;
}

void reset_before_draw() {
	cut_stroke_already_drawn = false;
	extrusion_base_already_drawn = false;
	(*base_mesh).V.resize(0, 3);
	(*base_mesh).F.resize(0, 3);
	(*base_mesh).vertex_is_fixed.resize(0);
	(*base_mesh).edge_boundary_markers.resize(0);
	(*base_mesh).new_mapped_indices.resize(0);
	(*base_mesh).sharp_edge.resize(0);
	(*base_mesh).mesh_to_patch_indices.resize(0);
	(*base_mesh).patches.clear(); //Request new patches
	(*base_mesh).face_patch_map.clear();
	viewer.data().clear();
	stroke_collection.clear();
	next_added_stroke_ID = 2;
}

void button_down(OculusVR::ButtonCombo pressed, Eigen::Vector3f& pos_left, Eigen::Vector3f& pos_right) {
	if (pressed == OculusVR::ButtonCombo::TRIG_RIGHT && (selected_tool_mode == ADD || selected_tool_mode == CUT || selected_tool_mode == EXTRUDE)) {
		if (stroke_collection.size() == 0) { //Don't go into these modes when there is no mesh yet
			prev_tool_mode = FAIL;
			return;
		}
	}
	else if (pressed == OculusVR::ButtonCombo::A && (selected_tool_mode == REMOVE || selected_tool_mode == CHANGE)) {
		if (stroke_collection.size() <= 1) {
			prev_tool_mode = FAIL;
			return;
		}
	}
	else if (pressed == OculusVR::ButtonCombo::THUMB_MOVE) {
		viewer.selected_data_index = pointer_mesh_index;
		if (prev_tool_mode != NAVIGATE) {
			prev_laser_show = viewer.data().show_laser;
		}
		viewer.data().show_laser = false;
		viewer.selected_data_index = base_mesh_index;
		V = viewer.data().V;
		for (int i = 0; i < stroke_collection.size(); i++) {
			stroke_collection[i].update_Positions(V, false);
		}
		draw_all_strokes();

		if (extrusion_base_already_drawn) {
			rotate_unpushed_extrusion();
		}
		if (cut_stroke_already_drawn) {
			rotate_unpushed_cut();
		}

		prev_tool_mode = NAVIGATE;
		return;
	}

	if (pressed == OculusVR::ButtonCombo::TRIG_RIGHT) {
		tool_mode = selected_tool_mode;
	}
	else if (pressed == OculusVR::ButtonCombo::NONE) {
		tool_mode = NONE;
	}
	else if (pressed == OculusVR::ButtonCombo::A) {
		if (selected_tool_mode == CHANGE || selected_tool_mode == REMOVE) {
			tool_mode = selected_tool_mode;
		}
		else {
			tool_mode = FAIL; //Pressing A is only allowed in CHANGE and REMOVE modes
		}
	}
	else if (pressed == OculusVR::ButtonCombo::X) {
		set_laser_points(pos_left, false);
		tool_mode = SMOOTH;
	}
	else if (pressed == OculusVR::ButtonCombo::Y) {
		//pos_left = viewer.oculusVR.get_left_hand_pos();
		set_laser_points(pos_left, false);
		if (!upsample_activated) {
			viewer.selected_data_index = pointer_mesh_index;
			prev_laser_show = viewer.data().show_laser;
			viewer.data().show_laser = true;
			viewer.oculusVR.right_hand_visible = false;
			viewer.selected_data_index = base_mesh_index;
		}
		if (tool_mode != UPSAMPLE) {
			tool_mode = UPSAMPLE;
		}
	}
	else if (pressed == OculusVR::ButtonCombo::GRIPTRIGBOTH) {
		tool_mode = TRANSLATE;
	}
	//Set laser points after button handling in case we have for example changed the position var
	set_laser_points(pos_right, !upsample_init); //When we're in UPSAMPLE mode (between first and second click), we want to show the laser of the left hand


	if (tool_mode == DRAW) {
		if (draw_should_block) { //User has been too close to first sample point (closing the stroke too much), so we're in blocked state till the buttons are released again
			return;
		}
		if (prev_tool_mode == NONE) {
			if (has_recentered) {
				viewer.oculusVR.set_start_action_view(viewer.core.get_view());
				added_stroke = new Stroke(&V, &F, 1);
				added_stroke->stroke_color = blue;
				added_stroke->addSegment(pos_right, viewer);
				prev_tool_mode = DRAW;
			}
			else {
				reset_before_draw();
				viewer.oculusVR.request_recenter();
				has_recentered = true;
			}
		}
		else if (prev_tool_mode == DRAW) {
			draw_should_block = added_stroke->addSegment(pos_right, viewer);
			return;
		}
	}
	else if (tool_mode == ADD) {
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) { //Adding a new control curve onto an existing mesh
			added_stroke = new Stroke(&V, &F, next_added_stroke_ID);
			added_stroke->stroke_color = blue;
			added_stroke->is_loop = false;
			added_stroke->addSegmentAdd(pos_right, viewer);
			prev_tool_mode = ADD;
		}
		else if (prev_tool_mode == ADD) {
			added_stroke->addSegmentAdd(pos_right, viewer);
		}
	}
	else if (tool_mode == REMOVE) {
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) {
			/*Eigen::Vector3f hit_pos;
			vector<igl::Hit> hits;

			if (igl::ray_mesh_intersect(pos, viewer.oculusVR.get_right_touch_direction(), V, F, hits)) { //Intersect the ray from the Touch controller with the mesh to get the 3D point
				hit_pos = (V.row(F(hits[0].id, 0))*(1.0 - hits[0].u - hits[0].v) + V.row(F(hits[0].id, 1))*hits[0].u + V.row(F(hits[0].id, 2))*hits[0].v).cast<float>();
			}
			else { //Hand ray did not intersect mesh
				prev_tool_mode = FAIL;
				return;
			}*/

			double closest_dist = INFINITY;
			double current_closest = closest_dist;
			int tmp_handleID, closest_stroke_idx;
			handleID = -1;
			for (int i = 1; i < stroke_collection.size(); i++) { //Skip stroke 0
				tmp_handleID = stroke_collection[i].selectClosestVertex(pos_right, closest_dist);

				//tmp_handleID = stroke_collection[i].selectClosestVertex(hit_pos, closest_dist);
				if ((closest_dist < current_closest) && (tmp_handleID != -1)) {
					current_closest = closest_dist;
					handleID = tmp_handleID;
					closest_stroke_ID = stroke_collection[i].get_ID();
					closest_stroke_idx = i;
				}
			}

			if (handleID == -1) {//User clicked too far from any of the stroke vertices
				prev_tool_mode = FAIL;
				return;
			}
			else if (closest_stroke_ID == prev_closest_stroke_ID) {
				remove_stroke_clicked++;
				draw_all_strokes();
			}
			else {
				remove_stroke_clicked = 1; //Start from 1
				if (prev_closest_stroke_idx != -1) {
					stroke_collection[prev_closest_stroke_idx].stroke_color = old_selected_stroke_color; //Reset the color of the previously selected stroke
				}
				//Redraw the original stroke and all added strokes, where the selected stroke is drawn in black.
				old_selected_stroke_color = stroke_collection[closest_stroke_idx].stroke_color;
				stroke_collection[closest_stroke_idx].stroke_color = black;
				prev_closest_stroke_idx = closest_stroke_idx;
				prev_closest_stroke_ID = closest_stroke_ID;
				draw_all_strokes();
			}

			if (remove_stroke_clicked == 2) { //Mechanism to force the user to click twice on the same stroke before removing it (safeguard)
				stroke_was_removed = true;
				stroke_collection[closest_stroke_idx].undo_stroke_add(edge_boundary_markers, sharp_edge, vertex_is_fixed);

				//We might have changed the patch structure (e.g. when removing a sharp cut stroke), so request new patches
				(*base_mesh).patches.clear();
				(*base_mesh).face_patch_map.clear();
				(*base_mesh).patches = Patch::init_patches(*base_mesh);

				stroke_collection.erase(stroke_collection.begin() + closest_stroke_idx);
				remove_stroke_clicked = 0; //Reset
				prev_closest_stroke_idx = -1;
			}

			prev_tool_mode = REMOVE;
		}
		else if (prev_tool_mode == REMOVE) {
			return; //For REMOVE we only take action upon button press and release 
		}
	}
	else if (tool_mode == PULL) { //Dragging an existing curve
		//FAIL below is needed when initialy hand position is too far away from the curves. Switch to FAIL instead of NONE to ensure proper button-release handling
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) {
			draw_all_strokes();
			select_dragging_handle(pos_right);

			if (handleID == -1) { //User clicked too far from any of the stroke vertices
				prev_tool_mode = FAIL;
				return;
			}
			CurveDeformation::startPullCurve(handleID, (*base_mesh).V, (*base_mesh).F, (*base_mesh).edge_boundary_markers);
			prev_tool_mode = PULL;

			//Record the mean of the minimum angle per triangle for each patch (value at start of pull)
			pull_path_starting_angles.resize((*base_mesh).patches.size());
			for (int i = 0; i < (*base_mesh).patches.size(); i++) {
				Eigen::MatrixXd angles;
				igl::internal_angles((*base_mesh).patches[i]->mesh.V, (*base_mesh).patches[i]->mesh.F, angles);
				pull_path_starting_angles[i] = angles.rowwise().minCoeff().mean();
			}
		}
		else if (prev_tool_mode == PULL) {
			CurveDeformation::pullCurve(pos_right.transpose().cast<double>(), (*base_mesh).V, (*base_mesh).edge_boundary_markers);

			for (int i = 0; i < (*base_mesh).patches.size(); i++) {
				(*base_mesh).patches[i]->update_patch_vertex_positions((*base_mesh).V);
			}

			int mesh_dependent_smoothing_iter = (int)(12500 / V.rows() + 1);
			for (int i = 0; i < mesh_dependent_smoothing_iter; i++) {
				SurfaceSmoothing::smooth(*base_mesh, dirty_boundary, false);
				for (int i = 0; i < stroke_collection.size(); i++) {
					stroke_collection[i].update_Positions(V, false);
				}
			}
			viewer.data().set_mesh(V, F);
			Eigen::MatrixXd N_corners;
			igl::per_corner_normals(V, F, 50, N_corners);
			viewer.data().set_normals(N_corners);
			draw_all_strokes();
		}
	}
	else if (tool_mode == CUT) {
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) {
			prev_tool_mode = CUT;
			if (cut_stroke_already_drawn) { //Clicked while cut stroke already drawn
				return;
			}

			draw_all_strokes(); //Will remove a possible old black (wrong) stroke
			added_stroke = new Stroke(&V, &F, next_added_stroke_ID); //For CUT, only increase the next added stroke ID upon button release (when we know it succeeded)
			added_stroke->stroke_color = red;
			viewer.oculusVR.set_start_action_view(viewer.core.get_view());
			added_stroke->addSegmentCut(pos_right, viewer);
		}
		else if (prev_tool_mode == CUT) {
			if (!cut_stroke_already_drawn) {
				added_stroke->addSegmentCut(pos_right, viewer);
			}
		}
	}
	else if (tool_mode == EXTRUDE) {
		if (extrude_base_should_block) {
			return;
		}
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) {
			prev_tool_mode = EXTRUDE;
			if (extrusion_base_already_drawn) {
				added_stroke = new Stroke(&V, &F, next_added_stroke_ID);
				added_stroke->stroke_color = blue;
				next_added_stroke_ID++;
				added_stroke->addSegmentExtrusionSilhouette(pos_right, viewer);
			}
			else {
				extrusion_base = new Stroke(&V, &F, next_added_stroke_ID);
				extrusion_base->stroke_color = red;
				next_added_stroke_ID++;
				extrude_base_should_block = extrusion_base->addSegmentExtrusionBase(pos_right, viewer);
			}
		}
		else if (prev_tool_mode == EXTRUDE) {
			if (extrusion_base_already_drawn) {
				added_stroke->addSegmentExtrusionSilhouette(pos_right, viewer);
			}
			else {
				extrude_base_should_block = extrusion_base->addSegmentExtrusionBase(pos_right, viewer);
			}
		}
	}
	else if (tool_mode == SMOOTH) {
		if (prev_tool_mode == NONE) {		
			prev_tool_mode = SMOOTH;

			//Turn lasers off
			viewer.selected_data_index = pointer_mesh_index;
			prev_laser_show = viewer.data().show_laser;
			viewer.data().show_laser = false;
			viewer.oculusVR.right_hand_visible = false;
			viewer.selected_data_index = base_mesh_index;

			select_dragging_handle(pos_left);
			if (handleID == -1) { //Perform general mesh smooth
				return;
			}
			rub_seamID = stroke_collection[closest_stroke_ID].stroke_ID;
			CurveRub::start_rubbing(handleID, V, F, edge_boundary_markers);
		}
		else if (prev_tool_mode == SMOOTH) {
			if (handleID == -1) {
				return;
			}
			Eigen::Vector3d double_pos = pos_left.cast<double>();
			Eigen::MatrixXd prevV = V;
			CurveRub::rubbing(double_pos, rub_seamID, V);
			for (int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V, false);
			}
			draw_all_strokes();
		}
	}
	else if (tool_mode == UPSAMPLE) {
		if (prev_tool_mode == NONE) {
			if (upsample_activated) {
				vector<igl::Hit> hits;
				if (igl::ray_mesh_intersect(pos_right, viewer.oculusVR.get_left_touch_direction(), V, F, hits)) { //Intersect the ray from the Touch controller with the mesh to get the 3D point
					Patch* hit_patch = ((*base_mesh).face_patch_map[hits[0].id]); //Find the patch that was hit

					for (int i = 0; i < (*base_mesh).patches.size(); i++) {
						(*base_mesh).patches[i]->update_patch_vertex_positions((*base_mesh).V);
					}

					double cur_avg_edge_length = igl::avg_edge_length((*hit_patch).mesh.V, (*hit_patch).mesh.F);
					if ((*hit_patch).mesh.F.rows() * 8 < Stroke::MAX_NR_TRIANGLES ) { //Make sure that upsampling doesn't result in too many faces
						(*hit_patch).upsample_patch((*base_mesh).V, (*base_mesh).F, (*base_mesh).face_patch_map, (*base_mesh).edge_boundary_markers, (*base_mesh).sharp_edge, (*base_mesh).vertex_is_fixed, replacing_vertex_bindings);
					}
					else {
						sound_error_beep();
						std::cerr << "Error: Could not upsample selected patch, as this would result in a patch with too many triangles (making the system unbearably slow). " << std::endl;
						prev_tool_mode = UPSAMPLE;
						return;
					}
					
					if (igl::is_edge_manifold((*base_mesh).F)) {
						(*base_mesh).patches.clear();
						(*base_mesh).face_patch_map.clear();
						(*base_mesh).patches = Patch::init_patches(*base_mesh);

						new_mapped_indices.setLinSpaced((*base_mesh).V.rows(), 0, (*base_mesh).V.rows() - 1);
						for (int i = 0; i < stroke_collection.size(); i++) {
							stroke_collection[i].update_vert_bindings(new_mapped_indices, replacing_vertex_bindings);
						}

						dirty_boundary = true;

						for (int i = 0; i < 8; i++) {
							SurfaceSmoothing::smooth(*base_mesh, dirty_boundary, false);
						}

						for (int i = 0; i < stroke_collection.size(); i++) {
							stroke_collection[i].update_Positions(V, true);
						}
					}
					else {						
						std::cerr << "Error: This shouldn't happen. Resampled mesh is no longer manifold" << std::endl;
					}
					
					viewer.data().clear();
					viewer.data().set_mesh(V, F);
					Eigen::MatrixXd N_corners;
					igl::per_corner_normals(V, F, 50, N_corners);
					viewer.data().set_normals(N_corners);

					draw_all_strokes();
				}
			}
			else {
				upsample_init = true;
			}
			prev_tool_mode = UPSAMPLE;
		}
		else if (prev_tool_mode == UPSAMPLE) {
			return; //Only take action on button press and release (not continuosly)
		}
	}
	else if (tool_mode == CHANGE) {
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) {
			double closest_dist = INFINITY;
			double current_closest = closest_dist;
			int tmp_handleID, closest_stroke_idx;
			handleID = -1;
			for (int i = 0; i < stroke_collection.size(); i++) {
				tmp_handleID = stroke_collection[i].selectClosestVertex(pos_right, closest_dist);

				if ((closest_dist < current_closest) && (tmp_handleID != -1)) {
					current_closest = closest_dist;
					handleID = tmp_handleID;
					closest_stroke_ID = stroke_collection[i].get_ID();
					closest_stroke_idx = i;
				}
			}

			if (handleID == -1) { //User clicked too far from any of the stroke vertices
				prev_tool_mode = FAIL;
				return;
			}
			else {
				stroke_collection[closest_stroke_idx].stroke_color = stroke_collection[closest_stroke_idx].stroke_color == red ? blue : red;
				stroke_collection[closest_stroke_idx].switch_stroke_edges_type(sharp_edge);
				draw_all_strokes();
				prev_tool_mode = CHANGE;
			}
		}
		else if (prev_tool_mode == CHANGE) {
			return; //Only take action upon button press and release
		}
	}
	else if (tool_mode == TRANSLATE) {
		if (V.rows() == 0) {
			return;
		}
		if (prev_tool_mode == NONE) {
			viewer.selected_data_index = pointer_mesh_index;
			prev_laser_show = viewer.data().show_laser;
			viewer.data().show_laser = false;
			viewer.selected_data_index = base_mesh_index;

			prev_translation_point = ((pos_right + pos_left) / 2.0).cast<double>();
			prev_scale_size = (pos_right - pos_left).norm();
			prev_tool_mode = TRANSLATE;
		}
		else if(prev_tool_mode == TRANSLATE) {
			Eigen::Vector3d cur_translation_point = ((pos_right + pos_left) / 2.0).cast<double>();
			Eigen::RowVector3d translation_vec = (cur_translation_point - prev_translation_point).transpose();
			
			double cur_scale_size = (pos_right - pos_left).norm();
			double scale_factor = cur_scale_size / prev_scale_size;

			Eigen::Vector3f mesh_translation = -V.colwise().mean().eval().cast<float>();
			Eigen::MatrixXf V_tmp(4, V.rows());
			V_tmp.block(0, 0, 3, V.rows()) = V.transpose().cast<float>();
			V_tmp.row(3) = Eigen::RowVectorXf::Constant(V.rows(), 1);

			Eigen::Matrix4f scaling = Eigen::Matrix4f::Identity();
			scaling.block(0, 0, 3, 3) *= scale_factor;

			Eigen::Matrix4f mesh_to_origin_translation = Eigen::Matrix4f::Identity();
			mesh_to_origin_translation.col(3).head(3) = mesh_translation;

			Eigen::Matrix4f move_translation = Eigen::Matrix4f::Identity();
			move_translation.col(3).head(3) = translation_vec.transpose().cast<float>();

			Eigen::Matrix4f place_back_translation = Eigen::Matrix4f::Identity();
			place_back_translation.col(3).head(3) = -mesh_translation;

			V = (place_back_translation*move_translation*scaling*mesh_to_origin_translation*V_tmp).topRows(3).transpose().cast<double>();

			prev_translation_point = cur_translation_point;
			prev_scale_size = cur_scale_size;

			viewer.data().set_mesh(V, F);
			Eigen::MatrixXd N_corners;
			igl::per_corner_normals(V, F, 50, N_corners);
			viewer.data().set_normals(N_corners);

			draw_all_strokes();
		}
	}
	else if (tool_mode == NONE) { //Have to finish up as if we're calling mouse_up()
		has_recentered = false;
		if (prev_tool_mode == NONE) {
			return;
		}
		else if (prev_tool_mode == DRAW) {
			prev_tool_mode = NONE;
			draw_should_block = false;

			if (added_stroke->toLoop()) { //Returns false if the stroke only consists of 1 point (user just clicked)
				viewer.data().add_edges(added_stroke->stroke3DPoints.bottomRows(1), added_stroke->stroke3DPoints.row(added_stroke->stroke3DPoints.rows()-2), Eigen::RowVector3d(0, 0, 1)); //Shows closing stroke

				added_stroke->generate3DMeshFromStroke(edge_boundary_markers, vertex_is_fixed, V, F, viewer);

				if (!igl::is_edge_manifold(F)) { //Check if the drawn stroke results in an edge-manifold mesh, otherwise sound a beep and revert
					sound_error_beep();
					added_stroke->update_Positions(V, false);
					Eigen::MatrixXd drawn_points = added_stroke->get3DPoints();
					vertex_is_fixed.resize(0);
					edge_boundary_markers.resize(0);
					dirty_boundary = true;
					viewer.data().set_edges(Eigen::MatrixXd(), Eigen::MatrixXi(), red); //Clear the non-original stroke edges
					viewer.data().add_edges(drawn_points.block(0, 0, drawn_points.rows() - 1, 3), drawn_points.block(1, 0, drawn_points.rows() - 1, 3), black); //Display the stroke in black to show that it went wrong
					prev_tool_mode = NONE;
					std::cerr << "This stroke results in a mesh that is not edge manifold. Will not be able to process, try again (don't use self-intersecting strokes)." << std::endl;
					viewer.update_screen_while_computing = false;
					return;
				}

				Eigen::MatrixXi EV, FE, EF;
				igl::edge_topology(V, F, EV, FE, EF);
				sharp_edge.resize(EV.rows());
				sharp_edge.setZero(); //Set all edges to smooth after initial draw

				stroke_collection.push_back(*added_stroke);
				(*base_mesh).patches = Patch::init_patches(*base_mesh);

				dirty_boundary = true;
				stroke_collection.back().update_Positions(V, true);

				for (int j = 0; j < 10; j++) {
					SurfaceSmoothing::smooth(*base_mesh, dirty_boundary, false);

				}
				viewer.data().clear(); //Don't clear before smoothing, otherwise the user won't see the drawn stroke will the shape is inflating

				viewer.data().set_mesh(V, F);
				Eigen::MatrixXd N_corners;
				igl::per_corner_normals(V, F, 50, N_corners);
				viewer.data().set_normals(N_corners);
				
				draw_all_strokes();

				std::cerr << "Created a mesh with " << V.rows() << " vertices, and " << F.rows() << " faces." << std::endl;
			}
		}
		else if (prev_tool_mode == ADD) {
			dirty_boundary = true;
			if (!added_stroke->has_points_on_mesh) {
				prev_tool_mode = NONE;
				sound_error_beep();
				std::cerr << "The stroke you are trying to add does not have any points that lie on the mesh surface. Will not be able to process, please try again." << std::endl;
				viewer.update_screen_while_computing = false;
				return;
			}
			else if (added_stroke->get3DPoints().rows() < 2) {
				prev_tool_mode = NONE;
				sound_error_beep();
				std::cerr << "The stroke you are trying to add does not have enough points. Please try again." << std::endl;
				Eigen::MatrixXd drawn_points = added_stroke->get3DPoints();
				int nr_edges = drawn_points.rows() - 1;
				viewer.data().add_edges(drawn_points.topRows(nr_edges), drawn_points.middleRows(1, nr_edges), black); //Display the stroke in black to show that it went wrong
				viewer.update_screen_while_computing = false;
				return;
			}

			bool success;
			if (added_stroke->starts_on_mesh && added_stroke->ends_on_mesh && !added_stroke->has_been_outside_mesh) {
				success = LaplacianRemesh::remesh_open_path(*base_mesh, *added_stroke, replacing_vertex_bindings, viewer);
			}
			else if (!added_stroke->starts_on_mesh && !added_stroke->ends_on_mesh) { //Stroke like a cut stroke (starts and ends off mesh to wrap around)
			   //Need to remesh like it's a cut but without removing the inside faces (we get 2 loops of boundary vertices that both need to be stitched, so can't use remesh_open_path)
				added_stroke->prepend_first_point(viewer);
				added_stroke->append_final_point(viewer);
				added_stroke->toLoop();
				success = LaplacianRemesh::remesh_cutting_path(*base_mesh, *added_stroke, replacing_vertex_bindings, viewer);
			}
			else {
				sound_error_beep();
				std::cerr << "An added stroke either needs both the end- & startpoint to be outside of the mesh, or all points to be on the mesh. Please try again. " << std::endl;
				success = false;
			}

			if (!success) {
				prev_tool_mode = NONE;
				sound_error_beep();
				draw_all_strokes();
				Eigen::MatrixXd drawn_points = added_stroke->get3DPoints();
				int nr_edges = drawn_points.rows() - 1;
				viewer.data().add_edges(drawn_points.topRows(nr_edges), drawn_points.middleRows(1, nr_edges), black); //Display the stroke in black to show that it went wrong
				viewer.update_screen_while_computing = false;
				return;
			}

			stroke_collection.push_back(*added_stroke);
			next_added_stroke_ID++;
			(*base_mesh).patches.clear();
			(*base_mesh).face_patch_map.clear();
			(*base_mesh).patches = Patch::init_patches(*base_mesh);

			int nr_removed = 0, original_collection_size = stroke_collection.size();
			for (int i = 0; i < original_collection_size - 1; i++) { //Don't update the added stroke, as this is done inside the remeshing already
				if (!stroke_collection[i - nr_removed].update_vert_bindings(new_mapped_indices, replacing_vertex_bindings)) { //Stroke dies, don't need to do stroke.undo_stroke_add, cause all its vertices also cease to exist
					stroke_collection.erase(stroke_collection.begin() + i - nr_removed);
					nr_removed++;
					dirty_boundary = true;
				}
			}

			for (int i = 0; i < initial_smooth_iter / 3; i++) {
				SurfaceSmoothing::smooth(*base_mesh, dirty_boundary, false);
			}

			//Update the stroke positions, will update with resampled stroke points
			for (int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V, true);
			}

			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			Eigen::MatrixXd N_corners;
			igl::per_corner_normals(V, F, 50, N_corners);
			viewer.data().set_normals(N_corners);

			draw_all_strokes();
		}
		else if (prev_tool_mode == REMOVE && stroke_was_removed) { //Only redraw if we actually removed a stroke (otherwise we draw unnecessary)
			dirty_boundary = true;
			for (int i = 0; i < initial_smooth_iter / 3; i++) {
				SurfaceSmoothing::smooth(*base_mesh, dirty_boundary, false);
			}

			for (int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V, false);
			}

			viewer.data().set_mesh(V, F);
			Eigen::MatrixXd N_corners;
			igl::per_corner_normals(V, F, 50, N_corners);
			viewer.data().set_normals(N_corners);

			stroke_was_removed = false;
			draw_all_strokes();

			//Turn left on
			viewer.selected_data_index = pointer_mesh_index;
			viewer.oculusVR.left_hand_visible = true;
			viewer.selected_data_index = base_mesh_index;
		}
		else if (prev_tool_mode == PULL && handleID != -1) {
			Eigen::MatrixXd angles;
			double mean_min_angle;
			bool stroke_structure_changed = false;
			for (int i = 0; i < (*base_mesh).patches.size(); i++) {
				igl::internal_angles((*base_mesh).patches[i]->mesh.V, (*base_mesh).patches[i]->mesh.F, angles);
				mean_min_angle = angles.rowwise().minCoeff().mean();

				if (((pull_path_starting_angles[i] - mean_min_angle) > 0.1) && mean_min_angle < (30.0 * (igl::PI/180.0)) && (*base_mesh).patches[i]->mesh.F.rows() * 4 < Stroke::MAX_NR_TRIANGLES) { //Only upsample a patch if its mean minimum triangle angle has significantly changed since the pull started, the new mean min angle is too big and upsampling won't lead to too many triangles 
					(*((*base_mesh).patches[i])).upsample_patch((*base_mesh).V, (*base_mesh).F, (*base_mesh).face_patch_map, (*base_mesh).edge_boundary_markers, (*base_mesh).sharp_edge, (*base_mesh).vertex_is_fixed, replacing_vertex_bindings);
				
					(*base_mesh).patches.clear();
					(*base_mesh).face_patch_map.clear();
					(*base_mesh).patches = Patch::init_patches(*base_mesh);
					
					new_mapped_indices.setLinSpaced((*base_mesh).V.rows(), 0, (*base_mesh).V.rows() - 1);
					for (int i = 0; i < stroke_collection.size(); i++) {
						stroke_collection[i].update_vert_bindings(new_mapped_indices, replacing_vertex_bindings);
					}
					stroke_structure_changed = true;
					dirty_boundary = true;
				}
			}

			for (int i = 0; i < initial_smooth_iter / 3; i++) {
				SurfaceSmoothing::smooth(*base_mesh, dirty_boundary, false);
			}

			for (int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V, stroke_structure_changed);
			}

			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			Eigen::MatrixXd N_corners;
			igl::per_corner_normals(V, F, 50, N_corners);
			viewer.data().set_normals(N_corners);
			draw_all_strokes();
		}
		else if (prev_tool_mode == CUT) {
			if (!added_stroke->has_points_on_mesh) {
				prev_tool_mode = NONE;
				viewer.update_screen_while_computing = false;
				return;
			}
			if (cut_stroke_already_drawn) { //User had already drawn the cut stroke and has now clicked/drawn the final stroke for removing the part. User should click on the part that should be removed
				dirty_boundary = true;
				vector<igl::Hit> hits;
				if (!igl::ray_mesh_intersect(pos_right, viewer.oculusVR.get_right_touch_direction(), V, F, hits)) { //Intersect the ray from the Touch controller with the mesh to get the 3D point
					viewer.update_screen_while_computing = false;
					prev_tool_mode = NONE;
					return;
				}
				int clicked_face = hits[0].id;
				int cut_success = MeshCut::cut(*base_mesh, *added_stroke, cut_surface_path, clicked_face, replacing_vertex_bindings, viewer);
				if (cut_success == 0) { //Catches the following cases: when the cut removes all mesh vertices/faces, when the first/last cut point aren't correct (face -1 in SurfacePath) or when sorting takes "infinite" time
					cut_stroke_already_drawn = false;
					handle_failed_cut();
					return;
				}
				else if (cut_success == -1) { //User clicked too close to the cutting line
					sound_error_beep();
					viewer.update_screen_while_computing = false;
					prev_tool_mode = NONE;
					return;
				}
				(*base_mesh).patches.clear();
				(*base_mesh).face_patch_map.clear();
				(*base_mesh).patches = Patch::init_patches(*base_mesh);

				stroke_collection.push_back(*added_stroke);
				next_added_stroke_ID++; //Stroke is successfully pushed back

				int nr_removed = 0, original_collection_size = stroke_collection.size();
				for (int i = 0; i < original_collection_size; i++) {
					if (!stroke_collection[i - nr_removed].update_vert_bindings(new_mapped_indices, replacing_vertex_bindings)) { //Stroke dies, don't need to do stroke.undo_stroke_add, cause all its vertices also cease to exist
						stroke_collection.erase(stroke_collection.begin() + i - nr_removed);
						nr_removed++;
						dirty_boundary = true;
					}
				}

				for (int i = 0; i < initial_smooth_iter / 2; i++) {
					SurfaceSmoothing::smooth(*base_mesh, dirty_boundary, false);
				}

				//Update the stroke positions after smoothing, will also add resampled points
				for (int i = 0; i < stroke_collection.size(); i++) {
					stroke_collection[i].update_Positions(V, true);
				}

				viewer.data().clear();
				viewer.data().set_mesh(V, F);
				Eigen::MatrixXd N_corners;
				igl::per_corner_normals(V, F, 50, N_corners);
				viewer.data().set_normals(N_corners);

				cut_stroke_already_drawn = false;
				draw_all_strokes();
			}
			else { //We're finished drawing the cut stroke, prepare for when user draws the final stroke to remove the part
				if (added_stroke->has_self_intersection(false)) {
					std::cerr << "Cut stroke contains a loop which is illegal. Please try again. " << std::endl;
					handle_failed_cut();
					return;
				}

				added_stroke->prepend_first_point(viewer);
				added_stroke->append_final_point(viewer);
				added_stroke->toLoop();

				bool success_cut_prepare = MeshCut::cut_prepare(*added_stroke, cut_surface_path);
				if (!success_cut_prepare) {
					handle_failed_cut();
					return;
				}
				draw_all_strokes();

				//Overlay the pink cut stroke
				Eigen::MatrixXd drawn_points = added_stroke->get3DPoints();
				int nr_edges = drawn_points.rows() - 1;
				viewer.data().add_edges(drawn_points.topRows(nr_edges), drawn_points.middleRows(1, nr_edges), Eigen::RowVector3d(1, 0, 1));
				cut_stroke_already_drawn = true;
			}
		}
		else if (prev_tool_mode == EXTRUDE) {
			if (extrusion_base_already_drawn) { //User has drawn the silhouette stroke for extrusion
				dirty_boundary = true;
				added_stroke->toLoop();
				added_stroke->is_loop = false; //Set to false manually, because we don't want curveDeformation to consider it as a loop (but we do need looped 3DPoints)

				bool success_extrude = false;
				if (added_stroke->get3DPoints().rows() > 6) {
					success_extrude = MeshExtrusion::extrude_main(*base_mesh, base_surface_path, *added_stroke, *extrusion_base, base_model, base_view, base_proj, base_viewport, replacing_vertex_bindings);
				}
				else {
					std::cerr << "Drawn stroke has too little sample points. Please try again." << std::endl;
				}

				if (!success_extrude) { //Catches the case that the extrusion base removes all faces/vertices or that the extrusion silhouette has too little samples/is just a click
					next_added_stroke_ID -= 2;
					extrusion_base_already_drawn = false;
					handle_failed_extrusion_base(true);
					Eigen::MatrixXd drawn_points = added_stroke->get3DPoints();
					int nr_edges = drawn_points.rows() - 1;
					viewer.data().add_edges(drawn_points.topRows(nr_edges - 1), drawn_points.middleRows(1, nr_edges - 1), black); //Display the stroke in black to show that it went wrong
					viewer.selected_data_index = pointer_mesh_index;
					viewer.data().show_laser = true;
					viewer.selected_data_index = base_mesh_index;
					return;
				}


				(*base_mesh).patches.clear();
				(*base_mesh).face_patch_map.clear();
				(*base_mesh).patches = Patch::init_patches(*base_mesh);

				stroke_collection.push_back(*extrusion_base);
				stroke_collection.push_back(*added_stroke);

				int nr_removed = 0, original_collection_size = stroke_collection.size();
				for (int i = 0; i < original_collection_size - 1; i++) { //Skip the newly added stroke since it is already updated inside extrude_main()
					if (!stroke_collection[i - nr_removed].update_vert_bindings(new_mapped_indices, replacing_vertex_bindings)) {
						//Stroke dies, don't need to do stroke.undo_stroke_add, cause all its vertices also cease to exist
						stroke_collection.erase(stroke_collection.begin() + i - nr_removed);
						nr_removed++;
						dirty_boundary = true;
					}
				}

				for (int i = 0; i < initial_smooth_iter / 2; i++) {
					SurfaceSmoothing::smooth(*base_mesh, dirty_boundary, false);
				}

				//Update the stroke positions after smoothing, also updates resampled stroke points
				for (int i = 0; i < stroke_collection.size(); i++) {
					stroke_collection[i].update_Positions(V, true);
				}

				viewer.data().clear();
				viewer.data().set_mesh(V, F);
				Eigen::MatrixXd N_corners;
				igl::per_corner_normals(V, F, 50, N_corners);
				viewer.data().set_normals(N_corners);

				draw_all_strokes();
				viewer.selected_data_index = pointer_mesh_index;
				viewer.data().show_laser = true;
				viewer.oculusVR.left_hand_visible = false; 			//Turn left off again
				viewer.selected_data_index = base_mesh_index;
				extrusion_base_already_drawn = false;
			}
			else { //mouse released after extrusion base drawn
				extrude_base_should_block = false;
				extrusion_base->toLoop();
				if (!extrusion_base->has_points_on_mesh) {
					next_added_stroke_ID--;
					viewer.update_screen_while_computing = false;
					prev_tool_mode = NONE;
					return;
				}
				else if (extrusion_base->has_been_outside_mesh || extrusion_base->has_self_intersection(true)) {
					next_added_stroke_ID--;
					if (extrusion_base->has_self_intersection(true)) {
						std::cerr << "Extrusion base contains a loop which is illegal. Please try again. " << std::endl;
					}
					else {
						std::cerr << "Extrusion base has parts which are drawn outside of the mesh surface, which is not allowed. Try again. " << std::endl;
					}
					handle_failed_extrusion_base(true);
					return;
				}

				dirty_boundary = true;

				base_model = viewer.core.get_model();
				base_view = viewer.oculusVR.get_start_action_view(); //Use the same view as has been used for the initial drawing of the base on the mesh, otherwise it is guaranteed to shift
				base_proj = viewer.core.get_proj();
				base_viewport = viewer.core.viewport;

				extrusion_base->resample_and_smooth_3DPoints(base_model, base_view, base_proj, base_viewport);

				bool succes_extrude_prepare = MeshExtrusion::extrude_prepare(*extrusion_base, base_surface_path); //Don't need to update all strokes here, since it didn't remove any vertices
				if (!succes_extrude_prepare) { //Catches the case that face == -1 in SurfacePath
					next_added_stroke_ID--;
					handle_failed_extrusion_base(true);
					return;
				}

				extrusion_base_already_drawn = true;

				

				draw_all_strokes();
				draw_extrusion_base(); //Need to draw the extrusion base separately, since it isn't added to the stroke_collection yet.
				viewer.selected_data_index = pointer_mesh_index;
				viewer.data().show_laser = false;
				viewer.oculusVR.left_hand_visible = true; 			//Turn left on
				viewer.selected_data_index = base_mesh_index;
			}
		}
		else if (prev_tool_mode == FAIL) {
			viewer.update_screen_while_computing = false;
		}
		else if (prev_tool_mode == SMOOTH) {
				if (V.rows() == 0) {
					prev_tool_mode = NONE;
					viewer.update_screen_while_computing = false;
					return;
				}

				for (int i = 0; i < (*base_mesh).patches.size(); i++) {
					(*base_mesh).patches[i]->update_patch_vertex_positions((*base_mesh).V);
				}

				for (int i = 0; i < 8; i++) {
					SurfaceSmoothing::smooth(*base_mesh, dirty_boundary, false);
				}
				for (int i = 0; i < stroke_collection.size(); i++) {
					stroke_collection[i].update_Positions(V, false);
				}

				viewer.data().set_mesh(V, F);
				Eigen::MatrixXd N_corners;
				igl::per_corner_normals(V, F, 50, N_corners);
				viewer.data().set_normals(N_corners);

				draw_all_strokes();
				viewer.oculusVR.right_hand_visible = true;
		}
		else if (prev_tool_mode == UPSAMPLE) {
			if (upsample_activated)  { //Done, button release after 2nd click
				upsample_activated = false; 
				upsample_init = false;
				viewer.selected_data_index = pointer_mesh_index;
				viewer.data().show_laser = prev_laser_show;
				viewer.oculusVR.right_hand_visible = true;
				viewer.selected_data_index = base_mesh_index;
			}
			else {
				upsample_activated = true;
			}
		}
		else if (prev_tool_mode == CHANGE) {
			//We might have changed the patch structure (e.g. when removing a sharp boundary stroke), so request new patches
			(*base_mesh).patches.clear();
			(*base_mesh).face_patch_map.clear();
			(*base_mesh).patches = Patch::init_patches(*base_mesh);
			dirty_boundary = true;

			SurfaceSmoothing::smooth(*base_mesh, dirty_boundary, true); //Force reset of the matrices. We might have changed a stroke type without changing the number of patches
			for (int i = 0; i < initial_smooth_iter / 2; i++) {
				SurfaceSmoothing::smooth(*base_mesh, dirty_boundary, false);
			}

			viewer.data().set_mesh(V, F);
			Eigen::MatrixXd N_corners;
			igl::per_corner_normals(V, F, 50, N_corners);
			viewer.data().set_normals(N_corners);

			draw_all_strokes();
		}
		else if (prev_tool_mode == TRANSLATE) {
			viewer.selected_data_index = pointer_mesh_index;
			viewer.data().show_laser = prev_laser_show;
			viewer.selected_data_index = base_mesh_index;
		}
		else if (prev_tool_mode == NAVIGATE) {
			viewer.selected_data_index = pointer_mesh_index;
			viewer.data().show_laser = prev_laser_show;
			viewer.selected_data_index = base_mesh_index;
		}

		prev_tool_mode = NONE;
		viewer.update_screen_while_computing = false;
		return;
	}
	else if (tool_mode == FAIL) {
		prev_tool_mode = FAIL;
	}

	return;
}

void set_floor_mesh(Viewer& viewer, Eigen::MatrixXd& V_floor, Eigen::MatrixXi& F_floor, std::string cur_dir) {
	viewer.data().set_mesh(V_floor, F_floor);
	Eigen::MatrixXd V_uv(V_floor.rows(), 2);
	V_uv.col(0) = V_floor.col(0);
	V_uv.col(1) = V_floor.col(2);
	V_uv.col(0) = V_uv.col(0).array() - V_uv.col(0).minCoeff();
	V_uv.col(0) = V_uv.col(0).array() / V_uv.col(0).maxCoeff();
	V_uv.col(1) = V_uv.col(1).array() - V_uv.col(1).minCoeff();
	V_uv.col(1) = V_uv.col(1).array() / V_uv.col(1).maxCoeff();
	V_uv = V_uv.array() * 2; //Amount of repetitions of the texture

	Eigen::Matrix<unsigned char, Dynamic, Dynamic> texR, texG, texB;
	std::string texture_file = cur_dir + "\\..\\data\\free\\floor6square.png";
	int width, height, n;
	unsigned char *data = stbi_load(texture_file.c_str(), &width, &height, &n, 4);

	if (!data) {
		std::cerr << "Could not load floor texture." << std::endl;
	}
	texR.resize(height, width);
	texG.resize(height, width);
	texB.resize(height, width);

	for (unsigned j = 0; j < height; ++j) {
		for (unsigned i = 0; i < width; ++i) {
			// used to flip with libPNG, but I'm not sure if
			// simply j*width + i wouldn't be better
			// stb_image uses horizontal scanline an starts top-left corner
			texR(i, j) = data[4 * ((width - 1 - i) + width * (height - 1 - j))];
			texG(i, j) = data[4 * ((width - 1 - i) + width * (height - 1 - j)) + 1];
			texB(i, j) = data[4 * ((width - 1 - i) + width * (height - 1 - j)) + 2];
		}
	}
	stbi_image_free(data);

	viewer.data().set_uv(V_uv);
	viewer.data().show_texture = true;
	viewer.data().show_lines = false;
	viewer.data().set_texture(texR, texG, texB);
}

void set_wall_mesh(Viewer& viewer, Eigen::MatrixXd& V_wall, Eigen::MatrixXi& F_wall, Eigen::Vector2d& cols_to_use, std::string cur_dir) {
	viewer.data().set_mesh(V_wall, F_wall);
	Eigen::MatrixXd V_uv(V_wall.rows(), 2);
	V_uv.col(0) = V_wall.col(cols_to_use(0));
	V_uv.col(1) = V_wall.col(cols_to_use(1));
	V_uv.col(0) = V_uv.col(0).array() - V_uv.col(0).minCoeff();
	V_uv.col(0) = V_uv.col(0).array() / V_uv.col(0).maxCoeff();
	V_uv.col(1) = V_uv.col(1).array() - V_uv.col(1).minCoeff();
	V_uv.col(1) = V_uv.col(1).array() / V_uv.col(1).maxCoeff();

	V_uv = V_uv.array() * 3; //Amount of repetitions of the texture

	Eigen::Matrix<unsigned char, Dynamic, Dynamic> texR, texG, texB;
	std::string texture_file = cur_dir + "\\..\\data\\free\\wall2square.jpg";
	int width, height, n;
	unsigned char *data = stbi_load(texture_file.c_str(), &width, &height, &n, 4);

	if (!data) {
		std::cerr << "Could not load wall texture." << std::endl;
	}
	texR.resize(height, width);
	texG.resize(height, width);
	texB.resize(height, width);

	for (unsigned j = 0; j < height; ++j) {
		for (unsigned i = 0; i < width; ++i) {
			// used to flip with libPNG, but I'm not sure if
			// simply j*width + i wouldn't be better
			// stb_image uses horizontal scanline an starts top-left corner
			texR(i, j) = data[4 * ((width - 1 - i) + width * (height - 1 - j))];
			texG(i, j) = data[4 * ((width - 1 - i) + width * (height - 1 - j)) + 1];
			texB(i, j) = data[4 * ((width - 1 - i) + width * (height - 1 - j)) + 2];
		}
	}
	stbi_image_free(data);

	viewer.data().set_uv(V_uv);
	viewer.data().show_texture = true;
	viewer.data().show_lines = false;
	viewer.data().set_texture(texR, texG, texB);
}

bool load_scene(std::string fname) {
	igl::deserialize(appState, "AppState", fname.c_str());

	stroke_collection.resize(appState.stroke_collection_size);
	for (int i = 0; i < appState.stroke_collection_size; i++) {
		igl::deserialize(stroke_collection[i], "obj" + std::to_string(i), fname.c_str());
		stroke_collection[i].setV(&V);
		stroke_collection[i].setF(&F);
	}

	edge_boundary_markers = appState.edge_boundary_markers;
	vertex_is_fixed = appState.vertex_is_fixed;
	sharp_edge = appState.sharp_edge;
	new_mapped_indices = appState.new_mapped_indices;
	replacing_vertex_bindings = appState.replacing_vertex_bindings;
	draw_all_strokes();

	return true;
}

bool save_scene(igl::opengl::glfw::Viewer& viewer, std::string fname) {
	viewer.save_scene(fname);
	igl::serialize(appState, "AppState", fname.c_str());

	for (int i = 0; i < stroke_collection.size(); i++) {
		igl::serialize(stroke_collection[i], "obj" + std::to_string(i), fname.c_str(), false);
	}
	return true;
}

bool callback_load_scene(Viewer& viewer, std::string cur_dir) {
	std::string fname = igl::file_dialog_open();
	if (fname.length() == 0) {
		return false;
	}
	viewer.data().clear();
	viewer.load_scene(fname);
	V = viewer.data().V;
	F = viewer.data().F;

	load_scene(fname);

	base_mesh = new Mesh(V, F, edge_boundary_markers, vertex_is_fixed, new_mapped_indices, sharp_edge, 0);
	(*base_mesh).patches.clear();
	(*base_mesh).face_patch_map.clear();
	(*base_mesh).patches = Patch::init_patches(*base_mesh);

	viewer.data().set_mesh(V, F);
	Eigen::MatrixXd N_corners;
	igl::per_corner_normals(V, F, 50, N_corners);
	viewer.data().set_normals(N_corners);
	return true;
}

bool callback_save_scene(Viewer& viewer, std::string cur_dir) {
	std::string fname = igl::file_dialog_save();
	if (fname.length() == 0) {
		return false;
	}
	save_scene(viewer, fname);
	return true;
}

void menu_opened() {
	menu.set_active();
	menu.callback_draw_viewer_window = menu_HUD;
	viewer.selected_data_index = pointer_mesh_index;
	prev_laser_show = viewer.data().show_laser;
	viewer.data().show_laser = true;
	viewer.selected_data_index = base_mesh_index;
	viewer.oculusVR.menu_active = true;
}

void menu_closed() {
	menu.set_inactive();
	menu.callback_draw_viewer_window = current_tool_HUD;
	viewer.selected_data_index = pointer_mesh_index;
	viewer.data().show_laser = prev_laser_show;
	viewer.selected_data_index = base_mesh_index;
	viewer.oculusVR.menu_active = false;
}

void reset_trackers() {
	if (extrusion_base_already_drawn) { //If we're switching from an unfinished extrusion
		next_added_stroke_ID--;
	}
	cut_stroke_already_drawn = false;
	extrusion_base_already_drawn = false;
	if (prev_closest_stroke_idx != -1) {
		stroke_collection[prev_closest_stroke_idx].stroke_color = old_selected_stroke_color; //Reset the color of the previously selected stroke
		prev_closest_stroke_ID = -1;
	}
	remove_stroke_clicked = 0;
	draw_all_strokes();
}

int main(int argc, char *argv[]) {
	//Init stroke selector
	added_stroke = new Stroke(&V, &F, 1);
	base_mesh = new Mesh(V, F, edge_boundary_markers, vertex_is_fixed, new_mapped_indices, sharp_edge, 0);

	viewer.core.light_position << 0.0f, -2.8f, -0.3f;

	Eigen::RowVector3d lfb(-2.5, 0, -2.5); //Left front bottom
	Eigen::RowVector3d lft(-2.5, 3, -2.5); //Left front top
	Eigen::RowVector3d rfb(2.5, 0, -2.5); //Right front bottom
	Eigen::RowVector3d rft(2.5, 3, -2.5); //Right front top
	Eigen::RowVector3d lbb(-2.5, 0, 2.5); //Left back bottom
	Eigen::RowVector3d lbt(-2.5, 3, 2.5); //Left back top
	Eigen::RowVector3d rbb(2.5, 0, 2.5); //Right back bottom
	Eigen::RowVector3d rbt(2.5, 3, 2.5); //Right back top


	V_floor.row(0) = lfb;
	V_floor.row(1) = rfb;
	V_floor.row(2) = rbb;
	V_floor.row(3) = lbb;

	F_floor.row(0) << 0, 3, 1;
	F_floor.row(1) << 3, 2, 1;

	V_roof.row(0) = lft;
	V_roof.row(1) = rft;
	V_roof.row(2) = rbt;
	V_roof.row(3) = lbt;

	F_roof.row(0) << 1, 3, 0;
	F_roof.row(1) << 1, 2, 3;

	V_left_wall.row(0) = lfb;
	V_left_wall.row(1) = lbb;
	V_left_wall.row(2) = lbt;
	V_left_wall.row(3) = lft;

	F_left_wall.row(0) << 1, 0, 3;
	F_left_wall.row(1) << 1, 3, 2;

	V_front_wall.row(0) = lfb;
	V_front_wall.row(1) = rfb;
	V_front_wall.row(2) = rft;
	V_front_wall.row(3) = lft;
	
	F_front_wall.row(0) << 1, 3, 0;
	F_front_wall.row(1) << 1, 2, 3;

	V_right_wall.row(0) = rfb;
	V_right_wall.row(1) = rbb;
	V_right_wall.row(2) = rbt;
	V_right_wall.row(3) = rft;

	F_right_wall.row(0) << 1, 3, 0;
	F_right_wall.row(1) << 1, 2, 3;

	V_back_wall.row(0) = rbb;
	V_back_wall.row(1) = lbb;
	V_back_wall.row(2) = lbt;
	V_back_wall.row(3) = rbt;

	F_back_wall.row(0) << 1, 3, 0;
	F_back_wall.row(1) << 1, 2, 3;


	char cur_dir[256];
	GetCurrentDirectoryA(256, cur_dir);

	if (argc == 2) {
		// Read mesh
		igl::readOFF(argv[1], V, F);
		callback_load_scene(viewer, std::string(cur_dir));
	}
	else {
		set_floor_mesh(viewer, V_floor, F_floor, std::string(cur_dir));
		viewer.data().uniform_colors(
			Eigen::Vector3d(igl::SILVER_AMBIENT[0], igl::SILVER_AMBIENT[1], igl::SILVER_AMBIENT[2]),
			Eigen::Vector3d(igl::SILVER_DIFFUSE[0], igl::SILVER_DIFFUSE[1], igl::SILVER_DIFFUSE[2]),
			Eigen::Vector3d(igl::SILVER_SPECULAR[0], igl::SILVER_SPECULAR[1], igl::SILVER_SPECULAR[2]));

		viewer.append_mesh();
		Eigen::Vector2d cols_to_use(0, 2);
		set_wall_mesh(viewer, V_roof, F_roof, cols_to_use, std::string(cur_dir)); //Roof
		viewer.data().uniform_colors(
			Eigen::Vector3d(igl::SILVER_AMBIENT[0], igl::SILVER_AMBIENT[1], igl::SILVER_AMBIENT[2]),
			Eigen::Vector3d(igl::SILVER_DIFFUSE[0], igl::SILVER_DIFFUSE[1], igl::SILVER_DIFFUSE[2]),
			Eigen::Vector3d(igl::SILVER_SPECULAR[0], igl::SILVER_SPECULAR[1], igl::SILVER_SPECULAR[2]));

		viewer.append_mesh();
		cols_to_use << 2, 1;
		set_wall_mesh(viewer, V_left_wall, F_left_wall, cols_to_use, std::string(cur_dir));
		viewer.data().uniform_colors(
			Eigen::Vector3d(igl::SILVER_AMBIENT[0], igl::SILVER_AMBIENT[1], igl::SILVER_AMBIENT[2]),
			Eigen::Vector3d(igl::SILVER_DIFFUSE[0], igl::SILVER_DIFFUSE[1], igl::SILVER_DIFFUSE[2]),
			Eigen::Vector3d(igl::SILVER_SPECULAR[0], igl::SILVER_SPECULAR[1], igl::SILVER_SPECULAR[2]));

		viewer.append_mesh();
		cols_to_use << 0, 1;
		set_wall_mesh(viewer, V_front_wall, F_front_wall, cols_to_use, std::string(cur_dir));
		viewer.data().uniform_colors(
			Eigen::Vector3d(igl::SILVER_AMBIENT[0], igl::SILVER_AMBIENT[1], igl::SILVER_AMBIENT[2]),
			Eigen::Vector3d(igl::SILVER_DIFFUSE[0], igl::SILVER_DIFFUSE[1], igl::SILVER_DIFFUSE[2]),
			Eigen::Vector3d(igl::SILVER_SPECULAR[0], igl::SILVER_SPECULAR[1], igl::SILVER_SPECULAR[2]));

		viewer.append_mesh();
		cols_to_use << 2, 1;
		set_wall_mesh(viewer, V_right_wall, F_right_wall, cols_to_use, std::string(cur_dir));
		viewer.data().uniform_colors(
			Eigen::Vector3d(igl::SILVER_AMBIENT[0], igl::SILVER_AMBIENT[1], igl::SILVER_AMBIENT[2]),
			Eigen::Vector3d(igl::SILVER_DIFFUSE[0], igl::SILVER_DIFFUSE[1], igl::SILVER_DIFFUSE[2]),
			Eigen::Vector3d(igl::SILVER_SPECULAR[0], igl::SILVER_SPECULAR[1], igl::SILVER_SPECULAR[2]));

		viewer.append_mesh();
		cols_to_use << 0, 1;
		set_wall_mesh(viewer, V_back_wall, F_back_wall, cols_to_use, std::string(cur_dir));
		viewer.data().uniform_colors(
			Eigen::Vector3d(igl::SILVER_AMBIENT[0], igl::SILVER_AMBIENT[1], igl::SILVER_AMBIENT[2]),
			Eigen::Vector3d(igl::SILVER_DIFFUSE[0], igl::SILVER_DIFFUSE[1], igl::SILVER_DIFFUSE[2]),
			Eigen::Vector3d(igl::SILVER_SPECULAR[0], igl::SILVER_SPECULAR[1], igl::SILVER_SPECULAR[2]));


		viewer.append_mesh();
		base_mesh_index = viewer.selected_data_index;
		viewer.data().set_mesh(V, F);
		Eigen::MatrixXd N_corners;
		igl::per_corner_normals(V, F, 50, N_corners);
		viewer.data().set_normals(N_corners);
	}

	viewer.append_mesh(); //For laser ray/point
	viewer.data().show_laser = false;
	pointer_mesh_index = viewer.selected_data_index;
	viewer.selected_data_index = base_mesh_index;
	viewer.plugins.push_back(&menu);

	viewer.init_oculus();

	GLuint img_texture = 0, img_texture1 = 0, img_texture2 = 0, img_texture3 = 0, img_texture4 = 0, img_texture5 = 0, img_texture6 = 0, img_texture7 = 0, img_texture8 = 0, img_texture9 = 0;
	int img_width, img_height, nrChannels;
	std::string filename = std::string(cur_dir) + "\\..\\data\\free\\draw.png";
	unsigned char *img_data = stbi_load(filename.c_str(), &img_width, &img_height, &nrChannels, 4);
	if (!img_data) {
		std::cerr << "Could not load image 1." << std::endl;
	}

	glGenTextures(1, &img_texture);
	GLenum err = glGetError();
	glBindTexture(GL_TEXTURE_2D, img_texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	void* im_texID_draw = (void *)(intptr_t)img_texture;

	filename = std::string(cur_dir) + "\\..\\data\\free\\pull.png";
	img_data = stbi_load(filename.c_str(), &img_width, &img_height, &nrChannels, 4);
	if (!img_data) {
		std::cerr << "Could not load image 2." << std::endl;
	}

	glGenTextures(1, &img_texture2);
	glBindTexture(GL_TEXTURE_2D, img_texture2);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	void* im_texID_pull = (void *)(intptr_t)img_texture2;

	filename = std::string(cur_dir) + "\\..\\data\\free\\add.png";
	img_data = stbi_load(filename.c_str(), &img_width, &img_height, &nrChannels, 4);
	if (!img_data) {
		std::cerr << "Could not load image 3." << std::endl;
	}

	glGenTextures(1, &img_texture3);
	glBindTexture(GL_TEXTURE_2D, img_texture3);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	void* im_texID_add = (void *)(intptr_t)img_texture3;

	filename = std::string(cur_dir) + "\\..\\data\\free\\cut.png";
	img_data = stbi_load(filename.c_str(), &img_width, &img_height, &nrChannels, 4);
	if (!img_data) {
		std::cerr << "Could not load image 4." << std::endl;
	}

	glGenTextures(1, &img_texture4);
	glBindTexture(GL_TEXTURE_2D, img_texture4);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	void* im_texID_cut = (void *)(intptr_t)img_texture4;

	filename = std::string(cur_dir) + "\\..\\data\\free\\bump.png";
	img_data = stbi_load(filename.c_str(), &img_width, &img_height, &nrChannels, 4);
	if (!img_data) {
		std::cerr << "Could not load image 5." << std::endl;
	}

	glGenTextures(1, &img_texture5);
	glBindTexture(GL_TEXTURE_2D, img_texture5);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	void* im_texID_extrude = (void *)(intptr_t)img_texture5;

	filename = std::string(cur_dir) + "\\..\\data\\free\\remove.png";
	img_data = stbi_load(filename.c_str(), &img_width, &img_height, &nrChannels, 4);
	if (!img_data) {
		std::cerr << "Could not load image 6." << std::endl;
	}

	glGenTextures(1, &img_texture6);
	glBindTexture(GL_TEXTURE_2D, img_texture6);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	void* im_texID_remove = (void *)(intptr_t)img_texture6;

	filename = std::string(cur_dir) + "\\..\\data\\free\\save.png";
	img_data = stbi_load(filename.c_str(), &img_width, &img_height, &nrChannels, 4);
	if (!img_data) {
		std::cerr << "Could not load image 7." << std::endl;
	}

	glGenTextures(1, &img_texture7);
	glBindTexture(GL_TEXTURE_2D, img_texture7);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	void* im_texID_save = (void *)(intptr_t)img_texture7;

	filename = std::string(cur_dir) + "\\..\\data\\free\\load.png";
	img_data = stbi_load(filename.c_str(), &img_width, &img_height, &nrChannels, 4);
	if (!img_data) {
		std::cerr << "Could not load image 8." << std::endl;
	}

	glGenTextures(1, &img_texture8);
	glBindTexture(GL_TEXTURE_2D, img_texture8);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	void* im_texID_load = (void *)(intptr_t)img_texture8;

	filename = std::string(cur_dir) + "\\..\\data\\free\\change.png";
	img_data = stbi_load(filename.c_str(), &img_width, &img_height, &nrChannels, 4);
	if (!img_data) {
		std::cerr << "Could not load image 9." << std::endl;
	}

	glGenTextures(1, &img_texture9);
	glBindTexture(GL_TEXTURE_2D, img_texture9);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	void* im_texID_change = (void *)(intptr_t)img_texture9;

	void* im_texID_cur = im_texID_draw;

	current_tool_HUD = [&]() {
		ImGui::SetNextWindowSize(ImVec2(200.0f, 200.0f), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
		ImGuiStyle& style = ImGui::GetStyle();
		style.WindowPadding = ImVec2(0, 0);
		style.WindowRounding = 0.0f;
		style.WindowTitleAlign = ImVec2(0.5f, 0.5f);
		style.ItemInnerSpacing = ImVec2(0, 0);
		style.Colors[ImGuiCol_WindowBg] = ImVec4(0.00f, 0.00f, 0.00f, 1.0f);
		ImGui::Begin("Current Tool", 0, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoScrollbar);
		ImGui::SameLine(0.5*(200 - 178)); //Take distance_to_vert around image and offset by half of it
		ImGui::Image(im_texID_cur, ImVec2(178.0f, 178.0f), ImVec2(0, 0), ImVec2(1, 1), ImColor(255, 255, 255, 255), ImColor(0, 0, 0, 0));
		ImGui::End();
	};

	menu_HUD = [=, &im_texID_cur]() {
		const ImVec2 texsize = ImGui_ImplGlfwGL3_GetTextureSize();
		ImColor icon_background_color(0, 0, 0, 255);
		ImVec2 uv_start(0, 0);
		ImVec2 uv_end(1, 1);
		ImVec2 icon_size(320, 320);
		ImGuiStyle& style = ImGui::GetStyle();
		style.WindowPadding = ImVec2(0, 0);
		style.WindowRounding = 0.0f;
		style.WindowTitleAlign = ImVec2(0.5f, 0.5f);
		style.ItemInnerSpacing = ImVec2(0, 0);
		style.Colors[ImGuiCol_WindowBg] = ImVec4(0.00f, 0.00f, 0.00f, 1.0f);

		ImGui::SetNextWindowSize(ImVec2(1008.0f, 1036.0f), ImGuiSetCond_FirstUseEver); //How much of the menu content is displayed
		ImGui::SetNextWindowPos(ImVec2(2.0f, 0.0f), ImGuiSetCond_FirstUseEver);

		ImGui::Begin("Selection Menu", 0, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoScrollbar);
		int frame_padding = 5;
		ImGuiIO& io = ImGui::GetIO();
		ImGui::SetWindowFontScale(2.f);

		ImGui::PushID(0);
		if (ImGui::ImageButton(im_texID_draw, icon_size, uv_start, uv_end, frame_padding, icon_background_color)) {
			reset_trackers();
			selected_tool_mode = DRAW;
			im_texID_cur = im_texID_draw;
			viewer.selected_data_index = pointer_mesh_index;
			viewer.data().show_laser = false;
			viewer.oculusVR.left_hand_visible = true; 			//Turn left on
			prev_laser_show = viewer.data().show_laser;
			viewer.selected_data_index = base_mesh_index;
			menu_closed();
		}
		ImGui::PopID();
		ImGui::SameLine();
		ImGui::PushID(1);
		if (ImGui::ImageButton(im_texID_pull, icon_size, uv_start, uv_end, frame_padding, icon_background_color)) {
			reset_trackers();
			selected_tool_mode = PULL;
			im_texID_cur = im_texID_pull;
			viewer.selected_data_index = pointer_mesh_index;
			viewer.data().show_laser = false;
			viewer.oculusVR.left_hand_visible = true; 			//Turn left on
			prev_laser_show = viewer.data().show_laser;
			viewer.selected_data_index = base_mesh_index;
			menu_closed();
		}
		ImGui::PopID();
		ImGui::SameLine();
		ImGui::PushID(2);
		if (ImGui::ImageButton(im_texID_add, icon_size, uv_start, uv_end, frame_padding, icon_background_color)) {
			reset_trackers();
			selected_tool_mode = ADD;
			im_texID_cur = im_texID_add;
			viewer.selected_data_index = pointer_mesh_index;
			viewer.data().show_laser = false;
			viewer.oculusVR.left_hand_visible = true; 			//Turn left on
			prev_laser_show = viewer.data().show_laser;
			viewer.selected_data_index = base_mesh_index;
			menu_closed();
		}
		ImGui::PopID();

		ImGui::PushID(3);
		if (ImGui::ImageButton(im_texID_cut, icon_size, uv_start, uv_end, frame_padding, icon_background_color)) {
			reset_trackers();
			selected_tool_mode = CUT;
			im_texID_cur = im_texID_cut;
			viewer.selected_data_index = pointer_mesh_index;
			viewer.data().show_laser = true;
			viewer.oculusVR.left_hand_visible = false; 			//Turn left off
			prev_laser_show = viewer.data().show_laser;
			viewer.selected_data_index = base_mesh_index;
			menu_closed();
		}
		ImGui::PopID();
		ImGui::SameLine();
		ImGui::PushID(4);
		if (ImGui::ImageButton(im_texID_extrude, icon_size, uv_start, uv_end, frame_padding, icon_background_color)) {
			reset_trackers();
			selected_tool_mode = EXTRUDE;
			im_texID_cur = im_texID_extrude;
			viewer.selected_data_index = pointer_mesh_index;
			viewer.data().show_laser = true;
			viewer.oculusVR.left_hand_visible = false; 			//Turn left off
			prev_laser_show = viewer.data().show_laser;
			viewer.selected_data_index = base_mesh_index;
			menu_closed();
		}
		ImGui::PopID();
		ImGui::SameLine();
		ImGui::PushID(5);
		if (ImGui::ImageButton(im_texID_remove, icon_size, uv_start, uv_end, frame_padding, icon_background_color)) {
			reset_trackers();
			selected_tool_mode = REMOVE;
			im_texID_cur = im_texID_remove;
			viewer.selected_data_index = pointer_mesh_index;
			viewer.data().show_laser = false;
			viewer.oculusVR.left_hand_visible = false; 			//Turn left off
			prev_laser_show = viewer.data().show_laser;
			viewer.selected_data_index = base_mesh_index;
			menu_closed();
		}
		ImGui::PopID();

		ImGui::PushID(6);
		if (ImGui::ImageButton(im_texID_save, icon_size, uv_start, uv_end, frame_padding, icon_background_color)) {
			reset_trackers();
			appState.init(vertex_is_fixed, edge_boundary_markers, sharp_edge, new_mapped_indices, replacing_vertex_bindings, stroke_collection.size());
			selected_tool_mode = FAIL;
			std::cerr << "Please specify how you want to save the scene. " << std::endl;
			callback_save_scene(viewer, cur_dir);
			im_texID_cur = im_texID_save;
			viewer.selected_data_index = pointer_mesh_index;
			viewer.data().show_laser = true;
			prev_laser_show = viewer.data().show_laser;
			viewer.selected_data_index = base_mesh_index;
			menu_closed();
		}
		ImGui::PopID();
		ImGui::SameLine();
		ImGui::PushID(7);
		if (ImGui::ImageButton(im_texID_load, icon_size, uv_start, uv_end, frame_padding, icon_background_color)) {
			viewer.oculusVR.request_recenter();
			reset_trackers();
			selected_tool_mode = FAIL;
			std::cerr << "Please select the file you want to load a scene from. " << std::endl;
			bool success_load = callback_load_scene(viewer, cur_dir);
			dirty_boundary = true;
			im_texID_cur = im_texID_load;
			viewer.selected_data_index = pointer_mesh_index;
			viewer.data().show_laser = true;
			prev_laser_show = viewer.data().show_laser;
			viewer.selected_data_index = base_mesh_index;
			menu_closed();
		}
		ImGui::PopID();
		ImGui::SameLine();
		ImGui::PushID(8);
		if (ImGui::ImageButton(im_texID_change, icon_size, uv_start, uv_end, frame_padding, icon_background_color)) {
			reset_trackers();
			selected_tool_mode = CHANGE;
			im_texID_cur = im_texID_change;
			viewer.selected_data_index = pointer_mesh_index;
			viewer.data().show_laser = false;
			viewer.oculusVR.left_hand_visible = false; 			//Turn left off
			prev_laser_show = viewer.data().show_laser;
			viewer.selected_data_index = base_mesh_index;
			menu_closed();
		}
		ImGui::PopID();
		ImGui::End();
	};

	menu.callback_draw_viewer_window = current_tool_HUD;
	viewer.oculusVR.callback_button_down = button_down;
	viewer.oculusVR.callback_menu_opened = menu_opened;
	viewer.oculusVR.callback_menu_closed = menu_closed;
	viewer.data().show_lines = true;
	viewer.launch_oculus();
}

