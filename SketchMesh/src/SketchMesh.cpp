#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#else
#include <unistd.h>
#endif
#include <iostream>
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/per_face_normals.h>
#include <igl/cat.h>
#include <igl/edge_topology.h>
#include <igl/read_triangle_mesh.h>
#include "SketchMesh.h"
#include "Stroke.h"
#include "SurfaceSmoothing.h"
#include "CurveDeformation.h"
#include "MeshCut.h"
#include "MeshExtrusion.h"
#include "SurfacePath.h"


using namespace std;
using Viewer = igl::viewer::Viewer;

// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Per face normals, #F x3
Eigen::MatrixXd N_Faces;

//Per vertex indicator of whether vertex is on boundary (on boundary if == 1)
Eigen::VectorXi vertex_boundary_markers;
//Per vertex indicator of whether vertex is on original stroke (outline of shape) (on OG stroke if ==1)
Eigen::VectorXi part_of_original_stroke;
//Per edge indicator of whether the edge is sharp (if == 1 then sharp, otherwise smooth)
Eigen::VectorXi sharp_edge;
//Takes care of index mapping from before a cut/extrusion action to after (since some vertices are removed)
Eigen::VectorXi new_mapped_indices;

//General
enum ToolMode { DRAW, ADD, CUT, EXTRUDE, PULL, REMOVE, CHANGE, SMOOTH, NAVIGATE, NONE };
ToolMode tool_mode = NAVIGATE;
Stroke* initial_stroke;
Stroke* added_stroke;
Stroke* extrusion_base;
vector<Stroke> stroke_collection;

//Mouse interaction
bool skip_standardcallback = false;
int down_mouse_x = -1, down_mouse_y = -1;
bool mouse_is_down = false; //We need this due to mouse_down not working in the nanogui menu, whilst mouse_up does work there
bool mouse_has_moved = false;

//For smoothing
int initial_smooth_iter = 8;

//For selecting vertices
int handleID = -1;

//Variables for pulling a curve (and removing added control curves)
int turnNr = 0;
bool dirty_boundary = false;
int closest_stroke_ID, prev_closest_stroke_ID;

//Keeps track of the stroke IDs
int next_added_stroke_ID = 2; //Start at 2 because marker 1 belongs to the original boundary

//Variables for adding control curves
bool last_add_on_mesh = false;
unordered_map<int, int> backside_vertex_map;

//Variables for removing a control curve
bool stroke_was_removed = false;
int remove_stroke_clicked = 0;

//Variables for cutting
bool cut_stroke_already_drawn = false;

//Variables for extrusion
bool extrusion_base_already_drawn = false;
SurfacePath base_surface_path;
Eigen::Matrix4f base_model, base_view, base_proj;
Eigen::Vector4f base_viewport;

void draw_all_strokes(Viewer& viewer) {
	Eigen::MatrixXd added_points;
	viewer.data.set_points((Eigen::MatrixXd) initial_stroke->get3DPoints().block(0, 0, initial_stroke->get3DPoints().rows() - 1, 3), Eigen::RowVector3d(1, 0, 0)); //Display the original stroke points and clear all the rest. Don't take the last point
	viewer.data.set_edges(Eigen::MatrixXd(), Eigen::MatrixXi(), Eigen::RowVector3d(0, 0, 1)); //Clear the non-original stroke edges
	
	if(initial_stroke->is_loop) {
		viewer.data.set_stroke_points(igl::cat(1, (Eigen::MatrixXd) initial_stroke->get3DPoints().block(0, 0, initial_stroke->get3DPoints().rows() - 1, 3), (Eigen::MatrixXd) initial_stroke->get3DPoints().row(0))); //Create a loop and draw edges
	} else { //set_stroke_points always makes a loop, so don't use that when our stroke ain't a loop (anymore)
		added_points = initial_stroke->get3DPoints();
		viewer.data.add_edges(added_points.block(0, 0, added_points.rows() - 2, 3), added_points.block(1, 0, added_points.rows() - 2, 3), Eigen::RowVector3d(1, 0, 0));
	}

	int points_to_hold_back;
	for(int i = 0; i < stroke_collection.size(); i++) {
		added_points = stroke_collection[i].get3DPoints();
		points_to_hold_back = 1 + !stroke_collection[i].is_loop;
		viewer.data.add_points(added_points, stroke_collection[i].stroke_color);
		viewer.data.add_edges(added_points.block(0, 0, added_points.rows() - points_to_hold_back, 3), added_points.block(1, 0, added_points.rows() - points_to_hold_back, 3), stroke_collection[i].stroke_color);
	}
}

void draw_extrusion_base(Viewer& viewer) {
	int points_to_hold_back;
	Eigen::MatrixXd added_points = extrusion_base->get3DPoints();
	points_to_hold_back = 1 + extrusion_base->is_loop;
	viewer.data.add_points(added_points, extrusion_base->stroke_color);
	viewer.data.add_edges(added_points.block(0, 0, added_points.rows() - points_to_hold_back, 3), added_points.block(1, 0, added_points.rows() - points_to_hold_back, 3), extrusion_base->stroke_color);
}

void select_dragging_handle(int down_mouse_x, int down_mouse_y) {
	double closest_dist = INFINITY;
	handleID = initial_stroke->selectClosestVertex(down_mouse_x, down_mouse_y, closest_dist);
	double current_closest = closest_dist;
	closest_stroke_ID = -1;
	int tmp_handleID;
	for(int i = 0; i < stroke_collection.size(); i++) { //Additional strokes that cross the original stroke will never be selected as the pulled curve when the user clicks a vertex that also belongs to the original boundary, since their vertex positions are the same and we check for SMALLER distances
		tmp_handleID = stroke_collection[i].selectClosestVertex(down_mouse_x, down_mouse_y, closest_dist); //Returns the index into stroke3DPoints of the closest point
		if((closest_dist < current_closest) && (tmp_handleID != -1)) {
			current_closest = closest_dist;
			handleID = tmp_handleID;
			closest_stroke_ID = i;
		}
	}
}

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
	if(key == '1') {
		viewer.data.clear();
	} else if(key == 'D') { //use capital letters
		//Draw initial curve/mesh
		tool_mode = DRAW;
	} else if(key == 'P') {
		if(initial_stroke->empty2D()) { //Don't go into "pull mode" if there is no mesh yet
			return true;
		}
		tool_mode = PULL;
	} else if(key == 'N') {
		//Use navigation
		tool_mode = NAVIGATE;
	} else if(key == 'A') {
		//Add an extra control curve to an existing mesh
		if(initial_stroke->empty2D()) { //Don't go into "additional curve mode" if there is no mesh yet
			return true;
		}
		tool_mode = ADD;
	} else if(key == 'R') {
		if(stroke_collection.size() == 0) { //Don't go into "remove curve mode" if there is no additional curves
			return true;
		}
		remove_stroke_clicked = 0; //Reset because we might be left with a single click from the last round
		tool_mode = REMOVE;
	} else if(key == 'C') {
		if(initial_stroke->empty2D()) { //Don't go into "cut mode" if there is no mesh yet
			return true;
		}
		cut_stroke_already_drawn = false; //Reset because we might have stopped before finishing the cut last time
		tool_mode = CUT;
	} else if(key == 'E') {
		if(initial_stroke->empty2D()) { //Don't go into "extrude mode" if there is no mesh yet
			return true;
		}
		//We need to switch to NAVIGATE in order to draw the silhouette stroke, so we cannot reset extrusion_base_already_drawn here
		tool_mode = EXTRUDE;
	}


	return true;
}

bool callback_mouse_down(Viewer& viewer, int button, int modifier) {
	if(button == (int)Viewer::MouseButton::Right) {
		return false;
	}
	mouse_is_down = true;

	down_mouse_x = viewer.current_mouse_x;
	down_mouse_y = viewer.current_mouse_y;

	if(tool_mode == DRAW) { //Creating the first curve/mesh
		viewer.data.clear();
		stroke_collection.clear();
		next_added_stroke_ID = 2;
		initial_stroke->strokeReset();
		initial_stroke->strokeAddSegment(down_mouse_x, down_mouse_y);
		skip_standardcallback = true;
	} 
	else if(tool_mode == ADD) { //Adding a new control curve onto an existing mesh
		added_stroke = new Stroke(V, F, viewer, next_added_stroke_ID);
		next_added_stroke_ID++;
		added_stroke->strokeAddSegmentAdd(down_mouse_x, down_mouse_y); //If the user starts outside of the mesh, consider the movement as navigation
		skip_standardcallback = true;
	} 
	else if(tool_mode == REMOVE) {
		double closest_dist = INFINITY;
		double current_closest = closest_dist;
		int tmp_handleID, closest_stroke_idx;
		handleID = -1;
		for(int i = 0; i < stroke_collection.size(); i++) {
			tmp_handleID = stroke_collection[i].selectClosestVertex(down_mouse_x, down_mouse_y, closest_dist);
			if((closest_dist < current_closest) && (tmp_handleID != -1)) {
				current_closest = closest_dist;
				handleID = tmp_handleID;
				closest_stroke_ID = stroke_collection[i].get_ID();
				closest_stroke_idx = i;
			}
		}

		if(handleID == -1) {//User clicked too far from any of the stroke vertices
			return false;
		}
		if(closest_stroke_ID == prev_closest_stroke_ID) {
			remove_stroke_clicked++;
		} else {
			remove_stroke_clicked = 1; //Start from 1
			prev_closest_stroke_ID = closest_stroke_ID;
		}

		//Redraw the original stroke and all added strokes, where the selected stroke is drawn in black.
		Eigen::MatrixXd init_points = initial_stroke->get3DPoints();
		viewer.data.set_points(init_points.topRows(init_points.rows() - 1), Eigen::RowVector3d(1, 0, 0));
		Eigen::MatrixXd added_points = stroke_collection[closest_stroke_idx].get3DPoints();
		viewer.data.add_points(added_points.topRows(added_points.rows() - 1), Eigen::RowVector3d(0, 0, 0));
		for(int i = 0; i < stroke_collection.size(); i++) {
			if(stroke_collection[i].get_ID() == closest_stroke_ID) {
				continue;
			}
			added_points = stroke_collection[i].get3DPoints();
			viewer.data.add_points(added_points.topRows(added_points.rows() - 1), stroke_collection[i].stroke_color);
		}


		if(remove_stroke_clicked == 2) { //Mechanism to force the user to click twice on the same stroke before removing it
			stroke_was_removed = true;
			stroke_collection[closest_stroke_idx].undo_stroke_add(vertex_boundary_markers); //Sets the vertex_boundary_markers for the vertices of this stroke to 0 again
			stroke_collection.erase(stroke_collection.begin() + closest_stroke_idx);
			remove_stroke_clicked = 0;
		}

		skip_standardcallback = true;
	} 
	else if(tool_mode == PULL) { //Dragging an existing curve
		select_dragging_handle(down_mouse_x, down_mouse_y);

		if(handleID == -1) {//User clicked too far from any of the stroke vertices
			return false;
		}
		if(closest_stroke_ID == -1) {
			CurveDeformation::startPullCurve(*initial_stroke, handleID);
		} else {
			CurveDeformation::startPullCurve(stroke_collection[closest_stroke_ID], handleID);
		}
		skip_standardcallback = true;
	} 
	else if(tool_mode == NAVIGATE) { //Navigate through the screen
		skip_standardcallback = false; //We do want to use the navigation functionality
	} 
	else if(tool_mode == CUT) {
		if(cut_stroke_already_drawn) { //clicked while cut stroke already drawn
			return true;
		}

		added_stroke = new Stroke(V, F, viewer, next_added_stroke_ID);
		next_added_stroke_ID++;
		added_stroke->strokeAddSegmentCut(down_mouse_x, down_mouse_y);
		skip_standardcallback = true;
	} 
	else if(tool_mode == EXTRUDE) {
		if(extrusion_base_already_drawn) { //clicked while the extrude base was already drawn
			added_stroke = new Stroke(V, F, viewer, next_added_stroke_ID);
			next_added_stroke_ID++;
			added_stroke->strokeAddSegmentExtrusionSilhouette(down_mouse_x, down_mouse_y);
		} else { //clicked with no extrude base yet
			extrusion_base = new Stroke(V, F, viewer, next_added_stroke_ID);
			next_added_stroke_ID++;
			extrusion_base->strokeAddSegmentExtrusionBase(down_mouse_x, down_mouse_y);
		}
		skip_standardcallback = true;
	}

	return skip_standardcallback; //Will make sure that we use standard navigation responses if we didn't do special actions and vice versa
}

bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y) {
	if(!skip_standardcallback) {
		return false;
	}

	if(viewer.down) { //Only consider it to be moving if the button was held down
		mouse_has_moved = true;
	}


	if(tool_mode == DRAW && viewer.down) { //If we're still holding the mouse down
		initial_stroke->strokeAddSegment(mouse_x, mouse_y);
		return true;
	} 
	else if(tool_mode == ADD && viewer.down) {
		last_add_on_mesh = added_stroke->strokeAddSegmentAdd(mouse_x, mouse_y);
		return true;
	} 
	else if(tool_mode == EXTRUDE && viewer.down) {
		if(extrusion_base_already_drawn) {
			added_stroke->strokeAddSegmentExtrusionSilhouette(mouse_x, mouse_y);
		} else {
			extrusion_base->strokeAddSegmentExtrusionBase(mouse_x, mouse_y);
		}
		return true;
	} 
	else if(tool_mode == PULL && viewer.down && handleID != -1) {
		double x = mouse_x;
		double y = viewer.core.viewport(3) - mouse_y;

		Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
		int global_handleID;
		if(closest_stroke_ID == -1) {
			global_handleID = initial_stroke->get_vertex_idx_for_point(handleID);
		} else {
			global_handleID = stroke_collection[closest_stroke_ID].get_vertex_idx_for_point(handleID);
		}
		Eigen::RowVector3f pt1(viewer.data.V(global_handleID, 0), viewer.data.V(global_handleID, 1), viewer.data.V(global_handleID, 2));
		Eigen::RowVector3f pr;
		igl::project(pt1, modelview, viewer.core.proj, viewer.core.viewport, pr);
		Eigen::RowVector3d pt = igl::unproject(Eigen::Vector3f(x, y, pr[2]), modelview, viewer.core.proj, viewer.core.viewport).transpose().cast<double>();

		if(turnNr == 0) { //increase the number to smooth less often
			CurveDeformation::pullCurve(pt, V, part_of_original_stroke);
			if(dirty_boundary) { //Smooth an extra time if the boundary is dirty, because smoothing once with a dirty boundary results in a flat mesh
				for(int i = 0; i < 2; i++) {
						SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
				}
			}
				SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);

			turnNr++;
		} else {
			turnNr++;
			if(turnNr == 4) {
				turnNr = 0;
			}
		}

		initial_stroke->update_Positions(V);
		for(int i = 0; i < stroke_collection.size(); i++) {
			stroke_collection[i].update_Positions(V);
		}

		viewer.data.set_mesh(V, F);
		viewer.data.compute_normals();
		draw_all_strokes(viewer);

		return true;
	} 
	else if(tool_mode == CUT && viewer.down) {
		if(!cut_stroke_already_drawn) {
			added_stroke->strokeAddSegmentCut(mouse_x, mouse_y);
		}
		return true;
	}

	return false;
}

bool callback_mouse_up(Viewer& viewer, int button, int modifier) {
	if(!mouse_is_down) {
		return true;
	}
	mouse_is_down = false;

	if(tool_mode == DRAW) {
		if(initial_stroke->toLoop()) {//Returns false if the stroke only consists of 1 point (user just clicked)
			//Give some time to show the stroke
#ifdef _WIN32
			Sleep(200);
#else
			usleep(200000);  /* sleep for 200 milliSeconds */
#endif
			backside_vertex_map = initial_stroke->generate3DMeshFromStroke(vertex_boundary_markers, part_of_original_stroke);
			F = viewer.data.F;
			V = viewer.data.V;

			if (!igl::is_edge_manifold(F)) { //Check if the drawn stroke results in an edge-manifold mesh, otherwise sound a beep and revert
#ifdef _WIN32
				Beep(500, 200);
#else
				cout << "error" << endl;
#endif				
				Eigen::MatrixXd drawn_points = initial_stroke->get3DPoints();
				initial_stroke->strokeReset();
				vertex_boundary_markers.resize(0);
				part_of_original_stroke.resize(0);
				dirty_boundary = true;

				viewer.data.clear();
				viewer.data.add_edges(drawn_points.block(0, 0, drawn_points.rows() - 1, 3), drawn_points.block(1, 0, drawn_points.rows() - 1, 3), Eigen::RowVector3d(0, 0, 0)); //Display the stroke in black to show that it went wrong
				return true;
			}

			Eigen::MatrixXi EV, FE, EF;
			igl::edge_topology(V, F, EV, FE, EF);
			sharp_edge.resize(EV.rows());
			sharp_edge.setZero(); //Set all edges to smooth after initial draw

			dirty_boundary = true;

			for(int i = 0; i < initial_smooth_iter; i++) {
				SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
			}

			viewer.data.set_mesh(V, F);
			viewer.data.compute_normals();

			//Overlay the drawn stroke
			int strokeSize = (vertex_boundary_markers.array() > 0).count();
			Eigen::MatrixXd strokePoints = V.block(0, 0, strokeSize, 3);
			viewer.data.set_points(strokePoints, Eigen::RowVector3d(1, 0, 0));
			viewer.data.set_stroke_points(igl::cat(1, strokePoints, (Eigen::MatrixXd) V.row(0)));
		}
		skip_standardcallback = false;
	} 
	else if(tool_mode == ADD) {
		dirty_boundary = true;
		if(!added_stroke->has_points_on_mesh) {
			mouse_has_moved = false;
			return true;
		}
		added_stroke->snap_to_vertices(vertex_boundary_markers);
		stroke_collection.push_back(*added_stroke);
		draw_all_strokes(viewer);
	} 
	else if(tool_mode == REMOVE && stroke_was_removed) { //Only redraw if we actually removed a stroke (otherwise we draw unnecessary)
		stroke_was_removed = false;
		dirty_boundary = true;

		draw_all_strokes(viewer);
	} 
	else if(tool_mode == PULL && handleID != -1 && mouse_has_moved) {
		for(int i = 0; i < 2; i++) {
			     SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
		}

		for(int i = 0; i < stroke_collection.size(); i++) {
			stroke_collection[i].update_Positions(V);
		}

		viewer.data.set_mesh(V, F);
		viewer.data.compute_normals(); 
		draw_all_strokes(viewer);
	} 
	else if(tool_mode == CUT) {
		if(!added_stroke->has_points_on_mesh) {
			mouse_has_moved = false;
			return true;
		}
		if(cut_stroke_already_drawn) { //User had already drawn the cut stroke and has now drawn the final stroke for removing the part
			dirty_boundary = true;
			added_stroke->append_final_point();
			added_stroke->toLoop();
			bool cut_success = MeshCut::cut(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, *added_stroke);
			if (!cut_success) { //Catches the cases that the cut removes all interior mesh vertices/faces and when the SurfacePath creation is stuck in an infinite loop
				cut_stroke_already_drawn = false;
				next_added_stroke_ID--; //Undo ID increment since stroke didn't actually get pushed back
#ifdef _WIN32
				Beep(500, 200);
#else
				cout << "error" << endl;
#endif	

				draw_all_strokes(viewer); //Will remove the pink cut stroke
				Eigen::MatrixXd drawn_points = added_stroke->get3DPoints();
				viewer.data.add_edges(drawn_points.block(0, 0, drawn_points.rows() - 1, 3), drawn_points.block(1, 0, drawn_points.rows() - 1, 3), Eigen::RowVector3d(0, 0, 0)); //Display the stroke in black to show that it went wrong
				return true;
			}
			stroke_collection.push_back(*added_stroke);

			initial_stroke->update_vert_bindings(new_mapped_indices, vertex_boundary_markers); //Don't test if the initial one dies
			
			int nr_removed = 0, original_collection_size = stroke_collection.size();
			for(int i = 0; i < original_collection_size; i++) {
				if(!stroke_collection[i - nr_removed].update_vert_bindings(new_mapped_indices, vertex_boundary_markers)) {
					//Stroke dies, don't need to do stroke.undo_stroke_add, cause all its vertices also cease to exist
					stroke_collection.erase(stroke_collection.begin() + i - nr_removed);
					nr_removed++;
					dirty_boundary = true;
					continue; //Go to the next stroke, don't update this ones' positions
				}
			}

			for(int i = 0; i < 8; i++) {
				SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
			}

			//Update the stroke positions after smoothing, in case their positions have changed (although they really shouldn't)
			initial_stroke->update_Positions(V);
			for(int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V);
			}

			viewer.data.clear();
			viewer.data.set_mesh(V, F);
			igl::per_face_normals(V, F, N_Faces);
			viewer.data.set_normals(N_Faces);
			viewer.core.align_camera_center(V, F);

			cut_stroke_already_drawn = false;
			draw_all_strokes(viewer);
		} else { //We're finished drawing the cut stroke, prepare for when user draws the final stroke to remove the part
			cut_stroke_already_drawn = true;
		}
	} 
	else if(tool_mode == EXTRUDE) {
		if(extrusion_base_already_drawn) { //User has drawn the silhouette stroke for extrusion
			dirty_boundary = true;
			added_stroke->toLoop();
			added_stroke->is_loop = false; //Set to false manually, because we don't want curveDeformation to consider it as a loop (but we do need looped 3DPoints)

			bool success_extrude = MeshExtrusion::extrude_main(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, base_surface_path, *added_stroke, *extrusion_base, base_model, base_view, base_proj, base_viewport);
			if (!success_extrude) { //Catches the case that the extrusion base removes all faces/vertices
				next_added_stroke_ID -= 2;
				extrusion_base_already_drawn = false;
#ifdef _WIN32
				Beep(700, 200);
#else
				cout << "error" << endl;
#endif	
				draw_all_strokes(viewer); //Will remove the drawn base & silhouette strokes
				Eigen::MatrixXd drawn_points = extrusion_base->get3DPoints();
				viewer.data.add_edges(drawn_points.block(0, 0, drawn_points.rows() - 1, 3), drawn_points.block(1, 0, drawn_points.rows() - 1, 3), Eigen::RowVector3d(0, 0, 0)); //Display the stroke in black to show that it went wrong
				drawn_points = added_stroke->get3DPoints();
				viewer.data.add_edges(drawn_points.block(0, 0, drawn_points.rows() - 2, 3), drawn_points.block(1, 0, drawn_points.rows() - 2, 3), Eigen::RowVector3d(0, 0, 0)); //Display the stroke in black to show that it went wrong
				return true;
			}
			stroke_collection.push_back(*extrusion_base);
			stroke_collection.push_back(*added_stroke);

			initial_stroke->update_vert_bindings(new_mapped_indices, vertex_boundary_markers);

			int nr_removed = 0, original_collection_size = stroke_collection.size();
			for(int i = 0; i < original_collection_size - 1; i++) { //Skip the newly added stroke since it is already updated inside extrude_main()
				if(!stroke_collection[i - nr_removed].update_vert_bindings(new_mapped_indices, vertex_boundary_markers)) {
					//Stroke dies, don't need to do stroke.undo_stroke_add, cause all its vertices also cease to exist (or in case of non-loop strokes that have a middle portion removed, the undo_stroke_add is done inside of the update_vert_bindings
					stroke_collection.erase(stroke_collection.begin() + i - nr_removed);
					nr_removed++;
					dirty_boundary = true;
					continue; //Go to the next stroke, don't update this ones' positions
				}
			}

			for(int i = 0; i < 3; i++) {
				SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
			}

			//Update the stroke positions after smoothing, in case their positions have changed (although they really shouldn't)
			initial_stroke->update_Positions(V);
			for(int i = 0; i < stroke_collection.size()-1; i++) {
				stroke_collection[i].update_Positions(V);
			}



			viewer.data.clear();
			viewer.data.set_mesh(V, F);
			igl::per_face_normals(V, F, N_Faces);
			viewer.data.set_normals(N_Faces);
			viewer.core.align_camera_center(V, F);

			extrusion_base_already_drawn = false;
			draw_all_strokes(viewer);
		} else { //mouse released after extrusion base drawn
			if(!extrusion_base->has_points_on_mesh) {
				mouse_has_moved = false;
				return true;
			}
			
			dirty_boundary = true;
			extrusion_base->toLoop();
			extrusion_base_already_drawn = true;
			bool success_extrude_prepare = MeshExtrusion::extrude_prepare(*extrusion_base, base_surface_path); //Don't need to update all strokes here, since it didn't remove any vertices
			if (!success_extrude_prepare) { //Catches the case that face == -1 in SurfacePath
#ifdef _WIN32
				Beep(900, 200);
#else
				cout << "error" << endl;
#endif	
				next_added_stroke_ID--;
				draw_all_strokes(viewer); //Removes the drawn base stroke
				Eigen::MatrixXd drawn_points = extrusion_base->get3DPoints();
				viewer.data.add_edges(drawn_points.block(0, 0, drawn_points.rows() - 1, 3), drawn_points.block(1, 0, drawn_points.rows() - 1, 3), Eigen::RowVector3d(0, 0, 0)); //Display the stroke in black to show that it went wrong
				extrusion_base_already_drawn = false;
				return true;
			}

			base_model = extrusion_base->viewer.core.model;
			base_view = extrusion_base->viewer.core.view;
			base_proj = extrusion_base->viewer.core.proj;
			base_viewport = extrusion_base->viewer.core.viewport;

			draw_all_strokes(viewer);
			draw_extrusion_base(viewer); //Need to draw the extrusion base separately, since it isn't added to the stroke_collection yet.
		}
	}

	mouse_has_moved = false;
	return skip_standardcallback;
}

bool callback_load_mesh(Viewer& viewer, string filename) {
	igl::readOFF(filename, V, F);
	viewer.data.clear();
	viewer.data.set_mesh(V, F);
	viewer.data.compute_normals();
	viewer.core.align_camera_center(viewer.data.V);

	std::cout << filename.substr(filename.find_last_of("/") + 1) << endl;
	return true;
}

int main(int argc, char *argv[]) {
	// Show the mesh
	Viewer viewer;
	viewer.callback_key_down = callback_key_down;
	viewer.callback_mouse_down = callback_mouse_down;
	viewer.callback_mouse_move = callback_mouse_move;
	viewer.callback_mouse_up = callback_mouse_up;
	viewer.core.point_size = 15;
	viewer.data.set_face_based(true);
	//viewer.callback_load_mesh = callback_load_mesh;

	viewer.callback_init = [&](igl::viewer::Viewer& viewer) {
		// Add new group
		viewer.ngui->addGroup("Inflation");

		// Expose a variable directly ...
		viewer.ngui->addVariable("Vertex Weights", SurfaceSmoothing::vertex_weight);
		viewer.ngui->addVariable("Edge Weights", SurfaceSmoothing::edge_weight);
		viewer.ngui->addVariable("factor", SurfaceSmoothing::factor);


		// Expose a variable directly ...
		viewer.ngui->addVariable("Initial smoothing iterations", initial_smooth_iter);


		// Add a button
		viewer.ngui->addButton("Perform 1 smoothing iteration", [&viewer]() {
			SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
			viewer.data.set_mesh(V, F);
			viewer.data.compute_normals();


			for(int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V);
			}

			//Overlay the updated stroke
			initial_stroke->update_Positions(V);
			Eigen::MatrixXd init_points = initial_stroke->get3DPoints();
			viewer.data.set_points(init_points.topRows(init_points.rows() - 1), Eigen::RowVector3d(1, 0, 0));
			viewer.data.set_stroke_points(init_points);

			viewer.data.set_edges(Eigen::MatrixXd(), Eigen::MatrixXi(), Eigen::RowVector3d(0, 0, 1));
			Eigen::MatrixXd added_points;
			int points_to_hold_back;
			for(int i = 0; i < stroke_collection.size(); i++) {
				added_points = stroke_collection[i].get3DPoints();
				points_to_hold_back = 1 + !stroke_collection[i].is_loop;
				viewer.data.add_points(added_points.topRows(added_points.rows() - 1), stroke_collection[i].stroke_color);
				viewer.data.add_edges(added_points.block(0, 0, added_points.rows() - points_to_hold_back, 3), added_points.block(1, 0, added_points.rows() - points_to_hold_back, 3), stroke_collection[i].stroke_color);
			}

		});

		viewer.ngui->addGroup("Curve deformation");
		viewer.ngui->addVariable<bool>("Smooth deformation", CurveDeformation::smooth_deform_mode);

		// call to generate menu
		viewer.screen->performLayout();
		return false;
	};

	//Init stroke selector
	initial_stroke = new Stroke(V, F, viewer, 0);
	if(argc == 2) {
		igl::readOFF(argv[1], V, F);
		//	callback_load_mesh(viewer, argv[1]);
	} else {
		//callback_load_mesh(viewer, "../data/cube.off");
	}

	callback_key_down(viewer, '1', 0);
	viewer.launch();
}


