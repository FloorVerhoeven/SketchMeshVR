#include "Stroke.h"
#include <igl/unproject_onto_mesh.h>
#include <igl/unproject.h>
#include <igl/triangle/triangulate.h>
#include <algorithm> 
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/dijkstra.h>
#include "Plane.h"
using namespace igl;
using namespace std;

std::chrono::steady_clock::time_point _time2, _time1;

Stroke::Stroke(const Eigen::MatrixXd &V_, const Eigen::MatrixXi &F_, igl::viewer::Viewer &v, int stroke_ID_) :
	V(V_),
	F(F_),
	viewer(v),
	stroke_ID(stroke_ID_) {
	stroke2DPoints = Eigen::MatrixXd::Zero(1, 2);
	stroke3DPoints = Eigen::MatrixX3d::Zero(1, 3);
	_time1 = std::chrono::high_resolution_clock::now();
	closest_vert_bindings.clear();
	has_points_on_mesh = false;
	stroke_color = Eigen::RowVector3d(0.8*(rand() / (double)RAND_MAX), 0.8*(rand() / (double)RAND_MAX), 0.8*(rand() / (double)RAND_MAX));
	has_been_reversed = false;
}

Stroke::Stroke(const Stroke& origin) :
	V(origin.V),
	F(origin.F),
	viewer(origin.viewer),
	stroke_ID(origin.stroke_ID),
	stroke2DPoints(origin.stroke2DPoints),
	stroke3DPoints(origin.stroke3DPoints),
	closest_vert_bindings(origin.closest_vert_bindings),
	has_points_on_mesh(origin.has_points_on_mesh),
	has_been_reversed(origin.has_been_reversed),
	stroke_color(origin.stroke_color),
	dep(origin.dep),
	is_loop(origin.is_loop) {
	_time1 = std::chrono::high_resolution_clock::now();
}

Stroke& Stroke::operator=(Stroke other) {
	swap(other);
	return *this;
}

void Stroke::swap(Stroke & tmp) {//The pointers to V and F will always be the same for all stroke instances, so no need to copy
	std::swap(this->viewer, tmp.viewer);
	std::swap(this->stroke_ID, tmp.stroke_ID);
	std::swap(this->stroke2DPoints, tmp.stroke2DPoints);
	std::swap(this->stroke3DPoints, tmp.stroke3DPoints);
	std::swap(this->closest_vert_bindings, tmp.closest_vert_bindings);
	std::swap(this->has_points_on_mesh, tmp.has_points_on_mesh);
	std::swap(this->has_been_reversed, tmp.has_been_reversed);
	std::swap(this->stroke_color, tmp.stroke_color);
	std::swap(this->dep, tmp.dep);
	std::swap(this->is_loop, tmp.is_loop);
}

Stroke::~Stroke() {}

/** Used for DRAW. Will add a new 2D point to the stroke (if it is new compared to the last point, and didn't follow up too soon) and will also add its unprojection as a 3D point with a z-value of 0. After adding the new point it will restart the timer. **/
void Stroke::strokeAddSegment(int mouse_x, int mouse_y) {
	//OpenGL has origin at left bottom, window(s) has origin at left top
	double x = mouse_x;
	double y = viewer.core.viewport(3) - mouse_y;
	if(!empty2D() && x == stroke2DPoints(stroke2DPoints.rows() - 1, 0) && y == stroke2DPoints(stroke2DPoints.rows() - 1, 1)) { //Check that the point is new compared to last time
		return;
	}

	if(!empty2D()) {
		_time2 = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_time2 - _time1).count();
		if(timePast < 20000000) {
			return;
		}
	}

	Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
	Eigen::RowVector3d pt(0, 0, 0);

	pt = igl::unproject(Eigen::Vector3f(x, y, 0.0f), modelview, viewer.core.proj, viewer.core.viewport).transpose().cast<double>();
	if(stroke2DPoints.rows() == 1 && empty2D()) {
		stroke2DPoints.row(0) << x, y;
		stroke3DPoints.row(0) << pt[0], pt[1], pt[2];
	} else {
		stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, Eigen::NoChange);
		stroke2DPoints.row(stroke2DPoints.rows() - 1) << x, y;

		stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
		stroke3DPoints.row(stroke3DPoints.rows() - 1) << pt[0], pt[1], pt[2];
	}
	closest_vert_bindings.push_back(stroke3DPoints.rows() - 1); //In the case of DRAW this will match the vertex indices, since we start from 0
	viewer.data.set_stroke_points(stroke3DPoints); //Will remove all previous points but is okay for draw since it's the only points there are

	_time1 = std::chrono::high_resolution_clock::now();
}

/** Used for ADD. Will add a new 2D point to the stroke (if it is new compared to the last point, and didn't follow up too soon) and will also add its unprojection onto the existing mesh as a 3D point. After adding the new point it will restart the timer. **/
bool Stroke::strokeAddSegmentAdd(int mouse_x, int mouse_y) {
	bool result = false;
	//OpenGL has origin at left bottom, window(s) has origin at left top
	double x = mouse_x;
	double y = viewer.core.viewport(3) - mouse_y;
	if(!empty2D() && x == stroke2DPoints(stroke2DPoints.rows() - 1, 0) && y == stroke2DPoints(stroke2DPoints.rows() - 1, 1)) { //Check that the point is new compared to last time
		return result;
	}

	if(!empty2D()) {
		_time2 = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_time2 - _time1).count();
		if(timePast < 10000000) {
			return result;
		}
	}

	Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
	Eigen::RowVector3d pt(0, 0, 0);
	int faceID = -1;

	Eigen::Vector3f bc;
	if(igl::unproject_onto_mesh(Eigen::Vector2f(x, y), modelview, viewer.core.proj, viewer.core.viewport, V, F, faceID, bc)) {
		pt = V.row(F(faceID, 0))*bc(0) + V.row(F(faceID, 1))*bc(1) + V.row(F(faceID, 2))*bc(2);
		double cur_dist, min_dist = INFINITY;
		int closest_vert_idx;
		for(int i = 0; i < 3; i++) {
			cur_dist = (pt - V.row(F(faceID, i))).norm();
			if(cur_dist < min_dist) {
				min_dist = cur_dist;
				closest_vert_idx = F(faceID, i);
			}
		}

		if(closest_vert_bindings.size() == 0 || closest_vert_idx != closest_vert_bindings.back()) {
			closest_vert_bindings.push_back(closest_vert_idx);
		}

		has_points_on_mesh = true;

		if(stroke2DPoints.rows() == 1 && empty2D()) {
            stroke2DPoints.row(0) << x, y;
			stroke3DPoints.row(0) << pt[0], pt[1], pt[2];
		} else {
			stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, Eigen::NoChange);
			stroke2DPoints.row(stroke2DPoints.rows() - 1) << x, y;

			stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
			stroke3DPoints.row(stroke3DPoints.rows() - 1) << pt[0], pt[1], pt[2];
		}

		if(stroke3DPoints.rows() > 1) {
			viewer.data.add_edges(stroke3DPoints.block(stroke3DPoints.rows() - 2, 0, 1, 3), stroke3DPoints.block(stroke3DPoints.rows() - 1, 0, 1, 3), Eigen::RowVector3d(0, 1, 0));
		}
		result = true;
	}

	_time1 = std::chrono::high_resolution_clock::now();
	return result;
}

/** Used for CUT. Will add a new 2D point to the stroke (if it is new compared to the last point, and didn't follow up too soon) and will also add its unprojection onto the existing as a 3D point. If the unprojection isn't on the mesh, it will unproject it and use it as the start or end point. After adding the new point it will restart the timer. **/
void Stroke::strokeAddSegmentCut(int mouse_x, int mouse_y) {
	//OpenGL has origin at left bottom, window(s) has origin at left top
	double x = mouse_x;
	double y = viewer.core.viewport(3) - mouse_y;
	if(!empty2D() && x == stroke2DPoints(stroke2DPoints.rows() - 1, 0) && y == stroke2DPoints(stroke2DPoints.rows() - 1, 1)) { //Check that the point is new compared to last time
		return;
	}

	if(!empty2D()) {
		_time2 = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_time2 - _time1).count();
		if(timePast < 10000000) {
			return;
		}
	}

	Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
	Eigen::RowVector3d pt(0, 0, 0);
	int faceID = -1;

	Eigen::Vector3f bc;
	if(igl::unproject_onto_mesh(Eigen::Vector2f(x, y), modelview, viewer.core.proj, viewer.core.viewport, V, F, faceID, bc)) {
		pt = V.row(F(faceID, 0))*bc(0) + V.row(F(faceID, 1))*bc(1) + V.row(F(faceID, 2))*bc(2);
		has_points_on_mesh = true;

		if(stroke2DPoints.rows() == 1 && empty2D()) {
			cout << "This shouldn't happen. Draw the first point outside of the mesh" << endl;
			return;
		} else {
			stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, Eigen::NoChange);
			stroke2DPoints.row(stroke2DPoints.rows() - 1) << x, y;

			stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
			stroke3DPoints.row(stroke3DPoints.rows() - 1) << pt[0], pt[1], pt[2];
		}

		viewer.data.add_points(stroke3DPoints.row(stroke3DPoints.rows() - 1), Eigen::RowVector3d(1, 0, 1));
		if(stroke3DPoints.rows() > 2) {
			viewer.data.add_edges(stroke3DPoints.block(stroke3DPoints.rows() - 2, 0, 1, 3), stroke3DPoints.block(stroke3DPoints.rows() - 1, 0, 1, 3), Eigen::RowVector3d(1, 0, 1));
		}

		just_came_from_mesh = true;
	} else if(stroke2DPoints.rows() > 1 && just_came_from_mesh) { //We need to add the final point of the stroke, even though it is off the mesh (needed in order to wrap around to the backside)
		cut_stroke_final_point = igl::unproject(Eigen::Vector3f(x, y, 0.95*dep), modelview, viewer.core.proj, viewer.core.viewport).transpose().cast<double>();
		cut_stroke_final_point_2D = Eigen::RowVectorXd::Zero(1, 2);
		cut_stroke_final_point_2D << x, y;
		just_came_from_mesh = false;
	} else if(stroke2DPoints.rows() == 1) { //Add the first point, which should be outside the mesh. Refresh the "first" point until we get to the last one before we enter the mesh
		Eigen::RowVector3d tmp = igl::unproject(Eigen::Vector3f(x, y, 0.95*dep), modelview, viewer.core.proj, viewer.core.viewport).cast<double>();
		stroke2DPoints.row(0) << x, y;
		stroke3DPoints.row(0) = tmp;
		stroke_color = Eigen::RowVector3d(1, 0, 0);
	}

	_time1 = std::chrono::high_resolution_clock::now();
	return;
}

/** Used for EXTRUDE. Extrusion base strokes need to be drawn entirely on the mesh (points outside of it will be ignored) and needs to surround at least one whole triangle. Will add a new 2D point to the stroke (if it is new compared to the last point, and didn't follow up too soon) and will also add its unprojection onto the existing mesh as a 3D point. After adding the new point it will restart the timer. The closest vertex bindings are handled in SurfacePath. **/
void Stroke::strokeAddSegmentExtrusionBase(int mouse_x, int mouse_y) {
	//OpenGL has origin at left bottom, window(s) has origin at left top
	double x = mouse_x;
	double y = viewer.core.viewport(3) - mouse_y;
	if(!empty2D() && x == stroke2DPoints(stroke2DPoints.rows() - 1, 0) && y == stroke2DPoints(stroke2DPoints.rows() - 1, 1)) { //Check that the point is new compared to last time
		return;
	}

	if(!empty2D()) {
		_time2 = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_time2 - _time1).count();
		if(timePast < 50000000) {
            return;
		}
	}

	Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
	Eigen::RowVector3d pt(0, 0, 0);
	int faceID = -1;

	Eigen::Vector3f bc;
	if(igl::unproject_onto_mesh(Eigen::Vector2f(x, y), modelview, viewer.core.proj, viewer.core.viewport, V, F, faceID, bc)) {
		pt = V.row(F(faceID, 0))*bc(0) + V.row(F(faceID, 1))*bc(1) + V.row(F(faceID, 2))*bc(2);

		has_points_on_mesh = true; 

		if(stroke2DPoints.rows() == 1 && empty2D()) {
			stroke2DPoints.row(0) << x, y;
			stroke3DPoints.row(0) << pt[0], pt[1], pt[2];
		} else {
			stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, Eigen::NoChange);
			stroke2DPoints.row(stroke2DPoints.rows() - 1) << x, y;

			stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
			stroke3DPoints.row(stroke3DPoints.rows() - 1) << pt[0], pt[1], pt[2];
		}
	} else {
		return;
	}

	if(stroke3DPoints.rows() > 1) {
		viewer.data.add_edges(stroke3DPoints.block(stroke3DPoints.rows() - 2, 0, 1, 3), stroke3DPoints.block(stroke3DPoints.rows() - 1, 0, 1, 3), Eigen::RowVector3d(0, 0, 1));
	}

	_time1 = std::chrono::high_resolution_clock::now();
}

/** Used for EXTRUDE. Extrusion silhouette strokes should be drawn in any of the side views. Will add a new 2D point to the stroke (if it is new compared to the last point, and didn't follow up too soon) and will also add its unprojection as a 3D point. These need to be projected to a different plane than the 2D points are in, before meshing (done in MeshExtrusion). After adding the new point it will restart the timer. The closest vertex bindings are handled in SurfacePath. **/
void Stroke::strokeAddSegmentExtrusionSilhouette(int mouse_x, int mouse_y) {
	//OpenGL has origin at left bottom, window(s) has origin at left top
	double x = mouse_x;
	double y = viewer.core.viewport(3) - mouse_y;
	if(!empty2D() && x == stroke2DPoints(stroke2DPoints.rows() - 1, 0) && y == stroke2DPoints(stroke2DPoints.rows() - 1, 1)) { //Check that the point is new compared to last time
		return;
	}

	if(!empty2D()) {
		_time2 = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_time2 - _time1).count();
		if(timePast < 50000000) {
            return;
		}
	}

	Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
	Eigen::RowVector3d pt(0, 0, 0);

	pt = igl::unproject(Eigen::Vector3f(x, y, 0.8f), modelview, viewer.core.proj, viewer.core.viewport).transpose().cast<double>();
	if(stroke2DPoints.rows() == 1 && empty2D()) {
		stroke2DPoints.row(0) << x, y;
		stroke3DPoints.row(0) << pt[0], pt[1], pt[2];
	} else {
		stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, Eigen::NoChange);
		stroke2DPoints.row(stroke2DPoints.rows() - 1) << x, y;

		stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
		stroke3DPoints.row(stroke3DPoints.rows() - 1) << pt[0], pt[1], pt[2];
	}

	if(stroke3DPoints.rows() > 1) {
		viewer.data.add_edges(stroke3DPoints.block(stroke3DPoints.rows() - 2, 0, 1, 3), stroke3DPoints.block(stroke3DPoints.rows() - 1, 0, 1, 3), Eigen::RowVector3d(0.2, 0.6, 1));
	}

	_time1 = std::chrono::high_resolution_clock::now();
}

/** Used for CUT. Adds the final point, which is the first point outside of the mesh. Doesn't get drawn. **/
void Stroke::append_final_point() {
	stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, Eigen::NoChange);
	stroke2DPoints.row(stroke2DPoints.rows() - 1) << cut_stroke_final_point_2D[0], cut_stroke_final_point_2D[1];

	stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
	stroke3DPoints.row(stroke3DPoints.rows() - 1) << cut_stroke_final_point[0], cut_stroke_final_point[1], cut_stroke_final_point[2];
}

void Stroke::strokeReset() {
	stroke2DPoints.resize(1, 2);
	stroke2DPoints.setZero();
	stroke3DPoints.resize(1, 3);
	stroke3DPoints.setZero();
	dep = -1;
	_time1 = std::chrono::high_resolution_clock::now();
	closest_vert_bindings.clear();
	has_been_reversed = false;
	is_loop = false;
}

bool Stroke::toLoop() {
	if(stroke2DPoints.rows() > 2) { //Don't do anything if we have only 1 line segment
		stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
		stroke3DPoints.row(stroke3DPoints.rows() - 1) << stroke3DPoints.row(0);
		if(closest_vert_bindings.size() > 0) {
			closest_vert_bindings.push_back(closest_vert_bindings[0]);
		}
		is_loop = true;
		return true;
	}
	is_loop = false;
	return false;
}

unordered_map<int, int> Stroke::generate3DMeshFromStroke(Eigen::VectorXi &vertex_boundary_markers, Eigen::VectorXi &part_of_original_stroke) {
	counter_clockwise(); //Ensure the stroke is counter-clockwise, handy later

	Eigen::MatrixXd original_stroke2DPoints = stroke2DPoints;
	stroke2DPoints = resample_stroke2D(original_stroke2DPoints);

	Eigen::MatrixXd V2_tmp, V2;
	Eigen::MatrixXi F2, F2_back, vertex_markers, edge_markers;
	Eigen::MatrixXi stroke_edges(stroke2DPoints.rows(), 2);
	for(int i = 0; i < stroke2DPoints.rows(); i++) {
		stroke_edges.row(i) << i, ((i + 1) % stroke2DPoints.rows());
	}
	igl::triangle::triangulate((Eigen::MatrixXd) stroke2DPoints, stroke_edges, Eigen::MatrixXd(0, 0), Eigen::MatrixXi::Constant(stroke2DPoints.rows(), 1, 1), Eigen::MatrixXi::Constant(stroke_edges.rows(), 1, 1), "YQq25", V2_tmp, F2, vertex_markers, edge_markers);
	V2 = Eigen::MatrixXd::Zero(V2_tmp.rows(), V2_tmp.cols() + 1);

	//zero mean in x and y
	V2_tmp.col(0) = V2_tmp.col(0).array() - V2_tmp.col(0).mean();
	V2_tmp.col(1) = V2_tmp.col(1).array() - V2_tmp.col(1).mean();

	V2.block(0, 0, V2_tmp.rows(), 2) = V2_tmp;


	generate_backfaces(F2, F2_back);
	int nr_boundary_vertices = (vertex_markers.array() == 1).count();
	int original_size = V2.rows();
	V2.conservativeResize(2 * original_size - nr_boundary_vertices, Eigen::NoChange); //Increase size, such that boundary vertices only get included one
	unordered_map<int, int> backside_vertex_map;
	int count = 0;
	for(int i = 0; i < original_size; i++) {
		if(vertex_markers(i) != 1) { //If the vertex is NOT on the boundary, duplicate it
			V2.row(original_size + count) << V2.row(i);
			backside_vertex_map.insert({i, original_size + count});
			backside_vertex_map.insert({original_size + count, i});
			count++;
		} else {
			backside_vertex_map.insert({i, i});
		}
	}

	V2.col(0) = V2.col(0).array() - V2.col(0).mean();
	V2.col(1) = V2.col(1).array() - V2.col(1).mean();

	//Add backside faces, using original vertices for boundary and copied vertices for interior
	original_size = F2_back.rows();
	F2.conservativeResize(original_size * 2, Eigen::NoChange);
	for(int i = 0; i < original_size; i++) {
		for(int j = 0; j < 3; j++) {
			auto got = backside_vertex_map.find(F2_back(i, j));
			F2(original_size + i, j) = got->second;
		}
	}

	Eigen::MatrixXd N_Faces, N_Vertices;
	igl::per_vertex_normals(V2, F2, PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM, N_Vertices);

	vertex_boundary_markers.resize(V2.rows());
	part_of_original_stroke.resize(V2.rows());
	for(int i = 0; i < V2.rows(); i++) {
		if(i >= vertex_markers.rows()) { //vertex can't be boundary (it's on backside)
			V2.row(i) = V2.row(i) + 0.1*N_Vertices.row(i);
			vertex_boundary_markers[i] = 0;
			part_of_original_stroke[i] = 0;
		} else {
			if(vertex_markers(i) == 1) { //Don't change boundary vertices
				vertex_boundary_markers[i] = 1;
				part_of_original_stroke[i] = 1;
				stroke3DPoints.row(i) = V2.row(i);
				continue;
			}
			V2.row(i) = V2.row(i) + 0.1*N_Vertices.row(i);
			vertex_boundary_markers[i] = 0;
			part_of_original_stroke[i] = 0;
		}
	}

	viewer.data.clear();
	viewer.data.set_mesh(V2, F2);
	
	igl::per_face_normals(V2, F2, N_Faces);
	viewer.data.set_normals(N_Faces);

	viewer.core.align_camera_center(viewer.data.V);
	return backside_vertex_map;
}

/** Makes the stroke counter clockwise. **/
void Stroke::counter_clockwise() {
	double total_area = 0;
	Eigen::Vector2d prev, next;
	prev = (Eigen::Vector2d) stroke2DPoints.row(stroke2DPoints.rows() - 1);
	for(int i = 0; i < stroke2DPoints.rows(); i++) {
		next = (Eigen::Vector2d) stroke2DPoints.row(i);
		total_area += (prev[1] + next[1]) * (next[0] - prev[0]);
		prev = next;
	}

	if(total_area > 0) { //reverse the vector
		stroke2DPoints = stroke2DPoints.colwise().reverse().eval();
		stroke3DPoints = stroke3DPoints.colwise().reverse().eval();
		has_been_reversed = true;
	}

}

/** Resamples the stroke by moving all vertices to the middle of their neighbors. Essentially has a smoothing effect. **/
Eigen::MatrixXd Stroke::resample_stroke2D(Eigen::MatrixXd & original_stroke2DPoints) {
	Eigen::MatrixXd new_stroke2DPoints = Eigen::MatrixXd::Zero(original_stroke2DPoints.rows(), 2);
	int nr_iterations = max(5.0, original_stroke2DPoints.rows() / 4.0);
	for(int i = 0; i < nr_iterations; i++) {
		move_to_middle(original_stroke2DPoints, new_stroke2DPoints);
		move_to_middle(new_stroke2DPoints, original_stroke2DPoints);
	}
	return new_stroke2DPoints;
}

void Stroke::move_to_middle(Eigen::MatrixXd &positions, Eigen::MatrixXd &new_positions) {
	int n = positions.rows();
	Eigen::Vector2d prev, cur, next;

	for(int i = 0; i < n; i++) {
		prev = positions.row(((i - 1) + n) % n);
		cur = positions.row(i % n);
		next = positions.row(((i + 1) + n) % n);

		new_positions(i, 0) = (cur[0] * 2 + prev[0] + next[0]) / 4;
		new_positions(i, 1) = (cur[1] * 2 + prev[1] + next[1]) / 4;
	}
}

void Stroke::generate_backfaces(Eigen::MatrixXi &faces, Eigen::MatrixXi &back_faces) {
	back_faces = faces.rowwise().reverse().eval();
}

/** Returns the ID of the stroke's 3D point that is closest to where the user clicked. Stores the distance to this point in closest_distance. If the user clicked too far away from any of the stroke's points, it will return -1. **/
int Stroke::selectClosestVertex(int mouse_x, int mouse_y, double& closest_distance) {
	double x = mouse_x;
	double y = viewer.core.viewport(3) - mouse_y;

	Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
	Eigen::RowVector3d pt, vert, clicked_pt(x, y, 0);
	double closest_dist = INFINITY, dist;
	int closest_ID;

	for(int i = 0; i < stroke3DPoints.rows(); i++) {
		vert = stroke3DPoints.row(i).transpose();
		igl::project(vert, modelview, viewer.core.proj, viewer.core.viewport, pt);
		dist = (pt - clicked_pt).squaredNorm();
		if(dist < closest_dist) {
			closest_dist = dist;
			closest_ID = i;
		}
	}

	closest_dist = sqrt(closest_dist);
	double stroke_diag = compute_stroke_diag();
	closest_distance = closest_dist;
	if(stroke3DPoints.rows() > 2 && closest_dist > 0.08*stroke_diag) { //Don't do pulling when the user clicks too far away from any curve (treat as navigation movement)
		return -1;
	} else if(stroke3DPoints.rows() == 2 && closest_dist > 0.08 * (stroke3DPoints.row(0).norm() * 2)) { //The stroke consist of a single point (that's "looped") so the diag will be 0. Use the point's norm instead
		return -1;
	}
	return closest_ID;
}

double Stroke::compute_stroke_diag() {
	Eigen::Vector3d maxBB = stroke3DPoints.colwise().maxCoeff();
	Eigen::Vector3d minBB = stroke3DPoints.colwise().minCoeff();
	return (maxBB - minBB).norm();
}

void Stroke::update_Positions(Eigen::MatrixXd V) {
	stroke3DPoints.resize(closest_vert_bindings.size(), Eigen::NoChange);
	for(int i = 0; i < closest_vert_bindings.size(); i++) { //Iterate over the (updated) closest_vert_bindings instead of over stroke3DPoints
		stroke3DPoints.row(i) = V.row(closest_vert_bindings[i]);
	}

	//In the case of extrusion silhouette strokes, closest_vert_bindings isn't looped. Don't make stroke3DPoints looped, because we already account for it not being a loop when drawing the curves
}

/** Remaps the stroke's vertex bindings after the mesh topology has changed. If a stroke becomes non-continous (can only happen to strokes that were non-looped at the start and that did not have their first and/or last point removed) the stroke will be removed and its boundary_vertex_markers unset. Looped strokes always stay continuous, since we can ony remove one continuous piece of them.
    Note that this method does not update the stroke3DPoints, use update_Positions() for that. **/
bool Stroke::update_vert_bindings(Eigen::VectorXi & new_mapped_indices, Eigen::VectorXi& vertex_boundary_markers) {
	vector<int> new_bindings;
	int first_included_after_remove = -1;
	bool points_were_removed = false, no_tracked_point_yet = true, stays_continuous = false, originally_is_loop = is_loop;
	
    for(int i = 0; i < closest_vert_bindings.size() - 1; i++) {	//Closest_vert_bindings is always a loop

		if(new_mapped_indices[closest_vert_bindings[i]] == -1) { //Removed vertex
			if(!originally_is_loop && (i==0 || i==closest_vert_bindings.size()-2)) {
				stays_continuous = true;
			}

			is_loop = false;
			points_were_removed = true;
			continue;
		}
		new_bindings.push_back(new_mapped_indices[closest_vert_bindings[i]]);
		if(points_were_removed && no_tracked_point_yet) {
			first_included_after_remove = new_bindings.size() - 1;
			no_tracked_point_yet = false;
		}
	}

	if(new_bindings.size() == 0) {
		return false; //Stroke ceases to exist.
	}

	if(first_included_after_remove != -1) {	//Check that we didn't remove ALL points in the stroke and that we actually need to reorder (because point 0 was removed)
		rotate(new_bindings.begin(), new_bindings.begin() + first_included_after_remove, new_bindings.end()); //Makes the stroke "continuous" again
	}

	new_bindings.push_back(new_bindings[0]); //Make the vertex bindings looped again
	set_closest_vert_bindings(new_bindings);

	if(!originally_is_loop && points_were_removed && !stays_continuous) { //For non-looped strokes that have a middle chunk removed (but not everything), unset their markers because we will remove the entire stroke
		undo_stroke_add(vertex_boundary_markers);
		return false;
	}


	return true;
}

/** Used for ADD. Snaps the drawn stroke3DPoints to the closest mesh vertices and sets their boundary_markers accordingly. **/
void Stroke::snap_to_vertices(Eigen::VectorXi &vertex_boundary_markers) {
	Eigen::VectorXd min_dist;
	Eigen::VectorXi previous;
	vector<int> result_path, added_stroke_final_vertices;
	int prev = -1;
	vector<vector<int>> adj_list;
	adjacency_list(F, adj_list);

	if((stroke3DPoints.row(0) - stroke3DPoints.row(stroke3DPoints.rows() - 1)).norm() < compute_stroke_diag() / 5.0) {
		is_loop = true;
	} else {
		is_loop = false;
	}

	stroke3DPoints.resize(0, 3);
	for(int i = 0; i < closest_vert_bindings.size() - 1; i++) {
		set<int> goal;
		goal.insert(closest_vert_bindings[i + 1]);
		dijkstra_compute_paths(closest_vert_bindings[i], goal, adj_list, min_dist, previous);
		dijkstra_get_shortest_path_to(closest_vert_bindings[i + 1], previous, result_path);

		for(int j = result_path.size() - 1; j >= 0; j--) {
			int idx = result_path[j];
			if(idx != prev) {
				added_stroke_final_vertices.push_back(idx);
				stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
				stroke3DPoints.row(stroke3DPoints.rows() - 1) << V.row(added_stroke_final_vertices.back());
				vertex_boundary_markers[added_stroke_final_vertices.back()] = stroke_ID;
				prev = idx;
			}
		}
	}

	closest_vert_bindings = added_stroke_final_vertices; //closest_vert_bindings will be updated when the topology changes, added_stroke_final_vertices not, so only use it as a temporary
	stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
	stroke3DPoints.row(stroke3DPoints.rows() - 1) << stroke3DPoints.row(0); // Also add a copy of the first point at the end for strokes that are later added, to keep consistency with the initial stroke
	if(closest_vert_bindings[0] != closest_vert_bindings[closest_vert_bindings.size() - 1]) {
		closest_vert_bindings.push_back(closest_vert_bindings[0]); //"Loop" the closest_vert_bindings because they are used to determine which 3DPoints to update, and we need 3DPoints to stay a loop for the sake of draw
	}
}

void Stroke::undo_stroke_add(Eigen::VectorXi& vertex_boundary_markers) {
	for(int i = 0; i < closest_vert_bindings.size(); i++) {
		vertex_boundary_markers[closest_vert_bindings[i]] = 0;
	}
}

int Stroke::get_ID() {
	return stroke_ID;
}

Eigen::MatrixXd Stroke::get_V() const {
	return V;
}

Eigen::MatrixXi Stroke::get_F() const {
	return F;
}

Eigen::MatrixXd Stroke::get_stroke2DPoints() const {
	return stroke2DPoints;
}

Eigen::MatrixX3d Stroke::get3DPoints() {
	return stroke3DPoints;
}

int Stroke::get_vertex_idx_for_point(int pt_idx) {
	return closest_vert_bindings[pt_idx];
}

vector<int> Stroke::get_closest_vert_bindings() {
	return closest_vert_bindings;
}

void Stroke::set3DPoints(Eigen::MatrixX3d new_3DPoints) {
	stroke3DPoints.resize(new_3DPoints.rows(), 3);
	stroke3DPoints = new_3DPoints;
}

void Stroke::set_closest_vert_bindings(vector<int> new_vert_bindings) {
	closest_vert_bindings.clear();
	closest_vert_bindings = new_vert_bindings;
}
