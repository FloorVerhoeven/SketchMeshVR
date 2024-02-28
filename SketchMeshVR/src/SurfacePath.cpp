#include <igl/edge_topology.h>
#include "SurfacePath.h"

Eigen::MatrixXi EV, FE, EF;
Stroke* origin_stroke;
Eigen::MatrixX3d looped_3DPoints;

using namespace std;

SurfacePath::SurfacePath() {
}

bool SurfacePath::create_from_open_path(const Stroke& stroke) {
	origin_stroke = new Stroke(stroke);
	path.clear();

	int prev_p = 0;
	int start_p = prev_p;
	int next_p;

	int faceID = origin_stroke->get_hit_faces()(prev_p, 0);
	int start_face = faceID;
	int n = 1, iter = 0;
	Eigen::RowVector3d pt(0, 0, 0);

	int nr3DPoints = origin_stroke->get3DPoints().rows();
	looped_3DPoints = origin_stroke->get3DPoints().topRows(nr3DPoints);

	igl::edge_topology(stroke.get_V(), stroke.get_F(), EV, FE, EF);

	while (true) {
		next_p = n;

		pt = looped_3DPoints.row(prev_p);
		PathElement newElement(faceID, PathElement::FACE, pt);
		path.push_back(newElement);

		if (next_p == start_p) { //Don't need to come back to the start face here since the open path ain't a loop
			break;
		}

		faceID = extend_path_extrude(prev_p, next_p, faceID);
		if (faceID == -1 || iter == 5000) {
			return false;
		}

		n = (next_p + 1) % nr3DPoints;
		prev_p = next_p;
		iter++;
	}

	return true;
}

/** Creates a SurfacePath that contains both the original points in stroke, and new points at the locations where stroke segments cross face edges. Won't wrap around to the backside of the mesh (because it will arrive at the first index again before having to switch direction). Used for extrusion **/
bool SurfacePath::create_from_stroke_extrude(const Stroke & stroke) {
	origin_stroke = new Stroke(stroke);
	path.clear();
	
	int prev_p = stroke.get_stroke2DPoints().rows() - 1;
	int start_p = prev_p;
	int next_p;

	int faceID = origin_stroke->get_hit_faces()(prev_p, 0);
	int start_face = faceID;
	int n = 0, iter = 0;
	Eigen::RowVector3d pt(0, 0, 0);

	int nr3DPoints = origin_stroke->get3DPoints().rows() - 1; //Don't take the loop duplicate
	looped_3DPoints = origin_stroke->get3DPoints().topRows(nr3DPoints);

	igl::edge_topology(stroke.get_V(), stroke.get_F(), EV, FE, EF);	

	while(true) {
		next_p = n;

		pt = looped_3DPoints.row(prev_p);
		PathElement newElement(faceID, PathElement::FACE, pt);
		path.push_back(newElement);

 		faceID = extend_path_extrude(prev_p, next_p, faceID);
		if (faceID == -1 || iter == 5000) {
			return false;
		}

		if (next_p == start_p && faceID == start_face) {
			break;
		}

		n = (next_p + 1) % nr3DPoints;
		prev_p = next_p;
		iter++;
	}
		pt = looped_3DPoints.row(start_p); //The origin_stroke's 3DPoints will be updated with the Path's vertex positions (in MeshExtrusion), so we need to keep it looped
		PathElement lastElement(faceID, PathElement::FACE, pt);
		path.push_back(lastElement);
	
	return true;
}

int SurfacePath::extend_path_extrude(int prev_p, int next_p, int faceID) {
	Eigen::Vector3d edge_cut_point;
	int edge = -1;
	int iter = 0;

	while (true) {
		if (origin_stroke->get_hit_faces()(next_p, 0) == faceID) { //next_p is in same face as prev_p
			return faceID;
		}

		edge = find_next_edge_extrude(next_p, prev_p, edge, faceID, edge_cut_point);
		iter++;
		if (edge == -1 || iter > 1000) { //Something is wrong with the stroke, exit gracefully
			return -1;
		}

		PathElement newElement(edge, PathElement::EDGE, edge_cut_point);
		path.push_back(newElement);

		faceID = (EF(edge, 0) == faceID) ? EF(edge, 1) : EF(edge, 0); //get the polygon on the other side of the edge
	}
}

//Finds the edge index of the edge to cross in order to get from polygon to (somewhere closer to) next_p
int SurfacePath::find_next_edge_extrude(int next_p, int prev_p, int prev_edge, int polygon, Eigen::Vector3d& edge_cut_point) {
	::Plane cutPlane(origin_stroke->get_hand_pos().row(prev_p), looped_3DPoints.row(prev_p), looped_3DPoints.row(next_p));
	int next_faceID = origin_stroke->get_hit_faces()(next_p, 0);

	for (int i = 0; i < 3; i++) {
		int edge = FE(polygon, i);
		if (edge != prev_edge) {
			int other_faceID = (EF(edge, 0) == polygon) ? EF(edge, 1) : EF(edge, 0);
			if (other_faceID == next_faceID) {
				edge_cut_point = cutPlane.cross_point(origin_stroke->get_V().row(EV(edge, 0)), origin_stroke->get_V().row(EV(edge, 1)));
				return edge;
			}
		}
	}

	//If we come here, it means that the next strokepoint (next_p) lies in a non-adjacent face (e.g. one or multiple faces away). Only then go for more elaborate checking 
	double t_val;
	for (int i = 0; i < 3; i++) {
		int edge = FE(polygon, i);
		if (edge != prev_edge) {
			Eigen::RowVector3d cut_point = cutPlane.cross_point(origin_stroke->get_V().row(EV(edge, 0)), origin_stroke->get_V().row(EV(edge, 1)), t_val);
			Eigen::RowVector3d seg_vec = looped_3DPoints.row(next_p) - looped_3DPoints.row(prev_p);
			seg_vec.normalize();
			Eigen::RowVector3d its_vec = cut_point - looped_3DPoints.row(prev_p);
			its_vec.normalize();
			if (its_vec.dot(seg_vec) >= 0 && t_val<= 0.9999999 && t_val>= 0.0000001) { //Else the projection points into the opposite direction of line segment
				edge_cut_point = cut_point;
				return edge;
			}
		}
	}

	return -1;
}

/** Creates a SurfacePath that contains both the original points in stroke, and new points at the locations where stroke segments cross face edges. Also wraps around to the backside of the mesh. Used for cutting **/
bool SurfacePath::create_from_stroke_cut(const Stroke & stroke) {
	origin_stroke = new Stroke(stroke);
	path.clear();

	int prev_p = 1; //Start at point 1, because point 0 is defined to be the last point at the beginning of the stroke to lie outside of the mesh
	int start_p = prev_p, next_p;
	int faceID = origin_stroke->get_hit_faces()(prev_p, 0), prev_faceID;
	int iter = 0;

	looped_3DPoints = create_loop_from_front_and_back(origin_stroke->get3DPoints(), origin_stroke->get3DPointsBack());
	int nr_looped_3DPoints = looped_3DPoints.rows();

	int start_face = faceID;
	int n = 2;
	Eigen::RowVector3d pt(0, 0, 0);

	igl::edge_topology(stroke.get_V(), stroke.get_F(), EV, FE, EF);
	bool on_front_side = true, first_iter = true;
	int edge = -1;
	while(true) {
		next_p = n;

		//Used to detect whether point prev_p is outside the mesh (when it's faceID is -1).
		//If it's outside of the mesh, we won't add a PathElement for it
		prev_faceID = faceID;
		int prev_p_tmp = (prev_p > nr_looped_3DPoints / 2) ? nr_looped_3DPoints - prev_p : prev_p;
		faceID = origin_stroke->get_hit_faces()(prev_p_tmp, !on_front_side);
		if (faceID == -1) { //We're dealing with one of the 2 off-mesh vertices
			faceID = prev_faceID; //Reset to the previous faceID
		} else {
			pt = looped_3DPoints.row(prev_p);
			PathElement newElement(faceID, PathElement::FACE, pt);
			path.push_back(newElement);
		}

		faceID = extend_path_cut(prev_p, next_p, faceID, on_front_side, edge, first_iter);
		if (faceID == -1 || iter == 5000) { //Either the cut stroke goes exactly over a vertex or we cannot find a closed path that starts and ends on the same face
			std::cerr << "Surface path could not be created for this cut stroke. Please try again. " << std::endl;
			return false;
		}

		if(next_p == start_p && faceID == start_face) {
			break;
		}

		n = (next_p + 1) % nr_looped_3DPoints;
		prev_p = next_p;
		iter++;
	}

	return true; //MeshCut will make the stroke3DPoints looped again, don't need to do it here (unlike for extrusion)
}

/** Determines the moving direction when going from prev_p to next_p and adds new vertices at mesh edges when the segment from prev_p to next_p crosses an edge. **/
int SurfacePath::extend_path_cut(int prev_p, int next_p, int faceID, bool& on_front_side, int& edge, bool& first_iter) {
	Eigen::RowVector3d plane_point;
	int nr_looped_3DPoints = looped_3DPoints.rows();
	int tmp_prev_p = (prev_p > nr_looped_3DPoints / 2) ? nr_looped_3DPoints - prev_p : prev_p;
	if (origin_stroke->get_hit_faces()(tmp_prev_p, 0) == -1) {
		int tmp_next_p = nr_looped_3DPoints - next_p; //Get the point that corresponds to next_p but is on the other side of the mesh
		plane_point = looped_3DPoints.row(tmp_next_p);
	}
	else {
		int tmp_prev_p = nr_looped_3DPoints - prev_p; 
		plane_point = looped_3DPoints.row(tmp_prev_p); //Get the point that corresponds to prev_p but is on the other side of the mesh
	}
	::Plane cutPlane(plane_point, looped_3DPoints.row(prev_p), looped_3DPoints.row(next_p));


	Eigen::RowVector3d start_pos = looped_3DPoints.row(prev_p);
	Eigen::RowVector3d end_pos = looped_3DPoints.row(next_p);
	Eigen::Vector3d edge_cut_point;
	int iter = 0;

	while(true) {
		//next_p and prev_P are used to index into hit_faces, which only has data stored for "the front half" and then has the info for the "back half" in the second column
		next_p = (next_p > nr_looped_3DPoints / 2) ? nr_looped_3DPoints - next_p : next_p; 
		prev_p = (prev_p > nr_looped_3DPoints / 2) ? nr_looped_3DPoints - prev_p : prev_p; 


		if(origin_stroke->get_hit_faces()(next_p, 0) == faceID || origin_stroke->get_hit_faces()(next_p, 1) == faceID){ //next_p is in same face as prev_p
			return faceID;
		}

		pair<int, int> strokeEdge(prev_p, next_p);
		edge = find_next_edge_cut(strokeEdge, edge, faceID, on_front_side, cutPlane, start_pos, end_pos, first_iter, edge_cut_point);
		first_iter = false;
		iter++;
		if(edge == -1 || iter > 1000) { //Something is wrong with the stroke, exit gracefully
			return -1;
		}

		PathElement newElement(edge, PathElement::EDGE, edge_cut_point);
		path.push_back(newElement);

		faceID = (EF(edge, 0) == faceID) ? EF(edge, 1) : EF(edge, 0); //Get the polygon on the other side of the edge

		if(origin_stroke->get_hit_faces()(prev_p, 0) == faceID || origin_stroke->get_hit_faces()(prev_p, 1) == faceID){
            //This means that the current point (prev_p) is projected into both "its own" polygon and into the polygon of next_p, while next_p's polygon is across an edge. Means that next_p's polygon is on the backside
			on_front_side = !on_front_side;
			return faceID;
		}
	} 
}

/** Finds and returns the edge ID of the edge that is being crossed by the segment from strokeEdge.start to strokeEdge.end. Returns -1 if no such edge exists. **/
int SurfacePath::find_next_edge_cut(pair<int, int> strokeEdge, int prev_edge, int polygon, bool on_front_side, ::Plane& cutPlane, Eigen::RowVector3d& start_pos, Eigen::RowVector3d& end_pos, bool first_iter, Eigen::Vector3d& edge_cut_point) {
	int next_faceID;
	if (origin_stroke->get_hit_faces()(strokeEdge.second, !on_front_side) == -1) { //Next point is a point outside of the mesh. Find the edge with prev_p's hit on the other side instead
		next_faceID = origin_stroke->get_hit_faces()(strokeEdge.first, on_front_side);
	}
	else {
		next_faceID = origin_stroke->get_hit_faces()(strokeEdge.second, !on_front_side);
	}
	for (int i = 0; i < 3; i++) {
		int edge = FE(polygon, i);
		if (edge != prev_edge) {
			int other_faceID = (EF(edge, 0) == polygon) ? EF(edge, 1) : EF(edge, 0);
			if (other_faceID == next_faceID) {
				edge_cut_point = cutPlane.cross_point(origin_stroke->get_V().row(EV(edge, 0)), origin_stroke->get_V().row(EV(edge, 1)));
				return edge;
			}
		}
	}

	
	//If we come here, it means that the next strokepoint (next_p) lies in a non-adjacent face (e.g. one or multiple faces away). Only then go for more elaborate checking 
	double t_val;
	for (int i = 0; i < 3; i++) {
		int edge = FE(polygon, i);
		if (edge != prev_edge) {
			Eigen::RowVector3d cut_point = cutPlane.cross_point(origin_stroke->get_V().row(EV(edge, 0)), origin_stroke->get_V().row(EV(edge, 1)), t_val);
			Eigen::RowVector3d seg_vec = end_pos - start_pos;
			seg_vec.normalize();
			Eigen::RowVector3d its_vec = cut_point - start_pos;
			its_vec.normalize();
			if ((t_val<=0.9999999 && t_val>=0.0000001) && (!first_iter || (first_iter && its_vec.dot(seg_vec)>=0))) { //Else we cross the "line" outside of its range or the first extension goes into the wrong direction
				edge_cut_point = cut_point;
				return edge;
			}
		}
	}

	return -1;
}

Eigen::MatrixX3d SurfacePath::create_loop_from_front_and_back(Eigen::MatrixX3d& front_3DPoints, Eigen::MatrixX3d& back_3DPoints) {
	Eigen::MatrixX3d result(front_3DPoints.rows() - 1 + back_3DPoints.rows(), 3); //Take all of 3DPointsBack, because we never add back-points for the points that are off the mesh
	result << front_3DPoints.topRows(front_3DPoints.rows() - 1), back_3DPoints.colwise().reverse().eval();
	return result;
}

//Used when the user navigates, to rotate the points that are inside the SurfacePath 
void SurfacePath::set_rotated_points(Eigen::MatrixXd& V) {
	for (int i = 0; i < path.size(); i++) {
		path[i].set_vertex(V.row(i));
	}
}

Eigen::MatrixXd SurfacePath::get3DPoints() {
	Eigen::MatrixXd points(path.size(), 3);
	for (int i = 0; i < path.size(); i++) {
		points.row(i) = path[i].get_vertex();
	}
	return points;
}

vector<PathElement> SurfacePath::get_path() {
	return path;
}

void SurfacePath::set_path(std::vector<PathElement> new_path){
	path = new_path;
}

int SurfacePath::get_origin_stroke_ID() {
	return origin_stroke->get_ID();
}

PathElement & SurfacePath::get_path_element(int i) {
	return path[i];
}