#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/project.h>
#include "SurfacePath.h"
#include "Plane.h"

Eigen::MatrixXi EV, FE, EF;
Stroke* origin_stroke;

using namespace std;

SurfacePath::SurfacePath() {

}

/** Creates a SurfacePath that contains both the original points in stroke, and new points at the locations where stroke segments cross face edges. Won't wrap around to the backside of the mesh (because it will arrive at the first index again before having to switch direction). Used for extrusion **/
bool SurfacePath::create_from_stroke_extrude(const Stroke & stroke) {
	origin_stroke = new Stroke(stroke);
	path.clear();
	Eigen::Matrix4f modelview = stroke.viewer.core.view * stroke.viewer.core.model;
	
	Eigen::Vector3f bc;
	int faceID = -1;
	int prev_p = stroke.get_stroke2DPoints().rows() - 1;
	int start_p = prev_p;
	int next_p;

	igl::unproject_onto_mesh(stroke.get_stroke2DPoints().row(prev_p).cast<float>(), modelview, stroke.viewer.core.proj, stroke.viewer.core.viewport, stroke.get_V(), stroke.get_F(), faceID, bc);
	int start_face = faceID;
	int n = 0;
	Eigen::RowVector3d pt(0, 0, 0);

	igl::edge_topology(stroke.get_V(), stroke.get_F(), EV, FE, EF);

	while(true) {
		next_p = n;

		pt = unproject_onto_polygon(stroke.get_stroke2DPoints().row(prev_p), faceID, modelview);
		PathElement newElement(faceID, PathElement::FACE, pt);
		path.push_back(newElement);

		bool forward;
		faceID = extend_path(prev_p, next_p, faceID, forward, modelview);
		if (faceID == -1) {
			return false;
		}

		if(!forward) {
			next_p = prev_p;
		}

		if(next_p == start_p && faceID == start_face) {
			break;
		}

		if(front_facing(faceID)) {
			n = next_p + 1;
		} else {
			n = next_p - 1;
		}
		prev_p = next_p;
	}
	pt = unproject_onto_polygon(stroke.get_stroke2DPoints().row(start_p), faceID, modelview);
	PathElement lastElement(faceID, PathElement::FACE, pt);
	path.push_back(lastElement);
	return true;
}


/** Creates a SurfacePath that contains both the original points in stroke, and new points at the locations where stroke segments cross face edges. Also wraps around to the backside of the mesh. Used for cutting **/
bool SurfacePath::create_from_stroke(const Stroke & stroke) {
	origin_stroke = new Stroke(stroke);
	path.clear();
	Eigen::Matrix4f modelview = stroke.viewer.core.view * stroke.viewer.core.model;

	Eigen::Vector3f bc;
	int faceID = -1;
	int prev_p = 1; //Start at point 1, because point 0 is defined to be the last point at the beginning of the stroke to lie outside of the mesh
	int start_p = prev_p;
	int next_p;
	igl::unproject_onto_mesh(stroke.get_stroke2DPoints().row(prev_p).cast<float>(), modelview, stroke.viewer.core.proj, stroke.viewer.core.viewport, stroke.get_V(), stroke.get_F(), faceID, bc);

	int start_face = faceID;
	int n = 2;
	Eigen::RowVector3d pt(0, 0, 0);

	igl::edge_topology(stroke.get_V(), stroke.get_F(), EV, FE, EF);

	while(true) {
		next_p = n;

		pt = unproject_onto_polygon(stroke.get_stroke2DPoints().row(prev_p), faceID, modelview);
		PathElement newElement(faceID, PathElement::FACE, pt);
		path.push_back(newElement);

		bool forward;
		faceID = extend_path(prev_p, next_p, faceID, forward, modelview);
		if (faceID == -1) {
			return false;
		}

		if(!forward) {
			next_p = prev_p;
		}

		if(next_p == start_p && faceID == start_face) {
			break;
		}

		if(front_facing(faceID)) {
			n = next_p + 1;
		} else {
			n = next_p - 1;
		}
		prev_p = next_p;
	}
	return true;
}

/** Determines the moving direction when going from prev_p to next_p and adds new vertices at mesh edges when the segment from prev_p to next_p crosses an edge. **/
int SurfacePath::extend_path(int prev_p, int next_p, int faceID, bool& forward, Eigen::Matrix4f modelview) {
	Eigen::Vector3d source, dir;
	Eigen::MatrixX2d stroke2DPoints = origin_stroke->get_stroke2DPoints();
	Eigen::MatrixX3d stroke3DPoints = origin_stroke->get3DPoints();
	int prev_p3D = origin_stroke->has_been_reversed ? prev_p + 1 : prev_p; //These two are needed because when a CW stroke is reversed to a CCW stroke, it generates a discrepancy between the points in the 2DPoints and 3DPoints. 3DPoints is a loop and has vertex 0 also at the end, while 2DPoints doesn't. Fixing this would break many other things
	int next_p3D = origin_stroke->has_been_reversed ? next_p + 1 : next_p;
	Eigen::Vector2d tmp = stroke2DPoints.row(prev_p);
	igl::unproject_ray(tmp, modelview, origin_stroke->viewer.core.proj, origin_stroke->viewer.core.viewport, source, dir);
	Plane cutPlane(source, stroke3DPoints.row(prev_p3D), stroke3DPoints.row(next_p3D));

	int edge = -1, proj_faceID = -1;
	pair<int, int> strokeEdge(prev_p, next_p);
	Eigen::Vector3f bc;
	int iter = 0;

	while(true) {
		if(is_projected_inside(stroke2DPoints.row(next_p), faceID, modelview)) {
			forward = true;
			return faceID;
		}

		edge = find_next_edge(strokeEdge, edge, faceID, modelview);
		iter++;
		if(edge == -1 || iter > 1000) { //Something is wrong with the stroke, exit gracefully
			return -1;
		}

		Eigen::Vector3d v = cutPlane.cross_point(origin_stroke->get_V().row(EV(edge, 0)), origin_stroke->get_V().row(EV(edge, 1)));
		PathElement newElement(edge, PathElement::EDGE, v);
		path.push_back(newElement);

		faceID = (EF(edge, 0) == faceID) ? EF(edge, 1) : EF(edge, 0); //get the polygon on the other side of the edge
		
		if(is_projected_inside(stroke2DPoints.row(prev_p), faceID, modelview)) {
            //This means that the current point (prev_p) is projected into both "its own" polygon and into the polygon of next_p, while next_p's polygon is across an edge. Means that next_p's polygon is on the backside
			forward = false;
			return faceID;
		}
	}
}

/** Finds and returns the edge ID of the edge that is being crossed by the segment from strokeEdge.start to strokeEdge.end. Returns -1 if no such edge exists. **/
int SurfacePath::find_next_edge(pair<int, int> strokeEdge, int prev_edge, int polygon, Eigen::Matrix4f modelview) {
	Eigen::MatrixX2d stroke2DPoints = origin_stroke->get_stroke2DPoints();
	Eigen::RowVector2d stroke_start = stroke2DPoints.row(strokeEdge.first).transpose();
	Eigen::RowVector2d stroke_end = stroke2DPoints.row(strokeEdge.second).transpose();
	Eigen::RowVector3d start, end;
	Eigen::RowVector3d tmp;
	for(int i = 0; i < 3; i++) {
		int edge = FE(polygon, i);
		if(edge != prev_edge) {
			tmp = origin_stroke->get_V().row(EV(edge, 0));
			igl::project(tmp, modelview, origin_stroke->viewer.core.proj, origin_stroke->viewer.core.viewport, start);

			tmp = origin_stroke->get_V().row(EV(edge, 1));
			igl::project(tmp, modelview, origin_stroke->viewer.core.proj, origin_stroke->viewer.core.viewport, end);

			if(edges2D_cross({stroke_start, stroke_end}, {start.block(0,0,1,2).transpose(), end.block(0,0,1,2).transpose()})) {
				return edge;
			}
		}
	}
	return -1;
}

/** Determines whether a pair of 2D segments crosses eachother. Follows principle from https://stackoverflow.com/questions/14176776/find-out-if-2-lines-intersect but slightly different **/
bool SurfacePath::edges2D_cross(pair<Eigen::Vector2d, Eigen::Vector2d> edge1, pair<Eigen::Vector2d, Eigen::Vector2d> edge2) {
	double a0, b0, c0, a1, b1, c1;
	a0 = edge1.first[1] - edge1.second[1];
    b0 = edge1.second[0] - edge1.first[0];
	c0 = edge1.second[1] * edge1.first[0] - edge1.second[0] * edge1.first[1];
	a1 = edge2.first[1] - edge2.second[1];
	b1 = edge2.second[0] - edge2.first[0];
	c1 = edge2.second[1] * edge2.first[0] - edge2.second[0] * edge2.first[1];

	if(((a0*edge2.first[0] + b0*edge2.first[1] + c0)*(a0*edge2.second[0] + b0*edge2.second[1] + c0) <= 0) &&
		((a1*edge1.first[0] + b1*edge1.first[1] + c1)*(a1*edge1.second[0] + b1*edge1.second[1] + c1) <= 0)) {
		return true;
	} else {
		return false;
	}
}

/** Returns whether a face is front facing in the current viewer. **/
bool SurfacePath::front_facing(int faceID) {
	return is_counter_clockwise(faceID);
}

/** Determines whether the vertices of an edge are in counter-clockwise order as seen from the current viewpoint. **/
bool SurfacePath::is_counter_clockwise(int faceID) {
	Eigen::Matrix4f modelview = origin_stroke->viewer.core.view * origin_stroke->viewer.core.model;
	double total_area = 0.0;
	Eigen::RowVector3d prev, next, pt, center(0, 0, 0);
	Eigen::Vector2d prev2D, next2D;
	for(int i = 0; i < 3; i++) {
		center += origin_stroke->get_V().row(origin_stroke->get_F()(faceID, i));
	}
	center /= 3;

	prev = origin_stroke->get_V().row(origin_stroke->get_F()(faceID, 2)) - center;
	igl::project(prev, modelview, origin_stroke->viewer.core.proj, origin_stroke->viewer.core.viewport, pt);
	prev2D = pt.leftCols(2).transpose();
	for(int i = 0; i < 3; i++) {
		next = origin_stroke->get_V().row(origin_stroke->get_F()(faceID, i)) - center;
		igl::project(next, modelview, origin_stroke->viewer.core.proj, origin_stroke->viewer.core.viewport, pt);
		next2D = pt.leftCols(2).transpose();
		total_area += (prev2D[1] + next2D[1]) * (next2D[0] - prev2D[0]);
		prev2D = next2D;
	}

	if(total_area > 0) {
		return false;
	}
	return true;

}

/** Checks if point v is projected in face as seen from the current viewer. Checks if the point is on the correct side of all face edges. **/
bool SurfacePath::is_projected_inside(Eigen::RowVector2d v, int face, Eigen::Matrix4f modelview) {
	int sign = -1;
	if(!front_facing(face)) {
		sign = 1;
	}

	Eigen::RowVector3d start, end, tmp;
	Eigen::Vector2d vec0, vec1;
	for(int i = 0; i < 3; i++) {
		tmp = origin_stroke->get_V().row((origin_stroke->get_F()(face, i)));
		igl::project(tmp, modelview, origin_stroke->viewer.core.proj, origin_stroke->viewer.core.viewport, start);
		tmp = origin_stroke->get_V().row((origin_stroke->get_F()(face, (i + 1) % 3)));
		igl::project(tmp, modelview, origin_stroke->viewer.core.proj, origin_stroke->viewer.core.viewport, end);

		vec0 = (v - start.block(0, 0, 1, 2)).transpose();
		vec1 = (end.block(0, 0, 1, 2) - start.block(0, 0, 1, 2)).transpose();

		if(cross_prod2D(vec0, vec1)*sign < 0) {
			return false;
		}

	}
	return true;
}

int SurfacePath::cross_prod2D(Eigen::Vector2d vec0, Eigen::Vector2d vec1) {
    return vec0[0] * vec1[1] - vec0[1] * vec1[0];
}

/** Takes a 2D point and unprojects this onto the face with ID faceID. This does not always have to be the first face in the line of sight to the 3D point (as compared to igl::unproject_onto_mesh, which unprojects to the first triangle that is hit). **/
Eigen::Vector3d SurfacePath::unproject_onto_polygon(Eigen::Vector2d point, int faceID, Eigen::Matrix4f modelview) {
	Eigen::Vector3d source, dir;
	Eigen::Vector3d tmp1 = origin_stroke->get_V().row(origin_stroke->get_F()(faceID, 0));
	Eigen::Vector3d tmp2 = origin_stroke->get_V().row(origin_stroke->get_F()(faceID, 1));
	Eigen::Vector3d tmp3 = origin_stroke->get_V().row(origin_stroke->get_F()(faceID, 2));

	double t, u, v;
	igl::unproject_ray(point, modelview, origin_stroke->viewer.core.proj, origin_stroke->viewer.core.viewport, source, dir);
	intersect_triangle1(source.data(), dir.data(), tmp1.data(), tmp2.data(), tmp3.data(), &t, &u, &v);

	return source + t*dir;
}

vector<PathElement> SurfacePath::get_path() {
	return path;
}

int SurfacePath::get_origin_stroke_ID() {
	return origin_stroke->get_ID();
}

PathElement & SurfacePath::get_path_element(int i) {
	return path[i];
}
