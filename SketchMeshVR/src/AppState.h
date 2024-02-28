#ifndef _APPSTATE_H_
#define _APPSTATE_H_

#include "Stroke.h"

class AppState {
public: 
	AppState() {};
	AppState& operator = (AppState& other) {
		vertex_is_fixed = other.vertex_is_fixed;
		edge_boundary_markers = other.edge_boundary_markers;
		sharp_edge = other.sharp_edge;
		new_mapped_indices = other.new_mapped_indices;
		replacing_vertex_bindings = other.replacing_vertex_bindings;
		stroke_collection_size = other.stroke_collection_size;
		return *this;
	}

	void init(Eigen::VectorXi vertex_is_fixed_, Eigen::VectorXi edge_boundary_markers_, Eigen::VectorXi sharp_edge_, Eigen::VectorXi new_mapped_indices_, Eigen::MatrixXi replacing_vertex_bindings_, int stroke_collection_size_) {
		vertex_is_fixed = vertex_is_fixed_;
		edge_boundary_markers = edge_boundary_markers_;
		sharp_edge = sharp_edge_;
		new_mapped_indices = new_mapped_indices_;
		replacing_vertex_bindings = replacing_vertex_bindings_;
		stroke_collection_size = stroke_collection_size_;
	};


	Eigen::VectorXi vertex_is_fixed;
	Eigen::VectorXi edge_boundary_markers;
	Eigen::VectorXi sharp_edge;
	Eigen::VectorXi new_mapped_indices;
	Eigen::MatrixXi replacing_vertex_bindings;
	int stroke_collection_size;
private:

};

#include <igl/serialize.h>
namespace igl {
	namespace serialization {

		inline void serialization(bool s, AppState& obj, std::vector<char>& buffer)
		{
			SERIALIZE_MEMBER(vertex_is_fixed);
			SERIALIZE_MEMBER(edge_boundary_markers);
			SERIALIZE_MEMBER(sharp_edge);
			SERIALIZE_MEMBER(new_mapped_indices);
			SERIALIZE_MEMBER(replacing_vertex_bindings);
			SERIALIZE_MEMBER(stroke_collection_size);
		}

		template<>
		inline void serialize(const AppState& obj, std::vector<char>& buffer)
		{
			serialization(true, const_cast<AppState&>(obj), buffer);
		}

		template<>
		inline void deserialize(AppState& obj, const std::vector<char>& buffer)
		{
			serialization(false, obj, const_cast<std::vector<char>&>(buffer));
		}
	}
}


#endif