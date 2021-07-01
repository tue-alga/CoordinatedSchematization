#pragma once
#include <vector>

template<class NT>
class DummyEdge
{
public:
	NT length;
	std::tuple<size_t, size_t> vertices;
	DummyEdge<NT>(NT in_length, size_t v1, size_t v2) { length = in_length; vertices = std::make_tuple(v1, v2); }
	DummyEdge<NT>() { 
	}
};

class DummyEulerianTrajectory
{
public:
	std::vector<size_t> EdgeIndices;
	std::vector<size_t>::iterator begin() { return EdgeIndices.begin(); }
	std::vector<size_t>::iterator end() { return EdgeIndices.end(); }
	DummyEulerianTrajectory(std::vector<size_t> in_indices) { EdgeIndices = in_indices; }
};

template<class NT>
class DummyVertex
{
public:
	NT lat, lon;
	std::vector<size_t> edges;
	size_t debug_index;
	DummyVertex(NT in_lat, NT in_lon, size_t in_debug_index = 0) { lat = in_lat; lon = in_lon; debug_index = in_debug_index; }
	DummyVertex() { lat = 0; lon = 0; debug_index = 0; }
};

template<class NT>
class SuperEdge
{
public:
	std::tuple<size_t, size_t> endpoint_super_vertices;
	std::vector<size_t> component_edges;
	NT length;
};

class SuperVertex
{
public:
	size_t base_vertex;
	std::vector<size_t> incident_superedges;
};