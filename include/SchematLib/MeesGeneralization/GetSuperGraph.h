#pragma once
#include "DS/DummyGraph.h"
#include <queue>
#include <map>
template <class NT>
class GetSuperGraph
{
private:
	enum EdgeColor{white,gray,black};
public:
	template <class EdgeIterator, class VertexIterator, class TakenIterator>
	void operator()(VertexIterator v_begin, VertexIterator v_beyond, EdgeIterator e_begin, EdgeIterator e_beyond, TakenIterator v_taken_begin, TakenIterator e_taken_begin, std::vector<SuperVertex> &supervertices,std::vector<SuperEdge<NT>> &superedges, std::vector<int> &superedge_ids) 
	{
		//Go to the first taken node
		auto starting_node = v_begin;
		auto taken_it = v_taken_begin;
		while (!(*taken_it))
		{
			starting_node++;
			taken_it++;
		}
		//Start a bfs from here, constructing the graph as we go
		std::queue<size_t> Q;
		std::map<size_t, size_t> vertex_conversion;
		std::vector<EdgeColor> edge_colors = std::vector<EdgeColor>(std::distance(e_begin, e_beyond), EdgeColor::white);
		SuperVertex first_supervertex;
		first_supervertex.base_vertex = std::distance(v_begin, starting_node);
		supervertices.push_back(first_supervertex);
		superedge_ids = std::vector<int>(e_beyond - e_begin, -1);
		vertex_conversion[first_supervertex.base_vertex] = supervertices.size() - 1;
		Q.push(first_supervertex.base_vertex);

		while (Q.size() > 0)
		{
			auto front = Q.front();
			Q.pop();

			//Get the node
			auto node = *(v_begin + front);
			//For each incident edge, start the contraction procedure turning it into a superedge.
			for (auto inc_edge_it = node.edges.begin(); inc_edge_it != node.edges.end(); inc_edge_it++)
			{
				//Check if the edge is taken and white, if not we skip it
				if (!(*(e_taken_begin + *inc_edge_it)) || edge_colors[*inc_edge_it] == EdgeColor::black)
					continue;
				//Start construction on a new superedge
				SuperEdge<NT> new_super_edge;
				new_super_edge.component_edges.push_back(*inc_edge_it);
				size_t curr_edge = *inc_edge_it;
				size_t curr_vert = front;
				bool done = false;
				while (!done)
				{
					//Get the edge
					auto edge = *(e_begin + curr_edge);
					//Get the other vertex
					size_t other_vertex_id;
					if (get<0>(edge.vertices) == curr_vert)
						other_vertex_id = get<1>(edge.vertices);
					else other_vertex_id = get<0>(edge.vertices);
					auto other_vertex = *(v_begin + other_vertex_id);
					//Check if the other vertex is degree 2
					size_t other_vertex_degree = 0;
					for (auto other_vert_inc_edge_it = other_vertex.edges.begin(); other_vert_inc_edge_it != other_vertex.edges.end(); other_vert_inc_edge_it++)
					{
						if (!(*(e_taken_begin + *other_vert_inc_edge_it)))
							continue;
						other_vertex_degree++;
					}
					if (other_vertex_degree == 2)
					{
						//The vertex should be collapsed
						//Get the new edge
						size_t other_edge_id;
						if (other_vertex.edges[0] == curr_edge)
							other_edge_id = other_vertex.edges[1];
						else other_edge_id = other_vertex.edges[0];
						//Add the new edge to the superedge
						new_super_edge.component_edges.push_back(other_edge_id);
						//Iterate the loop
						curr_edge = other_edge_id;
						curr_vert = other_vertex_id;
						continue;
					}
					//color the superedge black
					//new_super_edge.component_edges.push_back(other_edge_id);
					for (size_t comp_edge : new_super_edge.component_edges)
						edge_colors[comp_edge] = EdgeColor::black;

					//If the other endpoint has not been added yet; add it to the index
					if (vertex_conversion.find(other_vertex_id) == vertex_conversion.end())
					{
						//Make a new supervertex 
						SuperVertex new_supervertex;
						new_supervertex.base_vertex = other_vertex_id;
						supervertices.push_back(new_supervertex);
						vertex_conversion[other_vertex_id] = supervertices.size() - 1;
						//Push the point to the queue
						Q.push(other_vertex_id);
					}

					//The degree is either 1 or >2
					//The superedge is done

					new_super_edge.endpoint_super_vertices = make_tuple(vertex_conversion[front], vertex_conversion[other_vertex_id]);
					NT super_edge_length = 0;
					for (size_t edg : new_super_edge.component_edges)
					{
						super_edge_length += (e_begin + edg)->length;
						superedge_ids[edg] = superedges.size();
					}
					new_super_edge.length = super_edge_length;
					superedges.push_back(new_super_edge);
					supervertices[vertex_conversion[front]].incident_superedges.push_back(superedges.size() - 1);
					supervertices[vertex_conversion[other_vertex_id]].incident_superedges.push_back(superedges.size() - 1);
					done = true;
				}
			}
		}
	}
};