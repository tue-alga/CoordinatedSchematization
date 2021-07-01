#pragma once
#include <vector>

template<class NT, class Vertex>
class VertexWrapper
{ 
public: 
	NT path_cost;
	NT heuristical_cost;
	NT total_cost;
	shared_ptr<Vertex> vert;
	size_t parent_edge_index;
	size_t parent_vertex_index;
	VertexWrapper(Vertex v, Vertex t)
	{
		vert = std::make_shared<Vertex>(v);
		path_cost = -1; //Infinity
		total_cost = -1;
		heuristical_cost = (v.lat - t.lat) * (v.lat - t.lat) + (v.lon - t.lon) * (v.lon - t.lon); //Heuristic is square of euclidean distance to t
		parent_edge_index = ~0; //bitwise NOT to get the max value;
		parent_vertex_index = ~0;
	}
};

template<class NT, class Edge, class Vertex >
class AStar
{
public:

	template <class EdgeIterator, class VertexIterator>
	std::vector<size_t>operator()(EdgeIterator e_begin, EdgeIterator e_beyond, VertexIterator v_begin, VertexIterator v_beyond, size_t s, size_t t)
	{
		//Get the s and t vertices
		auto vert_s = *(std::next(v_begin, s));
		auto vert_t = *(std::next(v_begin, t));
		//Start by wrapping all vertices in the wrapper class;
		std::vector<std::shared_ptr<VertexWrapper<NT,DummyVertex<NT>>>> wrapped_vertices;
		for (auto v_it = v_begin; v_it != v_beyond; v_it++)
		{
			wrapped_vertices.push_back(std::make_shared<VertexWrapper<NT,DummyVertex<NT>>>(*v_it, vert_t));
		}
		// Using lambda to compare elements.
		auto cmp = [](std::tuple<size_t, NT, size_t> left, std::tuple<size_t, NT,size_t> right) { return std::get<1>(left) < std::get<1>(right); };
		std::multiset<std::tuple<size_t, NT,size_t>, decltype(cmp)> ms(cmp);
		//Start with node s;
		auto wrapped_s = *(std::next(wrapped_vertices.begin(), s));
		wrapped_s->path_cost = 0;
		wrapped_s->total_cost = wrapped_s->heuristical_cost; //cost = h(s)+0;
		ms.insert(make_tuple(s, wrapped_s->total_cost,~0));
		std::tuple<size_t, NT, size_t> front;
		while (!ms.empty())
		{
			front = (*(ms.begin()));
			//Get the node
			auto wrapped_front = *(std::next(wrapped_vertices.begin(), get<0>(front)));
			ms.erase(ms.begin()); //Remove the node from the queue;
			//Check if this entry in the queue was still relevant;
			if (wrapped_front->path_cost != -1 && get<1>(front) > wrapped_front->total_cost)
			{
				continue;
			}
			//If the queue entry is relevant, we update the vertex and enqueue all others;
			//Update goal vertex
			//Get the edge
			auto prev_edge = *(std::next(e_begin, get<2>(front)));
			//Get the other index;
			auto other_v = get<0>(prev_edge.vertices);
			if (other_v == get<0>(front))
				other_v = get<1>(prev_edge.vertices);
			wrapped_front->parent_edge_index = get<2>(front);
			wrapped_front->parent_vertex_index = other_v;
			wrapped_front->total_cost = get<1>(front);
			wrapped_front->path_cost = wrapped_front->total_cost - wrapped_front->heuristical_cost;
			if (get<0>(front) == t) //If we have found t, we are done.
			{
				break;
			}
			for (auto inc_edges_it = wrapped_front->vert->edges.begin(); inc_edges_it != wrapped_front->vert->edges.end(); inc_edges_it++)
			{
				//Skip the edge if it was the edge that was just handled
				if (*inc_edges_it == get<2>(front))
					continue;
				//Get the relevant edge
				auto inc_edge = *(std::next(e_begin, *inc_edges_it));
				//Get the other node;
				auto other_v = get<0>(inc_edge.vertices);
				if (other_v == get<0>(front))
					other_v = get<1>(inc_edge.vertices);
				auto wrapped_other_v = *std::next(wrapped_vertices.begin(), other_v);
				//Check if we are faster
				if (wrapped_other_v->path_cost == -1 || wrapped_front->path_cost + inc_edge.length < wrapped_other_v->path_cost)
				{
					//Add a tuple for this connection to the priority queue
					ms.insert(make_tuple(other_v, wrapped_front->path_cost + inc_edge.length + wrapped_other_v->heuristical_cost, *inc_edges_it));
				} //end if (wrapped_other_v->path_cost == -1 || wrapped_front->path_cost + inc_edge.length < wrapped_other_v->path_cost)
			} //end for (auto inc_edges_it = wrapped_front->vert->edges.begin(); inc_edges_it != wrapped_front->vert->edges.end(); inc_edges_it++)
		} //end while (!ms.empty())
		//Check if we found t
		if (get<0>(front) != t)
			throw "Did not find a path from s to t, they likely lie in different graph components";

		std::shared_ptr<VertexWrapper<NT,DummyVertex<NT>>> curr_vert = *(std::next(wrapped_vertices.begin(), get<0>(front)));
		std::vector<size_t> reverse_path;
		while (curr_vert->parent_edge_index != ~0)
		{
			reverse_path.push_back(curr_vert->parent_edge_index);
			curr_vert = *(std::next(wrapped_vertices.begin(), curr_vert->parent_vertex_index));
		}
		std::reverse(reverse_path.begin(), reverse_path.end());
		return reverse_path;
	}
};