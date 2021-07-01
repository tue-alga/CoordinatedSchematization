#pragma once
#include <vector>
#include <set>
#include <queue>
#include <optional>

template <class NT>
class ComputeCentrality
{
public:
	enum NodeColor{white, gray, black};

	ComputeCentrality() {}

	template <class SuperEdgeIterator, class SuperVertexIterator>
	void operator()(SuperVertexIterator s_v_begin, SuperVertexIterator s_v_beyond, SuperEdgeIterator s_e_begin, SuperEdgeIterator s_e_beyond, std::vector<int> stroke_ids, size_t num_strokes, std::vector<size_t> &covered_counts, std::vector<NT> &edge_cc)
	{
		//We assume the entire graph is a single connected component
		//Initialize our lists
		covered_counts = std::vector<size_t>(num_strokes, 0);
		std::vector<NT> vertex_cc;
		//std::queue<size_t> Q;
		std::vector<NT> distances;
		std::vector<std::optional<size_t>> prev_edge;
		std::vector<NodeColor> node_colors;
		std::priority_queue < std::tuple <NT, size_t, std::optional<size_t>>, std::vector<std::tuple<NT, size_t, std::optional<size_t>>>, std::greater<std::tuple<NT, size_t, std::optional<size_t>>> > Q;

		//Do a full Dijkstra pass starting at every single vertex
		for (auto bfs_vert_it = s_v_begin; bfs_vert_it != s_v_beyond; bfs_vert_it++)
		{

			//Start by clearing the values for this node.
			distances = std::vector<NT>(s_v_beyond - s_v_begin, -1);
			prev_edge = std::vector<std::optional<size_t>>(s_v_beyond - s_v_begin, std::optional<size_t>{});
			node_colors = std::vector<NodeColor>(s_v_beyond - s_v_begin, NodeColor::white);

			auto curr_bfs_vert_id = bfs_vert_it - s_v_begin;

			distances[curr_bfs_vert_id] = 0;
			node_colors[curr_bfs_vert_id] = NodeColor::gray;

			Q.push(make_tuple(0, curr_bfs_vert_id, std::optional<size_t>{}));
			while (Q.size() > 0)
			{
				auto dist = get<0>(Q.top());
				auto idx = get<1>(Q.top());
				auto prev_ref = get<2>(Q.top());
				Q.pop();
				//Check if the node still needs to be treated
				if (node_colors[idx] == NodeColor::black)
					continue;

				//Color the node black so it does not get treated again
				node_colors[idx] = NodeColor::black;
				//Lock in the data
				distances[idx] = dist;
				prev_edge[idx] = prev_ref;

				//Get the vertex
				auto vert = *(s_v_begin + idx);
				//Go over all incident edges
				for (auto inc_edge_it = vert.incident_superedges.begin(); inc_edge_it != vert.incident_superedges.end(); inc_edge_it++)
				{
					//Get the edge
					auto edge = *(s_e_begin + *inc_edge_it);
					//Get the other vertex
					size_t other_vertex_id;
					if (get<0>(edge.endpoint_super_vertices) == idx)
						other_vertex_id = get<1>(edge.endpoint_super_vertices);
					else if (get<1>(edge.endpoint_super_vertices) == idx)
						other_vertex_id = get<0>(edge.endpoint_super_vertices);
					else throw "bad edge.";

					//Check if the node is not black
					if (node_colors[other_vertex_id] == NodeColor::black)
						continue;
					//Color the node gray and push it
					node_colors[other_vertex_id] = NodeColor::gray;
					Q.push(make_tuple(distances[idx] + edge.length, other_vertex_id, *inc_edge_it));
				}

			}

			//After the queue is empty and thus the entire graph has been processed for this vertex, we can update the metrics
			//First the CC for the vertex
			NT sum_distances = 0;
			for (auto dist_it = distances.begin(); dist_it != distances.end(); dist_it++)
				if (*dist_it > 0) //Some distances may still be -1 due to the nodes not being part of the connected component
					sum_distances += *dist_it;
			vertex_cc.push_back(sum_distances);
			//For each vertex with index > the current (so we don't count every path twice), walk the shortest path and add to the counts
			for (auto end_vertex_it = bfs_vert_it + 1; end_vertex_it != s_v_beyond; end_vertex_it++)
			{
				std::set<size_t> touched_strokes;
				auto curr_vert_id = distance(s_v_begin, end_vertex_it);
				while (prev_edge[curr_vert_id])
				{
					touched_strokes.insert(stroke_ids[*(prev_edge[curr_vert_id])]);
					//covered_counts[covered_counts[*prev_edge]]++;
					auto edge = *(s_e_begin + *(prev_edge[curr_vert_id]));
					if (get<0>(edge.endpoint_super_vertices) == curr_vert_id)
						curr_vert_id = get<1>(edge.endpoint_super_vertices);
					else if (get<1>(edge.endpoint_super_vertices) == curr_vert_id)
						curr_vert_id = get<0>(edge.endpoint_super_vertices);
					else throw "bad edge";
				}
				if (curr_vert_id != distance(s_v_begin, bfs_vert_it))
					throw "failed during backtracing";
				for (size_t stroke_id : touched_strokes)
					covered_counts[stroke_id]++;
			}
		}

		//Compute the edge ccs
		for (auto s_e_it = s_e_begin; s_e_it != s_e_beyond; s_e_it++)
		{
			auto vert_ids = s_e_it->endpoint_super_vertices;
			edge_cc.push_back(2 / (vertex_cc[get<0>(vert_ids)] + vertex_cc[get<1>(vert_ids)]));
		}

	}
	
};