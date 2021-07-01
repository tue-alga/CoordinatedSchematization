#pragma once
#include <iostream>
#include <ostream>
#include <vector>
#include <queue>
namespace SchematLib::MeesGeneralization::SecondaryAlgs
{
 
template<class NT>
class ComponentCounter
{
    bool m_debugMessages = false;
public:
    struct Result
    {
        std::vector<std::size_t> componentids;
        std::size_t largest_component_idx = 0;
        std::size_t num_connected_components = 0;
        std::size_t num_leaf_nodes = 0;
        std::vector<size_t> leafnodes;
    };
    void setDebugMessages(bool value) { m_debugMessages = value; }
    bool debugMessages() const { return m_debugMessages; }

	//BFS to count the number of connected components and leaf vertices.
	template<typename Graph, class TakenIterator>
	Result operator()(const Graph& graph, //EdgeIterator begin, EdgeIterator beyond, 
        TakenIterator taken_begin, TakenIterator taken_end)
	{
        const char nl = '\n';
		if (m_debugMessages) std::cout << "[COMPONENT_COUNTER] Started." << nl;

        // Result object
        Result result;

		//ASSUMPTIONS ON INPUT:
		//Edges and vertices are sorted, their index is equal to their spot in the list
		//Edges are a tuple of two vertex indices
		//Vertices are a vector of incident edge indices
		result.componentids.resize(graph.num_edges(),0);

		//Create a color vector based on the taken vector
		std::vector<std::size_t> colors;
		for (auto bit = taken_begin; bit != taken_end; ++bit)
		{
            colors.push_back(*bit ? 1 : 0); //1= White, 0 = Black
		}
		//create the queue
		std::queue<size_t> Q;
		std::size_t e_it = 0;
		//Main loop, cycling through all edges
        std::size_t num_leaf_nodes = 0;
        std::size_t edge_count = 0;
        std::size_t largest_edge_count = 0;
		for (auto color : colors)
		{
            if (color != 1) // Not a white-colored edge
            {
                e_it++;
                continue;
            }
            //If the edge is white, it is the first edge we visit of a connected component
			++result.num_connected_components;
			edge_count = 0;
			Q.push(e_it);
			while (!Q.empty())
			{
				edge_count++;
				auto idx = Q.front();
				Q.pop();
				result.componentids[idx] = result.num_connected_components;
				//Get the edge so we can use the vertices
				*std::next(colors.begin(), idx) = 0; //color the edge black;
                auto edg = graph.edge(idx);
                auto v1 = graph.source(edg);
                auto v2 = graph.target(edg);
				std::vector<size_t> vertexlists{ v1,v2 };

                // Run over the vertices of the edge
				//for (auto vertex_list_it = vertexlists.begin(); vertex_list_it != vertexlists.end(); vertex_list_it++)
                for (auto v : vertexlists)
				{
					//auto vertex_list_sub_it = std::next(v_begin, *vertex_list_it);
					size_t incident_edge_count = 0;

					for (auto e : graph.out_edges(v))
					{
                        const auto eIndex = graph.getProperty(edge_index_t{}, e);
						auto edge_taken = *(std::next(taken_begin, eIndex));
                        if (!edge_taken) continue;
						
                        incident_edge_count++;
                        if (colors[eIndex] == 1) //If the edge hasn't been included in the BFS yet.
                        {
                            Q.push(eIndex);
                            colors[eIndex] = 2; //color the edge gray
                        }
					}
					if (incident_edge_count == 1)
					{
						num_leaf_nodes++;
						if (v1 == v)
						{
                            result.leafnodes.push_back(v1);
						}
                        else if (v2 == v)
                        {
                            result.leafnodes.push_back(v2);
                        }
                        else
                        {
                            throw std::runtime_error("vertex not found");
						}
					}
				}
			}
			if (edge_count > largest_edge_count)
			{
				largest_edge_count = edge_count;
				result.largest_component_idx = result.num_connected_components;
			}
            if (m_debugMessages) std::cout << "[COMPONENT_COUNTER] Found component with " << edge_count << " edges." << nl;
			e_it++;
		}
		if (m_debugMessages)
		{
			std::cout << "[COMPONENT_COUNTER] Finished." << nl;
            std::cout << "[COMPONENT_COUNTER] Connected components: " << result.num_connected_components << "." << nl;
            std::cout << "[COMPONENT_COUNTER] Leaf nodes: " << num_leaf_nodes << "." << std::endl;
		}
        result.num_leaf_nodes = num_leaf_nodes;
        return result;
	}
};
}