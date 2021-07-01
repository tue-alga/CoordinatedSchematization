#ifndef MOVETK_MOSTCOVEREDHEURISTIC_H
#define MOVETK_MOSTCOVEREDHEURISTIC_H
#include <numeric>
#include <set>


#include "SecondaryAlgs/RangeIterator.h"
#include "GraphInterface.h"

template <typename Graph, class NT>
class MostCoveredHeuristic
{
    template<typename GraphTrajectoryIterable>
    void computeEdgePriorities(const GraphTrajectoryIterable& trajectories, const Graph& graph, std::vector<std::size_t>& edgeIndices)
    {
        const std::size_t numberOfEdges = graph.num_edges();
        std::vector<size_t> edge_counts(numberOfEdges, 0);
        //Go over all trajectories, and count for each edge how often it occurs
        for (const auto& traj : trajectories)
        {
            for (const auto& e : traj) //Should be iterable
                edge_counts[e] += 1;
        }

        //Create a sorted list of edges, based on drive counts
        edgeIndices.resize(numberOfEdges);
        // Fill with indices for edges
        std::iota(edgeIndices.begin(), edgeIndices.end(), 0);
        // Rank the indices
        std::sort(edgeIndices.begin(), edgeIndices.end(), [&edge_counts, &graph](std::size_t i0, std::size_t i1)
        {
            // Prefer smaller edges when tied
            if (edge_counts[i0] == edge_counts[i1]) {
                return graph.getProperty(SchematLib::MeesGeneralization::edge_length_t{}, graph.edge(i0)) < graph.getProperty(SchematLib::MeesGeneralization::edge_length_t{}, graph.edge(i1));
            }
            return edge_counts[i0] > edge_counts[i1];
        });
    }

public:

	//We assume the Eulerian trajectories consist of a container of indices of the edge container
	template<class EulerianTrajectoryIterator>
	std::set<std::size_t> operator()(EulerianTrajectoryIterator trajectories_begin, EulerianTrajectoryIterator trajectories_beyond, const Graph& graph, NT budget)
	{
        using namespace SchematLib::MeesGeneralization::SecondaryAlgs;
        //Create a sorted list of edges, based on drive counts
        std::vector<std::size_t> edgeIndices;

        computeEdgePriorities(range(trajectories_begin, trajectories_beyond), graph, edgeIndices);

        // Find the selection that fits the budget
		std::set<std::size_t> selection;
		NT leftover_budget = budget;
        auto last = std::find_if(edgeIndices.begin(), edgeIndices.end(), [&leftover_budget, &graph](std::size_t i)
        {
            const auto budget_required = graph.getProperty(SchematLib::MeesGeneralization::edge_length_t{}, graph.edge(i));
            if(budget_required < leftover_budget)
            {
                leftover_budget -= budget_required;
                return false;
            }
            return true;
        });
        // First edge is already violating.
        if (last == edgeIndices.begin() && graph.getProperty(SchematLib::MeesGeneralization::edge_length_t{}, graph.edge(0)) > leftover_budget) return selection;

        selection.insert(edgeIndices.begin(), last);

		return selection;
	}
};

#endif //MOVETK_MOSTCOVEREDHEURISTIC_H
