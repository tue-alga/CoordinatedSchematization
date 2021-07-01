#ifndef TULIB_GREEDYTRAJECTORYPICKER
#define TULIB_GREEDYTRAJECTORYPICKER
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "DS/DummyGraph.h"
#include "GraphInterface.h"

namespace SchematLib::MeesGeneralization
{
	////Declare class for use in EdgeWrapper
	//template<class NT, class EulerianTrajectory,class Edge>
	//class EulerianTrajectoryWrapper;

	//template<class Edge, class NT, class EulerianTrajectory>
	//class EdgeWrapper
	//{
	//public:
	//	Edge edge;
	//	bool taken;
	//	std::vector<std::shared_ptr<EulerianTrajectoryWrapper<NT,EulerianTrajectory,Edge>>> covered_trajectories;

	//	EdgeWrapper(Edge e) :edge(e),covered_trajectories() { taken = false; }
	//};

	////Full definition of EulerianTrajectoryWrapper
	////NB: Only make instances  of this class at the start of the algorithm, before any edges are selected. 
	//template<class NT, class EulerianTrajectory, class Edge>
	//class EulerianTrajectoryWrapper
	//{
	//public:
	//	std::vector<std::shared_ptr<EdgeWrapper<Edge,NT,EulerianTrajectory>>> covered_edges;
	//	NT score;
	//	NT unselected_length;
	//	NT length; //stored as reciprocal for faster divisions
	//	NT score_over_length;

	//	EulerianTrajectoryWrapper(EulerianTrajectory trajectory, std::vector<std::shared_ptr<EdgeWrapper<Edge,NT,EulerianTrajectory>>> edges)
	//	{
	//		score = 0;
	//		unselected_length = 0;
	//		score_over_length = 0;
	//		for (auto e_it = trajectory.begin(); e_it != trajectory.end(); e_it++)
	//		{
	//			auto wrapped_edge = std::next(edges.begin(), *e_it);
	//			unselected_length += (*wrapped_edge)->edge.length;
	//			covered_edges.push_back(*wrapped_edge);
	//			//wrapped_edge->covered_trajectories.push_back(*this);
	//		}
	//		length = 1 / unselected_length;
	//	}
	//};

    template<typename NT>
    struct TrajectoryStats
    {
        NT score = 0;
        NT unselected_length = 0;
        NT one_over_length = 0;
        NT score_over_length = 0;

        void recomputeScoreOverLength()
        {
            score_over_length = score * one_over_length; // because reciprocal...
        }

        template<typename GraphTrajectory, typename Graph>
        void computeFromTrajectory(const GraphTrajectory& trajectory, const Graph& graph)
        {
            score = 0;
            unselected_length = 0;
            score_over_length = 0;
            for (const auto& edgeIndex : trajectory)
            {
                //auto wrapped_edge = std::next(edges.begin(), *e_it);
                unselected_length += graph.getProperty(SchematLib::MeesGeneralization::edge_length_t{}, graph.edge(edgeIndex));
                //covered_edges.push_back(*wrapped_edge);
                //wrapped_edge->covered_trajectories.push_back(*this);
            }
            one_over_length = NT{ 1 }  / unselected_length;
        }
    };

	template <typename Graph, typename NT, typename EulerianTrajectory>
	class GreedyTrajectoryPicker
	{
        bool m_debugMessages = false;
        bool m_rescore = false;
        bool m_removeWorst = false;

        // Random slack parameters
        NT m_length_epsilon = 1; //Make sure we are not ignoring captured trajectories by doing a double check on all trajectories with unselected length less than a slack parameter
        NT m_score_epsilon = 0.0005; //Make sure rounding errors do not get us below 0 values for length and score;

        // Trajectories that each edge covers
        std::vector<std::vector<std::size_t>> m_coveredTrajectoriesPerEdge = {};
        // Edge that each trajectory covers
        std::vector<std::vector<std::size_t>> m_coveredEdgePerTrajectory = {};
        // The graph
        const Graph* m_graph = nullptr;
        // Statistics for the trajectories
        std::vector<TrajectoryStats<NT>> m_trajectoryStats = {};

		//bool compare_score_over_length(EulerianTrajectoryWrapper<NT,EulerianTrajectory,Edge> t1, EulerianTrajectoryWrapper<NT, EulerianTrajectory, Edge> t2) { return t1.score_over_length < t2.score_over_length; }
        template<typename...Ts>
        void debugPrint(Ts...args)
		{
            if (m_debugMessages) {
                (std::cout << ... << args);
            }
		}
        void applyRemoveWorst(std::vector<std::size_t>& scoredTrajectories) 
		{
            auto trajectoryIt = scoredTrajectories.rbegin();
            bool changed = false;
            // Find worst to delete
            for (; trajectoryIt != scoredTrajectories.rend(); ++trajectoryIt)
            {
                const auto trajectoryIndex = *trajectoryIt;
                if (m_trajectoryStats[trajectoryIndex].unselected_length <= m_length_epsilon) //Trajectory is already picked
                    continue;

                changed = true;
                for (const auto& edgeIndex : m_coveredEdgePerTrajectory[trajectoryIndex])
                {
                    const auto edgeLength = m_graph->getProperty(SchematLib::MeesGeneralization::edge_length_t{}, m_graph->edge(edgeIndex));

                    NT old_score_component = edgeLength / m_coveredTrajectoriesPerEdge[edgeIndex].size();
                    auto it = std::find(m_coveredTrajectoriesPerEdge[edgeIndex].begin(), m_coveredTrajectoriesPerEdge[edgeIndex].end(), *trajectoryIt);
                    const bool found = it != m_coveredTrajectoriesPerEdge[edgeIndex].end();

                    if(!found) throw std::runtime_error("Trajectory not found");

                    m_coveredTrajectoriesPerEdge[edgeIndex].erase(it);

                    NT new_score_component = edgeLength / m_coveredTrajectoriesPerEdge[edgeIndex].size();
                    for (const auto& updateTrajectoryIndex : m_coveredTrajectoriesPerEdge[edgeIndex])
                    {
                        m_trajectoryStats[updateTrajectoryIndex].score += (new_score_component - old_score_component); //Update the scores
                        if (std::abs(m_trajectoryStats[updateTrajectoryIndex].score) < m_score_epsilon)
                            m_trajectoryStats[updateTrajectoryIndex].score = 0;
                    }
                }
            }
            if (changed)
                scoredTrajectories.erase(trajectoryIt.base());
		}

        void runIterations(std::vector<std::size_t>& scoredTrajectories, NT& budget, std::unordered_set<std::size_t>& takenEdges)
		{
            if (debugMessages()) std::cout << "[TRAJECTORY_PICKER] Running iterations on " << scoredTrajectories.size() << " trajectories " << std::endl;
            NT spent_budget = 0;
            size_t iteration = 0;
            bool done = false;
            while (!done)
            {
                iteration++;
                //Pick the trajectory with the lowest score.
                done = true;

                // Run in order of increasing score-over-length
                for (const std::size_t& trajectoryIndex : scoredTrajectories)
                {
                    if (m_trajectoryStats[trajectoryIndex].score <= 0 || m_trajectoryStats[trajectoryIndex].unselected_length > budget)//Trajectories with score 0 are already taken
                        continue;

                    bool rounding_error = true; //There is a chance due to rounding that a trajectory is selected even though all edges have already been selected.
                        //update all covered edges to have status taken, which causes them to update all covered trajectories
                    for (const auto& edgeIndex : m_coveredEdgePerTrajectory[trajectoryIndex])
                    {
                        // Already taken, so ignore
                        if (takenEdges.find(edgeIndex) != takenEdges.end()) continue;

                        rounding_error = false;

                        // Take the edge
                        takenEdges.insert(edgeIndex);

                        NT taken_length = m_graph->getProperty(SchematLib::MeesGeneralization::edge_length_t{}, m_graph->edge(edgeIndex));
                        NT score_component = taken_length / m_coveredTrajectoriesPerEdge[edgeIndex].size();
                        budget -= taken_length;
                        spent_budget += taken_length;

                        // Update 
                        for (const auto& coveredTrajectoryIndex : m_coveredTrajectoriesPerEdge[edgeIndex])
                        {
                            m_trajectoryStats[coveredTrajectoryIndex].unselected_length -= taken_length;
                            if (std::abs(m_trajectoryStats[coveredTrajectoryIndex].unselected_length) < m_score_epsilon)
                                m_trajectoryStats[coveredTrajectoryIndex].unselected_length = 0;

                            m_trajectoryStats[coveredTrajectoryIndex].score -= score_component; //The score is updated, but the score_over_length, that we use to select trajectories, is only updated when asked
                            if (std::abs(m_trajectoryStats[coveredTrajectoryIndex].score) < m_score_epsilon)
                                m_trajectoryStats[coveredTrajectoryIndex].score = 0;
                        } 
                    } 

                    if (rounding_error)
                    {
                        m_trajectoryStats[trajectoryIndex].score = 0;//fix the rounding for future iterations
                        m_trajectoryStats[trajectoryIndex].unselected_length = 0;
                        continue; //try again with the next trajectory
                    } //end if (rounding_error)
                    done = false;

                    break;
                } // end for (auto traj_it = wrapped_trajectories.begin(); traj_it != wrapped_trajectories.end(), traj_it++)

                //Removes the worst scoring trajectory. Only use in combination with rescoring
                if (m_removeWorst)
                {
                    applyRemoveWorst(scoredTrajectories);
                }

                if (m_rescore)
                {
                    for(const auto& trajectoryIndex: scoredTrajectories)
                    {
                        m_trajectoryStats[trajectoryIndex].recomputeScoreOverLength();
                    }
                    //resort the trajectories
                    std::sort(scoredTrajectories.begin(), scoredTrajectories.end(), [this](std::size_t i0, std::size_t i1) {
                        return m_trajectoryStats[i0].score_over_length < m_trajectoryStats[i1].score_over_length;
                    });
                } //end if (rescore)

            }//end while (!done)
		}
	public:
        // Parameters
        bool debugMessages() const
        {
            return m_debugMessages;
        }
        void setDebugMessage(bool value)
        {
            m_debugMessages = value;
        }
        void setRescore(const bool& rescore)
        {
            m_rescore = rescore;
        }
        bool rescore() const
        {
            return m_rescore;
        }
        void setRemoveWorst(const bool& removeWorst)
        {
            m_removeWorst = removeWorst;
        }
        bool removeWorst() const
        {
            return m_removeWorst;
        }


		//Scorer scorer;

		//ASSUMPTIONS ON INPUT:
		//Nonneggative edge weights
		//Each edge id is equal to its sorted place in the list of edges
		//Each edge has at least 1 trajectory running on it
		//Eulerian Trajectories are iterators over a list of edge ids
		 template<typename TrajectoriesIterable> 
		std::vector<bool> operator()(const TrajectoriesIterable& trajectories, 
            const Graph& graph, NT budget)
		{
            // Set and reset internal state
            m_graph = &graph;
            m_coveredTrajectoriesPerEdge.clear();
            m_coveredEdgePerTrajectory.clear();
            m_trajectoryStats.clear();

            // 
            m_coveredTrajectoriesPerEdge.resize(graph.num_edges(),{});
            std::unordered_set<std::size_t> takenEdges;
            const auto numberOfTrajectories = std::distance(trajectories.begin(), trajectories.end());
            m_coveredEdgePerTrajectory.resize(numberOfTrajectories, {});
            m_trajectoryStats.resize(numberOfTrajectories, {});
            const auto numberOfEdges = graph.num_edges();

            debugPrint("[TRAJECTORY_PICKER] Started.\n");

            // The budgets and spent budget
			NT original_budget = budget; //Store the original budget so we can easily compute how much length we end up taking

            // Go over all trajectories
            std::size_t index = 0;
            for(const auto& trajectory : trajectories)
            {
                // Compute statistics
                m_trajectoryStats[index].computeFromTrajectory(trajectory, graph);
                // Mark covered by edges and edges covered by trajectory
                for (const auto& eIndex : trajectory)
                {
                    m_coveredTrajectoriesPerEdge[eIndex].push_back(index);
                    m_coveredEdgePerTrajectory[index].push_back(eIndex);
                }
                ++index;
            }
            
			//Score the trajectories
            for (const auto& edge : graph.edges()) //Let each edge assign its score component 
			{
                const auto edgeIndex = graph.getProperty(edge_index_t{},edge);
                NT edge_length = graph.getProperty(SchematLib::MeesGeneralization::edge_length_t{}, edge);
				if (edge_length < 0) throw std::runtime_error("negative edge length");
                NT score_component = edge_length / m_coveredTrajectoriesPerEdge[edgeIndex].size(); // Potentially contains duplicates
                if (score_component < 0) throw std::runtime_error("negative score component");

                for(const std::size_t& trajectoryIndex: m_coveredTrajectoriesPerEdge[edgeIndex])
                {
                    m_trajectoryStats[trajectoryIndex].score += score_component;
                }
			}
            for(std::size_t trajectoryIndex = 0; trajectoryIndex < numberOfTrajectories; ++trajectoryIndex)
            {
                m_trajectoryStats[trajectoryIndex].recomputeScoreOverLength();
            }

            // Order of worst to best score-over-length trajectories
            std::vector<std::size_t> scoredTrajectories(numberOfTrajectories,0);
            std::iota(scoredTrajectories.begin(), scoredTrajectories.end(), 0);

            std::sort(scoredTrajectories.begin(), scoredTrajectories.end(), [this](std::size_t index0, std::size_t index1) {
                return m_trajectoryStats[index0].score_over_length < m_trajectoryStats[index1].score_over_length;
            });

            // Select edges
            runIterations(scoredTrajectories, budget, takenEdges);
			
			//Return the selected edges and record the statistics that we want to know.
			std::vector<bool> results(m_graph->num_edges(),false);
			for (const auto& edge: m_graph->edges())
			{
                const auto edgeIndex = m_graph->getProperty(edge_index_t{}, edge);
                results[edgeIndex] = takenEdges.find(edgeIndex) != takenEdges.end();
			}
			
			size_t k = 0;
			NT total_captured_length = 0;
			for (const auto& trajectory: scoredTrajectories)
			{
				if (m_trajectoryStats[trajectory].unselected_length > 0 && m_trajectoryStats[trajectory].unselected_length < m_length_epsilon)
				{
					bool rounding_error = true;
					for (const auto& edge : m_coveredEdgePerTrajectory[trajectory])
					{
						if (takenEdges.find(edge) == takenEdges.end())
							rounding_error = false;
					}//end for (auto cov_edge_it = wrapped_trajectory_it->covered_edges.begin(); cov_edge_it != wrapped_trajectory_it->covered_edges.end(); cov_edge_it++)
					if (rounding_error)
					{
                        m_trajectoryStats[trajectory].unselected_length = 0;//fix the error
					}
				}//end if (wrapped_trajectory_it->unselected_length > 0 && wrapped_trajectory_it->unselected_length < slack)

				if (m_trajectoryStats[trajectory].unselected_length == 0)
				{
					k++; //total trajectories captured
					total_captured_length += 1 / m_trajectoryStats[trajectory].one_over_length; //total length captured
				}
			}//end for (auto wrapped_trajectory_it = wrapped_trajectories.begin(); wrapped_trajectory_it != wrapped_trajectories.end(); wrapped_trajectory_it++)
			if (m_debugMessages)
			{
				std::cout << "[TRAJECTORY_PICKER] captured " << k << " trajectories." << std::endl;
                std::cout << "[TRAJECTORY_PICKER] spent " << (original_budget - budget) << " budget." << std::endl;
			}
			return results;
		 }// end void operator()(EdgeIterator e_begin, EdgeIterator e_beyond, EulerianTrajectoryIterator begin, EulerianTrajectoryIterator beyond, NT budget)

	};// end class GreedyTrajectoryPicker

} //end namespace tulib_algorithms

#endif