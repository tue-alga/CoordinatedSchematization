#ifndef SCHEMATLIB_TRAJECTORYSCHEMATIZATION_TRAJECTORYBUNDLER_H
#define SCHEMATLIB_TRAJECTORYSCHEMATIZATION_TRAJECTORYBUNDLER_H
#include <optional>
#include <queue>
#include <set>
#include <type_traits>
#include <tuple>
#include <vector>
#include <GCpp/DS/BoostEmbeddedGraph.h>
#include <GCpp/DS/Heap.h>
#include "TrajectoryBundleTraits.h"

namespace SchematLib::TrajectorySchematization
{
    template<typename EmbeddedGraph>
    class TrajectoryBundler
    {
        std::size_t m_minimum_bundle_element_complexity = 3;

        std::size_t m_max_bundle_complexity = 10;

        Models::NT m_minimum_bundle_length = 1000;//meters assumed

        double m_bundle_complexity_weight = 1.0;

        double m_bundle_weight_weight = 1.0;

        std::size_t m_min_matched_trajectories = 10;

        struct TrajectoryEdge
        {
            std::size_t trajectory;
            std::size_t edgeNum;
        };

        using Edge = typename EmbeddedGraph::edge_descriptor;
        using GT = boost::graph_traits<EmbeddedGraph>;
        using Vertex = typename EmbeddedGraph::vertex_descriptor;
        struct PathEdge
        {
            enum class Explore
            {
                Source, Target, Both
            };
            Edge edge;
            Explore toExplore;
            Vertex vertexToExplore(const EmbeddedGraph& graph) const
            {
                if (toExplore == Explore::Source) return boost::source(edge, graph);
                if (toExplore == Explore::Target) return boost::target(edge, graph);
                throw std::runtime_error("Vertex to explore only works for Source/Target");
            }

        };
        // Semantically tag the weight to avoid confusion with IDs of std::size_t
        struct EdgeWeight
        {
            std::size_t value;
            bool operator<(const EdgeWeight& other) const
            {
                return value < other.value;
            }
            bool operator>(const EdgeWeight& other) const
            {
                return value > other.value;
            }
            bool operator==(const EdgeWeight& other) const
            {
                return value == other.value;
            }
        };


        struct ComputationData
        {
            std::map<std::size_t, EdgeWeight> edgeWeights;
            std::map<std::size_t, std::vector<TrajectoryEdge>> edgeIdToTrajectoryEdges;
            std::map<std::size_t, std::set<std::size_t>> edgeIdToTrajectoryNum;
            std::map<std::size_t, typename GT::edge_descriptor> edgeIdToEdgeMap;
            const EmbeddedGraph& graph;
            ComputationData(const EmbeddedGraph& inputGraph) :graph(inputGraph) {}

            Models::NT edgeLength(std::size_t edgeId) const
            {
                return GCpp::DS::get_linear_edge_length(edgeIdToEdgeMap.at(edgeId), graph);
            }

            void precomputeComputationData(const std::vector<std::vector<std::size_t>>& trajectories,
                const std::unordered_set<std::size_t>& edgeIds
            )
            {
                std::size_t trajectoryNum = 0;
                for (const auto& traj : trajectories)
                {
                    std::size_t currentEdge = 0;
                    for (auto e : traj)
                    {
                        if (edgeIdToTrajectoryNum.find(e) == edgeIdToTrajectoryNum.end())
                        {
                            edgeIdToTrajectoryNum[e] = {};
                        }
                        edgeIdToTrajectoryNum[e].insert(trajectoryNum);
                        if (edgeIds.find(e) == edgeIds.end()) {
                            ++currentEdge;
                            continue;
                        }
                        if (edgeWeights.find(e) == edgeWeights.end()) edgeWeights[e] = EdgeWeight{ 0 };
                        ++edgeWeights[e].value;
                        if (edgeIdToTrajectoryEdges.find(e) == edgeIdToTrajectoryEdges.end()) edgeIdToTrajectoryEdges[e] = {};
                        edgeIdToTrajectoryEdges[e].push_back(TrajectoryEdge{ trajectoryNum, currentEdge });
                        ++currentEdge;
                    }
                    ++trajectoryNum;
                }
                GCpp::DS::computeEdgeIdToEdgeMap(graph, edgeIdToEdgeMap);
            }

            void updatePriority(
                const std::vector<std::size_t>& edges, const std::vector<std::size_t>& trajectories,
                GCpp::DS::Heap<EdgeWeight, std::size_t, std::greater<>>& priorityQueue)
            {
                for (auto e : edges)
                {
                    const auto newData = EdgeWeight{ priorityQueue.dataForId(e).value - trajectories.size() };
                    priorityQueue.updateKey(e, newData);
                    edgeWeights[e] = newData;
                }
            }

            /**
             * \brief Update edge ID to trajectories mapping by removing the given trajectories for the given edges.
             * \param edges The dges
             * \param trajectories The trajectories.
             */
            void updateEdgeIdToTrajectoryMapping(const std::vector<std::size_t>& edges, const std::vector<std::size_t>& trajectories)
            {
                for (auto edge : edges)
                {
                    auto& trajectorySet = edgeIdToTrajectoryNum[edge];
                    for (auto traj : trajectories)
                    {
                        trajectorySet.erase(traj);
                    }
                }
            }
        };

        /**
         * \brief Compute the score for the bundle, as defined by number of trajectories and complexity
         * \param matchedTrajectories 
         * \param complexity 
         * \return 
         */
        double computeScore(std::size_t matchedTrajectories, std::size_t complexity) const
        {
            return m_bundle_complexity_weight * static_cast<double>(complexity) + m_bundle_weight_weight * static_cast<double>(matchedTrajectories);
        }

       

        static bool isConnected(typename EmbeddedGraph::vertex_descriptor v, typename EmbeddedGraph::edge_descriptor e, const EmbeddedGraph& g)
        {
            return boost::source(e, g) == v || boost::target(e, g) == v;
        }
        static std::size_t getEdgeIdforConnection(typename EmbeddedGraph::vertex_descriptor v, typename EmbeddedGraph::edge_descriptor e, const EmbeddedGraph& g)
        {
            assert(isConnected(v, e, g));
            if(boost::source(e,g) == v)
            {
                auto[e, exists] = boost::lookup_edge(v, boost::source(e, g), g);
                if (!exists) std::tie(e, exists) = boost::lookup_edge(boost::source(e, g), v, g);
                assert(exists);
                return GCpp::DS::getEdgeId(g, e);
            }
            else
            {
                auto[e, exists] = boost::lookup_edge(v, boost::target(e, g), g);
                if (!exists) std::tie(e, exists) = boost::lookup_edge(boost::target(e, g), v, g);
                assert(exists);
                return GCpp::DS::getEdgeId(g, e);
            }
        }

        /**
         * \brief Search for a bundle starting with the edge, given as an ID. Outputs bundle, its score and matched trajectories. Returns whether
         * a bundle was found
         * \param edgeID ID of the edge to start at
         * \param data Computation data that is needed
         * \param graph The graph 
         * \param outputBundle The output bundle
         * \param outputScore The output score
         * \param matchedTrajectories The output matched trajectories.
         * \return Whether a bundle was found satisfying the criteria.
         */
        bool searchForBundle(std::size_t edgeID,
            const ComputationData& data,
            TrajectoryBundle& output,
            std::vector<std::size_t>& matchedTrajectories
            ) const
        {
            namespace it = GCpp::Helpers::Iterators;
            std::list<PathEdge> path;
            std::unordered_set<std::size_t> usedEdges;
            std::set<typename EmbeddedGraph::vertex_descriptor> exploredVertices;
            std::set<std::size_t> currentTrajectoriesMatched = data.edgeIdToTrajectoryNum.at(edgeID);
            
            // Heuristically define when to continue searching: we need atleast some minimum length path,
            // but we do not want to indefinitely increase it if the matched trajectories drastically decreases.
            if (data.edgeIdToTrajectoryNum.find(edgeID) == data.edgeIdToTrajectoryNum.end()) return false;

            path.push_back(PathEdge{data.edgeIdToEdgeMap.at(edgeID), PathEdge::Explore::Both});

            // Track the path length
            auto pathLength = data.edgeLength(edgeID);

            //Initial orientation does not matter
            usedEdges.insert(edgeID);
            double currentScore = computeScore(currentTrajectoriesMatched.size(), 1);

            // Initialize
            std::array<Vertex, 2> toExplore{
                boost::source(path.front().edge, data.graph),
                boost::target(path.front().edge, data.graph)
            };
            // Keep going while our score improves.
            while(true)
            {
                // Look at both ends of the path, going over the vertex edges and scoring them.
                double maxScore = 0;
                std::optional<Edge> bestEdge;
                std::optional<int> bestSide;
                int side = 0; //Side that we are expanding: 0=start of path, 1=end.
                for(auto v : toExplore)
                {
                    for(auto e : GCpp::Helpers::Iterators::range(boost::out_edges(v, data.graph)))
                    {
                        auto eID = GCpp::DS::getEdgeId(data.graph, e);
                        // Skip used edges
                        if (usedEdges.find(eID) != usedEdges.end()) continue;
                        // Skip edges without trajectories.
                        if (data.edgeIdToTrajectoryNum.find(eID) == data.edgeIdToTrajectoryNum.end()) continue;

                        // Determine improvement
                        auto set = data.edgeIdToTrajectoryNum.at(eID);
                        std::set<std::size_t> overlap;
                        std::set_intersection(set.begin(), set.end(), currentTrajectoriesMatched.begin(), currentTrajectoriesMatched.end(), std::inserter(overlap, std::begin(overlap)));

                        if (overlap.size() < m_min_matched_trajectories) continue;

                        auto score = computeScore(overlap.size(), path.size() + 1);
                        if(score > maxScore)
                        {
                            maxScore = score;
                            bestEdge = e;
                            bestSide = side;
                        }
                    }
                    ++side;
                }
                // Only stop when minimum size is reached and we are not improving.
                if(currentScore > maxScore && path.size() > m_minimum_bundle_element_complexity)
                {
                    break;
                }
                if (!bestEdge.has_value()) break;
                // Update score
                currentScore = maxScore;

                // Mark edge as used
                usedEdges.insert(GCpp::DS::getEdgeId(data.graph, bestEdge.value()));

                //Update path
                auto compareV = toExplore[bestSide.value()];
                PathEdge pathEdge;
                pathEdge.edge = bestEdge.value();
                pathEdge.toExplore = boost::source(bestEdge.value(), data.graph) == compareV ? PathEdge::Explore::Target : PathEdge::Explore::Source;

                // Update path and exploration vertices.
                if (bestSide.value() == 0) {
                    path.push_front(pathEdge);
                    toExplore[0] = pathEdge.vertexToExplore(data.graph);
                }
                else {
                    path.push_back(pathEdge);
                    toExplore[1] = pathEdge.vertexToExplore(data.graph);
                }

                // Update matched set of trajectories.
                auto set = data.edgeIdToTrajectoryNum.at(GCpp::DS::getEdgeId(data.graph, pathEdge.edge));
                std::set<std::size_t> overlap;
                std::set_intersection(set.begin(), set.end(), currentTrajectoriesMatched.begin(), currentTrajectoriesMatched.end(), std::inserter(overlap, std::begin(overlap)));
                currentTrajectoriesMatched = overlap;
            }
            output.matchedTrajectoriesCount= matchedTrajectories.size();
            output.edges.clear();
            for(auto el : path)
            {
                output.edges.push_back(GCpp::DS::getEdgeId(data.graph, el.edge));
            }
            output.matchedTrajectoriesCount = matchedTrajectories.size();
            // Assign the matched trajectories.
            matchedTrajectories.insert(matchedTrajectories.begin(), currentTrajectoriesMatched.begin(), currentTrajectoriesMatched.end());
            // Succesful path of required minimum length found.
            return path.size() >= m_minimum_bundle_element_complexity && matchedTrajectories.size() >= m_min_matched_trajectories;
        }
    public:
#pragma region Parameters
        void set_minimum_bundle_element_complexity(const std::size_t& minimum_bundle_element_complexity)
        {
            m_minimum_bundle_element_complexity = minimum_bundle_element_complexity;
        }
        std::size_t get_minimum_bundle_element_complexity() const
        {
            return m_minimum_bundle_element_complexity;
        }
        void set_max_bundle_complexity(const std::size_t& max_bundle_complexity)
        {
            m_max_bundle_complexity = max_bundle_complexity;
        }
        std::size_t get_max_bundle_complexity() const
        {
            return m_max_bundle_complexity;
        }
        void set_min_matched_trajectories(const std::size_t& min_match)
        {
            m_min_matched_trajectories = min_match;
        }
        std::size_t get_min_matched_trajectories() const
        {
            return m_min_matched_trajectories;
        }
#pragma endregion

        /**
         * \brief
         * \param graph Input simplified graph
         * \param trajectories Given as list of edge IDs.
         * \param output Vector of trajectory bundles, that contain score and bundle as list of edge IDs
         */
        void operator()(const EmbeddedGraph& graph, const std::vector<std::vector<std::size_t>>& trajectories, std::vector<TrajectoryBundle>& output) const
        {
            std::unordered_set<std::size_t> edgeIds;
            GCpp::DS::computeEdgeIdSet(graph, edgeIds);
            // Compute data
            ComputationData data(graph);
            data.precomputeComputationData(trajectories, edgeIds);
            
            GCpp::DS::Heap<EdgeWeight, std::size_t, std::greater<>> edgePriority;
            for (const auto& kv : data.edgeWeights)
            {
                edgePriority.insert(kv.second, kv.first);
            }
            // Perform a greedy search for bundle elements
            std::vector<typename decltype(edgePriority)::Node> poppedElements;
            while (output.size() < m_max_bundle_complexity)
            {
                auto topNode = edgePriority.extract();
                poppedElements.push_back(topNode);

                if (topNode.value.value < m_min_matched_trajectories) break;

                TrajectoryBundle bundle;
                std::vector<std::size_t> matchedTrajectories;
                if (searchForBundle(topNode.id,data, bundle, matchedTrajectories))
                {
                    output.push_back(bundle);
                    // Update weights
                    // TODO: delete all of trajectory?
                    data.updateEdgeIdToTrajectoryMapping(bundle.edges, matchedTrajectories);
                    for(auto el : poppedElements)
                    {
                        edgePriority.insert(el.value, el.id);
                    }
                    poppedElements.clear();
                    data.updatePriority(bundle.edges, matchedTrajectories, edgePriority);
                }
                else // Maybe it becomes better later?
                {
                    poppedElements.push_back(topNode);
                }
            }
            std::cout << "Found " << output.size() << "/" << m_max_bundle_complexity << " bundles \n";
        }
    };
}
#endif
