#ifndef SCHEMATLIB_TRAJECTORYSCHEMATIZATION_BIDIRECTIONALTRAJECTORYBUNDLER_H
#define SCHEMATLIB_TRAJECTORYSCHEMATIZATION_BIDIRECTIONALTRAJECTORYBUNDLER_H
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
    class BidirectionalTrajectoryBundler
    {
        std::size_t m_minimum_bundle_element_complexity = 3;

        std::size_t m_max_bundle_set_complexity = 10;

        Models::NT m_minimum_bundle_length = 1000;//meters assumed

        std::size_t m_min_matched_trajectories = 10;

        bool m_priotizie_by_length = false;

    private:
        struct TrajectoryEdge
        {
            std::size_t trajectory;
            std::size_t edgeNum;
            bool operator<(const TrajectoryEdge& other) const
            {
                return trajectory == other.trajectory ? edgeNum < other.edgeNum : trajectory < other.trajectory;
            }
            bool operator==(const TrajectoryEdge& other) const
            {
                return trajectory == other.trajectory && edgeNum == other.edgeNum;
            }
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
        enum class Direction
        {
            Forward = 0,
            Backward = 1
        };
        /**
         * \brief Matchre object for a trajectory, in a particular direction.
         */
        struct TrajectoryMatch
        {
            std::size_t trajectoryIndex;
            std::size_t startEdgeIndex;
            std::size_t endEdgeIndex;

            bool matches(Direction direction, std::size_t edgeId, const std::vector<std::vector<std::size_t>>& trajectories) const
            {
                if (direction == Direction::Backward)
                {
                    if (startEdgeIndex == 0) return false;
                    return trajectories[trajectoryIndex][startEdgeIndex - 1] == edgeId;
                }
                else
                {
                    if (endEdgeIndex == trajectories[trajectoryIndex].size() - 1) return false;
                    return trajectories[trajectoryIndex][endEdgeIndex + 1] == edgeId;
                }
            }
            void extend(Direction direction)
            {
                if (direction == Direction::Backward) --startEdgeIndex;
                else ++endEdgeIndex;
            }
        };

        using MatchingMap = std::map<std::size_t, TrajectoryMatch>;

        struct ComputationData
        {
            // Edge IDs of the given graph
            std::unordered_set<std::size_t> edgeIds;
            
            std::map<std::size_t, EdgeWeight> edgeWeights;
            // Indexed by edgeID, then trajectory index.
            std::map<std::size_t, std::map<std::size_t,TrajectoryEdge>> edgeIdToTrajectoryEdges;
            std::map<std::size_t, std::set<std::size_t>> edgeIdToTrajectoryNum;
            std::map<std::size_t, typename GT::edge_descriptor> edgeIdToEdgeMap;
            const EmbeddedGraph& graph;
            ComputationData(const EmbeddedGraph& inputGraph) :graph(inputGraph)
            {
                GCpp::DS::computeEdgeIdSet(graph, edgeIds);
            }

            void precomputeComputationData(const std::vector<std::vector<std::size_t>>& trajectories)
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
                        edgeIdToTrajectoryEdges[e][trajectoryNum] = TrajectoryEdge{ trajectoryNum, currentEdge };
                        ++currentEdge;
                    }
                    ++trajectoryNum;
                }
                GCpp::DS::computeEdgeIdToEdgeMap(graph, edgeIdToEdgeMap);
            }

            void updatePriority(
                const std::vector<std::size_t>& edges, const MatchingMap& trajectories,
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
            void updateEdgeIdToTrajectoryMapping(const std::vector<std::size_t>& edges, const MatchingMap& trajectories)
            {
                for (auto edge : edges)
                {
                    auto& trajectorySet = edgeIdToTrajectoryNum[edge];
                    auto& edgeToTrajectoryEdges = edgeIdToTrajectoryEdges[edge];
                    for (const auto& [trajId, trajMatching] : trajectories)
                    {
                        trajectorySet.erase(trajId);
                        edgeToTrajectoryEdges.erase(trajId);
                    }
                }
            }
        };
        
        static bool isConnected(typename EmbeddedGraph::vertex_descriptor v, typename EmbeddedGraph::edge_descriptor e, const EmbeddedGraph& g)
        {
            return boost::source(e, g) == v || boost::target(e, g) == v;
        }
        static std::size_t getEdgeIdforConnection(typename EmbeddedGraph::vertex_descriptor v, typename EmbeddedGraph::edge_descriptor e, const EmbeddedGraph& g)
        {
            assert(isConnected(v, e, g));
            if (boost::source(e, g) == v)
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
            bool flippedInitial,
            const std::vector<std::vector<std::size_t>>& mmTrajectories,
            TrajectoryBundle& bundle,
            MatchingMap& matchedTrajectories
        ) const
        {
            if (data.edgeIdToTrajectoryNum.find(edgeID) == data.edgeIdToTrajectoryNum.end()) {
                std::cout << "\tNo trajectories\n";
                return false;
            }

            namespace it = GCpp::Helpers::Iterators;
            // The subpath that matches
            std::list<PathEdge> path;
            std::unordered_set<std::size_t> usedEdges;
            Models::NT pathLength = 0;

            // Initial matched trajectories.
            for(auto el: data.edgeIdToTrajectoryNum.at(edgeID))
            {
                const auto edgeIndex = data.edgeIdToTrajectoryEdges.at(edgeID).at(el).edgeNum;
                matchedTrajectories[el] = TrajectoryMatch{ el, edgeIndex, edgeIndex };
            }

            // Heuristically define when to continue searching: we need atleast some minimum length path,
            // but we do not want to indefinitely increase it if the matched trajectories drastically decreases.
            path.push_back(PathEdge{ data.edgeIdToEdgeMap.at(edgeID), PathEdge::Explore::Both });
            pathLength += GCpp::DS::get_linear_edge_length(data.edgeIdToEdgeMap.at(edgeID), data.graph);

            //Initial orientation does not matter
            usedEdges.insert(edgeID);

            // Initialize
            std::array<Vertex, 2> toExplore{
                boost::source(path.front().edge, data.graph),
                boost::target(path.front().edge, data.graph)
            };

            if (flippedInitial) std::swap(toExplore[0], toExplore[1]);

            const Direction directions[2]{ Direction::Backward, Direction::Forward };

            // Keep going while our score improves.
            while (true)
            {
                // Look at both ends of the path, going over the vertex edges and scoring them.
                std::size_t maxMatchedTrajectoriesCount = 0;
                std::optional<Edge> bestEdge;
                std::optional<int> bestSide;
                Models::NT bestEdgeLength = 0;
                int side = 0; //Side that we are expanding: 0=start of path, 1=end.
                for (auto v : toExplore)
                {
                    for (auto e : GCpp::Helpers::Iterators::range(boost::out_edges(v, data.graph)))
                    {
                        auto eID = GCpp::DS::getEdgeId(data.graph, e);
                        // Skip used edges
                        if (usedEdges.find(eID) != usedEdges.end()) continue;
                        // Skip edges without trajectories.
                        if (data.edgeIdToTrajectoryNum.find(eID) == data.edgeIdToTrajectoryNum.end()) continue;

                        const auto eLength = GCpp::DS::get_linear_edge_length(data.edgeIdToEdgeMap.at(eID), data.graph);

                        // Count number of match when extending via this edge.
                        std::size_t totalSubMatchCount = 0;
                        for(const auto& kv: matchedTrajectories)
                        {
                            if (!kv.second.matches(directions[side], eID, mmTrajectories)) continue;

                            ++totalSubMatchCount;
                        }

                        if (totalSubMatchCount < m_min_matched_trajectories) continue;

                        if(m_priotizie_by_length)
                        {
                            if (eLength > bestEdgeLength)
                            {
                                maxMatchedTrajectoriesCount = totalSubMatchCount;
                                bestEdge = e;
                                bestSide = side;
                                bestEdgeLength = eLength;
                            }
                        }
                        else
                        {
                            if (totalSubMatchCount > maxMatchedTrajectoriesCount)
                            {
                                maxMatchedTrajectoriesCount = totalSubMatchCount;
                                bestEdge = e;
                                bestSide = side;
                                bestEdgeLength = eLength;
                            }
                        }
                        
                    }
                    ++side;
                }
                // Only stop when minimum size is reached and we are violating the minimum match count
                if ( pathLength > m_minimum_bundle_length && maxMatchedTrajectoriesCount < m_min_matched_trajectories)
                {
                    break;
                }
                if (!bestEdge.has_value()) break;

                pathLength += bestEdgeLength;

                const auto eID = GCpp::DS::getEdgeId(data.graph, bestEdge.value());

                // Mark edge as used
                usedEdges.insert(eID);
                
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
                std::set<std::size_t> toErase;
                for (auto& kv : matchedTrajectories)
                {
                    if (!kv.second.matches(directions[bestSide.value()], eID, mmTrajectories))
                    {
                        toErase.insert(kv.first);
                        continue;
                    }
                    auto& val = kv.second;
                    val.extend(directions[bestSide.value()]);
                }
                for (auto el : toErase) matchedTrajectories.erase(el);
            }
            bundle.matchedTrajectoriesCount = matchedTrajectories.size();
            bundle.edges.clear();
            bundle.bundleLength = pathLength;
            for (auto el : path)
            {
                bundle.edges.push_back(GCpp::DS::getEdgeId(data.graph, el.edge));
            }
            std::cout << "\tBundle with " << bundle.matchedTrajectoriesCount << " matched, length " << pathLength << '\n';
            // Succesful path of required minimum length found.
            return pathLength >= m_minimum_bundle_length && matchedTrajectories.size() >= m_min_matched_trajectories;
        }
    public:
#pragma region Parameters
        void set_priotizie_by_length(bool val)
        {
            m_priotizie_by_length = val;
        }
        bool get_priotizie_by_length() const
        {
            return m_priotizie_by_length;
        }
        void set_minimum_bundle_length(const Models::NT& minimum_bundle_length)
        {
            m_minimum_bundle_length = minimum_bundle_length;
        }
        Models::NT get_minimum_bundle_length() const
        {
            return m_minimum_bundle_length;
        }
        void set_max_bundle_set_complexity(const std::size_t& max_bundle_complexity)
        {
            m_max_bundle_set_complexity = max_bundle_complexity;
        }
        std::size_t get_max_bundle_set_complexity() const
        {
            return m_max_bundle_set_complexity;
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
            // Compute data
            ComputationData data(graph);
            data.precomputeComputationData(trajectories);

            GCpp::DS::Heap<EdgeWeight, std::size_t, std::greater<>> edgePriority;
            for (const auto& kv : data.edgeWeights)
            {
                edgePriority.insert(kv.second, kv.first);
            }
            // Perform a greedy search for bundle elements
            std::vector<typename decltype(edgePriority)::Node> poppedElements;
            while (output.size() < m_max_bundle_set_complexity)
            {
                auto topNode = edgePriority.extract();
                poppedElements.push_back(topNode);

                if (topNode.value.value < m_min_matched_trajectories) break;
                std::cout << "Edge "<< topNode.id << " with " << topNode.value.value << '\n';

                // Two directions
                TrajectoryBundle bundle, bundle2;
                MatchingMap matchedTrajectories, matchedTrajectories2;
                bool success = searchForBundle(topNode.id, data, false, trajectories, bundle, matchedTrajectories);
                bool success2 = searchForBundle(topNode.id, data, true, trajectories, bundle2, matchedTrajectories2);

                if (success || success2)
                {
                    std::cout << "Found bundle " << output.size() + 1 << "\n";
                    TrajectoryBundle bundleToAdd;
                    MatchingMap matchedTrajectoriesToAdd;
                    if(success && success2)
                    {
                        auto bundleBetter = bundle.matchedTrajectoriesCount > bundle2.matchedTrajectoriesCount;
                        if(m_priotizie_by_length)
                        {
                            bundleBetter = bundle.bundleLength > bundle2.bundleLength;
                        }
                        matchedTrajectoriesToAdd = bundleBetter ? (matchedTrajectories) : (matchedTrajectories2);
                        bundleToAdd = bundleBetter ? std::move(bundle) : std::move(bundle2);
                    }
                    else
                    {
                        matchedTrajectoriesToAdd = success? (matchedTrajectories) : (matchedTrajectories2);
                        bundleToAdd = success ? std::move(bundle) : std::move(bundle2);
                    }
                    output.push_back(bundleToAdd);

                    for (auto el : poppedElements)
                    {
                        edgePriority.insert(data.edgeWeights[el.id], el.id);
                    }
                    poppedElements.clear();
                    // Update weights
                    // TODO: delete of all of trajectory?
                    data.updateEdgeIdToTrajectoryMapping(bundleToAdd.edges, matchedTrajectoriesToAdd);
                    data.updatePriority(bundleToAdd.edges, matchedTrajectoriesToAdd, edgePriority);
                    
                }
            }
            std::cout << "Found " << output.size() << "/" << m_max_bundle_set_complexity << " bundles \n";
        }
    };
}
#endif
