#ifndef SCHEMATLIB_ALGORITHMS_MMTRAJECTORYPREPROCESS_H
#define SCHEMATLIB_ALGORITHMS_MMTRAJECTORYPREPROCESS_H
#include <GCpp/DS/Trajectory.h>

#include "ComputeThroughTrajectories.h"
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <GCpp/DS/Heap.h>

namespace SchematLib::Algorithms
{
    struct GraphEdgeLength
    {
        const Models::UndirectedEmbeddedGraph* m_graph = nullptr;
        using value_type = Models::NT;
        using reference = const Models::NT&;
        using key_type = Models::UndirectedEmbeddedGraph::edge_descriptor;
        using category = boost::readable_property_map_tag;

        std::map<std::size_t, Models::UndirectedEmbeddedGraph::edge_descriptor> edgeIdToEdge;

        GraphEdgeLength(const Models::UndirectedEmbeddedGraph* graph):m_graph(graph)
        {
            GCpp::DS::computeEdgeIdToEdgeMap(*m_graph, edgeIdToEdge);
        }
        Models::NT get(Models::UndirectedEmbeddedGraph::edge_descriptor edge) const
        {
            return GCpp::DS::get_linear_edge_length(edge, *m_graph);
        }
        Models::NT get(std::size_t edgeId)
        {
            return get(edgeIdToEdge.at(edgeId));
        }
    };
}
namespace boost
{
    inline SchematLib::Models::NT get(const SchematLib::Algorithms::GraphEdgeLength& gel, const SchematLib::Models::UndirectedEmbeddedGraph::edge_descriptor e)
    {
        return gel.get(e);
    }
}

namespace SchematLib::Algorithms
{
    template<typename Graph, typename EdgeWeightMap>
    class UpperboundedShortestPath
    {
    public:
        using Vertex = typename Graph::vertex_descriptor;
        using Edge = typename Graph::edge_descriptor;
        using Weight_t = std::decay_t<decltype(std::declval<EdgeWeightMap>().at(std::declval<Edge>()))>;
    private:
        const Graph& m_graph;
        const EdgeWeightMap& m_edgeWeights;
        Weight_t m_upperbound = std::numeric_limits<Weight_t>::max();
    public:
        UpperboundedShortestPath(const Graph& graph,const EdgeWeightMap& edgeWeights)
        : m_graph(graph),
        m_edgeWeights(edgeWeights){}

        void setUpperbound(Weight_t value)
        {
            m_upperbound = value;
        }
        Weight_t upperbound() const
        {
            return m_upperbound;
        }
        
        bool operator()(Vertex src, Vertex target, Weight_t& distance)
        {
            std::unordered_map<Vertex, Weight_t> distances;
            auto isInf = [&distances](auto vId)
            {
                return distances.find(vId) == distances.end();
            };
            
            GCpp::DS::Heap<Weight_t> prioQueue;

            std::unordered_set<std::size_t> visited;
            auto wasVisited = [&visited](auto vId)
            {
                return visited.find(vId) != visited.end();
            };
            auto markVisited = [&visited](auto vId)
            {
                visited.insert(vId);
            };

            // Insert value-id combo.
            prioQueue.insert(static_cast<Weight_t>(0), src);
            distances[src] = 0;
            while (!prioQueue.empty())
            {
                auto heapNode = prioQueue.extract();
                // Mark visited
                markVisited(heapNode.id);
                if (heapNode.id == target) break;

                for (const auto& e : GCpp::Helpers::Iterators::range(boost::out_edges(heapNode.id,m_graph)))
                {
                    // Skip vertices that were already visited
                    if (wasVisited(*e->m_sink)) continue;
                    const auto edgeSrc = boost::source(e, m_graph);
                    const auto edgeTarget = boost::target(e, m_graph);
                    const auto targetV = edgeSrc == heapNode.id ? edgeTarget : edgeSrc;

                    auto newDist = distances[heapNode.id] + m_edgeWeights.at(e);
                    if ((isInf(targetV) || newDist < distances[targetV]) && newDist < upperbound)
                    {
                        distance[targetV] = newDist;
                        if (prioQueue.containsId(targetV))
                            prioQueue.updateKey(targetV, newDist);
                        else
                            prioQueue.insert(newDist, targetV);
                    }
                }
            }
            if(distances.find(target) != distances.end())
            {
                distance = distances.at(target);
                return true;
            }
            return false;
        }

        bool operator()(Vertex src, Vertex target, Weight_t& distance,std::vector<Edge>& shortestPath)
        {
            std::unordered_map<Vertex, Weight_t> distances;
            std::unordered_map<Vertex, std::pair<Edge,Vertex>> predecessors;
            auto isInf = [&distances](auto vId)
            {
                return distances.find(vId) == distances.end();
            };

            GCpp::DS::Heap<Weight_t> prioQueue;

            std::unordered_set<std::size_t> visited;
            auto wasVisited = [&visited](auto vId)
            {
                return visited.find(vId) != visited.end();
            };
            auto markVisited = [&visited](auto vId)
            {
                visited.insert(vId);
            };

            // Insert value-id combo.
            prioQueue.insert(static_cast<Weight_t>(0), src);
            distances[src] = 0;
            while (!prioQueue.empty())
            {
                auto heapNode = prioQueue.extract();
                // Mark visited
                markVisited(heapNode.id);
                if (heapNode.id == target) break;

                for (const auto& e : GCpp::Helpers::Iterators::range(boost::out_edges(heapNode.id, m_graph)))
                {
                    // Skip vertices that were already visited
                    if (wasVisited(*e->m_sink)) continue;
                    const auto edgeSrc = boost::source(e, m_graph);
                    const auto edgeTarget = boost::target(e, m_graph);
                    const auto targetV = edgeSrc == heapNode.id ? edgeTarget : edgeSrc;

                    auto newDist = distances[heapNode.id] + m_edgeWeights.at(e);
                    if ((isInf(targetV) || newDist < distances[targetV]) && newDist < upperbound)
                    {
                        distance[targetV] = newDist;
                        if (prioQueue.containsId(targetV))
                            prioQueue.updateKey(targetV, newDist);
                        else
                            prioQueue.insert(newDist, targetV);
                        predecessors[targetV] = std::make_pair(e,heapNode.id);
                    }
                }
            }

            if (distances.find(target) != distances.end())
            {
                distance = distances.at(target);

                auto currV = target;
                while(true)
                {
                    auto [e,prevV] = predecessors[currV];
                    shortestPath.push_back(e);
                    currV = prevV;
                    if (currV == src) break;
                }
                std::reverse(shortestPath.begin(), shortestPath.end());
                return true;
            }
            return false;
        }
    };

    class MMTrajectoryPreprocess
    {
    public:

        using Trajectory = GCpp::DS::Trajectory<std::size_t, std::string, std::vector>;
        using path = std::filesystem::path;
        using Vertex = Models::UndirectedEmbeddedGraph::vertex_descriptor;

        struct MergeBackElement
        {
            OffroadElement srcElement;
            std::vector<std::size_t> rerouteEdges;
        };
        using SegmentationElement = std::variant<OffroadElement, Subtrajectory, MergeBackElement>;
    private:
        template<typename T>
        using Vector = std::vector<T>;
        Models::NT m_max_detour_length = 0;

        void determineShortestPathOwn(const Models::UndirectedEmbeddedGraph& graph,
            Models::UndirectedEmbeddedGraph::vertex_descriptor start,
            Models::UndirectedEmbeddedGraph::vertex_descriptor end,
            Models::NT& pathLength, std::vector<std::size_t>& pathEdgeIds)
        {
            if(start==end)
            {
                pathLength = 0;
                pathEdgeIds = {};
                return;
            }
            GraphEdgeLength gel(&graph);
            using Vertex = Models::UndirectedEmbeddedGraph::vertex_descriptor;
            // Run dijkstra for shortest path.
            std::map< Vertex, Vertex> predecessors;
            std::map< Vertex, Models::NT> distance;
            struct Prio
            {
                std::map< Vertex, Models::NT>* distance = nullptr;

                Models::NT getDist(const Vertex& v) const
                {
                    if (distance->find(v) == distance->end()) return std::numeric_limits<Models::NT>::max();
                    return distance->at(v);
                }

                bool operator()(const Vertex& v0, const Vertex& v1) const
                {
                    const auto d0 = getDist(v0);
                    const auto d1 = getDist(v1);
                    return d0==d1 ? v0 < v1 : getDist(v0) < getDist(v1);
                }
            };
            Prio prio;
            prio.distance = &distance;
            std::unordered_set<Vertex> processed;
            std::set<Vertex, Prio> prioQueue(prio);
            distance[start] = 0;
            prioQueue.insert(start);
            while(!prioQueue.empty())
            {
                auto curr = *prioQueue.begin();
                prioQueue.erase(prioQueue.begin());
                if (processed.find(curr) != processed.end()) continue;
                processed.insert(curr);
                if (curr == end) break;

                for(auto e : GCpp::Helpers::Iterators::range(boost::out_edges(curr,graph)))
                {
                    auto target = boost::source(e, graph) == curr ? boost::target(e, graph) : boost::source(e, graph);
                    if (processed.find(target) != processed.end()) continue;

                    auto edgeLen = gel.get(e);
                    const auto computedDist = prio.getDist(curr) + edgeLen;
                    if(computedDist < prio.getDist(target))
                    {
                        distance[target] = computedDist;
                        predecessors[target] = curr;
                        prioQueue.insert(target);
                    }
                }
            }
            // Find out the shortest path.
            if (distance.find(end) != distance.end() && predecessors.find(end) != predecessors.end())
            {
                pathLength = distance.at(end);
                // Reconstruct path via predecessors.
                auto curr = end;
                while (true)
                {
                    auto pred = predecessors[curr];
                    auto[edge, exists] = boost::lookup_edge(pred, curr, graph);
                    if (!exists) std::tie(edge, exists) = boost::lookup_edge(curr, pred, graph);
                    if (!exists) throw std::runtime_error("Invalid edge in shortest path");
                    pathEdgeIds.push_back(GCpp::DS::getEdgeId(graph, edge));
                    if (pred == start) break;
                    curr = pred;
                }
                std::reverse(pathEdgeIds.begin(), pathEdgeIds.end());
            }
            else
            {
                pathLength = std::numeric_limits<Models::NT>::max();
            }
        }

        void determineShortestPath(const Models::UndirectedEmbeddedGraph& graph,
            Models::UndirectedEmbeddedGraph::vertex_descriptor start, 
            Models::UndirectedEmbeddedGraph::vertex_descriptor end, 
            Models::NT& pathLength, std::vector<std::size_t>& pathEdgeIds)
        {
            GraphEdgeLength gel(&graph);
            // Run dijkstra for shortest path.
            std::map< Models::UndirectedEmbeddedGraph::vertex_descriptor, Models::UndirectedEmbeddedGraph::vertex_descriptor> predecessors;
            std::map< Models::UndirectedEmbeddedGraph::vertex_descriptor, Models::NT> distance;
            
            boost::dijkstra_shortest_paths(graph, start,
                boost::predecessor_map(boost::associative_property_map<decltype(predecessors)>(predecessors))
                .distance_map(boost::associative_property_map<decltype(distance)>(distance))
                .weight_map(gel)
                .distance_inf(std::numeric_limits<Models::NT>::max())
            );
            // Find out the shortest path.
            if(distance.find(end) != distance.end() && predecessors.at(end) != end)
            {
                pathLength = distance.at(end);
                // Reconstruct path via predecessors.
                auto curr = end;
                while(true)
                {
                    auto pred = predecessors[curr];
                    auto[edge, exists] = boost::lookup_edge(pred, curr, graph);
                    if (!exists) std::tie(edge, exists) = boost::lookup_edge(curr, pred, graph);
                    if (!exists) throw std::runtime_error("Invalid edge in shortest path");
                    pathEdgeIds.push_back(GCpp::DS::getEdgeId(graph, edge));
                    if (pred == start) break;
                    curr = pred;
                }
                std::reverse(pathEdgeIds.begin(), pathEdgeIds.end());
            }
            else
            {
                pathLength = std::numeric_limits<Models::NT>::max();
            }
        }

        /**
         * \brief Write the segmentations to the output stream.
         * Format is (per line):
         * <trajectoryStringId> <isOnRoad> <trajectoryEdgeCount> <edgeIds...>
         * Where offroad elements are separated, but still having the same trajectory id.
         * \param trajectoryName Name of the trajectory
         * \param trajectory Trajectory element
         * \param segmentation Segmentation for the trajectory
         * \param output The output stream
         */
        void writeSegmentation(const std::string& trajectoryName, const Trajectory& trajectory, const std::vector<SegmentationElement>& segmentation,
            const Models::UndirectedEmbeddedGraph& originalMap,
            std::ostream& onroadOutputStream, std::ostream& offroadOutputStream)
        {
            auto it = segmentation.begin();
            std::vector<std::size_t> trajectoryEdges;
            std::size_t trajectorySubSegment = 0;
            for(; it != segmentation.end(); ++it)
            {
                if(std::holds_alternative<Subtrajectory>(*it) || std::holds_alternative<MergeBackElement>(*it))
                {
                    onroadOutputStream << trajectoryName << "_" << trajectorySubSegment << " ";
                    ++trajectorySubSegment;
                    while(it != segmentation.end())
                    {
                        if (std::holds_alternative<Subtrajectory>(*it))
                        {
                            const auto& el = std::get<Subtrajectory>(*it);
                            for (auto i = el.startEdgeIndex; i < el.endEdgeIndex; ++i)
                            {
                                trajectoryEdges.push_back(trajectory.at(i));
                            }
                        }
                        else if(std::holds_alternative<MergeBackElement>(*it))
                        {
                            const auto& el = std::get<MergeBackElement>(*it);
                            for (auto e : el.rerouteEdges)
                            {
                                trajectoryEdges.push_back(e);
                            }
                        }
                        else
                        {
                            // Go back one, should be processed in next iteration
                            --it;
                            break;
                        }
                        ++it;
                    }
                    onroadOutputStream << trajectoryEdges.size() << " ";
                    for (auto el : trajectoryEdges)
                    {
                        onroadOutputStream << " " << el;
                    }
                    onroadOutputStream << '\n';
                    trajectoryEdges.clear();
                    // Jump out to avoid past end reading.
                    if (it == segmentation.end()) break;
                }
                // Must be offroad.
                else
                {
                    const auto& el = std::get<OffroadElement>(*it);
                    offroadOutputStream << trajectoryName << "_" << trajectorySubSegment << " " 
                        << GCpp::DS::getVertexId(originalMap,el.startV) << " " << GCpp::DS::getVertexId(originalMap, el.endV) << " "
                        << (el.endEdgeIndex-el.startEdgeIndex);
                    ++trajectorySubSegment;
                    for(auto i = el.startEdgeIndex; i < el.endEdgeIndex;++i)
                    {
                        offroadOutputStream << ' ' << trajectory.at(i);
                    }
                    offroadOutputStream << '\n';
                }
            }
        }
    public:
        void set_max_detour_length(Models::NT length)
        {
            m_max_detour_length = length;
        }
        Models::NT get_max_detour_length() const
        {
            return m_max_detour_length;
        }
        /**
         * \brief Preprocess trajectories for simplification visualization. Split trajectories in off/on-road, and remap
         * off-road elements to the subgraph that was selected if the shortest path detour is within the specified max detour length.
         * \param mapFile 
         * \param edgeSelection Subselection of the road graph that will be simplified
         * \param mmTrajectories 
         * \param output 
         */
        void operator()(const path& mapFile, const path& generalizationFile, const path& mmtrajectoriesFile,
            const path& onroadOutput, const path& offroadOutput)
        {
            ComputeThroughTrajectories throughComputer;
            std::map < std::string, std::vector < ComputeThroughTrajectories::SubtrajectoryType >> segmentations;
            std::map<std::string, ComputeThroughTrajectories::Trajectory> trajectories;
            throughComputer(mapFile, generalizationFile, mmtrajectoriesFile, segmentations, trajectories);

            // Read full map
            IO::ReadOsmMap undirectedReader;
            int epsg;
            Models::UndirectedEmbeddedGraph original, selection;
            undirectedReader.readCustomFormat(mapFile, original, epsg);

            IO::ReadGeneralization genReader;
            genReader(generalizationFile, selection);

            // Original
            GraphEdgeLength originalEdgeLengths(&original);

            // Mapping from vertex ID to vertex in selection graph.
            std::map<std::size_t, Vertex> vIdToVertMap;
            GCpp::DS::computeVertexIdToVertexMap(selection, vIdToVertMap);

            std::map < std::string, std::vector<SegmentationElement>> newSegmentations;
            
            // For each offroad element, determine if the detour length is sufficient
            for(const auto& trajSegmentation : segmentations)
            {
                newSegmentations[trajSegmentation.first] = {};
                auto& targetSegmentation = newSegmentations[trajSegmentation.first];
                for(const auto& el: trajSegmentation.second)
                {
                    if (!std::holds_alternative<OffroadElement>(el))
                    {
                        targetSegmentation.push_back(std::get<Subtrajectory>(el));
                        continue;
                    }

                    const auto& offroadEl = std::get<OffroadElement>(el);
                    // Check if we should reroute for a through element
                    if (offroadEl.isThrough())
                    {
                        // In original graph.
                        std::vector<std::size_t> pathEdgeIds;
                        Models::NT pathLength = 0;
                        if(offroadEl.startV == offroadEl.endV)
                        {
                            continue;
                            //throw std::runtime_error("Invalid offroad element");
                        }
                        std::cout << "StartV: " << offroadEl.startV << ", end V:L " << offroadEl.endV << "\n";
                        const auto& traj = trajectories.at(offroadEl.sourceTrajectory);
                        for(auto index = offroadEl.startEdgeIndex; index < offroadEl.endEdgeIndex; ++index)
                        {
                            pathEdgeIds.push_back(traj.at(index));
                            pathLength += originalEdgeLengths.get(pathEdgeIds.back());
                        }
                        // Rerouted shortest path in selection graph
                        auto reroutedStartV = vIdToVertMap.at(GCpp::DS::getVertexId(original, offroadEl.startV));
                        auto reroutedEndV = vIdToVertMap.at(GCpp::DS::getVertexId(original, offroadEl.endV));
                        std::vector<std::size_t> reroutedPathEdgeIds;
                        Models::NT reroutedPathLength;
                        // TODO: can do upperbounded s-t shortest path, depending on the max allowed detour length.
                        //determineShortestPath(selection, reroutedStartV, reroutedEndV, reroutedPathLength, reroutedPathEdgeIds);
                        determineShortestPathOwn(selection, reroutedStartV, reroutedEndV, reroutedPathLength, reroutedPathEdgeIds);
                        if (std::abs(pathLength - reroutedPathLength) < m_max_detour_length)
                        {
                            targetSegmentation.push_back(MergeBackElement{ offroadEl, reroutedPathEdgeIds });
                        }
                        else
                        {
                            targetSegmentation.push_back(offroadEl);
                        }
                    }
                    else
                    {
                        targetSegmentation.push_back(offroadEl);
                    }
                }
            }

            // Write output
            std::ofstream offroadOutputStream(offroadOutput.string());
            std::ofstream onroadOutputStream(onroadOutput.string());
            for(const auto& kv: newSegmentations)
            {
                writeSegmentation(kv.first, trajectories.at(kv.first), kv.second, original, onroadOutputStream, offroadOutputStream);
            }
        }
    };
}
#endif
