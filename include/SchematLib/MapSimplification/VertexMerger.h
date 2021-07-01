#ifndef SCHEMATLIB_MAPSIMPLIFICATION_VERTEXMERGER_H
#define SCHEMATLIB_MAPSIMPLIFICATION_VERTEXMERGER_H
#include "FaceMergeTraits.h"
#include <boost/graph/lookup_edge.hpp>
#include <GCpp/DS/BoostEmbeddedGraph.h>
#include <GCpp/Helpers/Timer.h>

#include "SimplficiationObserver.h"

namespace SchematLib::MapSimplification
{
    namespace MergeOrders
    {
        struct VertexOrder {};
        struct MostDenseVertexOrder {}; // Highest density first.
        struct LongestEdgeVertices {};
        struct ShortestEdgeVertices {};
    }

    template<typename EmbeddedGraph, typename NodeIndexStructure,
        typename SimplificationObserver = TrivialSimplificationObserver,
        typename OrderStrategy = MergeOrders::VertexOrder>
        class VertexMerger
    {
    public:
        using Vertex = typename boost::graph_traits<EmbeddedGraph>::vertex_descriptor;
        using Edge = typename boost::graph_traits<EmbeddedGraph>::edge_descriptor;
        using Point = decltype(GCpp::DS::get_vertex_location(std::declval<Vertex>(), std::declval<const EmbeddedGraph&>()));
        using NT = decltype(std::declval<const Point&>().x());

        /**
         * \brief Class to keep track of vertex mapping to other vertices for merging.
         */
        class MergeMapping
        {
            std::map<Vertex, Vertex> m_mapping;
        public:
            auto begin() const
            {
                return m_mapping.begin();
            }
            void clear() { m_mapping.clear(); }
            bool empty() const { return m_mapping.empty(); }
            auto end() const
            {
                return m_mapping.end();
            }

            bool isMapped(Vertex v) const
            {
                return m_mapping.find(v) != m_mapping.end();
            }

            Vertex getTarget(Vertex v) const
            {
                if (!isMapped(v)) throw std::runtime_error("Vertex does not have target! ");
                return m_mapping.at(v);
            }

            void mapVertex(Vertex src, Vertex target)
            {
                if (src == target) return;
                // TODO: some validity checking?
                assert(!isMapped(src));
                assert(!isMapped(target) || isMergeVertex(target));

                if (!isMapped(target)) m_mapping[target] = target;
                m_mapping[src] = target;
            }

            /**
             * \brief Returns whether the given vertex is a merge vertex: a vertex to which a group of other vertices will be merged.
             * \param v The vertex in question.
             * \return
             */
            bool isMergeVertex(Vertex v) const
            {
                return isMapped(v) && getTarget(v) == v;
            }
            bool hasMergeTarget(Vertex v) const
            {
                return isMapped(v) && getTargte(v) != v;
            }
        };
    private:
        // Threshold for merging vertices
        NT m_radius;
        // Strategy for the order of vertex merging
        OrderStrategy m_strategy;
        // Use mean position when merging vertices, instead of just using the merge vertex position.
        bool m_useMeanPosition = false;

        static_assert(std::is_same_v<OrderStrategy, MergeOrders::VertexOrder> ||
            std::is_same_v<OrderStrategy, MergeOrders::MostDenseVertexOrder> ||
            std::is_same_v<OrderStrategy, MergeOrders::LongestEdgeVertices> ||
            std::is_same_v<OrderStrategy, MergeOrders::ShortestEdgeVertices>);

        static std::size_t getVertexId(Vertex v, const EmbeddedGraph& graph)
        {
            return boost::get(GCpp::DS::vertex_id_t{}, graph, v);
        }
        static std::size_t getEdgeId(Edge e, const EmbeddedGraph& g)
        {
            return boost::get(GCpp::DS::edge_id_t{}, g, e);
        }
        static std::pair<Edge, bool> getEdge(Vertex v0, Vertex v1, const EmbeddedGraph& g)
        {
            auto result = boost::lookup_edge(v0, v1, g);
            /*if(!result.second)
            {
                return boost::lookup_edge(v1, v0, g);
            }*/
            return result;
        }

        /**
         * \brief Merge the given vertex with its neighbourhood if any vertices are near and the vertex is not already merged.
         * \param v The vertex
         * \param inputGraph The input graph
         * \param index The indexing structure for querying near vertices 
         * \param mapping The merge groups output structure. Check here if the vertex was already merged.
         */
        void computeNeighbourhoodMerge(Vertex v, const EmbeddedGraph& inputGraph, const NodeIndexStructure& index, MergeMapping& mapping, std::map<Vertex, Point>& positions) const
        {
            if (mapping.isMapped(v)) return;

            const auto vertexLocation = GCpp::DS::get_vertex_location(v, inputGraph);
            // Query for vertices within range
            std::set<Vertex> closeVertices;
            index.containedInDisk(vertexLocation.x(), vertexLocation.y(), m_radius, closeVertices);
            // Don't merge with anything
            if (closeVertices.empty()) return;

            // Count then number of vertices merged with the current vertex
            std::size_t mergeCount = 0;
            Point pos = vertexLocation;
            for (auto collapseV : closeVertices)
            {
                if (mapping.isMapped(collapseV)) continue;
                mapping.mapVertex(collapseV, v);
                ++mergeCount;
                if (m_useMeanPosition) pos += GCpp::DS::get_vertex_location(collapseV, inputGraph);
            }
            // We disallow the current vertex to be merged again if it merged with other vertices
            if (mergeCount > 0) {
                mapping.mapVertex(v, v);
                if (m_useMeanPosition)
                {
                    positions[v] = pos / (static_cast<NT>(mergeCount) + NT{ 1 });
                }
                else
                {
                    positions[v] = vertexLocation;
                }
            }
        }
        struct DensityAccessor
        {
            mutable std::vector<std::size_t> densities;
            const NodeIndexStructure* m_index;
            const EmbeddedGraph* m_graph;
            NT m_radius;
            DensityAccessor(std::size_t vertexCount, const NodeIndexStructure* index,
                const EmbeddedGraph* graph, NT radius) :densities(vertexCount, 0),
                m_index(index),
                m_graph(graph),
                m_radius(radius)
            {}

            void computeDensitiy(Vertex v0) const
            {
                if (densities[v0] == 0)
                {
                    const auto vertexLoc = GCpp::DS::get_vertex_location(v0, *m_graph);
                    std::set<Vertex> closeVertices;
                    m_index->containedInDisk(vertexLoc.x(), vertexLoc.y(), m_radius, closeVertices);
                    densities[v0] = closeVertices.size() + 1;
                }
            }

            bool operator()(const Vertex& v0, const Vertex& v1) const
            {
                computeDensitiy(v0);
                computeDensitiy(v1);
                return densities[v0] > densities[v1];
            }
        };

        /**
         * \brief Notify the observer of all structural changes we made
         * \param mergeMappings The mapping of vertex merging
         * \param sourceGraph The original graph
         * \param resultGraph The resulting graph after merging
         * \param oldToNewMap The mapping of old to new vertices.
         * \param observer The observer to notify of all events.
         */
        void triggerModificationEvents(const MergeMapping& mergeMappings, const EmbeddedGraph& sourceGraph, const EmbeddedGraph& resultGraph, 
            const std::map<Vertex,Vertex>& oldToNewMap, SimplificationObserver& observer) const
        {
            namespace it = GCpp::Helpers::Iterators;
            // Construct the groups, keyed by the merge vertex.
            std::map<Vertex, std::vector<Vertex>> groups;
            auto hasVert = [&groups](Vertex v) {return groups.find(v) != groups.end(); };
            for(const auto& mapping: mergeMappings)
            {
                if (!hasVert(mapping.second)) groups[mapping.second] = {};
                if(mapping.first != mapping.second) groups[mapping.second].push_back(mapping.first);
            }
            auto hasEdge = [](Vertex v0, Vertex v1, const EmbeddedGraph& graph, Edge& edge)
            {
                bool exists = false;
                std::tie(edge, exists) = boost::lookup_edge(v0, v1, graph);
                return exists;
            };

            for(const auto& group: groups)
            {
                // Vertex that all others will merge into.
                const auto mergeVertex = group.first;
                // Collect all vertices in the group.
                std::set<Vertex> affectedVertices;
                affectedVertices.insert(group.second.begin(), group.second.end());
                for(auto v: group.second)
                {
                    // Check all edges
                    for(auto e: it::range(boost::out_edges(v, sourceGraph)))
                    {
                        auto src = boost::source(e, sourceGraph);
                        auto target = boost::target(e, sourceGraph);
                        if (src != v) std::swap(src, target); //Target is now the other vertex.
                        // Intra-group edge collapses (not to merge vertex)
                        if(affectedVertices.find(target) != affectedVertices.end())
                        {
                            if(src < target) //Make sure to emit once only.
                            {
                                observer.handleEdgeCollapse(GCpp::DS::getEdgeId(sourceGraph,e), GCpp::DS::getVertexId(sourceGraph, mergeVertex));
                            }
                        }
                        // Connected to merge vertex.
                        else if(target == mergeVertex)
                        {
                            observer.handleEdgeCollapse(GCpp::DS::getEdgeId(sourceGraph, e), GCpp::DS::getVertexId(sourceGraph, mergeVertex));
                        }
                        // Edge between collapsing groups
                        else if(mergeMappings.isMapped(target))
                        {
                            Vertex otherMergeVertex;
                            if (!mergeMappings.isMergeVertex(target))
                            {
                                otherMergeVertex = mergeMappings.getTarget(target);
                                // Skip one side, other side will also try to reroute this edge.
                                if (mergeVertex > otherMergeVertex)
                                {
                                    continue;
                                }
                            }
                            else
                            {
                                otherMergeVertex = target;
                            }
                            // Check if the edge between the group merge vertices existed before, merge the edge to that one.
                            Edge mergerEdge;
                            if (hasEdge(mergeVertex, otherMergeVertex, sourceGraph, mergerEdge))
                            {
                                const auto eId = GCpp::DS::getEdgeId(sourceGraph, mergerEdge);
                                observer.handleEdgeEdgeMerge(getEdgeId(e, sourceGraph), eId);
                            }
                            // Otherwise, merge to the newly created edge that should be present.
                            else
                            {
                                auto mappedSrc = oldToNewMap.at(mergeVertex);
                                auto mappedTarget = oldToNewMap.at(otherMergeVertex);
                                auto[mergerEdge, exists] = getEdge(mappedSrc, mappedTarget, resultGraph);
                                if (!exists) {
                                    std::string msg = "Invalid state: mapped vertices do not have an edge: merge-v" + std::to_string(mergeVertex) + ": target " + std::to_string(target);
                                    throw std::runtime_error(msg);
                                }
                                observer.handleEdgeEdgeMerge(GCpp::DS::getEdgeId(sourceGraph, e), GCpp::DS::getEdgeId(resultGraph, mergerEdge));
                            }
                        }
                        // Edge between this group and a vertex that is not changed.
                        else 
                        {
                            Vertex otherMergeVertex = target;
                            // Check if the edge between the group merge vertices existed before, merge the edge to that one.
                            Edge mergerEdge;
                            if(hasEdge(mergeVertex, otherMergeVertex, sourceGraph, mergerEdge))
                            {
                                const auto eId = GCpp::DS::getEdgeId(sourceGraph, mergerEdge);
                                observer.handleEdgeEdgeMerge(GCpp::DS::getEdgeId(sourceGraph,e), eId);
                            }
                            // Otherwise, merge to the newly created edge that should be present.
                            else
                            {
                                auto mappedSrc = oldToNewMap.at(mergeVertex);
                                auto mappedTarget = oldToNewMap.at(otherMergeVertex);
                                auto[mergerEdge, exists] = getEdge(mappedSrc, mappedTarget, resultGraph);
                                if (!exists) {
                                    std::string msg = "Invalid state: mapped vertices do not have an edge: merge-v" + std::to_string(mergeVertex) + ": target " + std::to_string(target);
                                    throw std::runtime_error(msg);
                                }
                                observer.handleEdgeEdgeMerge(getEdgeId(e, sourceGraph), getEdgeId(mergerEdge, resultGraph));
                            }
                        }
                    }

                    // Emit vertex merging events.
                    observer.handleVertexMerge(getVertexId(v, sourceGraph), getVertexId(mergeVertex, sourceGraph));
                }
            }
        }



        void computeMergeGroups(const EmbeddedGraph& inputGraph, const NodeIndexStructure& index, MergeMapping& mapping, std::map<Vertex, Point>& newVertexLocations) const
        {
            using namespace GCpp::Helpers::Iterators;
            if constexpr (std::is_same_v<OrderStrategy, MergeOrders::VertexOrder>)
            {
                for (auto v : range(boost::vertices(inputGraph)))
                {
                    computeNeighbourhoodMerge(v, inputGraph, index, mapping, newVertexLocations);
                }
            }
            else if constexpr (std::is_same_v<OrderStrategy, MergeOrders::MostDenseVertexOrder>)
            {
                const auto vertexRange = boost::vertices(inputGraph);
                std::vector<Vertex> vertices(vertexRange.first, vertexRange.second);
                //std::vector<std::size_t> densities(vertices.size(), 0);
                // Sort by density,
                DensityAccessor densities(vertices.size(), &index, &inputGraph, m_radius);

                std::sort(vertices.begin(), vertices.end(), densities);
                for (auto v : vertices)
                {
                    computeNeighbourhoodMerge(v, inputGraph, index, mapping, newVertexLocations);
                }
            }
            else if constexpr (std::is_same_v<OrderStrategy, MergeOrders::LongestEdgeVertices> || std::is_same_v<OrderStrategy, MergeOrders::ShortestEdgeVertices>)
            {
                std::vector<Edge> edges;
                edges.reserve(boost::num_edges(inputGraph));
                auto edgeRange = boost::edges(inputGraph);
                std::copy(edgeRange.first, edgeRange.second, std::back_inserter(edges));
                // Sort such that largest edge is first.
                std::sort(edges.begin(), edges.end(), [&inputGraph](const auto& e0, const auto& e1)
                {
                    if constexpr (std::is_same_v<OrderStrategy, MergeOrders::LongestEdgeVertices>)
                    {
                        return GCpp::DS::get_linear_edge_length(e0, inputGraph) > GCpp::DS::get_linear_edge_length(e1, inputGraph);
                    }
                    else
                    {
                        return GCpp::DS::get_linear_edge_length(e0, inputGraph) < GCpp::DS::get_linear_edge_length(e1, inputGraph);
                    }
                });
                // Process
                for (auto e : edges)
                {
                    const auto src = boost::source(e, inputGraph);
                    const auto target = boost::target(e, inputGraph);
                    computeNeighbourhoodMerge(src, inputGraph, index, mapping, newVertexLocations);
                    computeNeighbourhoodMerge(target, inputGraph, index, mapping, newVertexLocations);
                }
            }
            else
            {
                static_assert(false);
            }
        }

        /**
         * \brief Given the merge mapping, construct a new embedded graph.
         * \param sourceGraph 
         * \param mapping 
         * \param newVertexLocations 
         * \param observer 
         * \param outputGraph 
         */
        void constructMerged(const EmbeddedGraph& sourceGraph, const MergeMapping& mapping, const std::map<Vertex, Point>& newVertexLocations, SimplificationObserver& observer,
            EmbeddedGraph& outputGraph) const
        {
            using namespace GCpp::Helpers::Iterators;
            std::map<Vertex, std::size_t> oldToNewVertex;
            std::map<std::size_t, Vertex> newToOldVertex;
            std::size_t totalOutVertices = 0;
            for (const auto& v : range(boost::vertices(sourceGraph)))
            {
                // If the vertex is not merged, or is only merged towards, record it
                if (!mapping.isMapped(v) || mapping.isMergeVertex(v))
                {
                    oldToNewVertex[v] = totalOutVertices;
                    newToOldVertex[totalOutVertices] = v;
                    ++totalOutVertices;
                }
            }
            outputGraph = EmbeddedGraph(totalOutVertices);
            outputGraph.m_use_canonical_edges_ids = true;
            // Copy vertex locations and ID's.
            for (const auto& v : range(boost::vertices(outputGraph)))
            {
                boost::put(GCpp::DS::vertex_location_t{}, outputGraph, v, newVertexLocations.at(newToOldVertex[v]));
                boost::put(GCpp::DS::vertex_id_t{}, outputGraph, v, boost::get(GCpp::DS::vertex_id_t{}, sourceGraph, newToOldVertex[v]));
            }
            // Create all edges in the output graph. Don't add duplicates.
            std::set<std::pair<Vertex, Vertex>> constructedEdges;
            auto edgeWasConstructed = [&constructedEdges](Vertex v0, Vertex v1)
            {
                return constructedEdges.find(std::make_pair(v0, v1)) != constructedEdges.end() &&
                    constructedEdges.find(std::make_pair(v1, v0)) != constructedEdges.end();
            };

            auto getNewVertex = [&oldToNewVertex, &mapping](const Vertex& oldVertex)
            {
                // The vertex is not a new vertex: it was merged
                if (oldToNewVertex.find(oldVertex) == oldToNewVertex.end())
                {
                    // Find where it mapped to, map that to new
                    return oldToNewVertex.at(mapping.getTarget(oldVertex));
                }
                return oldToNewVertex[oldVertex];
            };
            // Edge ID for newly constructed edges.
            std::size_t currentEdgeId = boost::get_property(sourceGraph, GCpp::DS::next_edge_id_t{});

            for (const auto& e : range(boost::edges(sourceGraph)))
            {
                auto src = boost::source(e, sourceGraph);
                auto target = boost::target(e, sourceGraph);
                // Translate to new graph
                const auto newSrc = getNewVertex(src);
                const auto newTarget = getNewVertex(target);
                // Check if edge not already added and edge not degenerate.
                if (newSrc == newTarget || edgeWasConstructed(newSrc, newTarget)) continue;

                auto[newEdge, wasCreated] = GCpp::DS::safe_add_edge(newSrc, newTarget, outputGraph);
                if (!wasCreated) continue;

                std::size_t edgeId = 0;
                // Find out which edge ID to use: if the mapped vertices of the edge formed an edge in the original,
                // copy the edge ID. Otherwise, a ''virtual'' (non-existent in original) edge is created, which will get
                // a new edge ID that is larger than the original IDs. 
                {
                    auto mappedSrc = mapping.isMapped(src) ? mapping.getTarget(src) : src;
                    auto mappedTarget = mapping.isMapped(target) ? mapping.getTarget(target) : target;
                    // Edge exists in original
                    auto[edgeDesc, exists] = getEdge(mappedSrc, mappedTarget, sourceGraph);
                    if (exists)
                    {
                        edgeId = boost::get(GCpp::DS::edge_id_t{}, sourceGraph, edgeDesc);
                    }
                    else
                    {
                        edgeId = currentEdgeId;
                        observer.handleNewEdge(edgeId, newSrc, newTarget);
                        ++currentEdgeId;
                    }
                }

                boost::put(GCpp::DS::edge_id_t{}, outputGraph, newEdge, edgeId);
                constructedEdges.insert(std::make_pair(newSrc, newTarget));
            }

            // Trigger all modifications
            triggerModificationEvents(mapping, sourceGraph, outputGraph, oldToNewVertex, observer);
            // Update the next edge ID on the graph.
            boost::get_property(outputGraph, GCpp::DS::next_edge_id_t{}) = currentEdgeId;
            boost::get_property(outputGraph, GCpp::DS::next_vertex_id_t{}) = boost::get_property(sourceGraph, GCpp::DS::next_vertex_id_t{});
        }
        static void correctOrThrow(bool condition, const std::string& correctMsg, const std::string& throwMsg)
        {
            if(condition) std::cout << "[VertexMerger] Check : " << correctMsg << '\n';
            else throw std::runtime_error(throwMsg);
        }
    public:
        
        VertexMerger(NT radius, OrderStrategy strategy = {}) :m_radius(radius), m_strategy(strategy) {}

        void setRadius(const NT& radius)
        {
            m_radius = radius;
        }
        NT radius() const
        {
            return m_radius;
        }
        void setUseMeanPosition(const bool& useMeanPosition)
        {
            m_useMeanPosition = useMeanPosition;
        }
        bool useMeanPosition() const
        {
            return m_useMeanPosition;
        }

        /**
         * \brief Merges vertices in the input graph, resulting in the output graph. Notifies the observer of structural changes that happen.
         * \param inputGraph The input graph
         * \param observer The observer that will be notified of changes 
         * \param outputGraph The output graph.
         */
        void operator()(const EmbeddedGraph& inputGraph, SimplificationObserver& observer, EmbeddedGraph& outputGraph) const
        {
            //{
            //    std::set<std::size_t> originalVertexIds;
            //    GCpp::DS::computeVertexIdSet(inputGraph, originalVertexIds);
            //    std::map<std::size_t, std::pair<std::size_t, std::size_t>> edgeMap;
            //    GCpp::DS::getEdgeIdToVerticesMapping(inputGraph, edgeMap);
            //    std::decay_t<SimplificationObserver> obs(originalVertexIds, edgeMap);
            //    // Time how long this takes.
            //    GCpp::Helpers::Timer timer;
            //    timer.start();
            //    FaceMergeTraits::SetupArrangement<EmbeddedGraph, SimplificationObserver> createArrangement;
            //    FaceMergeTraits::Arr arrangement;
            //    std::size_t nextV, nextE, nextF;
            //    createArrangement(inputGraph, arrangement, obs, nextV, nextE, nextF);
            //    timer.stop();
            //    std::cout << "[VM] CGAL arrangement in " << timer.elapsedMs() << " ms\n";

            //    // See if we find 'unreal' elements
            //    for(auto vIt = arrangement.vertices_begin(); vIt != arrangement.vertices_end(); ++vIt)
            //    {
            //        if (vIt->is_at_open_boundary()) std::cout << "Found vert at open boundary\n";
            //    }
            //    for(auto eIt = arrangement.edges_begin(); eIt != arrangement.edges_end(); ++eIt)
            //    {
            //        if (eIt->is_fictitious()) std::cout << "Fictitious edge found\n";
            //    }
            //    std::cout << "Done checking arrangement\n";
            //}
            if (!inputGraph.m_use_canonical_edges_ids) throw std::runtime_error("Graph needs to use canonical edge IDs");
            MergeMapping mergeGroups;
            std::map<Vertex, Point> vertexLocations;
            EmbeddedGraph temp = inputGraph;
            int iteration = 0;
            NodeIndexStructure currentIndex;
            currentIndex.construct(temp);
            while (true)
            {
                correctOrThrow(GCpp::DS::hasUniqueEdgeIds(temp),"Edge IDs unique","Not all edge IDs are unique!");
                //correctOrThrow(GCpp::DS::nextEdgeIdIsCorrect(temp), "Next edge ID correct", "Next edge ID is incorrect");
                {
                    auto nextEId = boost::get_property(temp, GCpp::DS::next_edge_id_t{});
                    bool success = true;
                    for(auto e : GCpp::Helpers::Iterators::range(boost::edges(temp)))
                    {
                        const auto eId = GCpp::DS::getEdgeId(temp, e);
                        if(SpecialEdgeIdFunctions::isCanonical(eId) && eId >= nextEId)
                        {
                            success = false;
                            break;
                        }
                    }
                    correctOrThrow(success, "Next edge ID correct", "Next edge ID is incorrect");
                }
                correctOrThrow(GCpp::DS::hasUniqueVertexIds(temp), "Vertex IDs unique","Not all vertex IDs are unique!");
                correctOrThrow(GCpp::DS::nextVertexIdIsCorrect(temp),"Next vert ID correct","Next vert ID is incorrect!");
                ++iteration;
                mergeGroups.clear();
                vertexLocations.clear();
                computeMergeGroups(temp, currentIndex, mergeGroups, vertexLocations);
                if (mergeGroups.empty() || iteration > 10) break;

                constructMerged(temp, mergeGroups, vertexLocations, observer, outputGraph);
                temp = outputGraph;
                currentIndex = {};
                currentIndex.construct(temp);
            }

        }
    };
}
#endif