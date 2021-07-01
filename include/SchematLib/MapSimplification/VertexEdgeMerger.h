#ifndef SCHEMATLIB_MAPSIMPLIFICATION_VERTEXEDGEMERGER_H
#define SCHEMATLIB_MAPSIMPLIFICATION_VERTEXEDGEMERGER_H
#include <iostream>
#include <GCpp/DS/BoostEmbeddedGraph.h>
#include <boost/graph/lookup_edge.hpp>

#include "SimplficiationObserver.h"

namespace SchematLib::MapSimplification
{
    //namespace MergeOrders
    //{
    //    struct VertexOrder {};
    //    struct MostDenseVertexOrder {}; // Highest density first.
    //    struct LongestEdgeVertices {};
    //    struct ShortestEdgeVertices {};
    //}

    template<typename EmbeddedGraph, typename EdgeIndexStructure, typename SimplificationObserver = TrivialSimplificationObserver>
    class VertexEdgeMerger
    {
    public:
        using Observer = SimplificationObserver;
        using Vertex = typename boost::graph_traits<EmbeddedGraph>::vertex_descriptor;
        using Edge = typename boost::graph_traits<EmbeddedGraph>::edge_descriptor;
        using Point = decltype(GCpp::DS::get_vertex_location(std::declval<Vertex>(), std::declval<const EmbeddedGraph&>()));
        using NT = decltype(std::declval<const Point&>().x());
    private:
        NT m_threshold;

        static auto lookupEdge(Vertex v0, Vertex v1, const EmbeddedGraph& graph)
        {
            auto result = boost::lookup_edge(v0, v1, graph);
            /*if(!result.second)
            {
                return boost::lookup_edge(v1, v0, graph);
            }*/
            return result;
        }

        static std::size_t getVertexIndex(Vertex vertex, const EmbeddedGraph& graph)
        {
            return boost::get(boost::vertex_index_t{}, graph, vertex);
        }
        static std::size_t getVertexId(Vertex vertex, const EmbeddedGraph& graph)
        {
            return boost::get(GCpp::DS::vertex_id_t{}, graph, vertex);
        }
        static std::size_t getEdgeId(Edge edge, const EmbeddedGraph& graph)
        {
            return boost::get(GCpp::DS::edge_id_t{}, graph, edge);
        }

        bool findIsInSlab(const Point& vertexPosition, Edge edge, const EmbeddedGraph& inputGraph, NT& edgeParameter, Point& location) const
        {
            // Compute if in slab
            auto edgeBegin = GCpp::DS::get_vertex_location(boost::source(edge, inputGraph), inputGraph);
            auto edgeEnd = GCpp::DS::get_vertex_location(boost::target(edge, inputGraph), inputGraph);

            auto diff = vertexPosition - edgeBegin;
            auto edgeDirection = edgeEnd - edgeBegin;//unnormalized
            edgeParameter = diff.dot(edgeDirection);
            // Get relative location on the edge (i.e. if in [0,1], then on the segment).
            edgeParameter /= edgeDirection.sqLength();

            // Projection of point is not on segment: we are not going to merge.
            if (edgeParameter < 0 || edgeParameter> 1) return false;

            // Return projection point
            location = edgeBegin + (edgeEnd - edgeBegin) * edgeParameter;
            return true;
        }

        bool mergeVertexToEdge(Vertex targetVert, Edge targetEdge, const Point& newVertexPoint, EmbeddedGraph& targetGraph, EdgeIndexStructure& index, SimplificationObserver& observer) const
        {
            {
                auto[edge, exists] = lookupEdge(boost::source(targetEdge, targetGraph), boost::target(targetEdge, targetGraph), targetGraph);
                if (!exists) throw std::runtime_error("Unexisting edge considered for merging!");
            }
            auto nextEdgeId = boost::get_property(targetGraph, GCpp::DS::next_edge_id_t{});
            // Merge into edge by adding extra edges if needed.
            auto startVert = boost::source(targetEdge, targetGraph);
            auto endVert = boost::target(targetEdge, targetGraph);
            /*std::cout << std::setprecision(std::numeric_limits<NT>::digits10);
            std::cout << "Merging vertex " << targetVert << "(deg="<< boost::degree(targetVert,targetGraph) <<") to "
            << boost::get(GCpp::DS::edge_id_t{}, targetGraph, targetEdge) << "("
            << boost::source(targetEdge, targetGraph) << "," << boost::target(targetEdge,targetGraph)
                <<")" << std::endl;
            std::cout << "\tLoc:" << GCpp::DS::get_vertex_location(targetVert, targetGraph)
                << " , projected " << newVertexPoint
                << "edge " << GCpp::DS::get_vertex_location(startVert, targetGraph) << " " << GCpp::DS::get_vertex_location(endVert, targetGraph) << '\n';*/

            // Assign the new vertex position. Needed for updates to the index!
            boost::put(GCpp::DS::vertex_location_t{}, targetGraph, targetVert, newVertexPoint);

            std::vector<std::size_t> replacementEdges;
            replacementEdges.reserve(2);
            {
                // Setup new edges first.
                auto[edge, edgeExists] = lookupEdge(startVert, targetVert, targetGraph);
                if (!edgeExists)
                {
                    auto[newEdge, wasInserted] = boost::add_edge(startVert, targetVert, targetGraph);
                    assert(wasInserted); //Edge should be inserted.
                    replacementEdges.push_back(nextEdgeId);
                    // Assign edge ID.
                    boost::put(GCpp::DS::edge_id_t{}, targetGraph, newEdge, nextEdgeId);
                    index.insert(newEdge, targetGraph);
                    observer.handleNewEdge(nextEdgeId, GCpp::DS::getVertexId(targetGraph, startVert), GCpp::DS::getVertexId(targetGraph, targetVert));
                    std::cout << "\tCreated new edge in the process: " << nextEdgeId << " from " << startVert<< " to " << targetVert << std::endl;
                    ++nextEdgeId;
                    
                }
                else
                {
                    replacementEdges.push_back(GCpp::DS::getEdgeId(targetGraph, edge));
                }
            }
            // Check end vert
            {
                auto[edge, edgeExists] = lookupEdge(targetVert,endVert, targetGraph);
                if (!edgeExists) {
                    auto[newEdge, wasInserted] = boost::add_edge(targetVert, endVert, targetGraph);
                    assert(wasInserted); //Edge should be inserted.
                    replacementEdges.push_back(nextEdgeId);
                    // Assign edge ID.
                    boost::put(GCpp::DS::edge_id_t{}, targetGraph, newEdge, nextEdgeId);
                    index.insert(newEdge, targetGraph);
                    observer.handleNewEdge(nextEdgeId, GCpp::DS::getVertexId(targetGraph, targetVert), GCpp::DS::getVertexId(targetGraph, endVert));
                    std::cout << "\tCreated new edge in the process: " << nextEdgeId << " from " << targetVert << " to " << endVert << std::endl;
                    ++nextEdgeId;
                }
                else
                {
                    replacementEdges.push_back(GCpp::DS::getEdgeId(targetGraph, edge));
                }
            }
            const auto numRemoved =index.remove(targetEdge, targetGraph);
            if (numRemoved != 1) { std::cout << "Expected 1 edge removal, got: " << numRemoved << '\n'; }

            observer.handleEdgeReplace(GCpp::DS::getEdgeId(targetGraph, targetEdge), replacementEdges);
            // Remove the old edge.
            boost::remove_edge(targetEdge, targetGraph);
            boost::get_property(targetGraph, GCpp::DS::next_edge_id_t{}) = nextEdgeId;
            return true;
        }

        struct EdgeDistanceCompare
        {
            Point m_reference;
            const EmbeddedGraph* m_graph = nullptr;
            std::map<Edge, NT> m_distMap;
            typename EdgeIndexStructure::SqDistance m_sqDist;

            NT distance(const Edge& e)
            {
                if (m_distMap.find(e) != m_distMap.end()) return m_distMap.at(e);
                // Compute and save distance
                auto start = GCpp::DS::get_vertex_location(boost::source(e, *m_graph), *m_graph);
                auto end = GCpp::DS::get_vertex_location(boost::target(e, *m_graph), *m_graph);
                //
                m_distMap[e] = m_sqDist(std::make_pair(start, end), m_reference);
                return m_distMap.at(e);
            }

            // 'Before' comparison operator
            bool operator()(const Edge& e0, const Edge& e1)
            {
                return distance(e0) < distance(e1);
            }
        };

        bool m_repeatUntilConvergence = false;

        void mergeVertices(Vertex vertex, Vertex vertexToMergeTo, EmbeddedGraph& graph, SimplificationObserver& observer) const
        {
            const auto mergeId = getVertexId(vertexToMergeTo, graph);
            auto nextEdgeId = boost::get_property(graph, GCpp::DS::next_edge_id_t{});
            observer.handleVertexMerge(getVertexId(vertex, graph), mergeId);
            // Edge between the vertices to merge: if it exists, collapse it to the merge vertex.
            auto[edge,edgeExists] = lookupEdge(vertex, vertexToMergeTo, graph);
            if(edgeExists)
            {
                observer.handleEdgeCollapse(getEdgeId(edge, graph), mergeId);
            }
            // Any edge that is connected to vertex will be connected to vertexToMerge. Merge to existing if possible,
            // otherwise, create a new edge.
            for(auto e : GCpp::Helpers::Iterators::range(boost::out_edges(vertex,graph)))
            {
                const auto src = boost::source(e, graph);
                const auto target = boost::target(e, graph);
                // The other vertex to which the merged vertex connects
                auto uncommon = src == vertex ? target : src;

                if (uncommon == vertexToMergeTo) continue;
                auto[mergeEdge, mergeEdgeExists] = lookupEdge(vertexToMergeTo, uncommon,graph);
                if(mergeEdgeExists)
                {
                    observer.handleEdgeEdgeMerge(GCpp::DS::getEdgeId(graph, e), GCpp::DS::getEdgeId(graph, mergeEdge));
                }
                else
                {
                    auto[newEdge, wasAdded] = boost::add_edge(vertexToMergeTo, uncommon, graph);
                    if (!wasAdded) throw std::runtime_error("Edge was not added while merging!");
                    boost::put(GCpp::DS::edge_id_t{}, graph, newEdge, nextEdgeId);
                    observer.handleNewEdge(nextEdgeId, GCpp::DS::getVertexId(graph, vertexToMergeTo), GCpp::DS::getVertexId(graph, uncommon));
                    observer.handleEdgeEdgeMerge(GCpp::DS::getEdgeId(graph, e), nextEdgeId);
                    ++nextEdgeId;
                }
            }
            // Clear the vertex, delete later
            boost::clear_vertex(vertex, graph);
            // Update the next edge ID.
            boost::get_property(graph, GCpp::DS::next_edge_id_t{}) = nextEdgeId;
        }
    public:
        void setThreshold(const NT& myVar)
        {
            m_threshold = myVar;
        }
        NT threshold() const
        {
            return m_threshold;
        }
        void setRepeatUntilConvergence(const bool& repeatUntilConvergence)
        {
            m_repeatUntilConvergence = repeatUntilConvergence;
        }
        bool repeatUntilConvergence() const
        {
            return m_repeatUntilConvergence;
        }


        void operator()(EmbeddedGraph& inputGraph, SimplificationObserver& observer) const
        {
            while (true)
            {
                // Did any vertrex merge
                bool didMerge = false;

                EdgeIndexStructure index;
                index.construct(inputGraph);

                // Vertices to delete due to vertex merging for fringe cases.
                std::set<Vertex, std::greater<>> verticesToDelete;

                namespace it = GCpp::Helpers::Iterators;
                // For each vertex, find closest edge.
                // If the closest is too far away, ignore
                // Determine if the edge slab contains the vertex, otherwise ignore.
                // Merge the vertex into the edge: what to do with other connections to the vertex?
                for (auto v : it::range(boost::vertices(inputGraph)))
                {
                    if(boost::out_degree(v,inputGraph) == 0)
                    {
                        observer.handleVertexDelete(boost::get(GCpp::DS::vertex_id_t{}, inputGraph, v));
                        verticesToDelete.insert(v);
                        continue;
                    }
                    const auto vertexPosition = GCpp::DS::get_vertex_location(v, inputGraph);
                    // Find the nearest edge
                    std::vector<Edge> nearestEdge;
                    index.intersectingDisk(vertexPosition, m_threshold, nearestEdge);

                    // Shouldn't happen, unless the index is empty
                    if (nearestEdge.empty())
                    {
                        continue;
                    }
                    // Sort by distance to vertex
                    EdgeDistanceCompare compare;
                    compare.m_graph = &inputGraph;
                    compare.m_reference = vertexPosition;
                    std::sort(nearestEdge.begin(), nearestEdge.end(), compare);
                    std::optional<Edge> selectedEdge;
                    for (auto e : nearestEdge)
                    {
                        // Don't collapse to an edge that v is part of.
                        if (boost::source(e, inputGraph) == v || boost::target(e, inputGraph) == v) continue;
                        selectedEdge = e;
                        break;
                    }

                    if (!selectedEdge.has_value()) continue;
                    auto edge = selectedEdge.value();
                    // Edge does not exist anymore (should not happen!
                    if (!lookupEdge(boost::source(edge, inputGraph), boost::target(edge, inputGraph), inputGraph).second)
                    {
                        continue;
                    }


                    NT projectionParameter = 0;
                    Point projectPoint;
                    if (!findIsInSlab(vertexPosition, edge, inputGraph, projectionParameter, projectPoint)) continue;

                    bool merged = false;
                    // Merge vertices that get too close together when projected.
                    if(projectionParameter < 0.01 || projectionParameter > 0.99)
                    {
                        // Merge to the appropriate vertex on the edge. The current vertex will disappear.
                        const auto mergeTarget = projectionParameter < 0.5 ? boost::source(edge, inputGraph)  : boost::target(edge, inputGraph);
                        verticesToDelete.insert(v);
                        mergeVertices(v, mergeTarget, inputGraph, observer);
                        merged = true;
                    }
                    else
                    {
                        // We are not deleting the vertex, so the iterators should remain valid 
                        merged = mergeVertexToEdge(v, edge, projectPoint, inputGraph, index, observer);
                    }
                    didMerge = merged || didMerge;
                }
                for(auto v : verticesToDelete)
                {
                    boost::remove_vertex(v, inputGraph);
                }
                if(m_repeatUntilConvergence && didMerge)
                {
                    std::cout << "[VEM] Going to next it =====================================\n";
                    continue;
                }
                break;
            }
        }
    };
}
#endif