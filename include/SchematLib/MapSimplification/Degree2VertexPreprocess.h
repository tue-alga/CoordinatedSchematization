#ifndef SCHEMATLIB_MAPSIMPLIFICATION_DEGREE2VERTEXPREPROCESS_H
#define SCHEMATLIB_MAPSIMPLIFICATION_DEGREE2VERTEXPREPROCESS_H
#include <boost/graph/lookup_edge.hpp>
#include <GCpp/DS/BoostEmbeddedGraph.h>

#include "SimplficiationObserver.h"

namespace SchematLib::MapSimplification
{
    template<typename EmbeddedGraph>
    class Degree2VertexPreprocess
    {
    public:
        using Vertex = typename boost::graph_traits<EmbeddedGraph>::vertex_descriptor;
        using Edge = typename boost::graph_traits<EmbeddedGraph>::edge_descriptor;
        using Point = decltype(GCpp::DS::get_vertex_location(std::declval<Vertex>(), std::declval<const EmbeddedGraph&>()));
        using NT = decltype(std::declval<const Point&>().x());
    private:
        NT m_maxTurningAngle;

        bool m_repeatUntilConvergence = false;

        struct Triple
        {
            Edge e0, e1;
            std::array<Vertex, 3> vertices;
            std::array<Vertex, 3> vIds;

            NT computeTurningAngle(const EmbeddedGraph& graph) const
            {
                std::array<Point, 3> points;
                for(auto i = 0; i < 3; ++i)
                {
                    points[i] = GCpp::DS::get_vertex_location(vertices[i], graph);
                }
                auto diff0 = points[1] - points[0];
                auto diff1 = points[2] - points[1];
                return std::acos(std::min<NT>(1.0, std::max<NT>(0.0, diff0.normalized().dot(diff1.normalized()))));
            }
            void setVerticesFromEdges(Vertex centralVertex, const EmbeddedGraph& graph)
            {
                vertices[1] = centralVertex;
                vertices[0] = boost::source(e0, graph) == centralVertex ? boost::target(e0, graph) : boost::source(e0, graph);
                vertices[2] = boost::source(e1, graph) == centralVertex ? boost::target(e1, graph) : boost::source(e1, graph);
                auto getVId = [&graph](auto v) {return boost::get(GCpp::DS::vertex_id_t{}, graph, v); };
                vIds = { getVId(vertices[0]),getVId(vertices[1]) ,getVId(vertices[2]) };
                if(vertices[0] != boost::source(e0,graph))
                {
                    e0 = boost::lookup_edge(vertices[0], vertices[1], graph).first;
                }
                if (vertices[2] != boost::target(e1, graph))
                {
                    e1 = boost::lookup_edge(vertices[1], vertices[2], graph).first;
                }
            }
        };

    public:
        using MergeGroupsContainer = std::map<std::size_t, std::size_t>;
        Degree2VertexPreprocess(NT maxTurningAngle) :m_maxTurningAngle(maxTurningAngle) {}

        /**
         * \brief Sets the maximum turning angle for which edges are collapsed to a single edge.s
         * \param maxTurningAngle 
         */
        void setMaxTurningAngle(const NT& maxTurningAngle)
        {
            m_maxTurningAngle = maxTurningAngle;
        }
        NT maxTurningAngle() const
        {
            return m_maxTurningAngle;
        }
        void setRepeatUntilConvergence(bool val)
        {
            m_repeatUntilConvergence = val;
        }
        bool repeatUntilConvergence() const
        {
            return m_repeatUntilConvergence;
        }

        template<typename SimplificationObserver=TrivialSimplificationObserver>
        void operator()(EmbeddedGraph& inputGraph, SimplificationObserver& observer)
        {
            while(true)
            {
                std::set<Vertex, std::greater<>> verticesToDelete;
                namespace it = GCpp::Helpers::Iterators;
                const auto vertexCount = boost::num_vertices(inputGraph);
                for (auto v : it::range(boost::vertices(inputGraph)))
                {
                    if (boost::out_degree(v, inputGraph) != 2) continue;
                    const auto vId = boost::get(GCpp::DS::vertex_id_t{}, inputGraph, v);
                    // Computation object for the triple of vertices
                    Triple triple;
                    auto edgeIt = boost::out_edges(v, inputGraph).first;
                    triple.e0 = *edgeIt;
                    triple.e1 = *std::next(edgeIt);
                    triple.setVerticesFromEdges(v, inputGraph);
                    auto getEdgeId = [&inputGraph](auto edge)
                    {
                        return boost::get(GCpp::DS::edge_id_t{}, inputGraph, edge);
                    };
                    // Compute turning angle
                    auto turningAngle = triple.computeTurningAngle(inputGraph);
                    if (turningAngle > m_maxTurningAngle) continue;

                    auto[existingEdge, doesExist] = boost::lookup_edge(triple.vertices[0], triple.vertices[2], inputGraph);
                    // Check if the edge we merge to exists: then merge the path to that. (should remove that at
                    // some point, because it implies a local triangle).
                    if (doesExist)
                    {
                        observer.handleEdgeReroute({ getEdgeId(triple.e0), getEdgeId(triple.e1) }, { getEdgeId(existingEdge) });
                        observer.handleVertexMerge(triple.vIds[1], triple.vIds[0]);
                        // Clear edges of the vertex
                        boost::clear_vertex(v, inputGraph);
                        verticesToDelete.insert(v);
                        continue;
                    }

                    // Create new edge
                    auto[edge, wasAdded] = boost::add_edge(triple.vertices[0], triple.vertices[2], inputGraph);
                    boost::put(GCpp::DS::edge_id_t{}, inputGraph, edge, getEdgeId(triple.e0));

                    // Notify observer
                    observer.handleEdgeEdgeMerge(getEdgeId(triple.e1), getEdgeId(triple.e0));
                    observer.handleVertexMerge(triple.vIds[1], triple.vIds[0]);

                    // Clear edges of the vertex
                    boost::clear_vertex(v, inputGraph);
                    verticesToDelete.insert(v);
                }

                std::cout << "Removing " << verticesToDelete.size() << " of " << vertexCount << "\n";

                // Optimization: removing vertices is relatively expensive for adjacency list boost graphs (currently using that).
                // So, if we delete more than 50% of the vertices, we are going to recreate the graph instead
                if(verticesToDelete.size() > vertexCount / 2)
                {
                    auto getVId = [](const auto& graph, auto v)
                    {
                        return boost::get(GCpp::DS::vertex_id_t{}, graph, v);
                    };
                    auto getEId = [](const auto& graph, auto e)
                    {
                        return boost::get(GCpp::DS::edge_id_t{}, graph, e);
                    };
                    EmbeddedGraph temp(boost::num_vertices(inputGraph) - verticesToDelete.size());
                    temp.m_use_canonical_edges_ids = true;
                    using VertIt = decltype(boost::vertices(inputGraph).first);
                    VertIt currentNewVert = boost::vertices(inputGraph).first;
                    std::unordered_map<std::size_t, VertIt> oldIdToNewVertMap;
                    auto vertHasMapping = [&oldIdToNewVertMap](auto vId)
                    {
                        return oldIdToNewVertMap.find(vId) != oldIdToNewVertMap.end();
                    };
                    for(auto v : it::range(boost::vertices(inputGraph)))
                    {
                        // Skip deleted
                        if (verticesToDelete.find(v) != verticesToDelete.end()) continue;
                        const auto vId = getVId(inputGraph, v);
                        oldIdToNewVertMap[vId] = currentNewVert;
                        boost::put(GCpp::DS::vertex_id_t{}, temp, *currentNewVert, vId);
                        boost::put(GCpp::DS::vertex_location_t{}, temp, *currentNewVert, GCpp::DS::get_vertex_location(v, inputGraph));
                        ++currentNewVert;
                    }
                    // These should all still map somewhere.
                    for(auto e: it::range(boost::edges(inputGraph)))
                    {
                        const auto srcId = getVId(inputGraph, boost::source(e, inputGraph));
                        const auto targetId = getVId(inputGraph, boost::target(e, inputGraph));
                        if (!vertHasMapping(srcId) || !vertHasMapping(targetId)) throw std::runtime_error("Invalid mapping!");
                        auto [createdEdge, wasCreated]= GCpp::DS::safe_add_edge(*oldIdToNewVertMap[srcId], *oldIdToNewVertMap[targetId], temp);
                        if (!wasCreated) {
                            if(GCpp::DS::getEdgeId(inputGraph,createdEdge) != getEId(inputGraph, e))
                            {
                                throw std::runtime_error("Edge ID incongruency");
                            }
                            continue;
                            //throw std::runtime_error("Duplicate edge?");
                        }
                        boost::put(GCpp::DS::edge_id_t{}, temp, createdEdge, getEId(inputGraph, e));
                    }
                    // Copy properties by hand...
                    const auto maxVId = boost::get_property(inputGraph, GCpp::DS::next_vertex_id_t{});
                    const auto maxEId = boost::get_property(inputGraph, GCpp::DS::next_edge_id_t{});
                    inputGraph = std::move(temp);
                    boost::get_property(inputGraph, GCpp::DS::next_vertex_id_t{}) = maxVId;
                    boost::get_property(inputGraph, GCpp::DS::next_edge_id_t{}) = maxEId;

                }
                else
                {
                    for (auto v : verticesToDelete)
                    {
                        boost::remove_vertex(v, inputGraph);
                    }
                }
                if((m_repeatUntilConvergence && verticesToDelete.empty()) || !m_repeatUntilConvergence)
                {
                    break;
                }
            }
        }
    };
}
#endif