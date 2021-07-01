#ifndef SCHEMATLIB_MAPSIMPLIFICATION_GCPPGRAPHINTERFACEIMPL_H
#define SCHEMATLIB_MAPSIMPLIFICATION_GCPPGRAPHINTERFACEIMPL_H
#include <tuple>
#include <boost/graph/lookup_edge.hpp>

#include "GraphInterface.h"
#include <GCpp/DS/BoostEmbeddedGraph.h>
namespace SchematLib::MapSimplification
{
    namespace detail
    {
        template<typename T, bool IsConst>
        struct MaybeConst;

        template<typename T>
        struct MaybeConst<const T,true>
        {
            static constexpr bool is_const = true;
            using type = std::decay_t<T>;
            const T& m_ref;
            MaybeConst(const T& ref):m_ref(ref){}
        };
        template<typename T>
        struct MaybeConst<T, false>
        {
            static constexpr bool is_const = true;
            using type = std::decay_t<T>;
            T& m_ref;
            MaybeConst(T& ref) :m_ref(ref) {}
        };
    }
    template<typename Point_>
    struct GraphTypes<GCpp::DS::BoostUndirectedEmbeddedGraph<Point_>,void>
    {
        template<typename It>
        using Range = GCpp::Helpers::Iterators::Range<It>;

        using Graph = GCpp::DS::BoostUndirectedEmbeddedGraph<Point_>;
        using GT = boost::graph_traits<Graph>;
        using Vertex = typename GT::vertex_descriptor;
        using VertexIterator = typename GT::vertex_iterator;
        using VertexIterable = Range<VertexIterator>;
        using EdgeIterator = typename GT::edge_iterator;
        using EdgeIterable = Range<EdgeIterator>;
        using Edge = typename GT::edge_descriptor;
        using VertexId = std::size_t;
        using EdgeId = std::size_t;
        using Point = Point_;
    };
    template<typename Point_>
    struct GraphInterface<GCpp::DS::BoostUndirectedEmbeddedGraph<Point_>, void>
    {
        using GraphType = GCpp::DS::BoostUndirectedEmbeddedGraph<Point_>;
        using SelfType = GraphInterface<GraphType, void>;
        GraphType& m_graph;
        using GT = GraphTypes<GraphType, void>;
        using Vertex = typename GT::Vertex;
        using VertexIterator = typename GT::VertexIterator;
        using VertexIterable = typename GT::VertexIterable;
        using EdgeIterator = typename GT::EdgeIterator;
        using EdgeIterable = typename GT::EdgeIterable;
        using Edge = typename GT::Edge;
        using VertexId = typename GT::VertexId;
        using EdgeId = typename GT::EdgeId;
        using Point = typename GT::Point;

        GraphInterface(GraphType& graph) :m_graph(graph) {}

        std::pair<Edge, bool> edgeExists(Vertex src, Vertex target) const
        {
            return boost::lookup_edge(src, target, m_graph);
        }

        EdgeId edgeId(Edge e) const
        {
            return boost::get(GCpp::DS::edge_id_t{}, m_graph, e);
        }

        VertexId vertexId(Vertex v) const
        {
            return boost::get(GCpp::DS::vertex_id_t{}, m_graph, v);
        }
        void setNewStartEdgeId(EdgeId value)
        {
            boost::put(m_graph, GCpp::DS::next_edge_id_t{}, value);
        }
        void setNewStartVertexId(VertexId value)
        {
            boost::put(m_graph, GCpp::DS::next_vertex_id_t{}, value);
        }
        Point vertexLocation(Vertex v) const
        {
            return GCpp::DS::get_vertex_location(v, m_graph);
        }
        std::size_t numberOfEdges() const
        {
            return boost::num_edges(m_graph);
        }
        std::size_t numberOfVertices() const
        {
            return boost::num_vertices(m_graph);
        }
        VertexIterable vertices() const
        {
            return GCpp::Helpers::Iterators::range(boost::vertices(m_graph));
        }
        EdgeIterable edges() const
        {
            return GCpp::Helpers::Iterators::range(boost::edges(m_graph));
        }

    };
}
#endif