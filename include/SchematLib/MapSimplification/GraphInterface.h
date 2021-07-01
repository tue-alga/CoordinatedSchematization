#ifndef SCHEMATLIB_MAPSIMPLIFICATION_GRAPHINTERFACE_H
#define SCHEMATLIB_MAPSIMPLIFICATION_GRAPHINTERFACE_H
#include <tuple>
namespace SchematLib::MapSimplification
{
    template<typename GraphType, typename Enabled = void>
    struct GraphTypes
    {
        
    };

    template<typename GraphType, typename Enabled=void>
    struct GraphInterface
    {
        using SelfType = GraphInterface<GraphType, Enabled>;
        const GraphType& m_graph;
        using GT = GraphTypes<GraphType, Enabled>;
        using Vertex = typename GT::Vertex;
        using VertexIterator = typename GT::VertexIterator;
        using VertexIterable = typename GT::VertexIterable;
        using EdgeIterator = typename GT::EdgeIterator;
        using EdgeIterable = typename GT::EdgeIterable;
        using Edge = typename GT::Edge;
        using VertexId = typename GT::VertexId;
        using EdgeId = typename GT::EdgeId;
        using Point = typename GT::Point;

        GraphInterface(const GraphType& graph):m_graph(graph){}

        std::pair<Edge, bool> edgeExists(Vertex src, Vertex target) const{}

        EdgeId edgeId(Edge e) const {}

        VertexId vertexId(Vertex v) const{}
        void setNewStartEdgeId(EdgeId value){}
        void setNewStartVertexId(VertexId value) {}
        Point vertexLocation(Vertex v) const{}
        std::size_t numberOfEdges() const{}
        std::size_t numberOfVertices() const {}
        VertexIterable vertices() const{}
        EdgeIterable edges() const {}

    };
}
#endif