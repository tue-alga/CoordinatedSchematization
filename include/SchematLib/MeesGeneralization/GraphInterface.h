#ifndef SCHEMATLIB_MEESGENERALIZATION_GRAPHINTERFACE_H
#define SCHEMATLIB_MEESGENERALIZATION_GRAPHINTERFACE_H
#include <unordered_map>

namespace SchematLib::MeesGeneralization
{
    /**
     * \brief Functor for acquiring a tagged property.
     * Override for the internal representation of the graph
     * \tparam WrappedGraph The internal representation of the graph
     * \tparam PropertyTag Tag of the property to acquire
     * \tparam KeyType Type of the key to use for acquiring the property
     */
    template<typename WrappedGraph, typename PropertyTag, typename KeyType>
    struct GetProperty
    {
        using Type = void;
        // Stub implementation, specialize with your own representation
        Type operator()(const WrappedGraph& graph, const KeyType& key)
        {
            return {};
        }
    };

    template<typename WrappedGraph, typename NT>
    struct GraphTypes
    {
        using Edge = void;
        using Vertex = void;
        using EdgeIterable = void;
        using OutEdgeIterable = void;
        using InEdgeIterable = void;
    };

    /**
     * \brief Archetype for graph needed in Generalization framework
     * \tparam WrappedGraph Internal representation of the graph. Graph is assumed to not modify this.
     * \tparam NT Number type to use for scalar properties
     */
    template<typename WrappedGraph, typename NT>
    struct Graph
    {
        // Reference to the internal representation
        const WrappedGraph& m_graph;

        // Types
        using GTs = GraphTypes<WrappedGraph, NT>;
        using Edge = typename GTs::Edge;
        using Vertex = typename GTs::Vertex;
        using count = std::size_t;


        Graph(const WrappedGraph& graph):m_graph(graph){}

        count num_edges()const;

        Vertex source(Edge e) const;
        Vertex target(Edge e) const;

        Edge edge(std::size_t index) const;

        typename GTs::OutEdgeIterable out_edges(Vertex v) const;

        typename GTs::InEdgeIterable in_edges(Vertex v) const;

        typename GTs::EdgeIterable edges() const;

        template<typename Tag, typename GraphElementDescriptor>
        auto getProperty(Tag, GraphElementDescriptor index)const
        {
            return GetProperty<WrappedGraph, Tag, GraphElementDescriptor>{}(m_graph, index);
        }
    };
    template<typename WrappedGraph>
    using SimpleEdgePath = std::vector <typename WrappedGraph::Edge>;

    // Tags for properties
    struct edge_length_t {};
    struct edge_index_t {};
}

#endif