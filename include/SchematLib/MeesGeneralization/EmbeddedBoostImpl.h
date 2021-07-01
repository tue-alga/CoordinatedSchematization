#ifndef SCHEMATLIB_MEESGENERALIZATION_EMBEDDEDBOOSTIMPL_H
#define SCHEMATLIB_MEESGENERALIZATION_EMBEDDEDBOOSTIMPL_H
#include <GCpp/DS/BoostEmbeddedGraph.h>
#include <GCpp/Helpers/Iterators.h>
#include <GCpp/Math/Vector.h>

#include "GraphInterface.h"
namespace SchematLib::MeesGeneralization
{
    template<typename NT, typename BoostEmbeddedGraphT>
    struct GraphWithEdgeLength
    {
        using Graph_t = BoostEmbeddedGraphT;
        using GT = typename boost::graph_traits<Graph_t>;
        using Edge = typename boost::graph_traits<Graph_t>::edge_descriptor;
        Graph_t graph;
        std::vector<NT> m_edgeLengths;
        // Fast lookup of edges, at the expense of duplicate data.
        std::map<std::size_t, Edge> m_edges; //Edge ID to edge
        std::vector<std::size_t> m_edgeIndexToId; //Indices are contiguous
        std::map<std::size_t, std::size_t> m_edgeIdToIndex;

        Graph_t& boostGraph() { return graph; }
        const Graph_t& boostGraph() const { return graph; }

        NT edgeLengthForEdgeId(std::size_t id) const
        {
            return m_edgeLengths[m_edgeIdToIndex.at(id)];
        }
        NT edgeLengthForEdgeIndex(std::size_t index) const
        {
            return m_edgeLengths[index];
        }
        Edge edge(std::size_t index) const
        {
            return m_edges.at(m_edgeIndexToId.at(index));
        }
        Edge edgeById(std::size_t id) const
        {
            return m_edges.at(id);
        }

        std::size_t edgeIdForIndex(std::size_t index) const
        {
            return m_edgeIndexToId.at(index);
        }
        std::size_t edgeIndexForId(std::size_t id) const
        {
            return m_edgeIdToIndex.at(id);
        }

        void computeEdgeMapping()
        {
            m_edges.clear();
            m_edgeIndexToId.clear();
            m_edgeIdToIndex.clear();

            m_edgeIndexToId.resize(boost::num_edges(graph));
            std::size_t currIndex = 0;
            for (auto e : GCpp::Helpers::Iterators::range(boost::edges(graph)))
            {
                const auto id = boost::get(GCpp::DS::edge_id_t{}, graph, e);
                m_edges[boost::get(GCpp::DS::edge_id_t{}, graph, e)] = e;
                m_edgeIndexToId[currIndex] = id;
                m_edgeIdToIndex[id] = currIndex;
                ++currIndex;
            }
        }

        void computeEdgeLengths()
        {
            m_edgeLengths.clear();
            m_edgeLengths.reserve(boost::num_edges(graph));
            for(const auto& e: GCpp::Helpers::Iterators::range(boost::edges(graph)))
            {
                const auto src = boost::source(e, graph);
                const auto target= boost::target(e, graph);
                m_edgeLengths.push_back(
                    (GCpp::DS::get_vertex_location(src, graph) - GCpp::DS::get_vertex_location(target, graph)).length()
                );
            }
        }
    };

    template<typename BoostGraph>
    inline bool verifyContiguousEdgeIndices(const BoostGraph& graph)
    {
        std::unordered_set<std::size_t> seenIndices;
        std::size_t minIndex = std::numeric_limits<std::size_t>::max();
        std::size_t maxIndex = 0;
        std::size_t edgeCount = 0;
        for(auto e : GCpp::Helpers::Iterators::range(boost::edges(graph)))
        {
            std::size_t index = boost::get(edge_index_t{}, graph, e);
            minIndex = std::min(minIndex, index);
            maxIndex = std::max(maxIndex, index);
            seenIndices.insert(index);
            ++edgeCount;
        }
        return maxIndex == edgeCount - 1 && minIndex == 0 && seenIndices.size() == edgeCount;
    }

    /**
     * \brief By default, return boost property.
     * \tparam WrappedGraph
     * \tparam PropertyTag
     */
    template<typename NT, typename BoostEmbeddedGraphT, typename Prop, typename KeyType>
    struct GetProperty<GraphWithEdgeLength<NT, BoostEmbeddedGraphT>, Prop, KeyType>
    {
        auto operator()(const GraphWithEdgeLength<NT, BoostEmbeddedGraphT>& graph, const KeyType& key)
        {
            return boost::get(Prop{}, graph.graph, key);
        }
    };

    /**
     * \brief Return edge_length
     * \tparam WrappedGraph
     * \tparam PropertyTag
     */
    template<typename NT, typename BoostEmbeddedGraphT>
    struct GetProperty<GraphWithEdgeLength<NT, BoostEmbeddedGraphT>, edge_length_t, typename GraphWithEdgeLength<NT, BoostEmbeddedGraphT>::Edge>
    {
        using Type = NT;
        using KeyType = typename GraphWithEdgeLength<NT, BoostEmbeddedGraphT>::Edge;
        Type operator()(const GraphWithEdgeLength<NT, BoostEmbeddedGraphT>& graph, const KeyType& key)
        {
            std::size_t idx = boost::get(GCpp::DS::edge_id_t{}, graph.graph, key);
            return graph.edgeLengthForEdgeId(idx);
        }
    };
    template<typename NT, typename BoostEmbeddedGraphT>
    struct GetProperty<GraphWithEdgeLength<NT, BoostEmbeddedGraphT>, edge_index_t, typename GraphWithEdgeLength<NT, BoostEmbeddedGraphT>::Edge>
    {
        using Type = std::size_t;
        using KeyType = typename GraphWithEdgeLength<NT, BoostEmbeddedGraphT>::Edge;
        Type operator()(const GraphWithEdgeLength<NT, BoostEmbeddedGraphT>& graph, const KeyType& key)
        {
            return graph.edgeIndexForId(boost::get(GCpp::DS::edge_id_t{}, graph.graph, key));
        }
    };

    template<typename NT, typename BoostEmbeddedGraphT>
    struct GraphTypes<GraphWithEdgeLength<NT,BoostEmbeddedGraphT>,NT>
    {
        // Graph traits of boost
        using GT = typename GraphWithEdgeLength<NT, BoostEmbeddedGraphT>::GT;
        using Edge = typename GraphWithEdgeLength<NT, BoostEmbeddedGraphT>::Edge;
        using Vertex = typename GT::vertex_descriptor;
        using EdgeIterable = GCpp::Helpers::Iterators::Range<typename GT::edge_iterator>;
    };

    template<typename NT, typename BoostEmbeddedGraphT>
    struct Graph<GraphWithEdgeLength<NT, BoostEmbeddedGraphT>,NT>
    {
        using WrappedGraph = GraphWithEdgeLength<NT, BoostEmbeddedGraphT>;
        const WrappedGraph& m_graph;

        using GTs = GraphTypes<WrappedGraph, NT>;
        using Edge = typename GTs::Edge;
        using Vertex = typename GTs::Vertex;


        Graph(const WrappedGraph& graph) :m_graph(graph)
        {
            const_cast<WrappedGraph*>(&graph)->computeEdgeMapping();
            const_cast<WrappedGraph*>(&graph)->computeEdgeLengths();
        }

        using count = std::size_t;

        count num_edges()const { return boost::num_edges(m_graph.graph); }

        Edge edge(std::size_t index) const
        {
            return m_graph.edge(index);
        }
        Edge edgeById(std::size_t id) const
        {
            return m_graph.edgeById(id);
        }

        std::size_t edgeIdForIndex(std::size_t index) const
        {
            return m_graph.edgeIdForIndex(index);
        }

        Vertex source(Edge e) const { return boost::source(e, m_graph.graph); }
        Vertex target(Edge e) const { return boost::target(e, m_graph.graph); }

        typename GTs::EdgeIterable edges() const
        {
            return GCpp::Helpers::Iterators::range(boost::edges(m_graph.graph));
        }
        auto out_edges(Vertex v) const
        {
            return GCpp::Helpers::Iterators::range(boost::out_edges(v, m_graph.graph));
        }
        auto in_edges(Vertex v) const
        {
            // TODO: this may break for undirected
            return GCpp::Helpers::Iterators::range(boost::in_edges(v, m_graph.graph));
        }

        template<typename Tag, typename GraphElementDescriptor>
        auto getProperty(Tag, GraphElementDescriptor index)const
        {
            return GetProperty<WrappedGraph, Tag, GraphElementDescriptor>{}(m_graph, index);
        }
    };
}
#endif