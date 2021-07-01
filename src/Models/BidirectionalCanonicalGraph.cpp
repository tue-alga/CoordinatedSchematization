#include "SchematLib/Models/BidirectionalCanonicalGraph.h"

namespace SchematLib::Models
{
    BidirectionalCanonicalGraph::custom_edge_iterator::custom_edge_iterator(const BaseGraph* graph,
                                                                            GT::edge_iterator baseIt,
                                                                            GT::edge_iterator endIt):
        GT::edge_iterator(baseIt),
        m_graph(graph), endIt(endIt)
    {
        toNextCanonical();
    }

    BidirectionalCanonicalGraph::custom_edge_iterator::custom_edge_iterator(const BaseGraph* graph,
                                                                            GT::edge_iterator baseIt,
                                                                            GT::edge_iterator endIt,
                                                                            bool uniqueEdges):
        GT::edge_iterator(baseIt),
        m_graph(graph), m_filterUniqueEdges(uniqueEdges), endIt(endIt)
    {
        toNextCanonical();
    }

    BidirectionalCanonicalGraph::custom_edge_iterator::custom_edge_iterator()
    {
    }

    BidirectionalCanonicalGraph::custom_edge_iterator& BidirectionalCanonicalGraph::custom_edge_iterator::operator++()
    {
        GT::edge_iterator::operator++();
        toNextCanonical();
        return *this;
    }

    void BidirectionalCanonicalGraph::custom_edge_iterator::toNextCanonical()
    {
        if (!m_filterUniqueEdges) return;
        if (*this == endIt) return;

        using SEF = MapSimplification::SpecialEdgeIdFunctions;
        auto eId = boost::get(GCpp::DS::edge_id_t{}, *m_graph, this->operator*());
        while(SEF::isNonCanonical(eId))
        {
            BaseIt::operator++();
            if (*this == endIt) break;
            eId = boost::get(GCpp::DS::edge_id_t{}, *m_graph, this->operator*());
        }
    }

    std::pair<BidirectionalCanonicalGraph::edge_iterator, BidirectionalCanonicalGraph::edge_iterator>
    BidirectionalCanonicalGraph::edges() const
    {
        auto its = boost::edges(m_graph);
        edge_iterator begin(&m_graph, its.first, its.second, !m_use_canonical_edges_ids);
        edge_iterator end(&m_graph, its.second, its.second, !m_use_canonical_edges_ids);
        return std::make_pair(begin, end);
    }

    BidirectionalCanonicalGraph::edge_descriptor BidirectionalCanonicalGraph::twin(edge_descriptor e) const
    {
        return boost::lookup_edge(boost::target(e, m_graph), boost::source(e, m_graph), m_graph).first;
    }

    std::size_t BidirectionalCanonicalGraph::number_of_edges() const
    {
        if (m_use_canonical_edges_ids) return boost::num_edges(m_graph);
        else
        {
            auto baseEdges = boost::num_edges(m_graph);
            assert(baseEdges % 2 == 0); // Require even number of edges
            return baseEdges/2;
        }
    }
}
