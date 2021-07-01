#ifndef SCHEMATLIB_MAPSIMPLIFICATION_REMOVETENDRILS_H
#define SCHEMATLIB_MAPSIMPLIFICATION_REMOVETENDRILS_H
#include <GCpp/DS/BoostEmbeddedGraph.h>
#include "SimplficiationObserver.h"
namespace SchematLib::MapSimplification
{
    template<typename EmbeddedGraph, typename SimplificationObserver= TrivialSimplificationObserver>
    class RemoveTendrils
    {
    public:
        using Vertex = typename boost::graph_traits<EmbeddedGraph>::vertex_descriptor;
        using Edge = typename boost::graph_traits<EmbeddedGraph>::edge_descriptor;
        using Point = decltype(GCpp::DS::get_vertex_location(std::declval<Vertex>(), std::declval<const EmbeddedGraph&>()));
        using NT = decltype(std::declval<const Point&>().x());
    private:
        NT m_maxLength;
    public:
        void setMaxLength(NT maxLength)
        {
            m_maxLength = maxLength;
        }
        NT maxLength() const
        {
            return m_maxLength;
        }

        /**
         * \brief
         * NOTE: we are assuming undirected graph here.
         * \param graph 
         */
        void operator()(EmbeddedGraph& graph, SimplificationObserver& observer) const
        {
            namespace it = GCpp::Helpers::Iterators;
            // Order in descending order to not invalidate descriptors when deleting.
            // TODO: this is very specific to boost, generalize this at some point.
            std::set<Vertex, std::greater<>> verticesToDelete;
            for (auto v : it::range(boost::vertices(graph)))
            {
                const auto vDegree = boost::out_degree(v, graph);
                if (verticesToDelete.find(v) != verticesToDelete.end()) continue;
                // Just delete if degree is zero (noise).
                if (vDegree == 0)
                {
                    verticesToDelete.insert(v);
                    continue;
                }
                // Not a leaf, continue
                if (vDegree != 1) continue;

                // Walk along the path starting from v until a degree > 2 vertex is encountered. Keep track of the length traverse,
                // if this is less than maxlength, delete all the vertices
                NT traversedLength = 0;
                std::vector<Vertex> traversedVertices;
                std::vector<Edge> traversedEdges;
                traversedVertices.push_back(v);
                Edge currentE = *boost::out_edges(v, graph).first;
                std::optional<std::size_t> endV;
                while (true)
                {
                    auto src = boost::source(currentE, graph);
                    auto target = boost::target(currentE, graph);
                    auto newV = src == traversedVertices.back() ? target : src;
                    NT length = GCpp::DS::get_linear_edge_length(currentE, graph);
                    traversedLength += length;
                    traversedEdges.push_back(currentE);
                    if (traversedLength > m_maxLength) {
                        traversedVertices.clear();
                        break;
                    }
                    if (boost::out_degree(newV, graph) > 2)
                    {
                        endV = boost::get(boost::vertex_index_t{}, graph, newV);
                        break;
                    }
                    if (boost::out_degree(newV, graph) <= 1)
                    {
                        traversedVertices.push_back(newV);
                        break;
                    }
                    // Find new edge that is not the edge that we just traversed.
                    for (auto edge : it::range(boost::out_edges(newV, graph)))
                    {
                        if (boost::source(edge, graph) != traversedVertices.back() && boost::target(edge, graph) != traversedVertices.back())
                        {
                            traversedVertices.push_back(newV);
                            currentE = edge;
                            break;
                        }
                    }
                    if (traversedVertices.back() != newV) throw std::runtime_error("Failed to find next edge in preprocessing");
                }
                if (!traversedVertices.empty())
                {
                    // Delete later to not invalidate the iterators we are using.
                    for (auto deleteV : traversedVertices) verticesToDelete.insert(deleteV);
                    std::optional<std::size_t> endVId;
                    if(endV.has_value())
                    {
                        endVId = boost::get(GCpp::DS::vertex_id_t{}, graph, endV.value());
                    }
                    // Notify observer that we are collapsing/deleting edges.
                    for(auto deleteE: traversedEdges)
                    {
                        auto edgeIndex = boost::get(GCpp::DS::edge_id_t{}, graph, deleteE);
                        
                        if(endV.has_value())
                        {
                            observer.handleEdgeCollapse(edgeIndex, endVId.value());
                        }
                        else
                        {
                            observer.handleEdgeDelete(edgeIndex);
                            
                        }
                    }
                    for (auto delV : traversedVertices)
                    {
                        const auto vId = boost::get(GCpp::DS::vertex_id_t{}, graph, delV);
                        if(!endVId.has_value() || endVId.value() != vId)
                        {
                            observer.handleVertexDelete(boost::get(GCpp::DS::vertex_id_t{}, graph, delV));
                        }
                    }
                }
            }
            // Erase all vertices to delete
            for (auto v : verticesToDelete)
            {
                boost::clear_vertex(v, graph);
                boost::remove_vertex(v, graph);
            }
        }
    };
}
#endif