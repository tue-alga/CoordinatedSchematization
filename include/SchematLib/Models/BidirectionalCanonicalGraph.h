#ifndef SCHEMATLIB_MODELS_BIDIRECTIONALCANONICALGRAPH_H
#define SCHEMATLIB_MODELS_BIDIRECTIONALCANONICALGRAPH_H
#include <GCpp/DS/BoostEmbeddedGraph.h>
#include <GCpp/Math/Vector.h>
#include "SchematLib/MapSimplification/SpecialEdgeIdFunctions.h"
#include "BaseTypes.h"
#include <iostream>
#define FORWARD_DEFINE(name) using name = SchematLib::Models::BidirectionalCanonicalGraph::GT::name


namespace GCpp::DS
{
    template<>
    inline auto compute_length<GCpp::Math::Vec2<SchematLib::Models::NT>>(const GCpp::Math::Vec2<SchematLib::Models::NT>& p0, const GCpp::Math::Vec2<SchematLib::Models::NT>& p1)->decltype(p0.x())
    {
        return (p0 - p1).length();
    }
}

namespace SchematLib::Models
{
    struct BidirectionalCanonicalGraph
    {
        using NT = long double;
        using Point = GCpp::Math::Vec2<NT>;
        using BaseGraph = GCpp::DS::BoostEmbeddedGraph<Point>;
        using GT = boost::graph_traits<BaseGraph>;
        FORWARD_DEFINE(vertex_descriptor);
        FORWARD_DEFINE(edge_descriptor);
        FORWARD_DEFINE(directed_category);
        FORWARD_DEFINE(edge_parallel_category);
        FORWARD_DEFINE(traversal_category);
        FORWARD_DEFINE(out_edge_iterator);
        FORWARD_DEFINE(degree_size_type);
        FORWARD_DEFINE(in_edge_iterator);
        FORWARD_DEFINE(vertex_iterator);
        FORWARD_DEFINE(vertices_size_type);
        FORWARD_DEFINE(vertices_size_type);
        //FORWARD_DEFINE(edge_iterator);
        FORWARD_DEFINE(edges_size_type);

        using graph_property_type = BaseGraph::graph_property_type;
        using edge_property_type = BaseGraph::edge_property_type;
        using vertex_property_type = BaseGraph::vertex_property_type;
        using vertex_bundled = BaseGraph::vertex_bundled;
        using edge_bundled = BaseGraph::edge_bundled;
        using graph_bundled = BaseGraph::graph_bundled;
        BaseGraph m_graph;
        // If true, boost interface functions will return all directed edges with canonical edge IDs.
        // Otherwise, only half of the edges will be returned.
        bool m_use_canonical_edges_ids = false;
        BidirectionalCanonicalGraph(){}
        explicit BidirectionalCanonicalGraph(std::size_t vertexCount):m_graph(vertexCount){}

        class custom_edge_iterator : public GT::edge_iterator
        {
        public:
            using BaseIt = GT::edge_iterator;
            custom_edge_iterator(const BaseGraph* graph, GT::edge_iterator baseIt, GT::edge_iterator endIt);

            custom_edge_iterator(const BaseGraph* graph, GT::edge_iterator baseIt, GT::edge_iterator endIt,
                                 bool uniqueEdges);
            custom_edge_iterator();

            custom_edge_iterator& operator++();
        private:
            void toNextCanonical();
            const BaseGraph* m_graph = nullptr;
            bool m_filterUniqueEdges = false;
            GT::edge_iterator endIt;
        };

        using edge_iterator = custom_edge_iterator;

        std::pair<edge_iterator, edge_iterator> edges() const;

        edge_descriptor twin(edge_descriptor e) const;

        std::size_t number_of_edges() const;
    };
}
// Define the boost operators
namespace boost
{
#define G SchematLib::Models::BidirectionalCanonicalGraph
    template<typename Property>
    struct property_map<G, Property>
    {
        using type = typename property_map<G::BaseGraph, Property>::type;
        using const_type = typename property_map<G::BaseGraph, Property>::const_type;
    };
    // To be safe
    template<>
    struct graph_traits<SchematLib::Models::BidirectionalCanonicalGraph>
    {
        FORWARD_DEFINE(vertex_descriptor);
        FORWARD_DEFINE(edge_descriptor);
        FORWARD_DEFINE(directed_category);
        FORWARD_DEFINE(edge_parallel_category);
        FORWARD_DEFINE(traversal_category);
        FORWARD_DEFINE(out_edge_iterator);
        FORWARD_DEFINE(degree_size_type);
        FORWARD_DEFINE(in_edge_iterator);
        FORWARD_DEFINE(vertex_iterator);
        FORWARD_DEFINE(vertices_size_type);
        FORWARD_DEFINE(vertices_size_type);
        FORWARD_DEFINE(edges_size_type);
        //FORWARD_DEFINE(edge_iterator);
        using edge_iterator = G::edge_iterator;
        static auto null_vertex() { return G::GT::null_vertex(); }
    };

    // 
    auto inline out_edges(G::vertex_descriptor v, const G& g)
    {
        return boost::out_edges(v,g.m_graph);
    }
    auto inline source(G::edge_descriptor e, const G& g)
    {
        return boost::source(e, g.m_graph);
    }
    auto inline target(G::edge_descriptor e, const G& g)
    {
        return boost::target(e, g.m_graph);
    }
    auto inline out_degree(G::vertex_descriptor v, const G& g)
    {
        return boost::out_degree(v, g.m_graph);
    }
    //
    auto inline in_edges(G::vertex_descriptor v, const G& g)
    {
        return boost::in_edges(v, g.m_graph);
    }
    auto inline in_degree(G::vertex_descriptor v, const G& g)
    {
        return boost::out_degree(v, g.m_graph);
    }
    auto inline degree(G::vertex_descriptor v, const G& g)
    {
        return boost::degree(v, g.m_graph);
    }
    //
    auto inline vertices(const G& g)
    {
        return boost::vertices(g.m_graph);
    }
    auto inline num_vertices(const G& g)
    {
        return boost::num_vertices(g.m_graph);
    }
    //
    auto inline edges(const G& g)
    {
        return g.edges();
    }
    auto inline num_edges(const G& g)
    {
        return g.number_of_edges();
    }
    //
    auto inline add_vertex(G& g)
    {
        return boost::add_vertex(g.m_graph);
    }
    auto inline clear_vertex(G::vertex_descriptor v, G& g)
    {
        boost::clear_vertex(v,g.m_graph);
    }
    auto inline remove_vertex(G::vertex_descriptor v,G& g)
    {
        boost::remove_vertex(v, g.m_graph);
    }
    auto inline add_edge(G::vertex_descriptor s, G::vertex_descriptor t, G& g)
    {
        boost::add_edge(t,s, g.m_graph);
        return boost::add_edge(s, t, g.m_graph);
    }
    auto inline remove_edge(G::vertex_descriptor s, G::vertex_descriptor t, G& g)
    {
        auto[otherE, exists] = boost::lookup_edge(t,s, g.m_graph);
        if (!exists) throw std::runtime_error("invalid edge");
        boost::remove_edge(otherE, g.m_graph);
        boost::remove_edge(s,t, g.m_graph);
    }
    auto inline remove_edge(G::edge_descriptor e,  G& g)
    {
        auto [otherE, exists] = boost::lookup_edge(boost::target(e, g.m_graph), boost::source(e, g.m_graph), g.m_graph);
        if (!exists) throw std::runtime_error("invalid edge");
        boost::remove_edge(otherE, g.m_graph);
        boost::remove_edge(e, g.m_graph);
    }
    template<typename Prop>
    auto inline get(Prop p, const G& g)
    {
        return boost::get(p, g.m_graph);
    }
    template<typename Prop>
    auto inline get(Prop p, G& g)
    {
        return boost::get(p, g.m_graph);
    }
    template<typename Prop,typename GraphEl>
    auto inline get(Prop p, const G& g, GraphEl el)
    {
        return boost::get(p, g.m_graph,el);
    }
    template<typename Prop, typename GraphEl>
    auto inline get(Prop p, G& g, GraphEl el)
    {
        return boost::get(p, g.m_graph,el);
    }
    template<typename Prop, typename GraphEl,typename Val>
    inline void put(Prop p, G& g, GraphEl el, const Val& val)
    {
        boost::put(p, g.m_graph, el,val);
    }
    template<>
    inline void put(GCpp::DS::edge_id_t, G& g, G::edge_descriptor e, const std::size_t& val)
    {
        using SEF = SchematLib::MapSimplification::SpecialEdgeIdFunctions;
        // Always use canonical ID's here. 
        boost::put(GCpp::DS::edge_id_t{}, g.m_graph, e, val);
        //std::cout << "Assigned " << val << " and " << SEF::flipCanonical(val) << " |" << SEF::NONCANONICAL_BIT<< "\n";
        boost::put(GCpp::DS::edge_id_t{}, g.m_graph, g.twin(e), SEF::flipCanonical(val));
    }
    template<>
    inline auto get(GCpp::DS::edge_id_t, const G& g, G::edge_descriptor e)
    {
        if(g.m_use_canonical_edges_ids)
        {
            return boost::get(GCpp::DS::edge_id_t{}, g.m_graph, e);
        }
        else
        {
            auto id = boost::get(GCpp::DS::edge_id_t{}, g.m_graph, e);
            return SchematLib::MapSimplification::SpecialEdgeIdFunctions::getEdgeId(id);
        }
    }
    // More properties...
    template < typename Tag >
    typename boost::graph_property<G::BaseGraph, Tag >::type& get_property(
        G& g, Tag)
    {
        return boost::get_property(g.m_graph, Tag{});
    }

    template < typename Tag >
    typename boost::graph_property<G::BaseGraph, Tag >::type const& get_property(
        G const& g, Tag)
    {
        return boost::get_property(g.m_graph, Tag{});
    }
}
#undef G
#undef FORWARD_DEFINE
#endif
