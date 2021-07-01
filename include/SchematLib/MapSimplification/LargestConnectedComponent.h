#ifndef SCHEMATLIB_MAPSIMPLIFICATION_LARGESTCONNECTEDCOMPONENT_H
#define SCHEMATLIB_MAPSIMPLIFICATION_LARGESTCONNECTEDCOMPONENT_H
#include <boost/graph/lookup_edge.hpp>
#include <GCpp/DS/BoostEmbeddedGraph.h>
#include <boost/graph/connected_components.hpp>
#include "SimplficiationObserver.h"

namespace SchematLib::MapSimplification
{
    template<typename EmbeddedGraph>
    class LargestConnectedComponent
    {
    public:
        using Vertex = typename boost::graph_traits<EmbeddedGraph>::vertex_descriptor;
        using Edge = typename boost::graph_traits<EmbeddedGraph>::edge_descriptor;
        using Point = decltype(GCpp::DS::get_vertex_location(std::declval<Vertex>(), std::declval<const EmbeddedGraph&>()));
        using NT = decltype(std::declval<const Point&>().x());
    public:
        LargestConnectedComponent() {}

        template<typename SimplificationObserver = TrivialSimplificationObserver>
        void operator()(EmbeddedGraph& inputGraph, SimplificationObserver& observer)
        {
            // Stored by vertex index
            std::vector< int > component(boost::num_vertices(inputGraph));
            int num = boost::connected_components(inputGraph, &component[0]);

            std::map<int, std::size_t> componentCounts;
            std::size_t maxCount = 0;
            int maxComp = -1;
            for(auto i = 0; i < component.size(); ++i)
            {
                const auto comp = component[i];
                if (componentCounts.find(comp) == componentCounts.end()) componentCounts[comp] = 0;
                ++componentCounts[comp];
                if(componentCounts[comp] > maxCount)
                {
                    maxCount = componentCounts[comp];
                    maxComp = comp;
                }
            }
            std::set<std::size_t, std::greater<>> verticesToDelete;
            for(std::size_t i = 0; i < component.size(); ++i)
            {
                if (component[i] == maxComp) continue;
                observer.handleVertexDelete(GCpp::DS::getVertexId(inputGraph, i));
                verticesToDelete.insert(i);
            }
            for(auto e : GCpp::Helpers::Iterators::range(boost::edges(inputGraph)))
            {
                if (component[boost::source(e, inputGraph)] == maxComp) continue;
                const auto eId = GCpp::DS::getEdgeId(inputGraph, e);
                if(SpecialEdgeIdFunctions::isCanonical(eId))
                {
                    observer.handleEdgeDelete(GCpp::DS::getEdgeId(inputGraph, e));
                }
            }
            while(!verticesToDelete.empty())
            {
                auto front = *verticesToDelete.begin();
                verticesToDelete.erase(verticesToDelete.begin());
                boost::clear_vertex(front, inputGraph);
                boost::remove_vertex(front, inputGraph);
            }
        }
    };
}
#endif