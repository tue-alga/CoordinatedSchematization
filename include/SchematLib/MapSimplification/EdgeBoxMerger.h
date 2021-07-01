#ifndef SCHEMATLIB_MAPSIMPLIFICATION_EDGEBOXMERGER_H
#define SCHEMATLIB_MAPSIMPLIFICATION_EDGEBOXMERGER_H
#include <GCpp/DS/BoostEmbeddedGraph.h>
#include <GCpp/Geometry/OrientedRect.h>
#include <GCpp/Geometry/Types.h>
#include <GCpp/Geometry/Intersect.h>
namespace SchematLib::MapSimplification
{
    template<typename EmbeddedGraph, typename EdgeIndexStructure>
    class EdgeBoxMerger
    {
    public:
        using Vertex = typename boost::graph_traits<EmbeddedGraph>::vertex_descriptor;
        using Edge = typename boost::graph_traits<EmbeddedGraph>::edge_descriptor;
        using Point = decltype(GCpp::DS::get_vertex_location(std::declval<Vertex>(), std::declval<const EmbeddedGraph&>()));
        using NT = decltype(std::declval<const Point&>().x());
    private:
        NT m_distance;
        EmbeddedGraph* m_graph = nullptr;

        EmbeddedGraph& graph() { return *m_graph; }

        void findSegmentsOfInterest(const Edge& edge, const EdgeIndexStructure& index,
            std::vector<std::size_t>& edgeOfInterest)
        {
            namespace geom = GCpp::Geometry;
            auto& inputGraph = this->graph();
            // Query index for segments passing through box.
            std::set<std::size_t> edges;

            const auto v0 = GCpp::DS::get_vertex_location(boost::source(edge, inputGraph), inputGraph);
            const auto v1 = GCpp::DS::get_vertex_location(boost::target(edge, inputGraph), inputGraph);

            geom::OrientedRect<NT> rect(v0, v1-v0, (v1-v0).length(), 0);
            rect.inflate(m_distance);

            // Find all segments intersecting the box
            index.intersectingSegments(rect.ccwVertices(), edges);

            geom::Segment seg0(rect.origin(), rect.origion() + rect.yDirection());
            geom::Segment seg1(rect.origin() + rect.xDirection(), rect.origion() + rect.yDirection() + rect.xDirection());

            // Find all segments intersecting the perpendicular boundaries of the box.
            auto curr = edges.begin();
            while (curr != edges.end())
            {
                curr = std::find_if(curr, edges.end(), [&inputGraph,&seg0, &seg1](auto edgeIndex)
                {
                    auto edge = getEdge(edgeIndex);//TODO
                    const auto segment = boost::get(GCpp::DS::edge_geometry_t{}, inputGraph, edge);
                    return geom::intersects(seg0, segment) || geom::intersects(seg1, segment);
                });
                if (curr == edges.end()) break;
                edgeOfInterest.push_back(*curr);

                // Advance one
                curr = std::next(curr);
            }
        }

        struct EdgeDescription
        {
            std::size_t index;
            NT length;
            Edge edge;
            bool operator<(const EdgeDescription& other) const { return length < other.length; }
        };

        void constructEdgeOrder(std::vector<EdgeDescription>& edgeOrder)
        {
            auto& inputGraph = this->graph();
            namespace it = GCpp::Helpers::Iterators;
            for (auto edge : it::range(boost::edges(inputGraph)))
            {
                const auto eIndex = boost::get(GCpp::DS::edge_id_t{}, inputGraph, edge);
                const auto eLength = GCpp::DS::get_linear_edge_length(edge, inputGraph);

                edgeOrder.push_back(EdgeDescription{ eIndex, eLength,edge });
            }
            // Sort according edge length (ascending).
            std::sort(edgeOrder.begin(), edgeOrder.end());
        }

        /**
         * \brief Merges one edge to another. 
         * \param from 
         * \param to 
         */
        void mergeEdges(Edge from, Edge to)
        {
            
        }

        void mergeWithClosestSegment(Edge edge, const std::vector<std::size_t>& targetEdges)
        {
            const auto edgeSegment = boost::get(GCpp::DS::edge_geometry_t{}, graph(), edge);
            // Compute closest distances
            std::vector<NT> distances;
            std::transform(targetEdges.begin(), targetEdges.end(), std::back_inserter(distances), [this,&edgeSegment](const auto& edgeIndex)
            {
                const auto targetEdge = graph().edge(edgeIndex);
                const auto targetEdgeSegment = boost::get(GCpp::DS::edge_geometry_t{}, graph(), targetEdge);
                return distance(edgeSegment, targetEdgeSegment);
            });
            // Find closest edge
            auto minElementIt = std::min_element(distances.begin(), distances.end());
            auto edgeToMerge = targetEdges[std::distance(distances.begin(), minElementIt)];
            // Merge
            mergeEdges(edge, graph().getEdge(edgeToMerge));
        }
    public:
        using MergeGroupsContainer = std::map<std::size_t, std::size_t>;
        EdgeBoxMerger(NT distance) :m_distance(distance){}

        void setDistance(const NT& distance)
        {
            m_distance = distance;
        }
        NT distance() const
        {
            return m_distance;
        }

        void operator()(const EmbeddedGraph& inputGraph, const EdgeIndexStructure& index, EmbeddedGraph& outputGraph)
        {
            outputGraph = inputGraph;
            this->operator()(outputGraph, index);
        }

        void operator()(EmbeddedGraph& inputGraph, const EdgeIndexStructure& index)
        {
            m_graph = &inputGraph;
            namespace it =GCpp::Helpers::Iterators;

            // Get the edge data in order of processing
            std::vector<EdgeDescription> edgeOrderIndices;
            constructEdgeOrder(inputGraph, edgeOrderIndices);
            
            std::unordered_set<std::size_t> processedEdges;
            auto wasProcessed = [&processedEdges](std::size_t edge) {return processedEdges.find(edge) != processedEdges.end(); };
            auto handleNewEdge = [](const auto& v0, const auto& v1)
            {
                
            };
            for(const auto& desc: edgeOrderIndices)
            {
                const auto eIndex = desc.index;
                if (wasProcessed(eIndex)) continue;
                const auto edge = desc.edge;

                // Acquire segments of interest that we potentially merge with this edge
                std::vector<std::size_t> targetEdges;
                findSegmentsOfInterest(inputGraph, desc.edge, index, targetEdges);
                if (targetEdges.empty()) continue;

                // Merge with closest segment. Note that this may generate new edges
                mergeWithClosestSegment(inputGraph, desc.edge, targetEdges);
            }
        }

    };
}
#endif