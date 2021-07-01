#ifndef SCHEMATLIBUI_EMBEDDEDGRAPH_H
#define SCHEMATLIBUI_EMBEDDEDGRAPH_H
#include <iostream>
#include "BaseTypes.h"
#include <GCpp/DS/BoostEmbeddedGraph.h>
#include <GCpp/Geometry/Index/GenericPointIndex.h>
#include <GCpp/Geometry/Index/GenericEdgeIndex.h>
#include <ogr_spatialref.h>
#include "BidirectionalCanonicalGraph.h"
namespace SchematLib::Models
{

    using EmbeddedGraph = GCpp::DS::BoostEmbeddedGraph<Point>;
    // Use this with special edge ID scheme to specify canonical direction which can be tracked.
    using BidirectionalEmbeddedGraph = GCpp::DS::BoostEmbeddedGraph<Point>;

    // Was this previously
    //using UndirectedEmbeddedGraph = GCpp::DS::BoostUndirectedEmbeddedGraph<Point>;
    using UndirectedEmbeddedGraph = BidirectionalCanonicalGraph;

    template<typename Graph>
    inline void testGraphEquality(const Graph& g1, const Graph& g2)
    {
        auto testStartEnd = boost::vertices(g1);
        auto graphStartEnd = boost::vertices(g2);
#define REQUIRE_EQ(val0,val1,msg) if((val0) != (val1)) {std::cout << "Failed because " << val0 << "!=" << val1; throw std::runtime_error(msg);}
        for (auto it = testStartEnd.first, it2 = graphStartEnd.first; it != testStartEnd.second; ++it, ++it2)
        {
            REQUIRE_EQ(GCpp::DS::getVertexId(g1, *it), GCpp::DS::getVertexId(g2, *it2), "Incongruent vertex IDs");
        }
        for (auto it = testStartEnd.first, it2 = graphStartEnd.first; it != testStartEnd.second; ++it, ++it2)
        {
            {
                std::set<Models::EmbeddedGraph::vertex_descriptor> verts;
                for (auto e : GCpp::Helpers::Iterators::range(boost::out_edges(*it2, g2)))
                {
                    const auto target = boost::target(e, g2);
                    if (verts.find(target) != verts.end())
                    {
                        std::cout << "Duplicate edge in base? Currv = " << GCpp::DS::getVertexId(g2, *it2) << '\n';
                        std::cout << "Duplicate:" << GCpp::DS::getVertexId(g2, target) << '\n';
                        std::cout << "Connected to ";
                        for (auto el : verts) std::cout << ' ' << GCpp::DS::getVertexId(g2, el);
                        throw std::runtime_error("Dup edge in base");
                    }
                    verts.insert(target);
                }
            }
            {
                std::set<Models::EmbeddedGraph::vertex_descriptor> verts;
                for (auto e : GCpp::Helpers::Iterators::range(boost::out_edges(*it, g1)))
                {
                    const auto target = boost::target(e, g1);
                    if (verts.find(target) != verts.end())
                    {
                        std::cout << "Duplicate edge in test? Currv = " << *it << '\n';
                        std::cout << "Connected to ";
                        for (auto el : verts) std::cout << ' ' << el;
                        throw std::runtime_error("Dup edge in test");
                    }
                    verts.insert(target);
                }
            }

            for (auto e : GCpp::Helpers::Iterators::range(boost::out_edges(*it, g1)))
            {
                auto[otherEdge, exists] = boost::lookup_edge(boost::source(e, g1), boost::target(e, g1), g2);
                REQUIRE_EQ(exists, true, "Edge not present in graph but is present in test");
                REQUIRE_EQ(GCpp::DS::getEdgeId(g1, e), GCpp::DS::getEdgeId(g2, otherEdge), "Incongruent edge IDs");
                REQUIRE_EQ(GCpp::DS::getVertexId(g1, *it), GCpp::DS::getVertexId(g2, boost::source(otherEdge, g2)), "Incongruent edge IDs");
            }
            for (auto e : GCpp::Helpers::Iterators::range(boost::out_edges(*it, g2)))
            {
                auto[otherEdge, exists] = boost::lookup_edge(boost::source(e, g2), boost::target(e, g2), g1);
                REQUIRE_EQ(exists, true, "Edge not present in graph but is present in test");
                REQUIRE_EQ(GCpp::DS::getEdgeId(g2, e), GCpp::DS::getEdgeId(g1, otherEdge), "Incongruent edge IDs");
                REQUIRE_EQ(GCpp::DS::getVertexId(g2, *it), GCpp::DS::getVertexId(g1, boost::source(otherEdge, g1)), "Incongruent edge IDs");
            }
        }
    }
#undef REQUIRE_EQ
    class ConvertedPointMaker
    {
    public:
        using TransformPtr = std::unique_ptr<OGRCoordinateTransformation, std::function<void(OGRCoordinateTransformation*)>>;
    private:
        OGRSpatialReference m_base;
        OGRSpatialReference m_target;
        TransformPtr m_transform;
        bool m_originalOrders = true;
        bool m_isIdentity = false;
        int srcEpsg=4326, targetEpsg = 4326;
    public:

        ConvertedPointMaker& operator=(const ConvertedPointMaker& other);

        void setOriginalGdalOrders(bool value);

        void setTargetFromEpsg(int epsg);

        void setBaseFromEpsg(int epsg);

        void setBaseToWGS84();

        bool hasTransform() const;

        void constructTransform();

        void constructTransform(int epsgSrc, int epsgTarget);

        void dumpData();

        Point operator()(std::initializer_list<NT> coordinates) const;
    };
    struct PointMaker
    {
        Point operator()(std::initializer_list<NT> coordinates) const;
    };

    struct SqDistance
    {
        NT operator()(const Point&  p0, const Point& p1) const
        {
            return (p1 - p0).sqLength();
        }

        NT operator()(const Point&  p0, const std::pair<Point, Point>& p1) const
        {
            auto diff = p0 - p1.first;
            auto dir = p1.second - p1.first;
            auto project = diff.dot(dir) / dir.sqLength();
            if (project < 0) return diff.sqLength();
            if (project > 1.0) return (p0 - p1.second).sqLength();
            return (p0 - (p1.first + project * dir)).sqLength();
        }

        NT operator()(const std::pair<Point, Point>& p1, const Point&  p0) const
        {
            return this->operator()(p0, p1);
        }
    };

    using GraphNodeIndex = GCpp::Geometry::Index::NetworkNodeIndex<EmbeddedGraph, NT, Point, PointMaker, SqDistance>;
    using UndirectedGraphNodeIndex = GCpp::Geometry::Index::NetworkNodeIndex<UndirectedEmbeddedGraph, NT, Point, PointMaker, SqDistance>;
    using UndirectedGraphEdgeIndex = GCpp::Geometry::Index::NetworkEdgeIndex<UndirectedEmbeddedGraph, NT, PointMaker, SqDistance>;

    struct OsmRanks
    {
        static const std::map<std::string, std::size_t> OSM_RANKS;

        static bool hasRank(const std::string& name);

        static std::size_t rank(const std::string& name);

        static std::size_t rank_or_default(const std::string& name, std::size_t defaultVal);
    };
}
//namespace GCpp::DS
//{
//    template<>
//    inline decltype(std::declval<SchematLib::Models::Point >().x()) compute_length<SchematLib::Models::Point>(const SchematLib::Models::Point& p0, const SchematLib::Models::Point& p1)
//    {
//        return (p0 - p1).length();
//    }
//}
#endif