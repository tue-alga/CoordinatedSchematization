#include "SchematLib/Models/EmbeddedGraph.h"

namespace SchematLib::Models {

    const std::map<std::string, std::size_t> OsmRanks::OSM_RANKS =
        {
            // Assign links same priority.
            {"motorway",0},
            {"motorway_link",0},
            {"trunk",1},
            {"trunk_link",1},
            {"primary",2},
            {"primary_link",2},
            {"secondary",3},
            {"secondary_link",3},
            {"tertiary",4},
            {"tertiary_link",4},
            {"residential",5}, //Flipped with unclassified
            {"unclassified",6},
            {"service",7},
            {"living_street",8}
        };

    ConvertedPointMaker& ConvertedPointMaker::operator=(const ConvertedPointMaker& other)
    {
        m_base = other.m_base;
        m_target = other.m_target;
        srcEpsg = other.srcEpsg;
        targetEpsg = other.targetEpsg;
        m_isIdentity = targetEpsg == srcEpsg;
        if (other.m_transform.get()) constructTransform();
        return *this;
    }

    void ConvertedPointMaker::setOriginalGdalOrders(bool value)
    {
        m_originalOrders = value;
    }

    void ConvertedPointMaker::setTargetFromEpsg(int epsg)
    {
        targetEpsg = epsg;
        m_target.importFromEPSG(epsg);
    }

    void ConvertedPointMaker::setBaseFromEpsg(int epsg)
    {
        srcEpsg = epsg;
        m_base.importFromEPSG(epsg);
    }

    void ConvertedPointMaker::setBaseToWGS84()
    {
        srcEpsg = 4326;
        m_base.importFromEPSG(4326);
    }

    bool ConvertedPointMaker::hasTransform() const
    {
        return m_transform.get() != nullptr;
    }

    void ConvertedPointMaker::constructTransform()
    {
        m_isIdentity = srcEpsg == targetEpsg;
        if(m_isIdentity)
        {
            return;
        }
        if (m_originalOrders)
        {
            m_base.SetAxisMappingStrategy(OSRAxisMappingStrategy::OAMS_TRADITIONAL_GIS_ORDER);
            m_target.SetAxisMappingStrategy(OSRAxisMappingStrategy::OAMS_TRADITIONAL_GIS_ORDER);
        }
        m_transform = TransformPtr(OGRCreateCoordinateTransformation(&m_base, &m_target), [](auto* crt)
        {
            OGRCoordinateTransformation::DestroyCT(crt);
        });
        if (!m_transform) throw std::runtime_error("Transform construction failed");
    }

    void ConvertedPointMaker::constructTransform(int epsgSrc, int epsgTarget)
    {
        targetEpsg = epsgTarget;
        srcEpsg = epsgSrc;
        if (epsgSrc == epsgTarget) { m_isIdentity = true; return; }
        m_base.importFromEPSG(epsgSrc);
        m_target.importFromEPSG(epsgTarget);
        constructTransform();
    }

    void ConvertedPointMaker::dumpData()
    {
        m_base.dumpReadable();
        m_target.dumpReadable();
    }

    Point ConvertedPointMaker::operator()(std::initializer_list<NT> coordinates) const
    {
        if (m_isIdentity) { return Point(*coordinates.begin(), *std::next(coordinates.begin())); }
        if (!m_transform) throw std::runtime_error("Coordinate transform not initialized");

        double x = *coordinates.begin();
        double y = *std::next(coordinates.begin());
        m_transform->Transform(1, &x, &y);
        return Point(x, y); //TODO concept Point
    }

    Point PointMaker::operator()(std::initializer_list<NT> coordinates) const
    {
        return Point(*coordinates.begin(), *std::next(coordinates.begin()));
    }

    bool OsmRanks::hasRank(const std::string& name)
    {
        return OSM_RANKS.find(name) != OSM_RANKS.end();
    }

    std::size_t OsmRanks::rank(const std::string& name)
    {
        return OSM_RANKS.at(name);
    }

    std::size_t OsmRanks::rank_or_default(const std::string& name, std::size_t defaultVal)
    {
        if (hasRank(name)) return rank(name);
        return defaultVal;
    }
}