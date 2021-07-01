#ifndef SCHEMATLIB_MAPSIMPLIFICATION_PATHMAPPER_H
#define SCHEMATLIB_MAPSIMPLIFICATION_PATHMAPPER_H
#include "EdgeTrackingObserver.h"

namespace SchematLib::MapSimplification{
struct PathMapper
{
    const EdgeTrackingObserver& m_data;
    using Data = EdgeTrackingObserver;

    using PartialEdgeDecodeResult = std::variant<Data::Vertex, Data::Edge, Data::Deleted, Data::ViaEdgePath, Data::Face>;
    using TerminalElement = std::variant<Data::Vertex, Data::Edge, Data::Deleted, Data::Face>;
    using DecodeVertexResult = std::variant<Data::Vertex, Data::Deleted, Data::Face>;

    bool m_verbose = true;

    PathMapper(const EdgeTrackingObserver& data) :m_data(data) {}

    void setVerbose(bool verbose)
    {
        m_verbose = verbose;
    }
    
    DecodeVertexResult decodeVertex(std::size_t vertex) const
    {
        std::size_t current = vertex;
        while (true)
        {
            auto result = m_data.m_vertexMapping.at(current);
            if (Data::Self::sameAs(result)) // Self
            {
                return { Data::Vertex{ vertex } };
            }
            if (Data::Deleted::sameAs(result)) //deleted
            {
                return { Data::Deleted{} };
            }
            if(Data::Face::sameAs(result))
            {
                return { Data::Face::get(result) };
            }
            if(Data::Vertex::sameAs(result))
            {
                current = std::get<Data::Vertex>(result).id;
            }
            else
            {
                throw std::runtime_error("Unknown vertex decode type");
            }
        }
    }
    PartialEdgeDecodeResult partialDecodeEdge(std::size_t edge) const
    {
        std::size_t current = edge;
        while (true)
        {
            auto result = m_data.m_edgeMapping.at(current);
            if (Data::Edge::sameAs(result))
            {
                current = std::get<Data::Edge>(result).id;
            }
            else if (Data::ViaEdgePath::sameAs(result))
            {
                return { std::get<Data::ViaEdgePath>(result) };
            }
            else if (Data::isDeleted(result))
            {
                return std::get<Data::Deleted>(result);
            }
            else if (Data::isVertex(result))
            {
                auto vertResult = decodeVertex(std::get<Data::Vertex>(result).id);
                if (Data::isDeleted(vertResult)) return Data::Deleted{};
                if (Data::Face::sameAs(vertResult)) return Data::Face::get(vertResult);
                return std::get<Data::Vertex>(vertResult);
            }
            else if (Data::Self::sameAs(result))
            {
                return Data::Edge{ Data::Self::get(result).id };
            }
            else if(Data::Face::sameAs(result))
            {
                return { Data::Face::get(result) };
            }
            else
            {
                throw std::runtime_error("Unknown edge decode type");
            }
        }
    }

#define IFVERBOSE(...) if(m_verbose){__VA_ARGS__}
    /**
     * \brief Maps an input edge ID path to a path in the simplified map.
     * Verify the paths if you use them!
     * \param edgeIdPath Path in the original graph, given as a sequence of edge ID's in the original graph.
     * \param mappedEdgeIdPath The path in the simplified graph, given as ID's of edges or vertices, wrapped in the appropriate type.
     */
    template<typename PathOutputIterator, typename = std::is_assignable<decltype(*std::declval<PathOutputIterator>()), TerminalElement>>
    void mapPath(const std::vector<std::size_t>& edgeIdPath, PathOutputIterator iterator) const
    {
        auto it = edgeIdPath.begin();
        auto end = edgeIdPath.end();
        std::size_t index = 0;
        for (; it != end; ++it,++index)
        {
            IFVERBOSE(std::cout << "[PM] Edge " << *it << '\n';)
            auto decodeResult = partialDecodeEdge(*it);

            if (Data::isViaEdgePath(decodeResult))
            {
                IFVERBOSE(std::cout << "========= Via edge path ======= \n";)
                const auto startPathId = std::get<Data::ViaEdgePath>(decodeResult).edgePathId;
                IFVERBOSE(std::cout << "Reconstruction via edgepath: " << startPathId << '\n';)
                IFVERBOSE(std::cout << "Current edge: " << *it << '\n';)

                // Find the edges of the input that match to the same edge path.
                auto innerIt = it + 1;
                for (; innerIt != end; ++innerIt)
                {
                    auto innerDecodeResult = partialDecodeEdge(*innerIt);
                    if (!Data::isViaEdgePath(innerDecodeResult)) break;
                    // Only capture edges going via the same subroute.
                    const auto innerPathId = std::get<Data::ViaEdgePath>(innerDecodeResult).edgePathId;
                    if (innerPathId != startPathId) break;
                }

                // This is the exclusion range.
                IFVERBOSE(std::cout << "Affected range size: " << std::distance(it, innerIt) << '\n';)
                // Range [it, innerIt] represents a subroute. Map it recursively.
                std::vector<std::size_t> ids;
                if (m_data.m_pathMapping.find(startPathId) == m_data.m_pathMapping.end())
                {
                    ids = m_data.m_edgePaths.at(startPathId).ids;
                }
                else
                {
                    bool continueIteration = false;
                    auto currentId = startPathId;
                    while (true)
                    {
                        if (m_data.m_pathMapping.find(currentId) == m_data.m_pathMapping.end())
                        {
                            ids = m_data.m_edgePaths.at(currentId).ids;
                            break;
                        }
                        auto mapping = m_data.m_pathMapping.at(currentId);
                        if (Data::isVertex(mapping))
                        {
                            it = innerIt;
                            *iterator = std::get<Data::Vertex>(mapping);
                            continueIteration = true;
                            break;
                        }
                        if (Data::isEdge(mapping))
                        {
                            it = innerIt;
                            *iterator = std::get<Data::Edge>(mapping);
                            continueIteration = true;
                            break;
                        }
                        if (Data::isDeleted(mapping))
                        {
                            it = innerIt;
                            *iterator = std::get<Data::Deleted>(mapping);
                            continueIteration = true;
                            break;
                        }
                        if (Data::isViaEdgePath(mapping))
                        {
                            currentId = std::get<Data::ViaEdgePath>(mapping).edgePathId;
                        }
                    }
                    // Continue outer iteration
                    if (continueIteration)
                    {
                        if (it == end) break;

                        IFVERBOSE(std::cout << "========= Mapped to terminal element, continuing ======= \n";)
                        continue;
                    }
                }
                IFVERBOSE(std::cout << "Ep ids: ";
                for (auto el : ids)
                {
                    std::cout << el << " ";
                }
                std::cout << "\n";
                std::cout << "New ids: ";
                for (auto it2 = innerIt; it2 != edgeIdPath.end(); ++it2)
                {
                    std::cout << *it2 << " ";
                }
                std::cout << "\n";)
                ids.insert(ids.cend(), innerIt, edgeIdPath.end());
                IFVERBOSE(std::cout << "New ids size: " << ids.size() << '\n';)
                // Recurse
                IFVERBOSE(std::cout << "========= Recurse ======= \n";)
                mapPath(ids, iterator);
                return;
            }
            if (Data::isVertex(decodeResult))
            {
                *iterator = TerminalElement{ std::get<Data::Vertex>(decodeResult) };
            }
            else if (Data::isEdge(decodeResult))
            {
                *iterator = TerminalElement{ std::get<Data::Edge>(decodeResult) };
            }
            else if (Data::isDeleted(decodeResult))
            {
                *iterator = TerminalElement{ std::get<Data::Deleted>(decodeResult) };
            }
        }
    }
};
#undef IFVERBOSE
}
#endif