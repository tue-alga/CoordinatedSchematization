#include "SchematLib/MapSimplification/EdgeTrackingObserver.h"
namespace SchematLib::MapSimplification
{
    void EdgeTrackingObserver::lockIncrementTimestamp()
    {
        if (m_timestepLocked) return;
        ++timeStamp;
        m_timestepLocked = true;
    }

    void EdgeTrackingObserver::unlockIncrementTimestamp()
    {
        m_timestepLocked = false;
    }

    bool EdgeTrackingObserver::hasEdge(std::size_t edgeID) const
    {
        return m_edgeMapping.find(edgeID) != m_edgeMapping.end();
    }

    bool EdgeTrackingObserver::hasVertex(std::size_t edgeID) const
    {
        return m_vertexMapping.find(edgeID) != m_vertexMapping.end();
    }

    std::pair<bool, std::size_t> EdgeTrackingObserver::edgePathExists(const EdgePath& path) const
    {
        if (m_edgePathToId.find(path) != m_edgePathToId.end())
        {
            return {true, m_edgePathToId.at(path)};
        }
        /*auto reverted = path.reverted();
        if (m_edgePathToId.find(reverted) != m_edgePathToId.end())
        {
            return { true, m_edgePathToId.at(reverted) };
        }*/
        return {false, 0};
    }

    std::size_t EdgeTrackingObserver::addEdgePath(EdgePath&& path)
    {
        auto localPath = std::move(path);
        localPath.id = m_currentEdgePathId;
        ++m_currentEdgePathId;
        m_edgePaths[localPath.id] = localPath;
        m_edgePathToId[localPath] = localPath.id;
        //m_edgePathToId[localPath.reverted()] = localPath.id;
        return localPath.id;
    }

    std::size_t EdgeTrackingObserver::addEdgePathIfNotExists(EdgePath&& path)
    {
        auto localPath = std::move(path);
        auto [exists, id] = edgePathExists(localPath);
        if (exists)
        {
            return id;
        }
        return addEdgePath(std::move(localPath));
    }

    void EdgeTrackingObserver::verifyMapsSelf(std::size_t edgeId) const
    {
        if (!isSelf(m_edgeMapping.at(edgeId)))
        {
            std::cout << "Does not map to self: " << edgeId << ", maps to index type " << m_edgeMapping.at(edgeId).
                index() << '\n';
            outputMapping(m_edgeMapping.at(edgeId), std::cout);
            std::cout << "\n";
            throw std::runtime_error("Expected edge to map to self");
        }
        if (m_useSpecialEdgeIds)
        {
            auto otherEId = SEF::flipCanonical(edgeId);
            if (!isSelf(m_edgeMapping.at(otherEId)))
            {
                std::cout << "Does not map to self: " << otherEId << ", maps to index type " << m_edgeMapping.
                    at(otherEId).index() << '\n';
                outputMapping(m_edgeMapping.at(otherEId), std::cout);
                std::cout << "\n";
                throw std::runtime_error("Expected edge to map to self");
            }
        }
    }

    void EdgeTrackingObserver::initialize(const std::set<std::size_t>& originalVertexIds,
                                          const std::map<std::size_t, std::pair<std::size_t, std::size_t>>&
                                          originalEdgeIds)
    {
        TimestampScope scope(this);
        // To make it easy for ourselves, trivially map originals to themselves.
        for (const auto& kv : originalEdgeIds)
        {
            if (SEF::isCanonical(kv.first))
            {
                handleNewEdge(kv.first, kv.second.first, kv.second.second);
                maxEId = std::max(maxEId, kv.first);
            }
        }
        for (auto v : originalVertexIds)
        {
            handleNewVertex(v);
            maxVId = std::max(maxVId, v);
        }
    }

    void EdgeTrackingObserver::selfMappingEdges(std::vector<std::size_t>& output) const
    {
        for (const auto& kv : m_edgeMapping)
        {
            if (isSelf(kv.second)) output.push_back(kv.first);
        }
    }

    void EdgeTrackingObserver::verifyEdgeIntegrity(const std::set<std::size_t>& edges) const
    {
        for (auto kv : m_edgeMapping)
        {
            if (m_useSpecialEdgeIds && SEF::isNonCanonical(kv.first))
            {
                if (isSelf(kv.second) && edges.find(SEF::flipCanonical(kv.first)) == edges.end())
                {
                    std::cout << kv.first << " -- " << SEF::flipCanonical(kv.first);
                    throw std::runtime_error("Missing non-canonical counter part");
                }
                continue;
            }
            if (isSelf(kv.second) && edges.find(kv.first) == edges.end())
            {
                std::cout << "Error: " << "Missing edge in graph " << kv.first << '\n';
                throw std::runtime_error("Missing edge in graph");
            }
        }
        // Reverse
        for (auto e : edges)
        {
            if (m_edgeMapping.find(e) == m_edgeMapping.end()) throw
                std::runtime_error("Edge in graph, but not tracker");
            if (!isSelf(m_edgeMapping.at(e))) throw std::runtime_error("Edge in graph maps to somehwere else!");
        }
    }

    void EdgeTrackingObserver::verifyVertexIntegrity(const std::set<std::size_t>& vertices) const
    {
        for (auto kv : m_vertexMapping)
        {
            if (isSelf(kv.second) && vertices.find(kv.first) == vertices.end())
            {
                std::cout << "Error: " << "Missing vertex in graph " << kv.first << '\n';
                //throw std::runtime_error("Missing vertex in graph");
            }
        }
        // Reverse
        for (auto v : vertices)
        {
            if (m_vertexMapping.find(v) == m_vertexMapping.end()) throw std::runtime_error(
                "Vertex in graph, but not tracker");
            if (!isSelf(m_vertexMapping.at(v)))
            {
                auto res = m_vertexMapping.at(v);
                std::cout << "Vertex in graph maps to somehwere else! v=" + std::to_string(v) + ", map type " +
                    std::to_string(res.index()) << " \n";
                outputMapping(m_vertexMapping.at(v), std::cout);
                throw std::runtime_error(
                    "Vertex in graph maps to somehwere else! v=" + std::to_string(v) + ", map type " + std::to_string(
                        res.index()));
            }
        }
    }

    void EdgeTrackingObserver::verifyFaceIntegrity(const std::set<std::size_t>& faces) const
    {
        for (auto kv : m_faceMapping)
        {
            if (isSelf(kv.second) && faces.find(kv.first) == faces.end())
            {
                std::cout << "Error: " << "Missing face in graph " << kv.first << '\n';
                throw std::runtime_error("Missing face in graph");
            }
        }
        // Reverse
        for (auto f : faces)
        {
            if (m_faceMapping.find(f) == m_faceMapping.end()) throw std::runtime_error(
                "Vertex in graph, but not tracker");
            if (!isSelf(m_faceMapping.at(f)))
            {
                auto res = m_faceMapping.at(f);
                std::cout << "Face in graph maps to somehwere else! v=" + std::to_string(f) + ", map type " +
                    std::to_string(res.index()) << " \n";
                outputMapping(res, std::cout);
                throw std::runtime_error(
                    "Face in graph maps to somehwere else! v=" + std::to_string(f) + ", map type " + std::to_string(
                        res.index()));
            }
        }
    }

    void EdgeTrackingObserver::verifyIntegrity(const std::set<std::size_t>& edges,
                                               const std::set<std::size_t>& vertices) const
    {
        verifyEdgeIntegrity(edges);
        verifyVertexIntegrity(vertices);
    }

    void EdgeTrackingObserver::verifyIntegrity(const std::set<std::size_t>& edges,
                                               const std::set<std::size_t>& vertices,
                                               const std::set<std::size_t>& faces) const
    {
        verifyEdgeIntegrity(edges);
        verifyVertexIntegrity(vertices);
        verifyFaceIntegrity(faces);
    }

    void EdgeTrackingObserver::handleNewVertex(std::size_t vertexId)
    {
        //TODO: reinstate
        //if (m_vertexMapping.find(vertexId) != m_vertexMapping.end()) throw std::runtime_error("Duplicate vertex id!");
        if (m_vertexMapping.find(vertexId) != m_vertexMapping.end()) return;
        m_vertexMapping[vertexId] = withCurrentTimestamp(Self{vertexId});
        ++m_eventCount;
        maxVId = std::max(maxVId, vertexId);
    }

    void EdgeTrackingObserver::handleVertexInsert(std::size_t edge, std::size_t newVertex, std::size_t newStartEdge,
                                                  std::size_t newEndEdge)
    {
        auto epId = addEdgePath(EdgePath{{newStartEdge, newEndEdge}, 0});
        TimestampScope scope(this);
        m_edgeMapping[edge] = withCurrentTimestamp(ViaEdgePath{epId, true});
        m_edgeMapping[newStartEdge] = withCurrentTimestamp(Self{newStartEdge});
        m_edgeMapping[newEndEdge] = withCurrentTimestamp(Self{newEndEdge});
        m_vertexMapping[newVertex] = withCurrentTimestamp(Self{newVertex});
        ++m_eventCount;
    }

    void EdgeTrackingObserver::handleVertexDelete(std::size_t vertex)
    {
        m_vertexMapping[vertex] = withCurrentTimestamp(Deleted{});
        ++m_eventCount;
    }

    void EdgeTrackingObserver::handleVertexMerge(std::size_t vertex, std::size_t mergeVertex)
    {
        if (vertex == mergeVertex) throw std::runtime_error("Invalid merge: merging to self!");
        //LogLine() << "Merging vertex  " << vertex << " to vertex  " << mergeVertex;
        m_vertexMapping[vertex] = withCurrentTimestamp(Vertex{mergeVertex});
        ++m_eventCount;
    }

    void EdgeTrackingObserver::handleJunctionReroute(const FaceTreeReroute<long double>& reroute)
    {
        const auto id = junctionId;
        m_junctions[id] = reroute;
        ++junctionId;
        TimestampScope scope(this);
        for (auto e : reroute.affected_edges())
        {
            m_edgeMapping[e] = withCurrentTimestamp(ViaSubgraph{id, e});
            if (m_useSpecialEdgeIds)
            {
                m_edgeMapping[SEF::flipCanonical(e)] = withCurrentTimestamp(ViaSubgraph{id, SEF::flipCanonical(e)});
            }
        }
        for (auto v : reroute.affected_vertices())
        {
            m_vertexMapping[v] = withCurrentTimestamp(ViaSubgraph{id, v});
        }
        for (auto v : reroute.vertices())
        {
            if (reroute.is_terminal_vertex(v)) continue;
            handleNewVertex(v);
        }
        for (const auto& e : reroute.edges())
        {
            if (SEF::isCanonical(e.id))
            {
                handleNewEdge(e.id, e.source, e.target);
            }
        }
        unlockIncrementTimestamp();
    }

    void EdgeTrackingObserver::handleNewEdge(std::size_t edgeId, std::size_t srcVert, std::size_t targetVert)
    {
        if (m_edges.find(edgeId) != m_edges.end()) throw std::runtime_error("Duplicate edge id!");
        auto [minV, maxV] = std::minmax(srcVert, targetVert);
        //LogLine() << "New edge: " << edgeId << " from " << srcVert << " to " << targetVert;
        m_edges[edgeId] = std::make_pair(minV, maxV);
        TimestampScope scope(this);
        m_edgeMapping[edgeId] = withCurrentTimestamp(Self{edgeId});
        ++m_eventCount;
        maxEId = std::max(maxEId, edgeId);
        if (m_useSpecialEdgeIds)
        {
            m_edges[SpecialEdgeIdFunctions::flipCanonical(edgeId)] = std::make_pair(maxV, minV);
            m_edgeMapping[SpecialEdgeIdFunctions::flipCanonical(edgeId)] = withCurrentTimestamp(Self{
                SpecialEdgeIdFunctions::flipCanonical(edgeId)
            });
        }
    }

    void EdgeTrackingObserver::handleEdgeReplace(std::size_t edgeId, const std::vector<std::size_t>& replacementEdges)
    {
        verifyMapsSelf(edgeId);
        handleEdgeReroute(std::vector<std::size_t>{edgeId}, replacementEdges);
    }

    void EdgeTrackingObserver::handleEdgeCollapse(std::size_t edge, std::size_t mergeVertex)
    {
        verifyMapsSelf(edge);
        TimestampScope scope(this);
        //LogLine() << "Edge collapse " << edge << " to vert " << mergeVertex;
        m_edgeMapping[edge] = withCurrentTimestamp(Vertex{mergeVertex});
        ++m_eventCount;
        if (m_useSpecialEdgeIds)
        {
            m_edgeMapping[SEF::flipCanonical(edge)] = withCurrentTimestamp(Vertex{mergeVertex});
        }
    }

    void EdgeTrackingObserver::handleEdgeEdgeMerge(std::size_t edge, std::size_t mergeEdge)
    {
        verifyMapsSelf(edge);
        verifyMapsSelf(mergeEdge);
        //LogLine() << "Edge merge " << edge << " to edge " << mergeEdge;
        if (isVertex(m_edgeMapping[mergeEdge])) throw std::runtime_error("Merging to edge that is merged to vertex!");
        TimestampScope scope(this);
        m_edgeMapping[edge] = withCurrentTimestamp(Edge{mergeEdge});
        ++m_eventCount;
        if (m_useSpecialEdgeIds)
        {
            m_edgeMapping[SEF::flipCanonical(edge)] = withCurrentTimestamp(Edge{SEF::flipCanonical(mergeEdge)});
        }
    }

    void EdgeTrackingObserver::handleEdgeDelete(std::size_t edge)
    {
        verifyMapsSelf(edge);
        TimestampScope scope(this);
        m_edgeMapping[edge] = withCurrentTimestamp(Deleted{});
        ++m_eventCount;
        if (m_useSpecialEdgeIds)
        {
            m_edgeMapping[SEF::flipCanonical(edge)] = withCurrentTimestamp(Deleted{});
        }
    }

    void EdgeTrackingObserver::handleEdgeReroute(const std::vector<std::size_t>& originalRoute,
                                                 const std::vector<std::size_t>& targetRoute)
    {
        if (targetRoute.size() == 0) throw std::runtime_error("Empty replacement path");
        TimestampScope scope(this);
        // Check
        {
            for (auto el : originalRoute) { verifyMapsSelf(el); }

            std::set<std::size_t> incoming;
            incoming.insert(originalRoute.begin(), originalRoute.end());
            for (auto e : targetRoute)
            {
                if (incoming.find(e) != incoming.end())
                {
                    std::cout << "Edge in input and output for reroute: " << e << '\n';
                    std::cout << "Original: ";
                    for (auto el : originalRoute)std::cout << " " << el;
                    std::cout << '\n';
                    std::cout << "Reroute: ";
                    for (auto el : targetRoute)std::cout << " " << el;
                    std::cout << '\n';
                    throw std::runtime_error("Edge in input and output for reroute!");
                }
            }
        }

        if (targetRoute.empty()) throw std::runtime_error("Empty route given to rerouting!");

        if (originalRoute.size() == 1 && targetRoute.size() == 1)
        {
            handleEdgeEdgeMerge(originalRoute.front(), targetRoute.front());
            return;
        }

        if (originalRoute.size() == 1)
        {
            auto targetEp = EdgePath{targetRoute, 0};
            auto [targetExists, targetId] = edgePathExists(targetEp);
            if (!targetExists)
            {
                targetId = addEdgePath(std::move(targetEp));
            }

            verifyMapsSelf(originalRoute[0]);
            m_edgeMapping[originalRoute[0]] = withCurrentTimestamp(ViaEdgePath{targetId, true,originalRoute[0] });
        }
        else if (targetRoute.size() == 1)
        {
            auto ep = EdgePath{originalRoute, 0};
            const auto id = addEdgePathIfNotExists(std::move(ep));
            m_pathMapping[id] = Edge{targetRoute[0]};
            for (auto e : originalRoute)
            {
                verifyMapsSelf(e);
                m_edgeMapping[e] = withCurrentTimestamp(ViaEdgePath{id, false,e});
            }
        }
        else
        {
            // Find if the original route is already an edge path.
            auto ep = EdgePath{originalRoute, 0};
            auto targetEp = EdgePath{targetRoute, 0};
            const auto id = addEdgePathIfNotExists(std::move(ep));
            const auto targetId = addEdgePathIfNotExists(std::move(targetEp));

            m_pathMapping[id] = withCurrentTimestamp(ViaEdgePath{targetId, false,id});
            // Findout if the route is a logical subsegment in all affected edges.
            for (auto e : originalRoute)
            {
                verifyMapsSelf(e);
                m_edgeMapping[e] = withCurrentTimestamp(ViaEdgePath{id, false,e});
            }
        }
        ++m_eventCount;
        if (m_useSpecialEdgeIds)
        {
            handleEdgeRerouteSpecialEdges(originalRoute, targetRoute);
        }
    }

    void EdgeTrackingObserver::handleEdgeRerouteSpecialEdges(const std::vector<std::size_t>& originalRouteIn,
                                                             const std::vector<std::size_t>& targetRouteIn)
    {
        std::vector<std::size_t> originalRoute = originalRouteIn;
        std::vector<std::size_t> targetRoute = targetRouteIn;
        SEF::reverseCanonicalComplement(originalRoute);
        SEF::reverseCanonicalComplement(targetRoute);
        if (targetRoute.size() == 0) throw std::runtime_error("Empty replacement path");
        // Check
        {
            std::set<std::size_t> incoming;
            incoming.insert(originalRoute.begin(), originalRoute.end());
            for (auto e : targetRoute)
            {
                if (incoming.find(e) != incoming.end())
                {
                    std::cout << "Edge in input and output for reroute: " << e << '\n';
                    throw std::runtime_error("Edge in input and output for reroute!");
                }
            }
        }

        if (targetRoute.empty()) throw std::runtime_error("Empty route given to rerouting!");
        if (originalRoute.size() == 1)
        {
            auto targetEp = EdgePath{targetRoute, 0};
            auto [targetExists, targetId] = edgePathExists(targetEp);
            if (!targetExists)
            {
                targetId = addEdgePath(std::move(targetEp));
            }
            m_edgeMapping[originalRoute[0]] = withCurrentTimestamp(ViaEdgePath{targetId, true,originalRoute[0] });
        }
        else if (targetRoute.size() == 1)
        {
            auto ep = EdgePath{originalRoute, 0};
            const auto id = addEdgePathIfNotExists(std::move(ep));
            m_pathMapping[id] = Edge{targetRoute[0]};
            for (auto e : originalRoute)
            {
                m_edgeMapping[e] = withCurrentTimestamp(ViaEdgePath{id, false,e });
            }
        }
        else
        {
            // Find if the original route is already an edge path.
            auto ep = EdgePath{originalRoute, 0};
            auto targetEp = EdgePath{targetRoute, 0};
            const auto id = addEdgePathIfNotExists(std::move(ep));
            const auto targetId = addEdgePathIfNotExists(std::move(targetEp));

            m_pathMapping[id] = withCurrentTimestamp(ViaEdgePath{targetId, false,id});
            // Findout if the route is a logical subsegment in all affected edges.
            for (auto e : originalRoute)
            {
                m_edgeMapping[e] = withCurrentTimestamp(ViaEdgePath{id, false,e});
            }
        }
        ++m_eventCount;
    }
}
