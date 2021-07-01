#ifndef SCHEMATLIB_MAPSIMPLIFICATION_EDGETRACKINGOBSERVER_H
#define SCHEMATLIB_MAPSIMPLIFICATION_EDGETRACKINGOBSERVER_H
#include <vector>
#include <iostream>
#include <map>
#include <set>
#include <variant>
#include <GCpp/DS/BoostEmbeddedGraph.h>
#include "SpecialEdgeIdFunctions.h"
namespace SchematLib::MapSimplification
{
    struct PathMapper;


    template<typename Val>
    struct BaseResult
    {
        std::size_t timeStamp = 0;
        BaseResult(){}
        BaseResult(std::size_t timeStamp):timeStamp(timeStamp){}
        template<typename ResultVariant>
        static bool sameAs(const ResultVariant& opt)
        {
            return std::holds_alternative<Val>(opt);
        }
        template<typename ResultVariant>
        static const Val& get(const ResultVariant& opt)
        {
            return std::get<Val>(opt);
        }
    };

    template<typename NT>
    class FaceTreeReroute
    {
        using SEF = SpecialEdgeIdFunctions;
        // Subgraph is assumed to be a tree with terminal vertices.
        std::set<std::size_t> terminalVertices;
        std::map<std::size_t, std::map<std::size_t, std::size_t>> vertToEdge;//Specified bidirectionally, with canonical ID's. Second map containing vert -> edge.
        std::map<std::size_t, NT> edgeLengths; // Edge length per edge, for rerouting purposes. Single sided, so flip canonical if not present.
        
        // Edges that are affected by the subgraph reroute, in the simplification before this reroute
        std::set<std::size_t> affectedEdges;
        std::set<std::size_t> affectedVertices;
        // Only edges pointing to terminal vertex. By canonicality, flipped canonical edge ID points away.
        std::map<std::size_t, std::size_t> affected_edge_to_endterminal;

        //TODO: ?
        std::map<std::size_t, std::size_t> incomingOutgoingEdges; //Map incoming/outgoing edge IDs to vertex ID in the subgraph.

        bool connects_to_terminal(std::size_t edgeId) const
        {
            return affected_edge_to_endterminal.find(edgeId) != affected_edge_to_endterminal.end() ||
                affected_edge_to_endterminal.find(SEF::flipCanonical(edgeId)) != affected_edge_to_endterminal.end();
        }
        bool target_is_terminal(std::size_t edgeId) const
        {
            return affected_edge_to_endterminal.find(edgeId) != affected_edge_to_endterminal.end();
        }
        bool source_is_terminal(std::size_t edgeId) const
        {
            return affected_edge_to_endterminal.find(SEF::flipCanonical(edgeId)) != affected_edge_to_endterminal.end();
        }
        std::optional<std::size_t> connected_terminal(std::size_t edgeId) const
        {
            if (affected_edge_to_endterminal.find(edgeId) != affected_edge_to_endterminal.end()) return affected_edge_to_endterminal.at(edgeId);
            if(affected_edge_to_endterminal.find(SEF::flipCanonical(edgeId)) != affected_edge_to_endterminal.end())return affected_edge_to_endterminal.at(SEF::flipCanonical(edgeId));
            return {};
        }
        std::optional<std::size_t> source_terminal(std::size_t edgeId) const
        {
            if (affected_edge_to_endterminal.find(SEF::flipCanonical(edgeId)) != affected_edge_to_endterminal.end())return affected_edge_to_endterminal.at(SEF::flipCanonical(edgeId));
            return {};
        }
        std::optional<std::size_t> target_terminal(std::size_t edgeId) const
        {
            if (affected_edge_to_endterminal.find(edgeId) != affected_edge_to_endterminal.end()) return affected_edge_to_endterminal.at(edgeId);
            return {};
        }

        void shortest_path(std::size_t startV, std::size_t endV, std::vector<std::size_t>& outEdges) const
        {
            if (startV == endV)
            {
                return;
            }
            using Vertex = std::size_t;
            using Edge = std::size_t;
            // Run dijkstra for shortest path.
            std::map< Vertex, std::pair<Vertex, Edge>> predecessors;
            std::map< Vertex, NT> distance;
            struct Prio
            {
                std::map< Vertex, NT>* distance = nullptr;

                NT getDist(const Vertex& v) const
                {
                    if (distance->find(v) == distance->end()) return std::numeric_limits<NT>::max();
                    return distance->at(v);
                }

                bool operator()(const Vertex& v0, const Vertex& v1) const
                {
                    const auto d0 = getDist(v0);
                    const auto d1 = getDist(v1);
                    return d0 == d1 ? v0 < v1 : getDist(v0) < getDist(v1);
                }
            };
            Prio prio;
            prio.distance = &distance;
            std::unordered_set<Vertex> processed;
            std::set<Vertex, Prio> prioQueue(prio);
            distance[startV] = 0;
            prioQueue.insert(startV);
            while (!prioQueue.empty())
            {
                auto curr = *prioQueue.begin();
                prioQueue.erase(prioQueue.begin());
                if (processed.find(curr) != processed.end()) continue;
                processed.insert(curr);
                if (curr == endV) break;

                for (const auto& kv: vertToEdge.at(curr))
                {
                    auto target = kv.first;
                    const auto edge = kv.second;
                    if (processed.find(target) != processed.end()) continue;

                    auto edgeLen = getLength(edge);
                    const auto computedDist = prio.getDist(curr) + edgeLen;
                    if (computedDist < prio.getDist(target))
                    {
                        distance[target] = computedDist;
                        predecessors[target] = std::make_pair(curr,edge);
                        prioQueue.insert(target);
                    }
                }
            }
            // Find out the shortest path.
            if (distance.find(endV) == distance.end() || predecessors.find(endV) == predecessors.end())
            {
                throw std::runtime_error("Expected a shortest path!");
            }

            // Reconstruct path via predecessors.
            auto curr = endV;
            while (true)
            {
                auto[predV, e] = predecessors[curr];
                outEdges.push_back(e);
                if (predV == startV) break;
                curr = predV;
            }
            std::reverse(outEdges.begin(), outEdges.end());
        }

        std::optional<std::size_t> edge_to(std::size_t startV, std::size_t endV) const
        {
            for(auto el : vertToEdge.at(startV))
            {
                if (el.first == endV) return el.second;
            }
            return {};
        }
    public:
        struct Edge { std::size_t id; std::size_t source; std::size_t target; };
        FaceTreeReroute(){}
        void set_terminal_vertices(const std::set<std::size_t>& terminalVerticesIn)
        {
            terminalVertices = terminalVerticesIn;
            for (auto el : terminalVertices) add_subgraph_vertex(el);
        }
        FaceTreeReroute(const std::set<std::size_t>& terminalVertices): terminalVertices(terminalVertices)
        {
            for (auto el : terminalVertices) add_subgraph_vertex(el);
        }

        template<typename It, typename ToVertexId>
        FaceTreeReroute(It terminalVerticesStart, It terminalVerticesEnd, ToVertexId&& toVertexId)
        {
            for(auto it = terminalVerticesStart; it != terminalVerticesEnd; ++it)
            {
                terminalVertices.insert(toVertexId(*it));
                add_subgraph_vertex(toVertexId(*it));
            }
        }

        std::set<std::size_t> vertices()const
        {
            std::set<std::size_t> vs;
            for (const auto& kv : vertToEdge) { vs.insert(kv.first); }
            return vs;
        }
        std::vector<Edge> edges() const
        {
            std::vector<Edge> es;
            for (const auto& [src,edges] : vertToEdge)
            {
                for(const auto& [target,eId] : edges)
                {
                    es.push_back(Edge{ eId,src,target });
                }
        }
            return es;
        }

        bool is_terminal_vertex(std::size_t v) const
        {
            return terminalVertices.find(v) != terminalVertices.end();
        }
        /**
         * \brief Affected vertices that are capture by the reroute. Excludes terminals
         * \return 
         */
        const std::set<std::size_t>& affected_vertices() const
        {
            return affectedVertices;
        }
        const std::set<std::size_t>& affected_edges() const
        {
            return affectedEdges;
        }

        NT getLength(std::size_t eId) const
        {
            if (edgeLengths.find(eId) == edgeLengths.end()) return edgeLengths.at(SpecialEdgeIdFunctions::flipCanonical(eId));
            return edgeLengths.at(eId);
        }

        void add_affected_edge(std::size_t edgeId, std::size_t sourceV, std::size_t targetV)
        {
            affectedEdges.insert(edgeId);
            affectedEdges.insert(SEF::flipCanonical(edgeId));
            if(terminalVertices.find(sourceV) != terminalVertices.end())
            {
                affected_edge_to_endterminal[SEF::flipCanonical(edgeId)] = sourceV;
            }
            else
            {
                affectedVertices.insert(sourceV);
            }
            if (terminalVertices.find(targetV) != terminalVertices.end())
            {
                affected_edge_to_endterminal[edgeId] = targetV;
            }
            else
            {
                affectedVertices.insert(targetV);
            }
        }

        void add_subgraph_vertex(std::size_t vId)
        {
            if (vertToEdge.find(vId) != vertToEdge.end()) throw std::runtime_error("DUplicate subgraph vert!");
            vertToEdge[vId] = {};
        }
        void add_subgraph_edge(std::size_t v0, std::size_t v1, std::size_t eId, NT length)
        {
            //if (vertToEdge.find(v0) == vertToEdge.end()) throw std::runtime_error("Add vertex first!");
            vertToEdge.at(v0)[v1] = eId;
            vertToEdge.at(v1)[v0] = SpecialEdgeIdFunctions::flipCanonical(eId);
            edgeLengths[eId] = length;
        }

        bool is_edge_affected(std::size_t edgeId) const
        {
            return affectedEdges.find(edgeId) != affectedEdges.end();
        }
        bool is_vertex_affected(std::size_t vertexId) const
        {
            return affectedEdges.find(vertexId) != affectedEdges.end();
        }

        template<typename HEHandle>
        void add_edge_from_handle(HEHandle he)
        {
            add_affected_edge(he->data().id.value(),he->source()->data().id, he->target()->data().id);
        }

        template<typename EdgeIdIterator>
        void decode(EdgeIdIterator affectedStart, EdgeIdIterator affectedEnd, std::vector<std::size_t>& outEdges) const
        {
            std::optional<EdgeIdIterator> firstEdgeWithTerminal;
            std::size_t startOffset = 0;
            for(auto it =affectedStart; it != affectedEnd; ++it,++startOffset)
            {
                if(source_is_terminal(*it))
                {
                    firstEdgeWithTerminal = it;
                    break;
                }
            }
            EdgeIdIterator start = firstEdgeWithTerminal.has_value() ? firstEdgeWithTerminal.value() : affectedStart;
            std::optional<EdgeIdIterator> lastEdgeWithTerminal;
            std::size_t endOffsetCounter = firstEdgeWithTerminal.has_value() ? startOffset : 0;
            std::size_t endOffset = 0;
            for(auto it = start; it != affectedEnd; ++it,++endOffsetCounter)
            {
                if (target_is_terminal(*it))
                {
                    lastEdgeWithTerminal = it;
                    endOffset = endOffsetCounter;
                }
            }
            if (!firstEdgeWithTerminal.has_value() && !lastEdgeWithTerminal.has_value()) return;
            if(firstEdgeWithTerminal.has_value())
            {
                std::size_t startVert = source_terminal(*firstEdgeWithTerminal.value()).value();
                if (lastEdgeWithTerminal.has_value())
                {
                    std::size_t endVert = target_terminal(*lastEdgeWithTerminal.value()).value();
                    std::cout << "[DecodeSubgraph] Through path, offsets:"<<startOffset << ',' << endOffset<<", totalSize " << endOffsetCounter<< "\n";
                    shortest_path(startVert, endVert, outEdges);
                }
                else
                {
                    std::cout << "[DecodeSubgraph] Dangling end, offsets:" << startOffset << ',' << endOffset << ", totalSize " << endOffsetCounter << "\n";
                    // Add single edge
                    outEdges = { vertToEdge.at(startVert).begin()->second };
                }
            }
            else if(lastEdgeWithTerminal.has_value())
            {
                std::cout << "[DecodeSubgraph] Dangling start, offsets:" << startOffset << ',' << endOffset << ", totalSize " << endOffsetCounter << "\n";
                std::size_t endVert = target_terminal(*lastEdgeWithTerminal.value()).value();
                outEdges = { edge_to(vertToEdge.at(endVert).begin()->first, endVert).value() };
            }
        }
    };

    /**
     * \brief Observer for structural changes of a map simplification.
     * Key assumption: initially, edge IDs and vertex IDs are in the range of the number of edges
     * resp. vertices.
     * Stub implementation that does nothing.
     */
    struct EdgeTrackingObserver
    {
        using SEF = SpecialEdgeIdFunctions;

        friend struct PathMapper;
        // Semantic tagging of mapping. Not necessarily efficient
        struct Vertex : public BaseResult<Vertex>
        {
            std::size_t id;
            Vertex(std::size_t id):id(id){}
            Vertex():id(0){}
        };
        struct Edge : public BaseResult<Edge>
        {
            std::size_t id;
            Edge(std::size_t id) :id(id) {}
            Edge() :id(0) {}
        };
        struct ViaSubgraph : public BaseResult<ViaSubgraph>
        {
            std::size_t id;
            std::size_t srcId;
            ViaSubgraph(std::size_t idIn): id(idIn){}
            ViaSubgraph(std::size_t idIn, std::size_t srcIdIn) : id(idIn),srcId(srcIdIn) {}
            ViaSubgraph():id(0){}
        };

        struct ViaEdgePath : public BaseResult<ViaEdgePath>
        {
            std::size_t edgePathId;
            bool singleToMany;
            std::size_t srcEdge = 0; //Edge being mapped
            ViaEdgePath(std::size_t id, bool singleToMany) :edgePathId(id),singleToMany(singleToMany) {}
            ViaEdgePath(std::size_t id, bool singleToMany,std::size_t srcEdge) :edgePathId(id), singleToMany(singleToMany), srcEdge(srcEdge) {}
            ViaEdgePath() :edgePathId(0),singleToMany(false) {}
        };
        struct EdgePath : public BaseResult<EdgePath>
        {
            std::vector<size_t> ids;
            std::size_t id;
            EdgePath(): id(0){}
            EdgePath(std::vector<size_t>&& idsIn, std::size_t idIn) : ids(std::move(idsIn)), id(idIn){}
            EdgePath(const std::vector<size_t>& idsIn, std::size_t idIn) : ids(idsIn), id(idIn) {}
            bool operator<(const EdgePath& other) const
            {
                const auto isSmallerSize = ids.size() < other.ids.size();
                const auto maxIndex = std::min(ids.size(), other.ids.size());
                for (std::size_t i = 0; i < maxIndex; ++i)
                {
                    if (ids[i] != other.ids[i]) return ids[i] < other.ids[i];
                }
                return isSmallerSize;
            }
            [[nodiscard]] EdgePath reverted() const
            {
                EdgePath copy = *this;
                std::reverse(copy.ids.begin(), copy.ids.end());
                return copy;
            }
        };
        struct Self : public BaseResult<Self>
        {
            std::size_t id;
            Self():id(0){}
            Self(std::size_t id):id(id){}
        };
        struct Deleted : public BaseResult<Deleted> {};
        struct Face : public BaseResult<Face>
        {
            std::size_t id; bool isOuter = false;
            Face(): id(0){}
            Face(std::size_t id, bool isOuter):id(id),isOuter(isOuter){}
        };
        static constexpr const char* prefix = "[EdgeTrackingObserver]";
        struct EdgeToFaceMapping { std::size_t leftF=0; std::size_t rightF = 0; };
        struct LogLine
        {
            LogLine()
            {
                std::cout << EdgeTrackingObserver::prefix;
            }
            template<typename T>
            LogLine& operator<<(const T& t)
            {
                std::cout << t;
                return *this;
            }
            template<typename T>
            LogLine& operator<<(const std::vector<T>& t)
            {
                for (auto el : t)
                {
                    std::cout << ' ';
                    std::cout << el;
                }
                return *this;
            }
            ~LogLine() {
                std::cout << '\n';
            }
        };

        template<typename ResultVariant>
        static bool isVertex(const ResultVariant& opt)
        {
            return std::holds_alternative<Vertex>(opt);
        }
        template<typename ResultVariant>
        static bool isEdge(const ResultVariant& opt)
        {
            return std::holds_alternative<Edge>(opt);
        }
        template<typename ResultVariant>
        static bool isEdgePath(const ResultVariant& opt)
        {
            return std::holds_alternative<EdgePath>(opt);
        }
        template<typename ResultVariant>
        static bool isViaEdgePath(const ResultVariant& opt)
        {
            return std::holds_alternative<ViaEdgePath>(opt);
        }
        template<typename ResultVariant>
        static bool isDeleted(const ResultVariant& opt)
        {
            return std::holds_alternative<Deleted>(opt);
        }
        template<typename ResultVariant>
        static bool isSelf(const ResultVariant& opt)
        {
            return std::holds_alternative<Self>(opt);
        }
        template<typename ResultVariant>
        static bool isMappedViaEdgePath(const ResultVariant& opt)
        {
            return std::holds_alternative<ViaEdgePath>(opt);
        }
        template<typename ResultVariant>
        static bool isFace(const ResultVariant& opt)
        {
            return std::holds_alternative<Face>(opt);
        }

        using MappingOptionsType = std::variant<Vertex, Edge, Deleted, Self, ViaEdgePath,Face, ViaSubgraph>;
        using ReconstructResult = std::variant<Vertex, Edge, Deleted, ViaEdgePath, Face, ViaSubgraph>;
        using EdgePathResult = std::variant<Vertex, Edge, EdgePath, Deleted>;
        using VertexMapping = std::variant<Self, Vertex, Deleted, Face, ViaSubgraph>;
        using PathMapping = std::variant<Vertex, Edge, Deleted, ViaEdgePath, Face>;
        using FaceMapping = std::variant<Face, Edge, Vertex, Deleted>;

        // Maps ID's !
        std::map<std::size_t, MappingOptionsType> m_edgeMapping;
        std::map<std::size_t, VertexMapping> m_vertexMapping; // Map to edge?
        std::map<std::size_t, std::pair<std::size_t, std::size_t>> m_edges; // EdgeID to src-sink vertex id (ordered low to high!)

        // Edge path mapping
        std::map<std::size_t, PathMapping> m_pathMapping;
        // Edge paths
        std::map<std::size_t, EdgePath> m_edgePaths;
        // Edge paths back map to IDs
        std::map<EdgePath, std::size_t> m_edgePathToId;

        std::map<std::size_t, Face> m_faces;
        // Where does a face go
        std::map<std::size_t, std::variant<Self,Face, Edge, Vertex, Deleted>> m_faceMapping;
        // To which face do edges belong
        std::map<std::size_t, EdgeToFaceMapping> m_edgeToFaceMapping;

        std::map<std::size_t, FaceTreeReroute<long double>> m_junctions;

        // Use special edge IDs encoding some canonical direction.
        bool m_useSpecialEdgeIds = false;

        std::size_t maxVId = 0;
        std::size_t maxEId = 0;
        std::size_t maxFId = 0;

        std::size_t junctionId = 0;

        // Stamp events by 'time' to define an order.
        std::size_t timeStamp = 1;

        // Next ID to use for edge paths
        std::size_t m_currentEdgePathId = 0;

        // Number of events seen.
        std::size_t m_eventCount = 0;

        bool m_timestepLocked = false;
        void lockIncrementTimestamp();

        void unlockIncrementTimestamp();

        struct TimestampScope
        {
            EdgeTrackingObserver* m_parent;
            TimestampScope(EdgeTrackingObserver* parent):m_parent(parent)
            {
                m_parent->lockIncrementTimestamp();
            }
            ~TimestampScope()
            {
                m_parent->unlockIncrementTimestamp();
            }
        };

        std::size_t nextTimeStep()
        {
            if (m_timestepLocked) return timeStamp;
            ++timeStamp;
            return timeStamp;
        }

        bool hasEdge(std::size_t edgeID) const;

        bool hasVertex(std::size_t edgeID) const;

        /**
         * \brief Determine if the edge path exists. Tries both forward and reversed variants (since undirected)
         * \param path The edge path
         * \return Pair of bool and ID, where the bool signals if the element exists.
         */
        std::pair<bool, std::size_t> edgePathExists(const EdgePath& path) const;

        std::size_t addEdgePath(EdgePath&& path);

        std::size_t addEdgePathIfNotExists(EdgePath&& path);

        static MappingOptionsType fromStream(std::istream& stream)
        {
            char c;
            std::size_t id;
            stream >> c;
            switch (c)
            {
            case 'v': {
                stream >> id;
                return Vertex{ id };
            }
            case 'e': {
                stream >> id;
                return Edge{ id };
            }
            case 'd':
            {
                return Deleted{};
            }
            case 's':
            {
                stream >> id;
                return Self{ id };
            }
            case 'g':
            {
                stream >> id;
                return ViaSubgraph{ id };
            }
            case 'f':
            {
                bool isOuter = false;
                stream >> id >> isOuter;
                return Face{ id, isOuter };
            }
            case 'p':
            {
                bool singleToMany = false;
                stream >> id >> singleToMany;
                return ViaEdgePath{ id,  singleToMany};
            }
            default:
                throw std::runtime_error(std::string{ "Unrecognized mapping option type " } +c);
            }
        }

        template<typename T>
        struct Holder
        {
            using type = T;
        };
        template<typename Variant, typename Input, typename T>
        static void assignIfType(Variant& variant, const Input& input, Holder<T>)
        {
            if(std::holds_alternative<T>(input))
            {
                variant = std::get<T>(input);
            }
        }

        template<typename Input, typename...VarTypes>
        static void getFromInput(const Input& input, std::variant<VarTypes...>& out)
        {
            (assignIfType(out, input, Holder<VarTypes>{}), ...);
        }

        struct VariantOutputter
        {
            std::ostream* m_output = nullptr;
            std::ostream& stream() { return *m_output; }
            VariantOutputter(std::ostream& output):m_output(&output){}

            void operator()(const Self& self)
            {
                stream() << "s " << self.id;
            }
            void operator()(const Vertex& v)
            {
                stream() << "v " << v.id;
            }
            void operator()(const Deleted& v)
            {
                stream() << "d";
            }
            void operator()(const Edge& e)
            {
                stream() << "e " << e.id;
            }
            void operator()(const ViaEdgePath& v)
            {
                stream() << "p " << v.edgePathId;
            }
            void operator()(const Face& f)
            {
                stream() << "f " << f.id << ' ' << f.isOuter;
            }
            void operator()(const ViaSubgraph& f)
            {
                stream() << "g " << f.id ;
            }
        };

        template<typename Opt>
        static void outputMapping(const Opt& opt, std::ostream& output)
        {
            VariantOutputter outputter(output);
            std::visit(outputter, opt);
        }

        template<typename El>
        void setCurrentTimestamp(El& el)
        {
            el.timeStamp = nextTimeStep();
        }
        template<typename El>
        El withCurrentTimestamp(const El& el)
        {
            El ret = el;
            ret.timeStamp = nextTimeStep();
            return ret;
        }

        void save(std::ostream& output) const
        {
            output << m_edgeMapping.size() << '\n';
            for(const auto& kv: m_edgeMapping)
            {
                output << kv.first << ' ';
                outputMapping(kv.second, output);
            }
            output << m_vertexMapping.size() << '\n';
            for (const auto& kv : m_vertexMapping)
            {
                output << kv.first << ' ';
                outputMapping(kv.second, output);
                output << '\n';
            }
            output << m_edges.size() << '\n';
            for (const auto& kv : m_edges)
            {
                output << kv.first << ' ';
                output << kv.second.first << ' ' << kv.second.second;
                output << '\n';
            }
            output << m_pathMapping.size() << '\n';
            for (const auto& kv : m_pathMapping)
            {
                output << kv.first << ' ';
                outputMapping(kv.second, output);
                output << '\n';
            }
            output << m_edgePaths.size() << '\n';
            for(const auto& kv: m_edgePaths)
            {
                output << kv.first << ' ' << kv.second.ids.size();
                for(auto el : kv.second.ids)
                {
                    output << ' ' << el;
                }
                output << '\n';
            }
            output << m_currentEdgePathId;
        }
        void load(std::istream& input)
        {
            std::size_t size;
            input >> size;
            for(std::size_t i = 0; i < size; ++i)
            {
                std::size_t edgeId;
                input >> edgeId;
                MappingOptionsType mapping = fromStream(input);
                m_edgeMapping[edgeId] = mapping;
            }
            input >> size;
            for (std::size_t i = 0; i < size; ++i)
            {
                std::size_t vId;
                input >> vId;
                VertexMapping mapping;
                getFromInput(fromStream(input), mapping);
                m_vertexMapping[vId] = mapping;
            }
            input >> size;
            for (std::size_t i = 0; i < size; ++i)
            {
                std::size_t eId,v0,v1;
                input >> eId >> v0 >> v1;
                m_edges[eId] = std::make_pair(v0, v1);
            }
            input >> size;
            for (std::size_t i = 0; i < size; ++i)
            {
                std::size_t pId;
                input >> pId;
                PathMapping mapping;
                getFromInput(fromStream(input), mapping);
                m_pathMapping[pId] = mapping;
            }
            input >> size;
            for (std::size_t i = 0; i < size; ++i)
            {
                std::size_t pId;
                std::size_t pathSize;
                input >> pId >> pathSize;
                std::vector<std::size_t> path;
                for(std::size_t j = 0;j < pathSize; ++j)
                {
                    std::size_t eId;
                    input >> eId;
                    path.push_back(eId);
                }
                m_edgePaths[pId] = EdgePath{ std::move(path), pId };
            }
            input >> m_currentEdgePathId;
        }

        EdgeTrackingObserver()
        {

        }

        void verifyMapsSelf(std::size_t edgeId) const;

        void initialize(const std::set<std::size_t>& originalVertexIds,
                        const std::map<std::size_t, std::pair<std::size_t, std::size_t>>& originalEdgeIds);

        /**
         * \brief Returns all edges mapping to themselves. Should be part of the current simplification.
         * \param output The edge ids of edges mapping to themselves.
         */
        void selfMappingEdges(std::vector<std::size_t>& output) const;
#pragma region Verification
        template<typename EmbeddedGraph>
        void verifyIntegrity(const EmbeddedGraph& graph) const
        {
            std::set<std::size_t> eIds, vIds;
            GCpp::DS::computeEdgeIdSet(graph, eIds);
            GCpp::DS::computeVertexIdSet(graph, vIds);
            verifyIntegrity(eIds, vIds);
        }

        void verifyEdgeIntegrity(const std::set<std::size_t>& edges) const;

        void verifyVertexIntegrity(const std::set<std::size_t>& vertices) const;

        void verifyFaceIntegrity(const std::set<std::size_t>& faces) const;

        void verifyIntegrity(const std::set<std::size_t>& edges, const std::set<std::size_t>& vertices) const;

        void verifyIntegrity(const std::set<std::size_t>& edges, const std::set<std::size_t>& vertices,
                             const std::set<std::size_t>& faces) const;
#pragma endregion
#pragma region Vertex event handlers
        void handleNewVertex(std::size_t vertexId);

        /**
         * \brief Handles the case where a vertex is inserted in an edge, splitting the edge in two sub edges
         * \param edge The edge to be split
         * \param newVertex
         * \param newStartEdge
         * \param newEndEdge
         */
        void handleVertexInsert(std::size_t edge, std::size_t newVertex, std::size_t newStartEdge,
                                std::size_t newEndEdge);

        /**
         * \brief Handles the case when an edge is deleted and has no representation anymore
         * \param edge
         */
        void handleVertexDelete(std::size_t vertex);

        /**
        * \brief Handles the case where a vertex is merged to another. Note that edges that will potentially be merged to the vertex emit
        * their own events.
        * \param vertex The original vertex
        * \param mergeVertex The merge vertex
        */
        void handleVertexMerge(std::size_t vertex, std::size_t mergeVertex);
#pragma endregion
#pragma region Face event handlers
        void handleNewFace(std::size_t faceId)
        {
            // TODO: when we do faces again, reinstate.
            //if (m_faces.find(faceId) != m_faces.end()) throw std::runtime_error("Duplicate face ID");
            TimestampScope scope(this);
            m_faces[faceId] = Face{ faceId ,false};
            m_faceMapping[faceId] = withCurrentTimestamp(Self{ faceId });
            maxFId = std::max(maxFId, faceId);
        }
        void handleEdgeToFaceConnection(std::size_t edgeId, std::size_t leftF, std::size_t rightF)
        {

        }
        void handleFaceMerge(std::size_t face, std::size_t targetFace)
        {
            m_faceMapping[face] = withCurrentTimestamp(Face{ targetFace,false });
        }
        void handleFaceCollapseToVertex(std::size_t face, std::size_t vertex)
        {
            m_faceMapping[face] = withCurrentTimestamp(Vertex{ vertex});
        }

        void handleEdgeToFaceMerge(std::size_t edge, std::size_t face)
        {
            TimestampScope scope(this);
            m_edgeMapping[edge] = withCurrentTimestamp(Face{ face,false });
            if(m_useSpecialEdgeIds)
            {
                m_edgeMapping[SpecialEdgeIdFunctions::flipCanonical(edge)] = withCurrentTimestamp(Face{ face,false });
            }
        }
        void handleVertexToFaceMerge(std::size_t vertex, std::size_t face)
        {
            m_vertexMapping[vertex] = withCurrentTimestamp(Face{ face,false });
        }

#pragma endregion
#pragma region Edge event handlers

        void handleInitialEdgeCount(std::size_t edgeCount) {}

        void handleInitialVertexCount(std::size_t vertexCount) {}

        void handleJunctionReroute(const FaceTreeReroute<long double>& reroute);

        void handleNewEdge(std::size_t edgeId, std::size_t srcVert, std::size_t targetVert);

        void handleEdgeReplace(std::size_t edgeId, const std::vector<std::size_t>& replacementEdges);

        /**
         * \brief Handles an edge collapse, where the given edge (ID) is completely merge to a vertex
         * \param edge
         * \param mergeVertex The target merge vertex
         */
        void handleEdgeCollapse(std::size_t edge, std::size_t mergeVertex);

        /**
         * \brief Handles an edge merge, where the edge is represented by the new edge.
         * \param edge
         * \param mergeEdge
         */
        void handleEdgeEdgeMerge(std::size_t edge, std::size_t mergeEdge);


        /**
         * \brief Handles the case when an edge is deleted and has no representation anymore
         * \param edge
         */
        void handleEdgeDelete(std::size_t edge);


        /**
         * \brief Handles the event when we reindex an edge.
         * \param edgeId The original edge ID
         * \param newEdgeId The new edge ID.
         */
        void handleEdgeReindex(std::size_t edgeId, std::size_t newEdgeId)
        {
            throw std::runtime_error("Don't reindex edges!");
            ++m_eventCount;
        }

        void handleEdgeReroute(const std::vector<std::size_t>& originalRoute,
                               const std::vector<std::size_t>& targetRoute);

        void handleEdgeRerouteSpecialEdges(const std::vector<std::size_t>& originalRouteIn,
                                           const std::vector<std::size_t>& targetRouteIn);

        void handleEdgeRerouteVertexRoute(const std::vector<std::size_t>& originalRout, const std::vector<std::size_t>& newRoute)
        {
            throw std::runtime_error("Don't reroute via vertices!");
            ++m_eventCount;
        }
#pragma endregion
    };
}
#endif