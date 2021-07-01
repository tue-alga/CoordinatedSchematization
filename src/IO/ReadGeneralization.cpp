#include "SchematLib/IO/ReadGeneralization.h"
namespace SchematLib::IO {
    std::string_view ReadGeneralization::strip(const std::string& string)
    {
        if (string.empty()) return string;
        std::size_t start = 0;
        std::size_t end = string.size() - 1;
        for (; start < string.size(); ++start)
        {
            if (!std::isspace(string[start])) break;
        }
        if (start == string.size()) return "";
        // Check end condition
        for (; end >= start; --end)
        {
            if (!std::isspace(string[end])) break;
        }
        return std::string_view(string.data() + start, (end - start) + 1);
    }

    bool ReadGeneralization::startsWith(const std::string& str, char c)
    {
        return !str.empty() && str[0] == c;
    }

    bool ReadGeneralization::startsWith(const std::string& str, const std::string& needle)
    {
        if (str.size() < needle.size()) return false;
        for (std::size_t i = 0; i < needle.size(); ++i)
        {
            if (str[i] != needle[i])return false;
        }
        return true;
    }

    void ReadGeneralization::operator()(const std::filesystem::path& generalizationFile,
                                        Models::UndirectedEmbeddedGraph& graph)
    {
        namespace fs = std::filesystem;
        if (!fs::exists(generalizationFile))
        {
            throw std::runtime_error("Generalization path does not exist");
        }
        // TODO: make this safe for other graphs at some point

        std::ifstream stream(generalizationFile.string());
        std::string line;
        ParsingState state = ParsingState::None;
        std::size_t expectedEdges = 0;
        std::size_t expectedVertices = 0;
        std::size_t currentVertex = 0;
        std::size_t nextEdgeId = 0;
        std::size_t maxVertId = 0;

        while (std::getline(stream, line))
        {
            // Ignore empty lines.
            if (line.empty()) continue;

            auto view = strip(line);
            std::string strippedLine(view);
            std::stringstream lineStream(strippedLine);

            if (startsWith(strippedLine, '#'))
            {
                if (state == ParsingState::None)
                {
                    char c;
                    lineStream >> c;
                    lineStream >> expectedVertices;
                    state = ParsingState::Vertices;
                    graph = Models::UndirectedEmbeddedGraph(expectedVertices);
                    graph.m_use_canonical_edges_ids = false;
                }
                else if (state == ParsingState::Vertices)
                {
                    state = ParsingState::Edges;
                    char c;
                    lineStream >> c;
                    lineStream >> expectedEdges;
                }
                else
                {
                    throw std::runtime_error("Unrecognized '#' line in input");
                }
                continue;
            }
            // Data
            if (state == ParsingState::Vertices)
            {
                Models::NT x, y;
                std::size_t id;
                lineStream >> x >> y >> id;
                boost::put(GCpp::DS::vertex_location_t{}, graph, currentVertex, Models::Point(x, y));

                if (lineStream.fail())
                {
                    boost::put(GCpp::DS::vertex_id_t{}, graph, currentVertex, currentVertex);
                    maxVertId = std::max(maxVertId, currentVertex);
                }
                else
                {
                    boost::put(GCpp::DS::vertex_id_t{}, graph, currentVertex, id);
                    maxVertId = std::max(maxVertId, id);
                }
                ++currentVertex;
            }
            else if (state == ParsingState::Edges)
            {
                std::size_t src, target, eId;
                lineStream >> src >> target >> eId;
                auto result = GCpp::DS::safe_add_edge(src, target, graph);
                if(!result.second)
                {
                    throw std::runtime_error("Did not add edge?");
                }
                assert(result.second);
                boost::put(GCpp::DS::edge_id_t{}, graph, result.first, eId);
                nextEdgeId = std::max(nextEdgeId, eId + 1);
            }
        }
        boost::get_property(graph, GCpp::DS::next_vertex_id_t{}) = maxVertId+1;
        boost::get_property(graph, GCpp::DS::next_edge_id_t{}) = nextEdgeId;
    }

    void WriteGeneralization::operator()(const std::filesystem::path& generalizationFile,
        const Models::UndirectedEmbeddedGraph& graph)
    {
        std::ofstream out(generalizationFile.string());
        if (!out.is_open()) throw std::runtime_error("Could not open output file " + generalizationFile.string());
        this->operator()(out, graph);
    }

    void WriteGeneralization::operator()(std::ostream& output, const Models::UndirectedEmbeddedGraph& graph)
    {
        namespace it = GCpp::Helpers::Iterators;
        output << "# " << boost::num_vertices(graph) << '\n';
        for(auto v : it::range(boost::vertices(graph)))
        {
            const auto loc = GCpp::DS::get_vertex_location(v, graph);
            output << loc.x() << " " << loc.y() << GCpp::DS::getVertexId(graph, v) << '\n';
        }
        for(auto e: it::range(boost::edges(graph)))
        {
            output << boost::source(e, graph) << " " << boost::target(e, graph) << " " << GCpp::DS::getEdgeId(graph, e) << '\n';
        }
    }

    void WriteGeneralization::operator()(const std::filesystem::path& generalizationFile,
        const MapSimplification::FaceMergeTraits::Arr& graph)
    {
        std::ofstream stream(generalizationFile.string());
        this->operator()(stream, graph);
    }

    void WriteGeneralization::operator()(std::ostream& output, const MapSimplification::FaceMergeTraits::Arr& graph)
    {
        namespace it = GCpp::Helpers::Iterators;
        output << "# " << graph.number_of_vertices() << '\n';
        for (auto v : graph.vertex_handles())
        {
            const auto loc = v->point();
            output << CGAL::to_double(loc.x()) << " " << CGAL::to_double(loc.y()) << v->data().id << '\n';
        }
        output << "# " << graph.number_of_edges();
        for (auto e : graph.edge_handles())
        {
            output << e->source()->data().id << " " << e->target()->data().id << " " << e->data().id.value() << '\n';
        }
        output << "# " << graph.number_of_faces();
        for (MapSimplification::FaceMergeTraits::Arr::Face_const_handle f : graph.face_handles())
        {
            if (f->is_unbounded() ) continue;
            output << f->data().id.value() << " ";
            std::vector<std::size_t> edgeIds;
            auto circ = f->outer_ccb();
            MapSimplification::FaceMergeTraits::for_each_ccb_edge(f, [&edgeIds](auto edge, auto& cmd)
            {
                edgeIds.push_back(edge->data().id.value());
            });
            output << edgeIds.size();
            for (auto e : edgeIds)
            {
                output << ' ' << e;
            }
        }
    }
}
