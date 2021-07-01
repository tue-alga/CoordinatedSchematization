#include "SchematLib/IO/ReadWriteOsmMap.h"


#include <iostream>
#include <GCpp/IO/BoostOsmReader.h>
#include <GCpp/String.h>

namespace iter = GCpp::Helpers::Iterators;
namespace fs = std::filesystem;

namespace SchematLib::IO
{
    struct ForwardRoadTag
    {
        std::string operator()(const std::string& tag) const
        {
            return tag;
        }
    };


    void ReadOsmMap::setMapEpsg(int epsg)
    {
        m_epsg = epsg;
    }

    int ReadOsmMap::epsg() const
    {
        return m_epsg;
    }

    void ReadOsmMap::operator()(const std::filesystem::path& generalizationFile, Models::EmbeddedGraph& graph,
                                std::unordered_map<std::size_t, std::string>& highwayTags)
    {
        using Graph = Models::EmbeddedGraph;
        auto mapFilePath = generalizationFile;
        if(mapFilePath.extension() != OSM_EXTENSION)
        {
            throw std::runtime_error("Not reading .osm file! Instead:" + mapFilePath.extension().string());
        }

        std::ifstream stream(mapFilePath);
        if (!stream.is_open()) throw std::runtime_error("Could not open map file stream for " + mapFilePath.string());
        // When including specific road types, we need to read priorities via the highway tag.
        using NodeCreator = GCpp::IO::VertexCreator< Models::EmbeddedGraph>;

        using EdgeCreator = GCpp::IO::ExternalStoreEdgeCreator<Graph, std::unordered_map, ForwardRoadTag>;
        using Reader = GCpp::IO::BoostOsmReaderWithTags<Graph, Models::ConvertedPointMaker, NodeCreator, EdgeCreator>;
        Reader reader;
        // Setup the transform
        reader.creator().setBaseToWGS84();
        reader.creator().setTargetFromEpsg(m_epsg);
        reader.creator().setOriginalGdalOrders(true);
        reader.creator().constructTransform();

        NodeCreator nodeCreator;
        ForwardRoadTag priorityComputer;
        EdgeCreator edgeCreator("highway", highwayTags, priorityComputer);
        reader.parse(stream, graph, nodeCreator, edgeCreator);

        namespace it = GCpp::Helpers::Iterators;
        // Print out tally
        std::map<std::string, std::size_t> tally;
        for (auto e : it::range(boost::edges(graph)))
        {
            const auto eId = boost::get(GCpp::DS::edge_id_t{}, graph, e);
            if (highwayTags.find(eId) != highwayTags.end())
            {
                const auto val = highwayTags.at(eId);
                if (tally.find(val) != tally.end()) tally[val] += 1;
                else tally[val] = 1;
            }
        }
        for (auto pair : tally)
        {
            std::cout << "Roads of type " << pair.first << ":" << pair.second << "\n";
        }
    }

    bool ReadOsmMap::hasHighwaysFile(const std::filesystem::path& mapInput)
    {
        auto path = mapInput;
        path.replace_extension(WriteMapFromOsm::HIGHWAYS_EXTENSION);
        return fs::exists(path);
    }

    bool ReadOsmMap::hasMappingFile(const std::filesystem::path& mapInput)
    {
        auto path = mapInput;
        path.replace_extension(WriteMapFromOsm::MAPPING_EXTENSION);
        return fs::exists(path);
    }

    bool ReadOsmMap::hasUndirectedFile(const std::filesystem::path& mapInput)
    {
        auto path = mapInput;
        path.replace_extension(WriteMapFromOsm::UNDIRECTED_GRAPH_EXTENSION);
        return fs::exists(path);
    }

    bool ReadOsmMap::hasDirectedFile(const std::filesystem::path& mapInput)
    {
        auto path = mapInput;
        path.replace_extension(WriteMapFromOsm::DIRECTED_GRAPH_EXTENSION);
        return fs::exists(path);
    }

    bool ReadOsmMap::readFileMissing(const std::filesystem::path& osmPath)
    {
        return !hasDirectedFile(osmPath) || !hasHighwaysFile(osmPath) || !hasUndirectedFile(osmPath) || !hasMappingFile(osmPath);
    }

    template<typename T>
    std::from_chars_result from_string_view(const std::string_view& view, T& output)
    {
        return std::from_chars(view.data(), view.data() + view.size(), output);
    }
    template<>
    std::from_chars_result from_string_view(const std::string_view& view, std::string& string)
    {
        string = std::string(view);
        return std::from_chars_result{nullptr, std::errc::no_message};
    }

    void ReadOsmMap::readCustomFormat(const std::filesystem::path& generalizationFile, Models::EmbeddedGraph& graph,
        std::unordered_map<std::size_t, std::string>& highwayTags, int& epsg)
    {
        auto highWayFile = generalizationFile;
        highWayFile.replace_extension(WriteMapFromOsm::HIGHWAYS_EXTENSION);
        if (!fs::exists(highWayFile)) throw std::runtime_error("Could not find .highways next to input file!");
        readCustomFormat(generalizationFile, graph,epsg);

        readHighwayTags(highWayFile, highwayTags);
    }

    void ReadOsmMap::readHighwayTags(const std::filesystem::path& highwayInput,
        std::unordered_map<std::size_t, std::string>& highwayTags)
    {
        if (!fs::exists(highwayInput)) throw std::runtime_error("Could not find .highways next to input file!");
        if (highwayInput.extension() != WriteMapFromOsm::HIGHWAYS_EXTENSION) throw std::runtime_error("Expected file with extension " + std::string{ WriteMapFromOsm::HIGHWAYS_EXTENSION });

        std::ifstream stream(highwayInput.string());
        std::string line;
        std::vector<std::string_view> parts;
        while (std::getline(stream, line))
        {
            parts.clear();
            GCpp::String::splitOnWS(line, parts);
            if (parts.size() != 2) throw std::runtime_error("Expected 2 elements on a line in highway!");
            std::size_t id;
            std::string highway;
            from_string_view(parts[0], id);
            from_string_view(parts[1], highway);
            highwayTags[id] = highway;
        }
    }

    void ReadOsmMap::readCustomFormat(const std::filesystem::path& mapInput, Models::EmbeddedGraph& graph, int& epsg)
    {
        if (!fs::exists(mapInput)) throw std::runtime_error("Map path does not exist");
        if (mapInput.extension() != WriteMapFromOsm::DIRECTED_GRAPH_EXTENSION) {
            throw std::runtime_error("Invalid extension, expected " + std::string{ WriteMapFromOsm::DIRECTED_GRAPH_EXTENSION });
        }
        std::ifstream stream(mapInput.string());
        std::string line;
        ParsingState state = ParsingState::None;
        std::size_t expectedEdges = 0;
        std::size_t expectedVertices = 0;
        std::size_t currentVertex = 0;
        std::size_t nextEdgeId = 0;
        std::size_t maxVertId = 0;

        // EPSG line
        std::getline(stream, line);
        line = line.substr(std::string_view{ WriteMapFromOsm::EPSG_LINE_START }.size());
        epsg = std::stoi(line);

        while (std::getline(stream, line))
        {
            // Ignore empty lines.
            if (line.empty()) continue;

            auto view = GCpp::String::strip(line);
            std::string strippedLine(view);
            std::stringstream lineStream(strippedLine);

            if (GCpp::String::startsWith(strippedLine, '#'))
            {
                if (state == ParsingState::None)
                {
                    char c;
                    lineStream >> c;
                    lineStream >> expectedVertices;
                    state = ParsingState::Vertices;
                    graph = Models::EmbeddedGraph(expectedVertices);
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

                if (lineStream.bad() || lineStream.eof())
                {
                    boost::put(GCpp::DS::vertex_id_t{}, graph, currentVertex, currentVertex);
                    maxVertId = std::max(maxVertId, currentVertex);
                    ++currentVertex;
                }
                else
                {
                    boost::put(GCpp::DS::vertex_id_t{}, graph, currentVertex, id);
                    maxVertId = std::max(maxVertId, id);
                }
            }
            else if (state == ParsingState::Edges)
            {
                std::size_t src, target, eId;
                lineStream >> src >> target >> eId;
                auto result = boost::add_edge(src, target, graph);
                assert(result.second);
                boost::put(GCpp::DS::edge_id_t{}, graph, result.first, eId);
                nextEdgeId = std::max(nextEdgeId, eId + 1);
            }
        }
        boost::get_property(graph, GCpp::DS::next_vertex_id_t{}) = maxVertId+1;
        boost::get_property(graph, GCpp::DS::next_edge_id_t{}) = nextEdgeId;
    }

    void ReadOsmMap::readCustomFormat(const std::filesystem::path& generalizationFile,
                                      Models::UndirectedEmbeddedGraph& graph, int& epsg)
    {
        namespace fs = std::filesystem;
        if (!fs::exists(generalizationFile)) throw std::runtime_error("Generalization path does not exist");
        if (generalizationFile.extension() != WriteMapFromOsm::UNDIRECTED_GRAPH_EXTENSION) {
            throw std::runtime_error("Invalid extension, expected " + std::string{ WriteMapFromOsm::UNDIRECTED_GRAPH_EXTENSION });
        }
        std::ifstream stream(generalizationFile.string());
        std::string line;
        ParsingState state = ParsingState::None;
        std::size_t expectedEdges = 0;
        std::size_t expectedVertices = 0;
        std::size_t currentVertex = 0;
        std::size_t nextEdgeId = 0;
        std::size_t maxVertId = 0;

        // Read EPSG line
        std::getline(stream, line);
        line = line.substr(std::string_view{ WriteMapFromOsm::EPSG_LINE_START }.size());
        epsg = std::stoi(line);

        while (std::getline(stream, line))
        {
            // Ignore empty lines.
            if (line.empty()) continue;

            auto view = GCpp::String::strip(line);
            std::string strippedLine(view);
            std::stringstream lineStream(strippedLine);

            if (GCpp::String::startsWith(strippedLine, '#'))
            {
                if (state == ParsingState::None)
                {
                    char c;
                    lineStream >> c;
                    lineStream >> expectedVertices;
                    state = ParsingState::Vertices;
                    graph = Models::UndirectedEmbeddedGraph(expectedVertices);
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

                if (lineStream.bad() || lineStream.eof())
                {
                    boost::put(GCpp::DS::vertex_id_t{}, graph, currentVertex, currentVertex);
                    maxVertId = std::max(maxVertId, currentVertex);
                    ++currentVertex;
                }
                else
                {
                    boost::put(GCpp::DS::vertex_id_t{}, graph, currentVertex, id);
                    maxVertId = std::max(maxVertId, id);
                }
            }
            else if (state == ParsingState::Edges)
            {
                std::size_t src, target, eId;
                lineStream >> src >> target >> eId;
                auto result = boost::add_edge(src, target, graph);
                assert(result.second);
                boost::put(GCpp::DS::edge_id_t{}, graph, result.first, eId);
                nextEdgeId = std::max(nextEdgeId, eId + 1);
            }
        }
        boost::get_property(graph, GCpp::DS::next_vertex_id_t{}) = maxVertId;
        boost::get_property(graph, GCpp::DS::next_edge_id_t{}) = nextEdgeId;
    }

    void ReadOsmMap::readCustomFormat(const std::filesystem::path& generalizationFile,
        Models::UndirectedEmbeddedGraph& graph,
        std::unordered_map<std::size_t, std::size_t>& directedToUndirectedEdgeIdMap, int& epsg)
    {
        auto mappingPath = generalizationFile;
        mappingPath.replace_extension(WriteMapFromOsm::MAPPING_EXTENSION);
        if (!fs::exists(mappingPath)) throw std::runtime_error("Expected .undirmapping file next to map for mapping reading!");
        readCustomFormat(generalizationFile, graph, epsg);
        readMapping(mappingPath, directedToUndirectedEdgeIdMap);
    }

    void ReadOsmMap::readMapping(const std::filesystem::path& mappingFile, 
        std::unordered_map<std::size_t, std::size_t>& directedToUndirectedEdgeIdMap)
    {
        if (!fs::exists(mappingFile)) throw std::runtime_error("Expected .undirmapping file next to map for mapping reading!");
        std::ifstream stream(mappingFile.string());
        std::string line;
        while (std::getline(stream, line))
        {
            std::stringstream lineStream(line);
            std::size_t from, to;
            lineStream >> from >> to;
            directedToUndirectedEdgeIdMap[from] = to;
        }
    }

    void ProcessUndirectedToDirected::operator()(const std::filesystem::path& osmInput, int epsg)
    {
        ReadOsmMap reader;
        reader.setMapEpsg(epsg); // Set the EPSG
        Models::EmbeddedGraph graph;
        std::unordered_map<std::size_t, std::string> highwayTags;
        reader(osmInput, graph, highwayTags);
        WriteMapFromOsm writer;
        // Write highway tags and directed graph
        writer(osmInput, graph, highwayTags, epsg);

        // Compute undirected graph and write.
        Models::UndirectedEmbeddedGraph undirGraph;
        std::map < std::size_t, std::pair < std::size_t, bool>> mapping;
        GCpp::DS::convertDirectedToUndirectedWithEdgeOrder(graph, undirGraph, mapping);
        std::unordered_map<std::size_t, std::size_t> savableMap;
        for(auto kv : mapping)
        {
            savableMap[kv.first] = kv.second.first;
        }
        writer(osmInput, undirGraph, savableMap, epsg);
    }

    void ProcessUndirectedToDirected::verifyMapping(const std::filesystem::path& basePathToFiles)
    {
        // Read directed, undirected and mapping.
    }

    void WriteMapFromOsm::operator()(const std::filesystem::path& generalizationFile,
                                     const Models::EmbeddedGraph& graph, int epsg)
    {
        auto outputFile = generalizationFile;
        outputFile.replace_extension(DIRECTED_GRAPH_EXTENSION);
        std::ofstream out(outputFile.string());
        out << std::setprecision(std::numeric_limits<Models::NT>::digits10);
        out << "# EPSG " << epsg << '\n';
        out << "# " << boost::num_vertices(graph) << '\n';
        for(auto v : iter::range(boost::vertices(graph)))
        {
            const auto pos = GCpp::DS::get_vertex_location(v, graph);
            out << pos.x() << ' ' << pos.y() << ' ' << GCpp::DS::getVertexId(graph, v) << '\n';
        }
        out << "# " << boost::num_edges(graph) << '\n';
        for (auto e : iter::range(boost::edges(graph)))
        {
            const auto src = boost::source(e, graph);
            const auto target = boost::target(e, graph);
            out << src << ' ' << target << ' ' << GCpp::DS::getEdgeId(graph, e) << '\n';
        }
    }

    void WriteMapFromOsm::operator()(const std::filesystem::path& mapOutput, const Models::EmbeddedGraph& graph,
        const std::unordered_map<std::size_t, std::string>& highwayTags, int epsg)
    {
        this->operator()(mapOutput, graph, epsg);
        auto highwayPath = mapOutput;
        highwayPath.replace_extension(HIGHWAYS_EXTENSION);
        std::ofstream stream(highwayPath.string());
        for(const auto& kv: highwayTags)
        {
            stream << kv.first << ' ' << kv.second << '\n';
        }
    }

    void WriteMapFromOsm::operator()(const std::filesystem::path& generalizationFile,
                                     const Models::UndirectedEmbeddedGraph& graph, int epsg)
    {
        auto outputFile = generalizationFile;
        outputFile.replace_extension(UNDIRECTED_GRAPH_EXTENSION);
        std::ofstream out(outputFile.string());
        out << std::setprecision(std::numeric_limits<Models::NT>::digits10);
        out << "# EPSG " << epsg << '\n';
        out << "# " << boost::num_vertices(graph) << '\n';
        for (auto v : iter::range(boost::vertices(graph)))
        {
            const auto pos = GCpp::DS::get_vertex_location(v, graph);
            out << pos.x() << ' ' << pos.y() << ' ' << GCpp::DS::getVertexId(graph, v) << '\n';
        }
        out << "# " << boost::num_edges(graph) << '\n';
        for (auto e : iter::range(boost::edges(graph)))
        {
            const auto src = boost::source(e, graph);
            const auto target = boost::target(e, graph);
            out << src << ' ' << target << ' ' << GCpp::DS::getEdgeId(graph, e) << '\n';
        }
    }

    void WriteMapFromOsm::operator()(const std::filesystem::path& mapOutput,
        const Models::UndirectedEmbeddedGraph& graph,
        const std::unordered_map<std::size_t, std::size_t>& directedToUndirectedEdgeIdMap, int epsg)
    {
        this->operator()(mapOutput, graph, epsg);
        auto mappingPath = mapOutput;
        mappingPath.replace_extension(MAPPING_EXTENSION);
        std::ofstream stream(mappingPath.string());
        for(const auto& kv: directedToUndirectedEdgeIdMap)
        {
            stream << kv.first << ' ' << kv.second << '\n';
        }
    }
}
