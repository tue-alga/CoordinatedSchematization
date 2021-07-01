#ifndef SCHEMATLIB_IO_READWRITEOSMMAP_H
#define SCHEMATLIB_IO_READWRITEOSMMAP_H
#include <filesystem>
#include <fstream>
#include "SchematLib/Models/EmbeddedGraph.h"
namespace SchematLib::IO
{
    class ReadOsmMap
    {
        int m_epsg = 4326;
        enum class ParsingState
        {
            None,
            Vertices,
            Edges
        };
    public:
        static constexpr char* OSM_EXTENSION = ".osm";
        void setMapEpsg(int epsg);
        int epsg() const;
        void operator()(const std::filesystem::path& mapInput, Models::EmbeddedGraph& graph, std::unordered_map<std::size_t, std::string>& highwayTags);

        static bool hasHighwaysFile(const std::filesystem::path& mapInput);
        static bool hasMappingFile(const std::filesystem::path& mapInput);
        static bool hasUndirectedFile(const std::filesystem::path& mapInput);
        static bool hasDirectedFile(const std::filesystem::path& mapInput);
        static bool readFileMissing(const std::filesystem::path& osmPath);

        void readCustomFormat(const std::filesystem::path& mapInput, Models::EmbeddedGraph& graph, std::unordered_map<std::size_t, std::string>& highwayTags, int& epsg);
        void readHighwayTags(const std::filesystem::path& highwayInput, std::unordered_map<std::size_t, std::string>& highwayTags);

        void readCustomFormat(const std::filesystem::path& mapInput, Models::EmbeddedGraph& graph, int& epsg);
        void readCustomFormat(const std::filesystem::path& mapInput, Models::UndirectedEmbeddedGraph& graph, int& epsg);
        void readCustomFormat(const std::filesystem::path& mapInput, Models::UndirectedEmbeddedGraph& graph, std::unordered_map<std::size_t, std::size_t>& directedToUndirectedEdgeIdMap, int& epsg);
        void readMapping(const std::filesystem::path& mappingFile, std::unordered_map<std::size_t, std::size_t>& directedToUndirectedEdgeIdMap);
    };

    /**
     * \brief Generates highway, undirected mapping, undirected and directed graph file from OSM file.
     */
    class ProcessUndirectedToDirected
    {
    public:
        void operator()(const std::filesystem::path& osmInput, int epsg);

        /**
         * \brief Verifies the directed to undirected mapping
         * \param basePathToFiles 
         */
        void verifyMapping(const std::filesystem::path& basePathToFiles);
    };

    /**
     * \brief Writes a custom format map file based on read OSM data.
     */
    class WriteMapFromOsm
    {
    public:
        static constexpr char* HIGHWAYS_EXTENSION = ".highways";
        static constexpr char* MAPPING_EXTENSION = ".undirmapping";
        static constexpr char* UNDIRECTED_GRAPH_EXTENSION = ".undigraph";
        static constexpr char* DIRECTED_GRAPH_EXTENSION = ".digraph";
        static constexpr char* EPSG_LINE_START = "# EPSG";
        void operator()(const std::filesystem::path& mapOutput, const Models::EmbeddedGraph& graph, int epsg);
        /**
         * \brief Outputs digraph and highways file
         * \param mapOutput 
         * \param graph 
         * \param highwayTags 
         */
        void operator()(const std::filesystem::path& mapOutput, const Models::EmbeddedGraph& graph,
            const std::unordered_map<std::size_t, std::string>& highwayTags, int epsg);

        void operator()(const std::filesystem::path& mapOutput, const Models::UndirectedEmbeddedGraph& graph, int epsg);

        /**
         * \brief Outputs undirected graph and mapping file
         * \param mapOutput 
         * \param graph 
         * \param directedToUndirectedEdgeIdMap 
         */
        void operator()(const std::filesystem::path& mapOutput, const Models::UndirectedEmbeddedGraph& graph, 
            const std::unordered_map<std::size_t, std::size_t>& directedToUndirectedEdgeIdMap, int epsg);
    };
}
#endif