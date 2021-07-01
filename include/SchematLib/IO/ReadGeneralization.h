#ifndef SCHEMATLIB_IO_READGENERALIZATION_H
#define SCHEMATLIB_IO_READGENERALIZATION_H
#include <filesystem>
#include <fstream>

#include "SchematLib/MapSimplification/FaceMergeTraits.h"
#include "SchematLib/Models/EmbeddedGraph.h"
namespace SchematLib::IO
{
    struct ReadGeneralization
    {
        static std::string_view strip(const std::string& string);

        static bool startsWith(const std::string& str, char c);

        static bool startsWith(const std::string& str, const std::string& needle);

        enum class ParsingState
        {
            None,
            Vertices,
            Edges
        };

        void operator()(const std::filesystem::path& generalizationFile, Models::UndirectedEmbeddedGraph& graph);
    };
    struct WriteGeneralization
    {
        void operator()(const std::filesystem::path& generalizationFile, const Models::UndirectedEmbeddedGraph& graph);
        void operator()(std::ostream& output, const Models::UndirectedEmbeddedGraph& graph);

        // CGAL
        void operator()(const std::filesystem::path& generalizationFile, const MapSimplification::FaceMergeTraits::Arr& graph);
        void operator()(std::ostream& output, const MapSimplification::FaceMergeTraits::Arr& graph);
    };
}
#endif