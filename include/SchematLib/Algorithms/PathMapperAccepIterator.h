#ifndef SCHEMATLIB_ALGORITHMS_PATHMAPPERACCEPTITERATOR_H
#define SCHEMATLIB_ALGORITHMS_PATHMAPPERACCEPTITERATOR_H
#include "SchematLib/Models/SchematizationStages.h"
#include "SchematLib/MapSimplification/PathMapper.h"

namespace SchematLib::Algorithms
{
    class PathMapperAcceptIterator
    {
    public:
        using Obs = SchematLib::MapSimplification::EdgeTrackingObserver;
        using Mapper = SchematLib::MapSimplification::PathMapper;
    private:
        // Reference to reconstructed path
        std::vector<std::size_t>* m_reconstructed = nullptr;
        std::optional<std::size_t> prevPrevEdge;
        std::optional<std::size_t> prevEdge;
        bool m_verbose = true;
    public:
        PathMapperAcceptIterator() {}
        PathMapperAcceptIterator(std::vector<std::size_t>* reconstructed) :m_reconstructed(reconstructed) {}

        void setVerbose(bool verbose)
        {
            m_verbose = verbose;
        }

        PathMapperAcceptIterator& operator*() { return *this; }
        PathMapperAcceptIterator& operator=(const Mapper::TerminalElement& result)
        {
            // Ignore Vertex
            if (Obs::isVertex(result)) {
                if (m_verbose) { std::cout << "[Output] Mapping to vertex " << std::get<Obs::Vertex>(result).id << '\n'; }
                return *this;
            }
            if (Obs::isEdge(result))
            {
                auto e = std::get<Obs::Edge>(result);
                const auto id = e.id;
                // Deduplicate: allow at most two consecutive edges of the same ID. Reconstruct properly later.
                if (prevPrevEdge.has_value() && prevEdge.has_value() && (prevPrevEdge.value() != id || prevEdge.value() != id))
                {
                    m_reconstructed->push_back(e.id);
                    if (m_verbose) {
                        std::cout << "[Output] Mapping to edge " << e.id << '\n';
                    }
                }
                else
                {
                    if (m_verbose) {
                        std::cout << "[Output] Skipping edge " << e.id << '\n';
                    }
                }
                prevPrevEdge = prevEdge;
                prevEdge = id;

            }
            else if (Obs::isDeleted(result))
            {
                if (m_verbose) {
                    std::cout << "[Output] Deleted element" << '\n';
                }
                // Ignore
                return *this;
            }
            else
            {
                std::cerr << "ERROR!!!!!!!!!! Unknown type in result variant!\n";
            }

            return *this;
        }
    };
}
#endif
