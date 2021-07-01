#ifndef SCHEMATLIB_MODELS_SCHEMATIZATIONSTAGES_H
#define SCHEMATLIB_MODELS_SCHEMATIZATIONSTAGES_H
#include <vector>
#include "SchematLib/MapSimplification/FaceMergeTraits.h"

#include "EmbeddedGraph.h"
#include "SchematLib/MapSimplification/EdgeTrackingObserver.h"
#include <libzippp/libzippp.h>
#include <filesystem>

#include "SchematLib/IO/ReadGeneralization.h"

namespace SchematLib::Models
{
    class SchematizationStages
    {
        using Graph = Models::UndirectedEmbeddedGraph;
        std::vector< std::shared_ptr<Models::UndirectedEmbeddedGraph>> m_simplificationStages;
        std::vector<MapSimplification::EdgeTrackingObserver> m_simplificationTrackerPerStage;
        std::vector<std::string> m_simplificationStageName;

        std::size_t m_currentStage = 0;


        bool m_validate = true;
    public:
        void reset();

        struct CanonicalAccessGuard
        {
            friend class SchematizationStages;
            SchematizationStages* m_parent;
            bool m_use_canonical;
            std::vector<bool> previousAccess;
            bool m_wasReset = false;
            CanonicalAccessGuard(SchematizationStages* parent, bool use_canonical);
        public:
            CanonicalAccessGuard(CanonicalAccessGuard&& guard) = default;
            CanonicalAccessGuard(const CanonicalAccessGuard&) = delete;
            CanonicalAccessGuard& operator=(CanonicalAccessGuard&& guard) = default;
            CanonicalAccessGuard& operator=(const CanonicalAccessGuard&) = delete;
            ~CanonicalAccessGuard();
        };

        CanonicalAccessGuard guardChangeAccess(bool useCanonical);

        void save(const std::filesystem::path& path)
        {
            libzippp::ZipArchive archive(path.string());
            archive.open(libzippp::ZipArchive::Write);
            IO::WriteGeneralization writer;
            // Write stage graphs
            for(auto i = 0; i < m_currentStage; ++i)
            {
                std::stringstream output;
                writer(output, *m_simplificationStages[i]);
                auto result = output.str();
                archive.addData("stage_" + std::to_string(i) + ".graph", result.c_str(), result.size());
            }
            // Write stagenames
            {
                std::stringstream namesData;
                for (auto stageName : m_simplificationStageName)
                {
                    namesData << stageName << '\n';
                }
                auto result = namesData.str();
                (void)archive.addData("stage_names", result.c_str(), result.size());
            }
            // Write simplification trackers.

            archive.close();
        }
        void load(const std::filesystem::path& path)
        {

        }

        std::size_t totalStages() const;

        const std::vector<std::string>& stageNames() const;

        void setNames(const std::vector<std::string>& names);

        void resize(std::size_t numberOfStages);

        void addNewBlankStage();

        void verifyNotFirst() const;

        std::shared_ptr<Graph> stageGraphAt(std::size_t stage) const;

        std::size_t lastAssignedStage() const
        {
            return m_currentStage - 1;
        }

        bool hasStages() const
        {
            return m_currentStage > 0;
        }

        const MapSimplification::EdgeTrackingObserver& observerAt(std::size_t stage) const;

        std::string nameAt(std::size_t stage) const;


        std::shared_ptr<Graph> prevStageGraph() const;

        std::shared_ptr<Graph> currentStageGraph() const;

        std::size_t currentStageToAssign() const;

        void assignNextStage(const std::string& stageName, const Graph& graph,
                             const MapSimplification::EdgeTrackingObserver& observer);
    };
}
#endif
