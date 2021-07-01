#include "SchematLib/Models/SchematizationStages.h"

namespace SchematLib::Models {
    void SchematizationStages::reset()
    {
        m_currentStage = 0;
    }

    SchematizationStages::CanonicalAccessGuard::CanonicalAccessGuard(SchematizationStages* parent, bool use_canonical):
        m_parent(parent), m_use_canonical(use_canonical)
    {
        for (std::size_t i = 0; i < m_parent->totalStages(); ++i)
        {
            auto graph = m_parent->stageGraphAt(i);
            previousAccess.push_back(graph->m_use_canonical_edges_ids);
            graph->m_use_canonical_edges_ids = m_use_canonical;
        }
    }

    SchematizationStages::CanonicalAccessGuard::~CanonicalAccessGuard()
    {
        if (m_wasReset) return;
        for (std::size_t i = 0; i < m_parent->totalStages(); ++i)
        {
            auto graph = m_parent->stageGraphAt(i);
            graph->m_use_canonical_edges_ids = previousAccess[i];
        }
    }

    SchematizationStages::CanonicalAccessGuard SchematizationStages::guardChangeAccess(bool useCanonical)
    {
        return CanonicalAccessGuard(this, useCanonical);
    }

    std::size_t SchematizationStages::totalStages() const
    {
        return m_currentStage;
    }

    const std::vector<std::string>& SchematizationStages::stageNames() const
    {
        return m_simplificationStageName;
    }

    std::shared_ptr<SchematizationStages::Graph> SchematizationStages::stageGraphAt(std::size_t stage) const
    {
        return m_simplificationStages[stage];
    }

    const MapSimplification::EdgeTrackingObserver& SchematizationStages::observerAt(std::size_t stage) const
    {
        return m_simplificationTrackerPerStage[stage];
    }

    std::string SchematizationStages::nameAt(std::size_t stage) const
    {
        return m_simplificationStageName[stage];
    }

    void SchematizationStages::assignNextStage(const std::string& stageName, const Graph& graph,
                                               const MapSimplification::EdgeTrackingObserver& observer)
    {
        if (m_currentStage >= m_simplificationStageName.size())
        {
            addNewBlankStage();
        }
        std::cout << "Assigning stage " << m_currentStage << ": " << stageName << "\n";
        m_simplificationStageName[m_currentStage] = stageName;
        if(m_validate)
        {
            std::set<std::size_t> currentEdgeIds;
            std::set<std::size_t> currentVertexIds;
            GCpp::DS::computeEdgeIdSet(graph, currentEdgeIds);
            GCpp::DS::computeVertexIdSet(graph, currentVertexIds);
            if (!GCpp::DS::hasUniqueVertexIds(graph)) throw std::runtime_error("Duplicate vertex ID!");
            if (currentEdgeIds.size() != boost::num_edges(graph)) throw std::runtime_error("Duplicate edge ID!");
            observer.verifyIntegrity(currentEdgeIds, currentVertexIds);
        }
        m_simplificationTrackerPerStage[m_currentStage] = observer;

        *m_simplificationStages[m_currentStage] = graph;
        ++m_currentStage;
    }

    std::size_t SchematizationStages::currentStageToAssign() const
    {
        return m_currentStage;
    }

    std::shared_ptr<SchematizationStages::Graph> SchematizationStages::currentStageGraph() const
    {
        return m_simplificationStages[m_currentStage];
    }

    std::shared_ptr<SchematizationStages::Graph> SchematizationStages::prevStageGraph() const
    {
        verifyNotFirst();
        return m_simplificationStages[m_currentStage - 1];
    }

    void SchematizationStages::verifyNotFirst() const
    {
        if (m_currentStage == 0) throw std::runtime_error("Cannot access to prev stage when at 0!");
    }

    void SchematizationStages::setNames(const std::vector<std::string>& names)
    {
        if (names.size() != m_simplificationStageName.size()) throw std::runtime_error("Incompatible sizes");
        m_simplificationStageName = names;
    }
    void SchematizationStages::resize(std::size_t numberOfStages)
    {
        m_simplificationStages.clear();
        for (auto i = 0; i < numberOfStages; ++i) m_simplificationStages.push_back(std::make_shared<Models::UndirectedEmbeddedGraph>());
        m_simplificationTrackerPerStage.resize(numberOfStages, { });
        m_simplificationStageName.clear();
        m_simplificationStageName.resize(numberOfStages, "");
    }
    void SchematizationStages::addNewBlankStage()
    {
        m_simplificationStages.push_back(std::make_shared<Models::UndirectedEmbeddedGraph>());
        m_simplificationTrackerPerStage.push_back({ });
        m_simplificationStageName.push_back("");
    }
}
