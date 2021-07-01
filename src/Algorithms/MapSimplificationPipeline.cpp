#include "SchematLib/Algorithms/MapSimplificationPipeline.h"

#include <boost/property_tree/ptree.hpp>

#include "SchematLib/MapSimplification/BboxFaceMerge.h"
#include "SchematLib/MapSimplification/RemoveFaceTendrils.h"

namespace SchematLib::Algorithms {
    Models::NT MapSimplificationPipeline::maxBboxAxis() const
    {
        return std::max(m_bboxWidth, m_bboxHeight);
    }

    std::string MapSimplificationPipeline::stageName(const std::string& name) const
    {
        if (m_currentSettings.m_totalIterations == 1 && !m_currentSettings.m_untillConvergence)
        {
            return name;
        }
        return std::to_string(m_currentSettings.m_currentIteration) + ":" + name;
    }

    void MapSimplificationPipeline::applyRemoveTendrils(Graph& target, SimpObs& simplificationObserver)
    {
        MapSimplification::RemoveTendrils<Models::UndirectedEmbeddedGraph, SimpObs> removeTendrils;
        removeTendrils.setMaxLength(m_currentSettings.m_maxPruneLength);

        removeTendrils(target, simplificationObserver);
    }

    void MapSimplificationPipeline::applyDegree2Smoothing(Graph& target, SimpObs& simplificationObserver)
    {
        MapSimplification::Degree2VertexPreprocess<Models::UndirectedEmbeddedGraph> degree2Smoothing(
            m_currentSettings.m_maxTurningAngle / 180. * M_PI);
        degree2Smoothing.setRepeatUntilConvergence(m_currentSettings.m_edgeSmoothUntilConvergence);
        degree2Smoothing(target, simplificationObserver);
    }

    void MapSimplificationPipeline::applyDegree2SmoothingPostFM(Graph& target, SimpObs& simplificationObserver)
    {
        MapSimplification::Degree2VertexPreprocess<Models::UndirectedEmbeddedGraph> degree2Smoothing(
            m_currentSettings.m_maxTurningAnglePostFM / 180. * M_PI);
        degree2Smoothing.setRepeatUntilConvergence(m_currentSettings.m_edgeSmoothUntilConvergence);
        degree2Smoothing(target, simplificationObserver);
    }

    void MapSimplificationPipeline::applyVertexCollapse(Graph& target, SimpObs& simplificationObserver)
    {
        Graph tempCopy = target;
        const auto radius = m_currentSettings.m_radius * std::max(m_bboxHeight, m_bboxWidth);
        if (m_preprocessApproach == VerexOrderApproach::DenseVertexOrder)
        {
            //// Setup merger algorithm
            MapSimplification::VertexMerger<Models::UndirectedEmbeddedGraph, Models::UndirectedGraphNodeIndex, SimpObs,
                                            MapSimplification::MergeOrders::MostDenseVertexOrder> merger(radius);
            merger.setUseMeanPosition(m_currentSettings.m_useMeanPosition);
            //// Compute merge groups
            merger(tempCopy, simplificationObserver, target);
        }
        else
        {
            //// Setup merger algorithm
            MapSimplification::VertexMerger<Models::UndirectedEmbeddedGraph, Models::UndirectedGraphNodeIndex, SimpObs>
                merger(radius);
            merger.setUseMeanPosition(m_currentSettings.m_useMeanPosition);
            //// Compute merge groups
            merger(tempCopy, simplificationObserver, target);
        }
    }

    void MapSimplificationPipeline::applyVertexToEdgeMerge(Graph& target, SimpObs& simplificationObserver)
    {
        VertexEdgeMerger vertexEdgeMerger;
        vertexEdgeMerger.setThreshold(m_currentSettings.m_radius * maxBboxAxis());
        vertexEdgeMerger.setRepeatUntilConvergence(m_currentSettings.m_edgeMergeUntilConverge);
        vertexEdgeMerger(target, simplificationObserver);
    }

    void MapSimplificationPipeline::applyFaceMerge(Graph& target, SimpObs& simplificationObserver)
    {
        namespace ms = SchematLib::MapSimplification;
        // Alternative with bounding box.
        std::array<Models::Point, 2> minMaxBBPoints;
        GCpp::DS::computeBoundingBox(target, minMaxBBPoints);

        Models::NT area = (minMaxBBPoints.at(1).x() - minMaxBBPoints.at(0).x()) * (minMaxBBPoints.at(1).y() -
            minMaxBBPoints.at(0).y());
        auto areaThreshold = m_currentSettings.m_faceCollapseAreaFraction * area;
        using Selector_t = ms::FaceSelectApproach::SelectWithSmallArea<Models::NT, ms::FaceMergeTraits::Kernel>;
        Selector_t selector;
        selector.setAreaThreshold(areaThreshold);

        // using Selector_t = ms::FaceSelectApproach::SelectWithBoundedComplexity;
        // Selector_t selector(m_maxFaceComplexity);

        MapSimplification::FaceMerger<Models::UndirectedEmbeddedGraph,
                                      Selector_t, ms::FaceMergeOrders::LowestFaceComplexityOrder,
                                      ms::FaceMergeApproach::DeleteSmallFacesToNeighbourViaContraction,
                                      SimpObs> faceMerger(0, selector);
        faceMerger.setNormalizeCoordinates(false);
        faceMerger(target, simplificationObserver);
    }

    void MapSimplificationPipeline::applyContinuationFaceMerge(Graph& target, SimpObs& simplificationObserver)
    {
        namespace ms = SchematLib::MapSimplification;
        // Alternative with bounding box.
        std::array<Models::Point, 2> minMaxBBPoints;
        GCpp::DS::computeBoundingBox(target, minMaxBBPoints);

        Models::NT area = (minMaxBBPoints.at(1).x() - minMaxBBPoints.at(0).x()) * (minMaxBBPoints.at(1).y() -
            minMaxBBPoints.at(0).y());
        auto areaThreshold = m_currentSettings.m_faceCollapseAreaFraction * area;
        using Selector_t = ms::FaceSelectApproach::SelectWithSmallArea<Models::NT, ms::FaceMergeTraits::Kernel>;
        Selector_t selector;
        selector.setAreaThreshold(areaThreshold);

        // using Selector_t = ms::FaceSelectApproach::SelectWithBoundedComplexity;
        // Selector_t selector(m_maxFaceComplexity);

        MapSimplification::FaceMerger<Models::UndirectedEmbeddedGraph,
            Selector_t, ms::FaceMergeOrders::LowestFaceComplexityOrder,
            ms::FaceMergeApproach::CollapseGoodContinuityFacesByEdge,
            SimpObs> faceMerger(0, selector);
        faceMerger.setNormalizeCoordinates(false);
        faceMerger(target, simplificationObserver);
    }

    void MapSimplificationPipeline::applyBboxFaceMerge(Graph& target, SimpObs& simplificationObserver)
    {
        namespace ms = SchematLib::MapSimplification;
        // Alternative with bounding box.
        std::array<Models::Point, 2> minMaxBBPoints;
        GCpp::DS::computeBoundingBox(target, minMaxBBPoints);

        Models::NT area = (minMaxBBPoints.at(1).x() - minMaxBBPoints.at(0).x()) * (minMaxBBPoints.at(1).y() -
            minMaxBBPoints.at(0).y());
        auto areaThreshold = m_currentSettings.m_faceCollapseAreaFraction * area;
        using Selector_t = ms::FaceSelectApproach::SelectWithSmallArea<Models::NT, ms::FaceMergeTraits::Kernel>;
        Selector_t selector;
        selector.setAreaThreshold(areaThreshold);

        MapSimplification::FaceMerger<Models::UndirectedEmbeddedGraph,
            Selector_t, ms::FaceMergeOrders::LowestFaceComplexityOrder,
            ms::FaceMergeApproach::BboxFaceMerge,
            SimpObs> faceMerger(0, selector);
        faceMerger(target, simplificationObserver);
    }

    void MapSimplificationPipeline::applyLargestComponentSelect(Graph& target, SimpObs& simplificationObserver)
    {
        MapSimplification::LargestConnectedComponent<Graph> lcc;
        lcc(target, simplificationObserver);
    }

    void MapSimplificationPipeline::applyVertexCollapsePostFM(Graph& target, SimpObs& simplificationObserver)
    {
        Graph tempCopy = target;
        const auto radius = m_currentSettings.m_postRadius * std::max(m_bboxHeight, m_bboxWidth);
        if (m_preprocessApproach == VerexOrderApproach::DenseVertexOrder)
        {
            //// Setup merger algorithm
            MapSimplification::VertexMerger<Models::UndirectedEmbeddedGraph, Models::UndirectedGraphNodeIndex, SimpObs,
                                            MapSimplification::MergeOrders::MostDenseVertexOrder> merger(radius);
            merger.setUseMeanPosition(m_currentSettings.m_useMeanPosition);
            //// Compute merge groups
            merger(tempCopy, simplificationObserver, target);
        }
        else
        {
            //// Setup merger algorithm
            MapSimplification::VertexMerger<Models::UndirectedEmbeddedGraph, Models::UndirectedGraphNodeIndex, SimpObs>
                merger(radius);
            merger.setUseMeanPosition(m_currentSettings.m_useMeanPosition);
            //// Compute merge groups
            merger(tempCopy, simplificationObserver, target);
        }
    }

    void MapSimplificationPipeline::applyRemoveTendrilsPostFM(Graph& target, SimpObs& simplificationObserver)
    {
        MapSimplification::RemoveTendrils<Models::UndirectedEmbeddedGraph, SimpObs> removeTendrils;
        removeTendrils.setMaxLength(m_currentSettings.m_maxPruneLengthPostFM);

        removeTendrils(target, simplificationObserver);
    }

    void MapSimplificationPipeline::applyFaceTendrilRemove(Graph& target, SimpObs& simplificationObserver)
    {
        MapSimplification::RemoveFaceTendrils<Graph, SimpObs> remover;
        remover(target, simplificationObserver);
    }

    void MapSimplificationPipeline::Settings::loadFromTree(const boost::property_tree::ptree& tree)
    {
#define FROM_TREE(name) m_##name = tree.get<std::decay_t<decltype(m_##name)>>(#name)
        FROM_TREE(radius);
        FROM_TREE(postRadius);
        FROM_TREE(useMeanPosition);
        FROM_TREE(maxPruneLength);
        FROM_TREE(applyFaceMerge);
        FROM_TREE(maxFaceComplexity);
        FROM_TREE(faceCollapseAreaFraction);
        FROM_TREE(edgeMergeUntilConverge);
        FROM_TREE(edgeSmoothUntilConvergence);
        FROM_TREE(maxTurningAngle);
        FROM_TREE(maxPruneLengthPostFM);
        FROM_TREE(maxTurningAnglePostFM);
#undef FROM_TREE
    }

    void MapSimplificationPipeline::Settings::saveToTree(boost::property_tree::ptree& tree)
    {
#define TO_TREE(name) tree.put(#name, m_##name)
        TO_TREE(radius);
        TO_TREE(postRadius);
        TO_TREE(useMeanPosition);
        TO_TREE(maxPruneLength);
        TO_TREE(applyFaceMerge);
        TO_TREE(maxFaceComplexity);
        TO_TREE(faceCollapseAreaFraction);
        TO_TREE(edgeMergeUntilConverge);
        TO_TREE(edgeSmoothUntilConvergence);
        TO_TREE(maxTurningAngle);
        TO_TREE(maxPruneLengthPostFM);
        TO_TREE(maxTurningAnglePostFM);
#undef TO_TREE
    }

    void MapSimplificationPipeline::setSettings(const Settings& settings)
    {
        m_currentSettings = settings;
    }

    const MapSimplificationPipeline::Settings& MapSimplificationPipeline::settings() const
    {
        return m_currentSettings;
    }

    void MapSimplificationPipeline::setPreprocessApproach(
        const VerexOrderApproach& preprocessApproach)
    {
        m_preprocessApproach = preprocessApproach;
    }

    MapSimplificationPipeline::VerexOrderApproach MapSimplificationPipeline
    ::preprocessApproach() const
    {
        return m_preprocessApproach;
    }

    void MapSimplificationPipeline::setRadius(const Models::NT& radius)
    {
        m_currentSettings.m_radius = radius;
    }

    Models::NT MapSimplificationPipeline::radius() const
    {
        return m_currentSettings.m_radius;
    }

    void MapSimplificationPipeline::setPostRadius(const Models::NT& radius)
    {
        m_currentSettings.m_postRadius = radius;
    }

    Models::NT MapSimplificationPipeline::postRadius() const
    {
        return m_currentSettings.m_postRadius;
    }

    void MapSimplificationPipeline::setUseMeanPosition(const bool& useMeanPosition)
    {
        m_currentSettings.m_useMeanPosition = useMeanPosition;
    }

    bool MapSimplificationPipeline::useMeanPosition() const
    {
        return m_currentSettings.m_useMeanPosition;
    }

    void MapSimplificationPipeline::setMaxPruneLength(const Models::NT& maxPruneLength)
    {
        m_currentSettings.m_maxPruneLength = maxPruneLength;
    }

    Models::NT MapSimplificationPipeline::maxPruneLength() const
    {
        return m_currentSettings.m_maxPruneLength;
    }

    void MapSimplificationPipeline::setMaxPruneLengthPostFM(const Models::NT& maxPruneLength)
    {
        m_currentSettings.m_maxPruneLengthPostFM = maxPruneLength;
    }

    Models::NT MapSimplificationPipeline::maxPruneLengthPostFM() const
    {
        return m_currentSettings.m_maxPruneLengthPostFM;
    }

    void MapSimplificationPipeline::setApplyFaceMerge(const bool& applyFaceMerge)
    {
        m_currentSettings.m_applyFaceMerge = applyFaceMerge;
    }

    bool MapSimplificationPipeline::applyFaceMerge() const
    {
        return m_currentSettings.m_applyFaceMerge;
    }

    void MapSimplificationPipeline::setMaxFaceComplexity(const int& maxFaceComplexity)
    {
        m_currentSettings.m_maxFaceComplexity = maxFaceComplexity;
    }

    int MapSimplificationPipeline::maxFaceComplexity() const
    {
        return m_currentSettings.m_maxFaceComplexity;
    }

    void MapSimplificationPipeline::setFaceCollapseAreaFraction(
        const Models::NT& faceCollapseAreaFraction)
    {
        m_currentSettings.m_faceCollapseAreaFraction = faceCollapseAreaFraction;
    }

    Models::NT MapSimplificationPipeline::faceCollapseAreaFraction() const
    {
        return m_currentSettings.m_faceCollapseAreaFraction;
    }

    void MapSimplificationPipeline::setEdgeMergeUntilConverge(const bool& edgeMergeUntilConverge)
    {
        m_currentSettings.m_edgeMergeUntilConverge = edgeMergeUntilConverge;
    }

    bool MapSimplificationPipeline::edgeMergeUntilConverge() const
    {
        return m_currentSettings.m_edgeMergeUntilConverge;
    }

    void MapSimplificationPipeline::setEdgeSmoothUntilConvergence(
        const bool& edgeSmoothUntilConvergence)
    {
        m_currentSettings.m_edgeSmoothUntilConvergence = edgeSmoothUntilConvergence;
    }

    bool MapSimplificationPipeline::edgeSmoothUntilConvergence() const
    {
        return m_currentSettings.m_edgeSmoothUntilConvergence;
    }

    void MapSimplificationPipeline::setMaxTurningAngle(const Models::NT& maxTurningAngle)
    {
        m_currentSettings.m_maxTurningAngle = maxTurningAngle;
    }

    Models::NT MapSimplificationPipeline::maxTurningAngle() const
    {
        return m_currentSettings.m_maxTurningAngle;
    }

    void MapSimplificationPipeline::setMaxTurningAnglePostFM(const Models::NT& maxTurningAngle)
    {
        m_currentSettings.m_maxTurningAnglePostFM = maxTurningAngle;
    }

    Models::NT MapSimplificationPipeline::maxTurningAnglePostFM() const
    {
        return m_currentSettings.m_maxTurningAnglePostFM;
    }

    void MapSimplificationPipeline::apply(const Models::UndirectedEmbeddedGraph& inputGraph)
    {
        std::array<Point, 2> bboxPoints;
        GCpp::DS::computeBoundingBox(inputGraph, bboxPoints);
        m_bboxWidth = bboxPoints[1].x() - bboxPoints[0].x();
        m_bboxHeight = bboxPoints[1].y() - bboxPoints[0].y();

        namespace ms = ::SchematLib::MapSimplification;
        // Reset all stages.
        m_stages->reset();
        m_currentSettings.m_currentIteration = 0;

        SimpObs simplificationObserver;
        // Prepare observer
        // Copy the input graph
        Models::UndirectedEmbeddedGraph target = inputGraph;
        target.m_use_canonical_edges_ids = true; //Don't return canonical edge IDs during simplification.
        //
        // Setup the observer with all knowledge of the initial graph.
        {
            std::set<std::size_t> currentVertexIds;
            GCpp::DS::computeVertexIdSet(target, currentVertexIds);
            std::map<std::size_t, std::pair<std::size_t, std::size_t>> edgeMapping;
            GCpp::DS::getEdgeIdToVerticesMapping(target, edgeMapping);
            // Mark special IDs usage.
            simplificationObserver.m_useSpecialEdgeIds = true;
            simplificationObserver.initialize(currentVertexIds, edgeMapping);
        }


        // Step 0) original
        m_stages->assignNextStage(stageName("Original"), target, simplificationObserver);

        while (m_currentSettings.m_currentIteration < m_currentSettings.m_totalIterations || m_currentSettings.m_untillConvergence)
        {
            const std::size_t initialEventCount = simplificationObserver.m_eventCount;

            // 1) Preprocess (remove tendrils)
            applyRemoveTendrils(target, simplificationObserver);
            // Preprocessed
            m_stages->assignNextStage(stageName("Preprocessed"), target, simplificationObserver);

            applyLargestComponentSelect(target, simplificationObserver);
            // Largest connected component selected
            m_stages->assignNextStage(stageName("LCC"), target, simplificationObserver);

            // 1.5) Degree 2 smoothing
            applyDegree2Smoothing(target, simplificationObserver);
            m_stages->assignNextStage(stageName("Degree2 smoothed"), target, simplificationObserver);

            // 2) Vertex collapsing
            applyVertexCollapse(target, simplificationObserver);
            m_stages->assignNextStage(stageName("Vertex collapsed"), target, simplificationObserver);

            // 3) Vertex to edge collapsing
            applyVertexToEdgeMerge(target, simplificationObserver);
            m_stages->assignNextStage(stageName("Vertex->edge merged"), target, simplificationObserver);

            // 4) Face merging
            //Variants

            /*applyFaceMerge(target, simplificationObserver);
            m_stages->assignNextStage(stageName("Face merged"), target, simplificationObserver);*/

            /*applyContinuationFaceMerge(target, simplificationObserver);
            m_stages->assignNextStage(stageName("Face merged-continuation"), target, simplificationObserver);*/

            applyBboxFaceMerge(target, simplificationObserver);
            m_stages->assignNextStage(stageName("Face merged-bbox"), target, simplificationObserver);

            // 4.1) Remove tendrils afterwards
            //applyRemoveTendrils(target, simplificationObserver);
            applyRemoveTendrilsPostFM(target, simplificationObserver);
            m_stages->assignNextStage(stageName("Tendrils removed after FM"), target, simplificationObserver);

            // 4.2) Collapse vertices afterwards
            //applyVertexCollapse(target, simplificationObserver);
            applyVertexCollapsePostFM(target, simplificationObserver);
            m_stages->assignNextStage(stageName("Vertex collapsed after FM"), target, simplificationObserver);

            applyDegree2SmoothingPostFM(target, simplificationObserver);
            m_stages->assignNextStage(stageName("Degree 2 smoothing after FM"), target, simplificationObserver);

            applyFaceTendrilRemove(target, simplificationObserver);
            m_stages->assignNextStage(stageName("Interior tendril remove after FM"), target, simplificationObserver);

            // Update iteration
            ++m_currentSettings.m_currentIteration;
            // No new events: converged
            if (m_currentSettings.m_untillConvergence && initialEventCount == simplificationObserver.m_eventCount)
            {
                break;
            }
        }
    }
}
