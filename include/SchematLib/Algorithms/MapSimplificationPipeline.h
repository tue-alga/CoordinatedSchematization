#ifndef SCHEMATLIB_ALGORITHMS_MAPSIMPLIFICATIONPIPELINE_H
#define SCHEMATLIB_ALGORITHMS_MAPSIMPLIFICATIONPIPELINE_H
#include <boost/property_tree/ptree_fwd.hpp>

#include "SchematLib/Models/SchematizationStages.h"
#include "SchematLib/MapSimplification/Degree2VertexPreprocess.h"
#include "SchematLib/MapSimplification/FaceMerger.h"
#include "SchematLib/MapSimplification/LargestConnectedComponent.h"
#include "SchematLib/MapSimplification/RemoveTendrils.h"
#include "SchematLib/MapSimplification/VertexEdgeMerger.h"
#include "SchematLib/MapSimplification/VertexMerger.h"
#include <GCpp/DS/BoostEmbeddedGraph.h>

namespace SchematLib::Algorithms
{
    class MapSimplificationPipeline
    {
        std::shared_ptr<Models::SchematizationStages> m_stages;
    public:
        enum class VerexOrderApproach
        {
            DenseVertexOrder,
            VertexOrder
        };
    private:
        // Observer to use for changes
        using SimpObs = SchematLib::MapSimplification::EdgeTrackingObserver;
        using Graph = Models::UndirectedEmbeddedGraph;

        // Algorithms
        using VertexEdgeMerger = MapSimplification::VertexEdgeMerger<Models::UndirectedEmbeddedGraph, Models::UndirectedGraphEdgeIndex, SimpObs>;
        using Point = Models::Point;

        VerexOrderApproach m_preprocessApproach = VerexOrderApproach::DenseVertexOrder;
        // Parameters
        // Simplification parameters
        /*Models::NT m_radius = 0.01;
        Models::NT m_postRadius = 0.01;
        bool m_useMeanPosition = true;
        Models::NT m_maxPruneLength = 100.;
        bool m_applyFaceMerge = true;
        int m_maxFaceComplexity = 5;
        Models::NT m_faceCollapseAreaFraction = 0.01;
        bool m_edgeMergeUntilConverge = true;
        bool m_edgeSmoothUntilConvergence = true;
        Models::NT m_maxTurningAngle = 30;

        Models::NT m_maxPruneLengthPostFM = 100.;*/

        // Bounding box
        Models::NT m_bboxWidth = 0, m_bboxHeight = 0;

        Models::NT maxBboxAxis() const;

        // Allow multiple iterations?
        /*int m_currentIteration = 0;
        int m_totalIterations = 1;
        bool m_untillConvergence = false;*/

        //Models::NT m_maxTurningAnglePostFM = 30;

        std::string stageName(const std::string& name) const;

        void applyRemoveTendrils(Graph& target, SimpObs& simplificationObserver);

        void applyDegree2Smoothing(Graph& target, SimpObs& simplificationObserver);

        void applyDegree2SmoothingPostFM(Graph& target, SimpObs& simplificationObserver);

        void applyVertexCollapse(Graph& target, SimpObs& simplificationObserver);

        void applyVertexToEdgeMerge(Graph& target, SimpObs& simplificationObserver);

        void applyFaceMerge(Graph& target, SimpObs& simplificationObserver);

        void applyContinuationFaceMerge(Graph& target, SimpObs& simplificationObserver);

        void applyBboxFaceMerge(Graph& target, SimpObs& simplificationObserver);

        void applyLargestComponentSelect(Graph& target, SimpObs& simplificationObserver);

        void applyVertexCollapsePostFM(Graph& target, SimpObs& simplificationObserver);
        void applyRemoveTendrilsPostFM(Graph& target, SimpObs& simplificationObserver);

        void applyFaceTendrilRemove(Graph& target, SimpObs& simplificationObserver);
    public:
        struct Settings
        {
            // Parameters
            // Simplification parameters
            Models::NT m_radius = 0.01;
            Models::NT m_postRadius = 0.01;
            bool m_useMeanPosition = true;
            Models::NT m_maxPruneLength = 100.;
            bool m_applyFaceMerge = true;
            int m_maxFaceComplexity = 5;
            Models::NT m_faceCollapseAreaFraction = 0.01;
            bool m_edgeMergeUntilConverge = true;
            bool m_edgeSmoothUntilConvergence = true;
            Models::NT m_maxTurningAngle = 30;

            Models::NT m_maxPruneLengthPostFM = 100.;

            // Allow multiple iterations?
            int m_currentIteration = 0;
            int m_totalIterations = 1;
            bool m_untillConvergence = false;

            Models::NT m_maxTurningAnglePostFM = 30;
            void loadFromTree(const boost::property_tree::ptree& tree);
            void saveToTree(boost::property_tree::ptree& tree);
        };
    private:
        Settings m_currentSettings;
    public:
        void setSettings(const Settings& settings);
        const Settings& settings() const;
#pragma region Parameter setters/getters
        void setPreprocessApproach(const VerexOrderApproach& preprocessApproach);

        VerexOrderApproach preprocessApproach() const;

        void setRadius(const Models::NT& radius);

        Models::NT radius() const;

        void setPostRadius(const Models::NT& radius);

        Models::NT postRadius() const;

        void setUseMeanPosition(const bool& useMeanPosition);

        bool useMeanPosition() const;

        void setMaxPruneLength(const Models::NT& maxPruneLength);

        Models::NT maxPruneLength() const;

        void setMaxPruneLengthPostFM(const Models::NT& maxPruneLength);

        Models::NT maxPruneLengthPostFM() const;

        void setApplyFaceMerge(const bool& applyFaceMerge);

        bool applyFaceMerge() const;

        void setMaxFaceComplexity(const int& maxFaceComplexity);

        int maxFaceComplexity() const;

        void setFaceCollapseAreaFraction(const Models::NT& faceCollapseAreaFraction);

        Models::NT faceCollapseAreaFraction() const;

        void setEdgeMergeUntilConverge(const bool& edgeMergeUntilConverge);

        bool edgeMergeUntilConverge() const;

        void setEdgeSmoothUntilConvergence(const bool& edgeSmoothUntilConvergence);

        bool edgeSmoothUntilConvergence() const;

        void setMaxTurningAngle(const Models::NT& maxTurningAngle);

        Models::NT maxTurningAngle() const;

        void setMaxTurningAnglePostFM(const Models::NT& maxTurningAngle);

        Models::NT maxTurningAnglePostFM() const;
#pragma endregion

        MapSimplificationPipeline(std::shared_ptr<Models::SchematizationStages> stages) :m_stages(std::move(stages)) {}

        void apply(const Models::UndirectedEmbeddedGraph& inputGraph);
    };
}
#endif
