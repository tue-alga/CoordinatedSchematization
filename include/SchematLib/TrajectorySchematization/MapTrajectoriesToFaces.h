#ifndef SCHEMATLIB_TRAJECTORYSCHEMATIZAATION_H
#define SCHEMATLIB_TRAJECTORYSCHEMATIZAATION_H
#include "OffroadTraits.h"
#include "SchematLib/MapSimplification/FaceMergeTraits.h"
#include "SchematLib/MapSimplification/PathMapper.h"

namespace SchematLib::TrajectorySchematization
{
    template<typename EmbeddedGraph>
    class MapTrajectoriesToFaces
    {
    public:
        // By IDs.
        using FaceToTrajectoriesMap = std::map<std::size_t, std::set<std::size_t>>;
        using Traits = MapSimplification::FaceMergeTraits;
        using OffroadTraits = OffroadTraits<EmbeddedGraph>;
        using Trajectory = typename OffroadTraits::OffroadTrajectory; // List of edge IDs.

        struct FaceTrajectories
        {
            using VertexId = std::size_t;
            using TrajectoryId = std::size_t;
            std::vector<std::pair<VertexId, VertexId>> through;
            std::vector<VertexId> enter;
            std::vector<VertexId> exit;
            std::vector<TrajectoryId> inside;
        };

    private:
        struct ArrangementMapping
        {
            std::map<std::size_t, Traits::HE_handle> edgeIdToEdge;
            std::map<std::size_t, Traits::Vertex_handle> vertexIdToVertex;

            void compute(const Traits::Arr& arrangement)
            {
                for(auto he : arrangement.halfedge_handles())
                {
                    if(he->data().id.has_value())
                    {
                        edgeIdToEdge[he->data().id.value()] = he;
                    }
                }
                for(auto v: arrangement.vertex_handles())
                {
                    vertexIdToVertex[v->data().id] = v;
                }
            }
        };

        template<typename Observer>
        void handleTrajectory(const EmbeddedGraph& sourceGraph, const CGAL::Arr_trapezoid_ric_point_location<Traits::Arr>& pointLocation,
            const Trajectory& trajectory, std::size_t trajectoryNum, const Observer& observer,
            std::map<std::size_t, FaceTrajectories>& toFaceMapping)
        {
            if (trajectory.empty()) return;

            std::map<std::size_t, typename EmbeddedGraph::edge_descriptor> edgeIdtoEdge;
            GCpp::DS::computeEdgeIdToEdgeMap(sourceGraph, edgeIdtoEdge);

            // Find face that contains most of the trajectory.
            std::vector<Traits::Point_2> points;
            Traits traits;
            for(std::size_t i = 0; i < trajectory.size()+1; ++i)
            {
                typename EmbeddedGraph::vertex_descriptor v;
                if(i == 0)
                {
                    v = boost::source(edgeIdtoEdge[trajectory.at(0)], sourceGraph);
                }
                else
                {
                    v= boost::target(edgeIdtoEdge[trajectory.at(i-1)], sourceGraph);
                }
                points.push_back(traits.gcppToCgalPosition(GCpp::DS::get_vertex_location(v,sourceGraph)));
            }
            // Determine most hit face
            std::map<std::size_t, std::size_t> countPerFace;
            std::size_t maxCountFid = 0;
            std::size_t maxVal = 0;
            for(const auto& p : points)
            {
                auto result= pointLocation.locate(p);
                const Traits::Arr::Face_const_handle* f;
                if ((f = boost::get<Traits::Arr::Face_const_handle>(&result)))
                {
                    if ((*f)->is_unbounded()) continue;
                    if((*f)->data().id.has_value())
                    {
                        const auto fId = (*f)->data().id.value();
                        if (countPerFace.find(fId) == countPerFace.end()) countPerFace[fId] = 0;
                        countPerFace[fId]++;
                        if(countPerFace[fId] > maxVal)
                        {
                            maxVal = countPerFace[fId];
                            maxCountFid = fId;
                        }
                    }
                }
            }
            // TODO: may mean assigned to outer face, or we somehow only hit vertices.
            if (maxVal == 0) return;

            MapSimplification::PathMapper mapper(observer);
            // Determine type
            const auto& data = trajectory.trajectoryData();
            if(data.startOnBoundary)
            {
                if(data.endOnBoundary)
                {
                    // Map vertices
                    auto startResult = mapper.decodeVertex(data.startVId);
                    if (std::holds_alternative<MapSimplification::EdgeTrackingObserver::Deleted>(startResult)) throw std::runtime_error("Cannot delete vertices!");
                    auto startV = std::get < MapSimplification::EdgeTrackingObserver::Vertex>(startResult).id;
                    auto endV = std::get < MapSimplification::EdgeTrackingObserver::Vertex>(mapper.decodeVertex(data.endVId)).id;

                    toFaceMapping[maxCountFid].through.push_back(std::make_pair(startV, endV));
                }
                else
                {
                    auto vert = std::get < MapSimplification::EdgeTrackingObserver::Vertex>(mapper.decodeVertex(data.startVId)).id;
                    toFaceMapping[maxCountFid].enter.push_back(vert);
                }
            }
            else
            {
                if(data.endOnBoundary)
                {
                    auto vert = std::get < MapSimplification::EdgeTrackingObserver::Vertex>(mapper.decodeVertex(data.endVId)).id;
                    toFaceMapping[maxCountFid].exit.push_back(vert);
                }
                else
                {
                    toFaceMapping[maxCountFid].inside.push_back(trajectoryNum);
                }
            }
        }
    public:
        template<typename Observer>
        void operator()(const MapSimplification::FaceMergeTraits::Arr& arrangement,
                const EmbeddedGraph& sourceGraph,
                const std::vector<typename OffroadTraits::OffroadTrajectory>& offroadTrajectories,
                const Observer& observer,
                        std::map<std::size_t, FaceTrajectories>& toFaceMapping)
        {
            // Map trajectories to a face by majority vote for now.
            // Find out their original enter/exit vertices or through vertices.

            CGAL::Arr_trapezoid_ric_point_location<Traits::Arr> pointLocation(arrangement);

            // Compute mapping for edges
            std::map<std::size_t, Traits::HE_handle> edgeIdToEdge;
            for(auto he : arrangement.halfedge_handles())
            {
                if (!he->data().id.has_value()) throw std::runtime_error("Missing ID on edge");
                edgeIdToEdge[he->data().id.value()] = he;
            }
            // Potentially compute vertices for trajectory to determine directionality.
            // TODO
            std::size_t trajectoryNum = 0;
            for(const auto& traj: offroadTrajectories)
            {
                handleTrajectory(sourceGraph, pointLocation, traj, trajectoryNum, observer, toFaceMapping);
                ++trajectoryNum;
            }
        }
    };
}
#endif
