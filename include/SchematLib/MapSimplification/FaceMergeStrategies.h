#ifndef SCHEMATLIB_MAPSIMPLIFICATION_FACEMERGESTRATEGIES_H
#define SCHEMATLIB_MAPSIMPLIFICATION_FACEMERGESTRATEGIES_H
#include "FaceMergeTraits.h"

namespace SchematLib::MapSimplification
{
    namespace detail
    {
        template<typename CircularEgeIterator>
        struct ComputeFaceComplexity
        {
            std::size_t operator()(CircularEgeIterator edgeIterator)
            {
                auto it = edgeIterator;
                ++it;
                std::size_t complexity = 1;
                while (it != edgeIterator)
                {
                    ++it;
                    ++complexity;
                }
                return complexity;
            }
        };
    }
    namespace Predicates
    {
        // Functor to check if a vertex belongs to a path going into the interior of the given face.
        struct VertexBelongsToFaceInterior
        {
            using Traits = FaceMergeTraits;
            bool operator()(FaceMergeTraits::Face_handle face, Traits::Vertex_handle vert) const
            {
                auto it = vert->incident_halfedges();
                do
                {
                    if (it->face() != face || it->twin()->face() != face)
                    {
                        return false;
                    }
                    ++it;
                } while (it != vert->incident_halfedges());
                return true;
            }
        };
    }
    namespace Operators
    {
        // Removes faces that are nested in the given face (aka 'holes').
        struct RemoveHolesInFace
        {
            using Traits = FaceMergeTraits;

            void operator()(Traits::Arr& arr, typename Traits::Face_handle face, std::set<std::size_t>& affectedVertices, std::set<std::size_t>& affectedEdges) const
            {
                using FaceHandle = typename Traits::Face_handle;
                using HeHandle = typename Traits::HE_handle;
                std::map<std::size_t, typename Traits::Arr::Ccb_halfedge_circulator> halfedgesToDelete;
                typename Traits::Arr::Hole_iterator holeIt;
                for (holeIt = face->holes_begin(); holeIt != face->holes_end(); ++holeIt)
                {
                    std::map<std::size_t, FaceHandle> interiorFaces;
                    typename Traits::Arr::Ccb_halfedge_circulator eIt = *holeIt;
                    do
                    {
                        auto heHandle = eIt;
                        affectedVertices.insert(heHandle->source()->data().id);
                        affectedVertices.insert(heHandle->target()->data().id);
                        halfedgesToDelete[heHandle->data().id.value()] = eIt;
                        if (heHandle->face() != heHandle->twin()->face())
                        {
                            auto subIt = heHandle->twin();
                            do
                            {
                                affectedVertices.insert(subIt->source()->data().id);
                                affectedVertices.insert(subIt->target()->data().id);
                                halfedgesToDelete[subIt->data().id.value()] = subIt;
                            } while (subIt != heHandle->twin());
                        }
                        ++eIt;
                    } while (eIt != *holeIt);
                }
                for (const auto& kv : halfedgesToDelete)
                {
                    arr.remove_edge(kv.second);
                }
            }
        };

        struct RemoveInteriorFaceTendrils
        {
            using Traits = FaceMergeTraits;
            RemoveHolesInFace removeHoles;
            Predicates::VertexBelongsToFaceInterior vertBelongsToInterior;

            void operator()(Traits::Arr& arr, Traits::Face_handle face, std::set<std::size_t>& affectedVertices, std::set<std::size_t>& affectedEdges) const
            {
                removeHoles(arr, face, affectedVertices, affectedEdges);

                std::map<std::size_t, Traits::HE_handle> halfedgesToDelete;

                auto handle = face->outer_ccb();
                do
                {
                    if (handle->face() == face && handle->twin()->face() == face)
                    {
                        if (halfedgesToDelete.find(handle->twin()->data().id.value()) != halfedgesToDelete.end())
                        {
                            ++handle;
                            continue;
                        }
                        affectedEdges.insert(handle->data().id.value());
                        halfedgesToDelete[handle->data().id.value()] = handle;
                        if (vertBelongsToInterior(face, handle->source()))
                        {
                            affectedVertices.insert(handle->source()->data().id);
                        }
                        if (vertBelongsToInterior(face, handle->target()))
                        {
                            affectedVertices.insert(handle->target()->data().id);
                        }
                    }
                    ++handle;
                } while (handle != face->outer_ccb());
                if (halfedgesToDelete.size() > 0) std::cout << "[RemoveInteriorTendrils] Removing " << halfedgesToDelete.size() << " halfedges\n";
                for (const auto& kv : halfedgesToDelete)
                {
                    arr.remove_edge(kv.second);
                }
            }
        };

        template<typename Traits=FaceMergeTraits>
        struct ComputeFaceCenter
        {
            typename Traits::Point_2 operator()(typename Traits::Face_handle face)
            {
                std::size_t complexity = 0;
                auto edgeIt = face->outer_ccb();
                // Compute center of face.
                auto center = edgeIt->source()->point() - CGAL::ORIGIN;
                ++edgeIt;
                ++complexity;
                while (edgeIt != face->outer_ccb())
                {
                    auto prevIt = edgeIt;
                    center += edgeIt->source()->point() - CGAL::ORIGIN;
                    ++edgeIt;
                    ++complexity;
                    if (prevIt->target() != edgeIt->source()) throw std::runtime_error("Unexpected direction!");
                }
                center /= (double)complexity;
                typename Traits::Point_2 preOutput = CGAL::ORIGIN + center;
                typename Traits::Point_2 output;

                // Force into not too terrible resolution
                Traits::reducePrecision(preOutput, output);
                return output;
            }
        };

        template<typename Traits=FaceMergeTraits>
        struct ComputeMeanPosition
        {
            template<typename IterableWithCgalPoint, typename IteratorToPoint>
            typename Traits::Point_2 operator()(IterableWithCgalPoint begin, IterableWithCgalPoint end, IteratorToPoint&& func, const std::vector<typename Traits::Point_2>& extraPoints={})
            {
                // Compute center of face.
                auto center = func(begin) - CGAL::ORIGIN;
                auto curr = std::next(begin);
                std::size_t complexity = 1;
                for(; curr != end; ++curr)
                {
                    center += func(curr) - CGAL::ORIGIN;
                    ++complexity;
                }
                for(const auto& pnt: extraPoints)
                {
                    center += pnt - CGAL::ORIGIN;
                }
                complexity += extraPoints.size();
                center /= typename Traits::CGAL_NT((double)complexity);
                typename Traits::Point_2 preOutput = CGAL::ORIGIN + center;
                typename Traits::Point_2 output;

                // Force into not too terrible resolution
                Traits::reducePrecision(preOutput, output);
                return output;
            }
        };

        /**
         * \brief Splits a (simple) face into edges between degree >2 vertices, and subpaths between degree >2 vertices.
         * \tparam Traits Trait class specifying handle types
         */
        template<typename Traits>
        struct ComputeBoundaryTypes
        {
            void operator()(typename Traits::Face_handle face, std::vector<typename Traits::Vertex_handle>& verticesToKeep,
                std::vector<std::vector<typename Traits::Halfedge_handle>>& degree2Subpaths,
                std::vector<typename Traits::Halfedge_handle>& simpleRouteEdges
            )
            {
                // Find degree-2 subpaths on the face boundary.
                // Collect vertices
                using HE = typename Traits::Halfedge_handle;
                auto heHandle = face->outer_ccb();
                while (heHandle->source()->degree() <= 2)
                {
                    --heHandle;
                    // TODO: check equality to outer_ccb() to avoid potential infinite loop.
                    if (heHandle == face->outer_ccb()) throw std::runtime_error("Got a degree 2 only face!");

                    // Todo: delete and assign to face.
                }
                auto start = heHandle;
                do
                {
                    const auto srcDegree = heHandle->source()->degree();
                    const auto targetDegree = heHandle->target()->degree();
                    if (srcDegree > 2)
                    {
                        std::cout << "[FH] Marked vertex for edge\n";
                        verticesToKeep.push_back(heHandle->source());
                    }
                    if (srcDegree > 2 && targetDegree > 2)
                    {
                        std::cout << "[FH] Simple edge" << heHandle->data() << ", vs=" << heHandle->source()->data().id << ',' << heHandle->target()->data().id << '\n';
                        simpleRouteEdges.push_back(heHandle);
                    }
                    else
                    {
                        std::cout << "[FH] Subpath edge" << heHandle->data() << ", vs=" << heHandle->source()->data().id << ',' << heHandle->target()->data().id << '\n';
                        if (degree2Subpaths.empty()) degree2Subpaths.push_back(std::vector<HE>{});
                        // New path
                        // Capture circularity.
                        if (!degree2Subpaths.front().empty() && degree2Subpaths.front().front()->source() == heHandle->target() && heHandle->target()->degree() <= 2)
                        {
                            degree2Subpaths.front().insert(degree2Subpaths.front().begin(), heHandle);
                        }
                        // 
                        else
                        {
                            if (!degree2Subpaths.back().empty() &&
                                (degree2Subpaths.back().back()->target() != heHandle->source() || heHandle->source()->degree() > 2)) {
                                degree2Subpaths.push_back(std::vector<HE>{});
                            }
                            degree2Subpaths.back().push_back(heHandle);
                        }
                    }
                    ++heHandle;
                } while (heHandle != start);
            }

            void operator()(typename Traits::Face_handle face, 
                std::vector<std::vector<typename Traits::Halfedge_handle>>& ccwSubpaths
                )
            {
                using Subpath = std::vector<typename Traits::Halfedge_handle>;
                // Find degree-2 subpaths on the face boundary.
                // Collect vertices
                using HE = typename Traits::Halfedge_handle;
                auto heHandle = face->outer_ccb();
                while (heHandle->source()->degree() <= 2)
                {
                    --heHandle;
                    // TODO: check equality to outer_ccb() to avoid potential infinite loop.
                    if (heHandle == face->outer_ccb()) throw std::runtime_error("Got a degree 2 only face!");

                    // Todo: delete and assign to face.
                }
                auto start = heHandle;
                do
                {
                    const auto srcDegree = heHandle->source()->degree();
                    const auto targetDegree = heHandle->target()->degree();
                    if (srcDegree > 2 && targetDegree > 2)
                    {
                        ccwSubpaths.push_back(Subpath{ heHandle });
                    }
                    else
                    {
                        if (ccwSubpaths.empty()) ccwSubpaths.push_back(std::vector<HE>{});
                        // New path
                        if (!ccwSubpaths.back().empty() && srcDegree > 2) {
                            ccwSubpaths.push_back(std::vector<HE>{});
                        }
                        ccwSubpaths.back().push_back(heHandle);
                    }
                    ++heHandle;
                } while (heHandle != start);
            }
        };

        template<typename Traits=FaceMergeTraits>
        struct ComputeTurningAngle
        {
            using VH = typename Traits::Vertex_handle;
            typename Traits::NT operator()(VH v0, VH v1, VH v2) const
            {
                using NT_t = typename Traits::NT;
                auto diff0 = Traits::cgalPositionToGcpp(v1->point()) - Traits::cgalPositionToGcpp(v0->point());
                auto diff1 = Traits::cgalPositionToGcpp(v2->point()) - Traits::cgalPositionToGcpp(v1->point());
                return std::acos(std::clamp<NT_t>(diff0.dot(diff1) / diff0.length() / diff1.length(),0,1));
            }
        };
    }

    namespace FaceMergeApproach
    {
        class DeleteSmallFacesToNeighbourViaContraction
        {
            using Traits = FaceMergeTraits;

            // The factor above which we deem the effect minor when removing this face.
            double m_faceSizeFactor = 2;

            /**
             * \brief Computes face complexity
             * \tparam FaceIterator The face iterator type
             * \param it The face iterator
             * \return The complexity
             */
            static std::size_t computeComplexity(FaceMergeTraits::Face_handle it)
            {
                detail::ComputeFaceComplexity<decltype(it->outer_ccb())> computer;
                return computer(it->outer_ccb());
            }

            std::optional<Traits::HE_handle> findEdge(Traits::Vertex_handle source, Traits::Vertex_handle target)
            {
                auto circulator = target->incident_halfedges();
                auto it = circulator;
                do
                {
                    if (it->source() == source) return it;
                    ++it;
                } while (it != circulator);
                return {};
            }
            static std::size_t getEdgeId(Traits::HE_handle handle)
            {
                return handle->data().id.value();
            }
            bool vertBelongsToInterior(FaceMergeTraits::Face_handle face, Traits::Vertex_handle vert)
            {
                auto it = vert->incident_halfedges();
                do
                {
                    if (it->face() != face || it->twin()->face() != face)
                    {
                        return false;
                    }
                    ++it;
                } while (it != vert->incident_halfedges());
                return true;
            }

            void removeHoles(Traits::Arr& arr, typename Traits::Face_handle face, std::set<std::size_t>& affectedVertices, std::set<std::size_t>& affectedEdges)
            {
                using FaceHandle = typename Traits::Face_handle;
                using HeHandle = typename Traits::HE_handle;
                std::map<std::size_t, typename Traits::Arr::Ccb_halfedge_circulator> halfedgesToDelete;
                typename Traits::Arr::Hole_iterator holeIt;
                for (holeIt = face->holes_begin(); holeIt != face->holes_end(); ++holeIt)
                {
                    std::map<std::size_t, FaceHandle> interiorFaces;
                    typename Traits::Arr::Ccb_halfedge_circulator eIt = *holeIt;
                    do
                    {
                        auto heHandle = eIt;
                        affectedVertices.insert(heHandle->source()->data().id);
                        affectedVertices.insert(heHandle->target()->data().id);
                        halfedgesToDelete[heHandle->data().id.value()] = eIt;
                        if (heHandle->face() != heHandle->twin()->face())
                        {
                            auto subIt = heHandle->twin();
                            do
                            {
                                affectedVertices.insert(subIt->source()->data().id);
                                affectedVertices.insert(subIt->target()->data().id);
                                halfedgesToDelete[subIt->data().id.value()] = subIt;
                            } while (subIt != heHandle->twin());
                        }
                        ++eIt;
                    } while (eIt != *holeIt);
                }
                for (const auto& kv : halfedgesToDelete)
                {
                    arr.remove_edge(kv.second);
                }
            }

            void removeInteriorTendrils(Traits::Arr& arr, Traits::Face_handle face, std::set<std::size_t>& affectedVertices, std::set<std::size_t>& affectedEdges)
            {
                removeHoles(arr, face, affectedVertices, affectedEdges);

                std::map<std::size_t, Traits::HE_handle> halfedgesToDelete;

                auto handle = face->outer_ccb();
                do
                {
                    if (handle->face() == face && handle->twin()->face() == face)
                    {
                        // Only one of the two halfedges needs to be deleted.
                        if (halfedgesToDelete.find(handle->twin()->data().id.value()) != halfedgesToDelete.end())
                        {
                            ++handle;
                            continue;
                        }
                        affectedEdges.insert(handle->data().id.value());
                        halfedgesToDelete[handle->data().id.value()] = handle;
                        if (vertBelongsToInterior(face, handle->source()))
                        {
                            affectedVertices.insert(handle->source()->data().id);
                        }
                        if (vertBelongsToInterior(face, handle->target()))
                        {
                            affectedVertices.insert(handle->target()->data().id);
                        }
                    }
                    ++handle;
                } while (handle != face->outer_ccb());
                for (const auto& kv : halfedgesToDelete)
                {
                    arr.remove_edge(kv.second);
                }
            }
            static bool replacementIsFeasible(Traits::Arr& arrangement, Traits::Face_handle face, const Traits::Arr::Point_2& center)
            {
                using Segment = typename Traits::Arr::X_monotone_curve_2;
                Traits::Make_segment_2 makeSeg;

                auto ccb = face->outer_ccb();

                const auto holeCount = std::distance(face->holes_begin(), face->holes_end());
                std::cout << "[FM] Number of holes " << holeCount << ", total complexity:" << computeComplexity(face) << "\n";
                std::cout << "[FM] Arrangement complexity: |E|=" << arrangement.number_of_edges() << ", |V|:" << arrangement.number_of_vertices() << ", |F|=" << arrangement.number_of_faces() << "\n";
                std::size_t checkCount = 0;
                bool involvedUnbounded = false;
                // Setup lookup structure (point location)
                std::cout << "[FM] Building pointlocation\n";
                CGAL::Arr_trapezoid_ric_point_location<Traits::Arr> pointLocation(arrangement);
                std::cout << "[FM] Starting search\n";
                do
                {
                    if (ccb->twin()->face()->is_unbounded() && !involvedUnbounded)
                    {
                        involvedUnbounded = true;
                        std::cout << "[FM] Unbounded is involved\n";
                    }
                    auto next = ccb;
                    ++next;

                    // Internally interects, then reject for now.
                    Segment seg = makeSeg(center, ccb->target()->point());
                    std::vector< boost::variant<typename Traits::Vertex_handle, typename Traits::HE_handle, typename Traits::Face_handle>> elements;

                    CGAL::zone(arrangement, seg, std::back_inserter(elements), pointLocation);
                    std::cout << "[FM] Zone size: " << elements.size() << "\n";
                    if (elements.size() >= 2) // It should always be atleast one, since we intersect at least one vertex. (the source of the HE).
                    {
                        //TODO: for now only care about edge intersections. This breaks when we just intersect an unknown vertex.
                        bool intersectsEdge = false;
                        for (const auto& el : elements)
                        {
                            if (el.which() != 1) continue;

                            auto he = boost::get<typename Traits::HE_handle>(el);
                            if (he == ccb || he == next)
                            {
                                continue;
                            }
                            if (he->face() == face || he->twin()->face() == face)
                            {
                                continue;
                            }
                            intersectsEdge = true;
                            break;
                        }
                        if (intersectsEdge) {
                            return false;
                        }
                    }
                    ++checkCount;
                    ++ccb;
                } while (ccb != face->outer_ccb());
                std::cout << "[FM] Checked " << checkCount << " segments for feasibility \n";
                return true;
            }
            void findNonInvolvedHalfEdgesPerFace(Traits::Face_handle faceToCollapse, std::map <std::size_t, Traits::HE_handle > & handles)
            {
                auto handle = faceToCollapse->outer_ccb();
                do
                {
                    if (handle->twin()->face() == faceToCollapse)
                    {
                        ++handle;
                        continue;
                    }
                    auto twinFace = handle->twin()->face();
                    if (twinFace->data().id.has_value() && !twinFace->is_unbounded())
                    {
                        const auto fId = twinFace->data().id.value();
                        if (handles.find(fId) != handles.end())
                        {
                            ++handle;
                            continue;
                        }
                        // Find an edge that is not involved in the merge, so we can restore from that.
                        auto subIt = handle->twin()->ccb();
                        do
                        {
                            // Not aligned to the face we are collapsing.
                            if (subIt->twin()->face() != faceToCollapse)
                            {
                                break;
                            }
                            ++subIt;
                        } while (true);
                        handles[fId] = subIt;
                    }
                    ++handle;
                } while (handle != faceToCollapse->outer_ccb());
            }

            void computeFaceBoundaryPathTypes(Traits::Face_handle face, std::vector<Traits::Vertex_handle>& verticesToKeep,
                std::vector<std::vector<Traits::HE_handle>>& degree2Subpaths,
                std::vector<Traits::HE_handle>& simpleRouteEdges
            )
            {
                // Find degree-2 subpaths on the face boundary.
                // Collect vertices
                using HE = Traits::HE_handle;
                auto heHandle = face->outer_ccb();
                while (heHandle->source()->degree() <= 2)
                {
                    --heHandle;
                    // TODO: check equality to outer_ccb() to avoid potential infinite loop.
                    if (heHandle == face->outer_ccb()) throw std::runtime_error("Got a degree 2 only face!");

                    // Todo: delete and assign to face.
                }
                auto start = heHandle;
                do
                {
                    const auto srcDegree = heHandle->source()->degree();
                    const auto targetDegree = heHandle->target()->degree();
                    if (srcDegree > 2)
                    {
                        std::cout << "[FH] Marked vertex for edge\n";
                        verticesToKeep.push_back(heHandle->source());
                    }
                    if (srcDegree > 2 && targetDegree > 2)
                    {
                        std::cout << "[FH] Simple edge" << heHandle->data() << ", vs=" << heHandle->source()->data().id << ',' << heHandle->target()->data().id << '\n';
                        simpleRouteEdges.push_back(heHandle);
                    }
                    else
                    {
                        std::cout << "[FH] Subpath edge" << heHandle->data() << ", vs=" << heHandle->source()->data().id << ',' << heHandle->target()->data().id << '\n';
                        if (degree2Subpaths.empty()) degree2Subpaths.push_back(std::vector<HE>{});
                        // New path
                        // Capture circularity.
                        if (!degree2Subpaths.front().empty() && degree2Subpaths.front().front()->source() == heHandle->target() && heHandle->target()->degree() <= 2)
                        {
                            degree2Subpaths.front().insert(degree2Subpaths.front().begin(), heHandle);
                        }
                        // 
                        else
                        {
                            if (!degree2Subpaths.back().empty() &&
                                (degree2Subpaths.back().back()->target() != heHandle->source() || heHandle->source()->degree() > 2)) {
                                degree2Subpaths.push_back(std::vector<HE>{});
                            }
                            degree2Subpaths.back().push_back(heHandle);
                        }
                    }
                    ++heHandle;
                } while (heHandle != start);
            }

            void computeCenter(Traits::Face_handle face, Traits::Arr::Point_2& output)
            {
                std::size_t complexity = 0;
                auto edgeIt = face->outer_ccb();
                // Compute center of face.
                Traits::Kernel::Vector_2 center = edgeIt->source()->point() - CGAL::ORIGIN;
                ++edgeIt;
                ++complexity;
                while (edgeIt != face->outer_ccb())
                {
                    auto prevIt = edgeIt;
                    center += edgeIt->source()->point() - CGAL::ORIGIN;
                    ++edgeIt;
                    ++complexity;
                    if (prevIt->target() != edgeIt->source()) throw std::runtime_error("Unexpected direction!");
                }
                center /= (double)complexity;
                Traits::Arr::Point_2 preOutput = CGAL::ORIGIN + center;

                // Force into not too terrible resolution
                output = Traits::Arr::Point_2(Traits::CGAL_NT(CGAL::to_double(preOutput[0])), Traits::CGAL_NT(CGAL::to_double(preOutput[1])));
            }
        public:
            // Flag that we are going to try to remap faces. 
            bool recomputeAfterMerge() const
            {
                return false;
            }
            /**
             * \brief
             * \tparam Arrangement
             * \tparam FaceIterator
             * \tparam LookupStructure
             * \tparam Observer
             * \param arrangement The CGAL arrangement to use.
             * \param face
             * \param lookup The lookup structure to use. Should be attached to the arrangement!
             * \param remappedFaces
             * \param observer
             */
            template<typename SelectPredicate, typename Observer = EdgeTrackingObserver>
            void operator()(Traits::Arr& arrangement, Traits::Face_handle face, std::map<std::size_t, Traits::Face_handle>& remappedFaces, Observer& observer, const SelectPredicate& pred)
            {
                using HE = Traits::HE_handle;
                using Segment = typename Traits::Arr::X_monotone_curve_2;
                using FaceHandle = typename Traits::Arr::Face_handle;
                using Vector = decltype(face->outer_ccb()->source()->point() - CGAL::ORIGIN);
                using VH = Traits::Vertex_handle;
                struct SubpathDeletion
                {
                    std::vector<std::size_t> edgeIds;
                    std::vector<std::size_t> collapsingVerts;
                    VH startVert;
                    VH endVert;
                };

                const auto complexity = computeComplexity(face);
                face->data().beingDeleted = true;

                if (complexity == 2)
                {
                    arrangement.remove_edge(face->outer_ccb(), false, false);
                    observer.handleEdgeDelete(getEdgeId(face->outer_ccb())); //Lookup
                    return;
                }

                typename Traits::Arr::Point_2 center;
                computeCenter(face, center);

                // Contract face.
                // Verfy that the replacement is feasible for some constraints (planar).
                std::cout << "[FaceMergeHandler] Checking feasibility\n";
                if (!replacementIsFeasible(arrangement, face, center)) {
                    std::cout << "[FaceMergeHandler] Infeasible, continuing\n";
                    return;
                }

                // Cleanup the interior before inserting the vertex.
                std::set<std::size_t> interiorCleanAffectedVerts, interiorCleanAffectedEdges;
                std::cout << "[FaceMergeHandler] Removing tendrils\n";
                removeInteriorTendrils(arrangement, face, interiorCleanAffectedVerts, interiorCleanAffectedEdges);


                // Insert new vertex first
                // TODO: hope to not split an edge.
                auto newVertexHandle = arrangement.insert_in_face_interior(center, face);
                std::cout << "[FaceMergeHandler] Insert v with " << newVertexHandle->data().id << "\n";

                // Notify of internal removals
                if (interiorCleanAffectedVerts.size() > 0)
                {
                    std::cout << "[FaceMergeHandler] Delete tendril verts: ";
                    for (auto v : interiorCleanAffectedVerts)
                    {
                        observer.handleVertexMerge(v, newVertexHandle->data().id);
                        std::cout << ' ' << v;
                    }
                    std::cout << '\n';
                }
                if (interiorCleanAffectedEdges.size() > 0)
                {
                    std::cout << "[FaceMergeHandler] Delete tendril edges: ";
                    for (auto e : interiorCleanAffectedEdges)
                    {
                        observer.handleEdgeCollapse(e, newVertexHandle->data().id);
                        std::cout << ' ' << e;
                    }
                    std::cout << '\n';
                }

                // Store representative halfedges for faces, to restore the face data afterwards.
                std::map<std::size_t, HE> fIdToHE;
                findNonInvolvedHalfEdgesPerFace(face, fIdToHE);
                std::map<std::size_t, Traits::FaceData> faceDataMap;
                for (const auto& kv : fIdToHE)
                {
                    faceDataMap[kv.first] = kv.second->face()->data();
                }

                // Find degree-2 subpaths on the face boundary.
                // Collect vertices
                std::vector<VH> vertices;
                // Subpaths that route via the new center by 2 edges
                std::vector<std::vector<HE>> degree2Subpaths;
                // Edges that just route via the new center
                std::vector<HE> simpleRouteEdges;

                // Compute face boundary types, where subpaths are oriented in CCW manner in the face.
                computeFaceBoundaryPathTypes(face, vertices, degree2Subpaths, simpleRouteEdges);

                // Setup reroute with terminals
                FaceTreeReroute<Traits::NT> rerouter(vertices.begin(), vertices.end(), [](auto vh) {return vh->data().id; });

                rerouter.add_subgraph_vertex(newVertexHandle->data().id);

                // Setup affected edges
                Traits::for_each_ccb_edge(face, [&rerouter](auto he,auto& cmd)
                {
                    rerouter.add_affected_edge(he->data().id.value(), he->source()->data().id, he->target()->data().id);
                });

                const auto subpathTotal = std::accumulate(degree2Subpaths.begin(), degree2Subpaths.end(), std::size_t{ 0 }, [](auto total, const auto& el)
                {
                    return total + el.size();
                });

                // Special case: degree 2 path collapsing to single edge.
                if (simpleRouteEdges.size() == 1 && degree2Subpaths.size() == 1)
                {
                    std::vector<std::size_t> pathIds;
                    for (auto he : degree2Subpaths.front())
                    {
                        pathIds.push_back(getEdgeId(he));
                    }
                    // Reroute via the single edge, we don't need to add new vertices.
                    // NOTE: boundary routes are CCW, we need to align the reroute properly
                    const auto targetEdge = simpleRouteEdges[0]->twin();
                    observer.handleEdgeReroute(pathIds, { getEdgeId(targetEdge) });


                    // Fixed vertices of the edge that will remain
                    const auto firstFixedVert = targetEdge->source();
                    const auto secondFixedVert = targetEdge->target();

                    std::cout << "[FM] ### Special case : merging to edge " << getEdgeId(targetEdge) << ", with v=" << firstFixedVert->data().id << ","
                        << secondFixedVert->data().id << "\n";
                    std::vector<std::size_t> deletedVerts;
                    for (auto he : degree2Subpaths.front())
                    {
                        if (he->source() != firstFixedVert && he->source() != secondFixedVert)
                        {
                            // TODO: replace with merge to edge.
                            // Notify
                            observer.handleVertexDelete(he->source()->data().id);
                            deletedVerts.push_back(he->source()->data().id);
                        }
                        arrangement.remove_edge(he, true, true);
                    }
                    std::cout << "[FM]\tDeleted ";
                    for (auto v : deletedVerts) std::cout << v << ' ';
                    std::cout << '\n';
                    return;
                }

                // Handle reroutes of degree2 paths
                for (auto subPath : degree2Subpaths)
                {
                    // Delete all edges on the subpath
                    for (auto he : subPath)
                    {
                        arrangement.remove_edge(he, true, true);
                    }
                }

                // Remove old simple edges
                for (auto heHandle : simpleRouteEdges)
                {
                    arrangement.remove_edge(heHandle, false, false);
                }

                Traits::Make_segment_2 makeSeg;
                // Create the new edges
                for (auto vHandle : vertices)
                {
                    Segment seg = makeSeg(center, vHandle->point());
                    //NOTE: insert_at_vertices cannot be called with a segment that intersects something (apart from start and end vertex)
                    auto newHeHandle = arrangement.insert_at_vertices(seg, vHandle, newVertexHandle);
                    auto length = Traits::approximate_length(newHeHandle->source()->point(), newHeHandle->target()->point());
                    rerouter.add_subgraph_edge(newHeHandle->source()->data().id, newHeHandle->target()->data().id, newHeHandle->data().id.value(), length);
                }
                {
                    std::set<FaceHandle> handles;
                    for (auto& kv : fIdToHE)
                    {
                        if (handles.find(kv.second->face()) != handles.end()) throw std::runtime_error("Duplicate face after reconstruction");
                        handles.insert(kv.second->face());
                    }
                }
                // Restore face data
                for (const auto& kv : fIdToHE)
                {
                    kv.second->face()->set_data(faceDataMap[kv.first]);
                    remappedFaces[kv.first] = kv.second->face();
                    if (face->is_unbounded()) throw std::runtime_error("Cannot assign new unbounded face!");
                }
                observer.handleJunctionReroute(rerouter);
            }
        };


        class CollapseGoodContinuityFacesByEdge
        {
        public:
            using Traits = FaceMergeTraits;

            template<typename SelectPredicate, typename Observer = TrivialSimplificationObserver>
            void operator()(Traits::Arr& arrangement, Traits::Face_handle face,
                std::map<std::size_t, Traits::Face_handle>& remappedFaces, Observer& observer, const SelectPredicate& select)
            {
                std::vector<Traits::Face_handle> faces;
                std::vector<std::vector<Traits::Halfedge_handle>> interior;
                std::vector< FaceSegmentation> segmentations;
                find_face_continuation(face, select, faces, interior, segmentations);
                if (faces.size() == 1)
                {
                    //TODO: temporarily skip, but fix this later.
                    //split_single_face(face,segmentations[0],observer);
                    std::cout << "[ContinuityFM] Found single, not handled yet\n";
                    DeleteSmallFacesToNeighbourViaContraction deleter;
                    deleter(arrangement, face, remappedFaces, observer, select);
                    return;
                }
                if (faces.empty()) throw std::runtime_error("Expected non-empty face list");
                std::cout << "[ContinuityFM] Merging " << faces.size() << " faces\n";
                merge_continuation(arrangement, faces, interior, segmentations, observer);

            }

            // Flag that we are not going to try to remap faces. 
            bool recomputeAfterMerge() const
            {
                return true;
            }
            // Boundary path should span at most this angle.
            void set_max_turning_angle(Traits::NT angle)
            {
                m_max_turning_angle = Traits::CGAL_NT(angle);
            }
        private:
            Operators::ComputeBoundaryTypes<Traits> compute_boundary_types;
            Operators::ComputeTurningAngle<Traits> compute_turning_angle;
            struct FaceSegmentation
            {
                using BoundaryPath = std::vector<Traits::Halfedge_handle>;
                // Top paths
                std::vector<BoundaryPath> top;
                // Bottoms paths
                std::vector<BoundaryPath> bottom;

                /**
                 * \brief Return separating vertices for the given path, excluding the start
                 * \param paths
                 * \return
                 */
                static std::vector<Traits::Vertex_handle> separatingVertices(const std::vector<BoundaryPath>& paths)
                {
                    if (paths.empty()) return { };
                    std::vector<Traits::Vertex_handle> returnVal;
                    for (const auto& path : paths)
                    {
                        returnVal.push_back(path.back()->target());
                    }
                    return returnVal;
                }
            };
            template<typename Observer = EdgeTrackingObserver>
            void split_single_face(Traits::Face_handle face, const FaceSegmentation& seg, Observer& observer)
            {

            }
            struct FaceVertices
            {
                std::vector<Traits::Vertex_handle> top;
                std::vector<Traits::Vertex_handle> bottom;
                void setFromSegmentation(const FaceSegmentation& seg)
                {
                    top = FaceSegmentation::separatingVertices(seg.top);
                    bottom = FaceSegmentation::separatingVertices(seg.bottom);
                }
            };

            template<typename Observer = EdgeTrackingObserver>
            void merge_continuation(Traits::Arr& arrangement, std::vector<Traits::Face_handle> faces, const std::vector<std::vector<Traits::Halfedge_handle>>& interior,
                const std::vector< FaceSegmentation>& segs, Observer& observer)
            {

                

                // Compute start and end positions for first and last face
                Operators::ComputeFaceCenter<Traits> compute_center;
                auto startPos = compute_center(faces.front());
                auto endPos = compute_center(faces.back());

                // Compute midpoints for interior paths
                Operators::ComputeMeanPosition<Traits> compute_mean_pos;
                std::vector<Traits::Point_2> midPoints;
                for (const auto& path : interior)
                {
                    midPoints.push_back(compute_mean_pos(path.begin(), path.end(), [](auto pathIt)
                    {
                        return (*pathIt)->source()->point();
                    }, std::vector<Traits::Point_2>{path.back()->target()->point()}));
                }

                // Compute separating vertices that should be retained
                std::vector<FaceVertices> faceVerts;
                for (const auto& seg : segs)
                {
                    faceVerts.push_back({});
                    auto& currFaceVerts = faceVerts.back();
                    currFaceVerts.setFromSegmentation(seg);
                    if (!currFaceVerts.top.empty())
                    {
                        if (currFaceVerts.top.back()->degree() <= 3) currFaceVerts.top.pop_back();
                    }
                    if (!currFaceVerts.bottom.empty())
                    {
                        if (currFaceVerts.bottom.back()->degree() <= 3) currFaceVerts.bottom.pop_back();
                    }
                }
                // Add as terminal vertices
                std::set<std::size_t> terminalVerts;
                for (const auto& fv : faceVerts)
                {
                    for (auto vh : fv.top) terminalVerts.insert(vh->data().id);
                    for (auto vh : fv.bottom) terminalVerts.insert(vh->data().id);
                }

                FaceTreeReroute<Traits::NT> rerouter(terminalVerts);
                // Add all affected edges to our rerouter.
                for (const auto& path : interior)
                {
                    for (auto he : path)
                    {
                        rerouter.add_edge_from_handle(he);
                    }
                }
                // Delete all involved edges, collect vertices to reconnect.
                // Add all affected edges to our rerouter.
                for (const auto& path : interior)
                {
                    for (auto he : path)
                    {
                        arrangement.remove_edge(he);
                    }
                }
                for (const auto& seg : segs)
                {
                    for (const auto& path : seg.top)
                    {
                        for (auto he : path)
                        {
                            arrangement.remove_edge(he);
                        }
                    }
                    for (const auto& path : seg.bottom)
                    {
                        for (auto he : path)
                        {
                            arrangement.remove_edge(he);
                        }
                    }
                }

                for (const auto& seg : segs)
                {
                    for (const auto& path : seg.top)
                    {
                        for (auto he : path)
                        {
                            rerouter.add_edge_from_handle(he);
                        }
                    }
                    for (const auto& path : seg.bottom)
                    {
                        for (auto he : path)
                        {
                            rerouter.add_edge_from_handle(he);
                        }
                    }
                }

                // Track everything that we add.
                RecordAddedVerticesAndEdges recorder;
                recorder.attach(arrangement);

                // Reconnect first and last face
                auto v0 = recorder.insert_vertex(startPos);
                auto v1 = recorder.insert_vertex(endPos);
                for (std::size_t i = 0; i < faceVerts.front().top.size() - 1; ++i)
                {
                    CGAL::insert(arrangement, segmentBetweenVertices(v0, faceVerts.front().top[i]));
                }
                for (std::size_t i = 0; i < faceVerts.back().bottom.size() - 1; ++i)
                {
                    CGAL::insert(arrangement, segmentBetweenVertices(v1, faceVerts.back().bottom[i]));
                }

                // Insert midpoints
                std::vector<Traits::Vertex_handle> midPointVertices;
                for (const auto& mp : midPoints)
                {
                    midPointVertices.push_back(recorder.insert_vertex(mp));
                }

                // Reconnect intermediate faces
                for (std::size_t i = 1; i < faceVerts.size() - 1; ++i)
                {
                    project_segmentation(arrangement, faceVerts[i], i - 1, midPointVertices);
                }
                for (const auto& v : recorder.added_vertices())
                {
                    rerouter.add_subgraph_vertex(v->data().id);
                }
                for (auto he : recorder.added_edges())
                {
                    auto len = (Traits::cgalPositionToGcpp(he->source()->point()) - Traits::cgalPositionToGcpp(he->target()->point())).length();
                    rerouter.add_subgraph_edge(he->source()->data().id, he->target()->data().id, he->data().id.value(),len);
                }
                // Notify our observer
                observer.handleJunctionReroute(rerouter);
            }

            struct BoundaryInsertPoint
            {
                bool isFromTop = false;
                std::size_t subpathThatTargetsVert = 0; //Subpath with end the vertex that should be inserted.
                // Projection point on interior edge
                Traits::Point_2 point;
                // Vertex that was projected
                Traits::Vertex_handle srcVertex;
            };

            static std::vector<std::size_t> computePathIds(const std::vector<Traits::Halfedge_handle>& path, bool reverse = false)
            {
                std::vector<std::size_t> ids;
                if (reverse)
                {
                    for (auto it = path.rbegin(); it != path.rend(); ++it) ids.push_back((*it)->twin()->data().id.value());
                }
                else
                {
                    for (auto he : path) ids.push_back(he->data().id.value());
                }
                return ids;
            }

            static Traits::Segment_2 segmentBetweenVertices(Traits::Vertex_handle v0, Traits::Vertex_handle v1)
            {
                Traits::Make_segment_2 makeSeg;
                return makeSeg(v0->point(), v1->point());
            }

            /**
             * \brief Project elements of middle face
             * \param arrangement
             * \param faceVertices
             * \param index
             * \param newVertices
             */
            void project_segmentation(Traits::Arr& arrangement, const FaceVertices& faceVertices, std::size_t index, const std::vector<Traits::Vertex_handle>& newVertices)
            {
                // Vertices of the interior paths that were deleted
                const auto newV0 = newVertices[index];
                const auto newV1 = newVertices[index + 1];

                // Handle simple and nice case.
                if (faceVertices.bottom.size() <= 1 && faceVertices.top.size() <= 1)
                {
                    CGAL::insert(arrangement, segmentBetweenVertices(newV0, newV1));
                    if (!faceVertices.bottom.empty()) CGAL::insert(arrangement, segmentBetweenVertices(faceVertices.bottom.front(), newV1));
                    if (!faceVertices.top.empty()) CGAL::insert(arrangement, segmentBetweenVertices(faceVertices.top.front(), newV0));
                }
                // We need to split in the new edge for the degree>2 vertices on the boundaries.
                else
                {
                    Traits::Kernel::Line_2 support_line(newV0->point(), newV1->point());

                    // Project intermediate vertices on the line.
                    std::vector<BoundaryInsertPoint> projections;
                    for (std::size_t j = 0; j < faceVertices.top.size() - 1; ++j)
                    {
                        projections.push_back(BoundaryInsertPoint{
                            true,
                            j,
                            support_line.projection(faceVertices.top[j]->point()),
                            faceVertices.top[j]
                            });
                    }
                    for (std::size_t j = 0; j < faceVertices.bottom.size() - 1; ++j)
                    {
                        projections.push_back(BoundaryInsertPoint{
                            false,
                            j,
                            support_line.projection(faceVertices.bottom[j]->point()),
                            faceVertices.bottom[j]
                            });
                    }
                    // Somehow sort along edge.
                    auto distanceLine = support_line.perpendicular(newV0->point()).opposite();
                    struct Sorter
                    {
                        Traits::Kernel::Line_2 distanceLine; //Compute distances relative to this line.

                        bool operator()(const BoundaryInsertPoint& p0, const BoundaryInsertPoint& p1) const
                        {
                            return CGAL::has_smaller_signed_distance_to_line(distanceLine, p0.point, p1.point);
                        }
                    };
                    Sorter sorter{ distanceLine };
                    std::sort(projections.begin(), projections.end(), sorter);
                    auto start = projections.begin();
                    auto end = projections.end();
                    // TODO: make sure all points map to the correct range.
                    /*for(;start != end; ++start)
                    {
                        if (!distanceLine.has_on_negative_side(start->point)) break;
                    }*/
                    // Leave vertex pruning for later step
                    // Create the edges

                    // Keep track of insertions. Note that we expect some other algorithm to assign proper ID's!
                    RecordAddedVerticesAndEdges addedEls;
                    addedEls.attach(arrangement);
                    for (; start != end; ++start)
                    {
                        auto prevV = addedEls.added_vertices().size() == 0 ? newV0 : addedEls.added_vertices().back();
                        CGAL::insert_point(arrangement, start->point);
                        CGAL::insert(arrangement, segmentBetweenVertices(prevV, addedEls.added_vertices().back()));
                    }
                    CGAL::insert(arrangement, segmentBetweenVertices(addedEls.added_vertices().back(), newV1));

                    // Reconnect vertices at the boundary.
                    for (std::size_t j = 0; j < projections.size(); ++j)
                    {
                        CGAL::insert(arrangement, segmentBetweenVertices(addedEls.added_vertices().at(j), projections[j].srcVertex));
                    }
                    addedEls.detach();
                }
            }

            void cleanup_interior(const std::vector<std::vector<Traits::Halfedge_handle>>& interior, Traits::Arr& arrangement, std::vector<Traits::Vertex_handle>& newVertices)
            {
                // For every path, we insert a point at its mean, and remove all edges associated to the path
                for (std::size_t i = 0; i < interior.size(); ++i)
                {
                    const auto& path = interior[i];
                    // Compute mean
                    auto initial = path.front()->source()->point() - CGAL::ORIGIN;
                    for (std::size_t i = 0; i < path.size(); ++i)
                    {
                        initial += path[i]->target()->point() - CGAL::ORIGIN;
                    }
                    initial /= Traits::CGAL_NT(static_cast<Traits::NT>(path.size() + 1));
                    Traits::Kernel::Point_2 point;
                    Traits::reducePrecision(CGAL::ORIGIN + initial, point);


                    // This will contain the merged face after removing an edge.
                    Traits::Face_handle targetFace;
                    auto startVert = path.front()->source();
                    auto endVert = path.back()->target();
                    for (auto he : path) {
                        targetFace = arrangement.remove_edge(he);
                    }
                    // TODO: check point is in face? On average will be okay...?
                    Traits::Vertex_handle newVert = arrangement.insert_in_face_interior(point, targetFace);
                    newVertices.push_back(newVert);
                    CGAL::insert(arrangement, segmentBetweenVertices(startVert, newVert));
                    CGAL::insert(arrangement, segmentBetweenVertices(newVert, endVert));
                }
            }

            class RecordAddedVerticesAndEdges : public CGAL::Arr_observer<Traits::Arr>
            {
            public:
                void after_create_vertex(Vertex_handle vh) override
                {
                    m_added_vertices.push_back(vh);
                }
                const std::vector<Traits::Vertex_handle>& added_vertices()
                {
                    return m_added_vertices;
                }
                const std::vector<Traits::Halfedge_handle>& added_edges()
                {
                    return m_added_edges;
                }

                Traits::Vertex_handle insert_vertex(const Traits::Point_2& point)
                {
                    const auto curr_size = m_added_vertices.size();
                    auto vh = CGAL::insert_point(*arrangement(), point);
                    const auto new_size = m_added_vertices.size();
                    if (new_size != curr_size + 1)
                    {
                        m_added_vertices.push_back(vh);
                    }
                    return m_added_vertices.back();
                }
                Traits::Halfedge_handle insert_segment(const Traits::Segment_2& point)
                {
                    const auto curr_size = m_added_edges.size();
                    CGAL::insert(*arrangement(), point);
                    const auto new_size = m_added_edges.size();
                    if (new_size != curr_size + 1) throw std::runtime_error("Something is off with edge insertion");
                    return m_added_edges.back();
                }

                void after_create_edge(Halfedge_handle he) override
                {
                    m_added_edges.push_back(he);
                }
                void after_split_edge(Halfedge_handle, Halfedge_handle) override
                {
                    throw std::runtime_error("Not expecting splits!");
                }
            private:
                std::vector<Traits::Vertex_handle> m_added_vertices;
                std::vector<Traits::Halfedge_handle> m_added_edges;
            };

            struct Wedge
            {
                Traits::Kernel::Direction_2 lower, upper;
                bool m_exceedsPi = false;
                Wedge() {}
                Wedge(Traits::Kernel::Direction_2 startDir) :lower(startDir), upper(startDir) {}

                bool isValid() const
                {
                    return !m_exceedsPi;
                }

                // Prefers update to lower, 
                void expand(const Traits::Kernel::Direction_2& dir)
                {
                    if (m_exceedsPi) return;
                    if ( lower == dir || upper == dir || (lower != upper && dir.counterclockwise_in_between(lower, upper))) return;

                    if (dir.counterclockwise_in_between(upper, -lower)) upper = dir;
                    else if (dir.counterclockwise_in_between(-upper, lower)) lower = dir;
                    else
                    {
                        m_exceedsPi = true;
                    }
                }

                bool overlaps(const Wedge& other) const
                {
                    return containsDirection(other.lower) || containsDirection(other.upper) ||
                        other.containsDirection(lower) || other.containsDirection(upper);
                }
                bool contains(const Wedge& other) const
                {
                    return containsDirection(other.lower) && containsDirection(other.upper);
                }
                Wedge join(const Wedge& other) const
                {
                    Wedge w = *this;
                    w.expand(other.lower);
                    w.expand(other.upper);
                    return w;
                }

                bool containsDirection(const Traits::Kernel::Direction_2& dir) const
                {
                    return dir.counterclockwise_in_between(lower, upper) || dir == lower || dir == upper;
                }

                bool expand_with_max(const Traits::Halfedge_handle& dirEdge, Traits::CGAL_NT angle)
                {
                    return expand_with_max(dirFromHalfedge(dirEdge), angle);
                }
                bool expand_with_max(const Traits::Kernel::Direction_2& dir, Traits::CGAL_NT angle)
                {
                    if (would_exceed_max_within_pi(dir, angle)) return false;

                    expand(dir);
                    return isValid();
                }

                bool would_exceed_max_within_pi(const Traits::Kernel::Direction_2& dir, Traits::CGAL_NT angle) const
                {
                    Wedge w = *this;
                    w.expand(dir);
                    return w.exceeds_max_within_pi(angle);
                }

                bool exceeds_max_within_pi(Traits::CGAL_NT angle) const
                {
                    if ((!upper.counterclockwise_in_between(lower, -lower)&&upper != lower) || m_exceedsPi)
                    {
                        return true;
                    }
                    return Traits::CGAL_NT(Traits::approximate_angle(lower.vector(), upper.vector())) > angle;
                }
                Traits::NT approximate_angle() const
                {
                    if (m_exceedsPi) return M_PI;
                    return Traits::approximate_angle(lower.vector(), upper.vector());
                }
            };

            void compute_boundary(const std::vector<Traits::Face_handle>& faces, const std::vector<std::vector<Traits::Halfedge_handle>>& interior,
                std::list<Traits::Halfedge_handle>& boundary)
            {
                // Expecting that all faces are not unbounded!

                // TODO: what would be logical here.
                if (interior.size() <= 1) return;

                { // Start with first face
                    const auto& path = interior[0];
                    Traits::Halfedge_handle curr = path.back();
                    while (true)
                    {
                        curr = curr->next();
                        if (curr == path.front()) break;
                        boundary.push_back(curr);
                    }
                }

                for (std::size_t i = 1; i < interior.size() - 1; ++i)
                {
                    Traits::Halfedge_handle curr = interior[i].front()->twin();
                    //'Lower' connection
                    while (true)
                    {
                        curr = curr->next();
                        if (curr == interior[i + 1].front()) break;
                        boundary.push_back(curr);
                    }
                    //'Upper' connection
                    curr = interior[i].back()->twin();
                    while (true)
                    {
                        curr = curr->prev();
                        if (curr == interior[i + 1].back()) break;
                        boundary.push_front(curr);
                    }
                }
                { // Last face
                    const auto& path = interior.back();
                    Traits::Halfedge_handle curr = path.front()->twin();
                    while (true)
                    {
                        curr = curr->next();
                        if (curr == path.back()->twin()) break;
                        boundary.push_back(curr);
                    }
                }
            }

            template<typename SelectPredicate>
            void find_face_continuation(Traits::Face_handle face, const SelectPredicate& pred, std::vector<Traits::Face_handle>& faces, std::vector<std::vector<Traits::Halfedge_handle>>& interior,
                std::vector<FaceSegmentation>& segmentations)
            {
                faces.push_back(face);

                std::cout << "[ContinuityFM] Searching with face " << face->data().id.value() << " face\n";

                std::vector<std::vector<Traits::Halfedge_handle>> subpaths;
                compute_boundary_types(face, subpaths);

                std::cout << "[ContinuityFM] Boundary paths: " << subpaths.size() << "\n";

                std::optional<Wedge> bestOverAll;
                std::vector<Traits::Face_handle> bestFaces;
                for (const auto& subpath : subpaths)
                {
                    if (subpath.front()->twin()->face()->is_unbounded()) continue;
                    auto targetFace = subpath.front()->twin()->face();
                    // Construct wedge by directions of edges pointing 'into' the next face.
                    Wedge initialWedge(dirFromHalfedge(subpath.back()->next()->twin())); //TODO: setup
                    std::cout << "[ContinuityFM] Handling subpath\n";
                    std::cout << "[ContinuityFM] Wedge angle: " << initialWedge.approximate_angle() << "\n";
                    if (!initialWedge.expand_with_max(dirFromHalfedge(subpath.front()->prev()), m_max_turning_angle))
                    {
                        initialWedge.expand(dirFromHalfedge(subpath.front()->prev()));
                        std::cout << "[ContinuityFM] Wedge angle combined: " << initialWedge.approximate_angle() << ", which is too much\n";
                        continue;
                    }
                    Wedge bestWedge;
                    std::vector<Traits::Face_handle> outFaces = faces;
                    // If fail, continue with next
                    if (!find_face_continuation_with_target(face, targetFace, subpath, initialWedge, pred, outFaces, bestWedge))continue;
                    std::cout << "[ContinuityFM] Computed continuation with wedge angle: " << bestWedge.approximate_angle() << "\n";
                    if (!bestOverAll.has_value() || bestOverAll.value().approximate_angle() > bestWedge.approximate_angle())
                    {
                        bestOverAll = bestWedge;
                        bestFaces = outFaces;
                    }
                }
                if (bestOverAll.has_value())
                {
                    faces = bestFaces;
                }
                // Compute interior and segmentation
                if(faces.size() > 0)
                {
                    compute_segmentation(faces, interior, segmentations);
                }
            }

            void compute_segmentation(const std::vector<Traits::Face_handle>& faces, std::vector<std::vector<Traits::Halfedge_handle>>& interior,
                std::vector<FaceSegmentation>& segmentations)
            {
                for(std::size_t i = 0; i < faces.size(); ++i)
                {
                    segmentations.push_back(FaceSegmentation{});
                    auto& segmentation = segmentations.back();
                    std::vector<std::vector<Traits::Halfedge_handle>> ccwSubpaths;
                    compute_boundary_types(faces[i], ccwSubpaths);

                    std::optional<std::size_t> forwardPath;
                    std::optional<std::size_t> backwardPath;
                    if(i != 0)
                    {
                        for(std::size_t j = 0; j < ccwSubpaths.size(); ++j)
                        {
                            if(ccwSubpaths[j].front()->twin()->face() == faces[i-1])
                            {
                                backwardPath = j;
                                break;
                            }
                        }
                        if (!backwardPath.has_value()) throw std::runtime_error("Could not find backward path?");
                    }
                    if(i != faces.size() -1)
                    {
                        for (std::size_t j = 0; j < ccwSubpaths.size(); ++j)
                        {
                            if (ccwSubpaths[j].front()->twin()->face() == faces[i + 1])
                            {
                                forwardPath = j;
                                break;
                            }
                        }
                        if (!forwardPath.has_value()) throw std::runtime_error("Could not find forward path?");
                    }

                    if(forwardPath.has_value())
                    {
                        interior.push_back(ccwSubpaths[forwardPath.value()]);
                    }
                    auto nextIndex = [&ccwSubpaths](auto ind)
                    {
                        return (ind + 1) % ccwSubpaths.size();
                    };
                    if(forwardPath.has_value())
                    {
                        std::size_t j = nextIndex(forwardPath.value());
                        while(true)
                        {
                            if (j == forwardPath.value())break;
                            if (backwardPath.has_value() && j == backwardPath.value())break;
                            segmentation.top.push_back(ccwSubpaths[j]);
                            j = nextIndex(j);
                        }
                    }
                    if (backwardPath.has_value())
                    {
                        std::size_t j = nextIndex(backwardPath.value());
                        while (true)
                        {
                            if (j == backwardPath.value())break;
                            if (forwardPath.has_value() && j == forwardPath.value())break;
                            segmentation.bottom.push_back(ccwSubpaths[j]);
                            j = nextIndex(j);
                        }
                    }
                }
            }

            static Traits::Kernel::Direction_2 dirFromHalfedge(Traits::Halfedge_handle e)
            {
                return Traits::Kernel::Direction_2(e->target()->point() - e->source()->point());
            }


            template<typename SelectPredicate>
            bool find_face_continuation_with_target(Traits::Face_handle face, Traits::Face_handle targetFace, const std::vector<Traits::Halfedge_handle>& separatorPath,
                const Wedge& initialWedge, const SelectPredicate& pred, std::vector<Traits::Face_handle>& faces, Wedge& bestWedge)
            {
                if (!targetFace->has_outer_ccb() || !pred(targetFace, targetFace->outer_ccb()))
                {
                    std::cout << "[ContinuityFM] Face is not selected: has outer ccb?" << targetFace->has_outer_ccb() << '\n';
                    return false;
                }
                Wedge wedge = initialWedge;
                Wedge topWedge = wedge;
                Wedge bottomWedge = wedge;
                // Running top in backwards direction to stay in face.
                Traits::Halfedge_handle nextTop = separatorPath.back()->twin()->prev();
                Traits::Halfedge_handle nextBottom = separatorPath.front()->twin()->next();

                if (!bottomWedge.expand_with_max(dirFromHalfedge(nextBottom), m_max_turning_angle))
                {
                    std::cout << "[ContinuityFM] Bottom continuation does not fit" << '\n';
                    return false;
                }
                // Expand with top edge
                if (!topWedge.expand_with_max(dirFromHalfedge(nextTop->twin()), m_max_turning_angle))
                {
                    std::cout << "[ContinuityFM] Top continuation does not fit" << '\n';
                    return false;
                }
                std::vector<std::vector<Traits::Halfedge_handle>> subpaths;
                auto moduloInc = [&subpaths](auto curr)
                {
                    return (curr + 1) % subpaths.size();
                };
                auto moduloDec = [&subpaths](auto curr)
                {
                    return curr == 0 ? subpaths.size() - 1 : curr + 1;
                };
                compute_boundary_types(face, subpaths);
                std::set<std::size_t> ignorePaths;
                std::size_t startPath = 0;
                for (std::size_t i = 0; i < subpaths.size(); ++i)
                {
                    auto it = std::find(subpaths[i].begin(), subpaths[i].end(), nextTop->next());
                    if (it != subpaths[i].end()) {
                        startPath = i;
                        ignorePaths.insert(i);
                    }
                    if (subpaths[i].front()->twin()->face()->is_unbounded())ignorePaths.insert(i);
                }
                // For each boundary, determine which subpaths can de reached within acceptable angles
                std::set<std::size_t> reachable;
                std::size_t topReachable = moduloDec(startPath);
                std::size_t bottomReachable = moduloInc(startPath);
                std::map<std::size_t, Wedge> topReachableWedges;
                std::map<std::size_t, Wedge> bottomReachableWedges;
                while (true)
                {
                    nextTop = nextTop->prev();
                    if (topWedge.expand_with_max(nextTop->twin(), m_max_turning_angle))
                    {
                        std::cout << "- " << topWedge.approximate_angle();
                        if (nextTop->source()->degree() > 2) {
                            topReachable = moduloDec(topReachable);
                            topReachableWedges[topReachable] = topWedge;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
                std::cout << '\n';
                while (true)
                {
                    nextBottom = nextBottom->next();
                    if (topWedge.expand_with_max(nextBottom, m_max_turning_angle))
                    {
                        std::cout << "- " << topWedge.approximate_angle();
                        if (nextBottom->target()->degree() > 2) {
                            bottomReachable = moduloInc(bottomReachable);
                            if (bottomReachable == topReachable)
                            {
                                reachable.insert(bottomReachable);
                                bottomReachableWedges[bottomReachable] = bottomWedge;
                            }
                            else if (!reachable.empty())
                            {
                                reachable.insert(bottomReachable);
                                bottomReachableWedges[bottomReachable] = bottomWedge;
                            }
                        }
                    }
                    else
                    {
                        break;
                    }
                }

                std::cout << '\n';
                std::cout << "Reachable subpaths:";
                for (auto el : reachable) std::cout << " " << el;
                std::cout << '\n';

                std::optional<Wedge> recurseBestWedge;
                std::vector< Traits::Face_handle> bestFaces;
                for (auto reachablePath : reachable)
                {
                    if (ignorePaths.find(reachablePath) != ignorePaths.end()) continue;
                    std::vector<Traits::Face_handle> facesNext = faces;
                    facesNext.push_back(targetFace);
                    Wedge currentWedge;//TODO: combine top and bottom
                    currentWedge = topReachableWedges[reachablePath].join(bottomReachableWedges[reachablePath]);
                    std::cout << "[ContinuityFM] Continuing with wedge angle: " << currentWedge.approximate_angle() << "\n";
                    Wedge localWedge;
                    bool result = find_face_continuation_with_target(targetFace, subpaths[reachablePath].front()->twin()->face(), subpaths[reachablePath], currentWedge, pred, facesNext, localWedge);
                    if (!result)continue;
                    if (!recurseBestWedge.has_value() || recurseBestWedge.value().approximate_angle() > localWedge.approximate_angle())
                    {
                        recurseBestWedge = wedge;
                        bestFaces = facesNext;
                    }
                }
                if (!recurseBestWedge.has_value())
                {
                    Wedge w = initialWedge;
                    w.expand_with_max(dirFromHalfedge(nextBottom), m_max_turning_angle);
                    w.expand_with_max(dirFromHalfedge(nextTop->twin()), m_max_turning_angle);
                    bestWedge = w;
                    faces.push_back(targetFace);
                }
                else
                {
                    bestWedge = recurseBestWedge.value();
                    faces = bestFaces;
                }
                return true;
            }

            /**
             * \brief collects adjacent faces and associated edge between faces. Ignore unbounded.
             * \param face
             * \param adjacentFaces
             */
            void compute_adjacent_face(Traits::Face_handle face, std::map<Traits::Face_handle, Traits::HE_handle>& adjacentFaces)
            {
                Traits::for_each_ccb_edge(face, [&adjacentFaces](auto edge, auto& cmd)
                {
                    if (!edge->twin()->face()->is_unbounded())
                    {
                        adjacentFaces[edge->twin()->face()] = edge;
                    }
                });
            }
            Traits::CGAL_NT m_max_turning_angle = Traits::CGAL_NT(30. / 180. * M_PI);

            Traits::NT m_max_dist_for_vertex_merge = 1.;
        };

        /*template<typename...Approaches>
        class MultiApproach
        {
            std::tuple<Approaches...> m_approaches;
            template<typename Arrangement, typename FaceIterator, std::size_t...Is>
            void apply(Arrangement& arrangement, FaceIterator face, std::index_sequence<Is...>)
            {
                (std::get<Is>(m_approaches)(arrangement, face), ...);
            }
        public:
            template<typename Arrangement, typename FaceIterator>
            void operator()(Arrangement& arrangement, FaceIterator face)
            {
                apply(arrangement, face, std::index_sequence_for<Approaches...>{});
            }
        };*/
    }
}

#endif