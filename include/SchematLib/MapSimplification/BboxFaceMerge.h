#ifndef SCHEMATLIB_MAPSIMPLIFICATION_BBOXFACEMERGE_H
#define SCHEMATLIB_MAPSIMPLIFICATION_BBOXFACEMERGE_H

#include "FaceMergeStrategies.h"
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/convex_hull_traits_2.h>
#include <CGAL/min_quadrilateral_2.h>
namespace SchematLib::MapSimplification::FaceMergeApproach
{

    class BboxFaceMerge
    {
        using Traits = FaceMergeTraits;

        // The factor above which we deem the effect minor when removing this face.
        double m_faceSizeFactor = 2;

        template<typename T>
        static void apply_permutation(const std::vector<T>& source, std::vector<T>& target, const std::vector<std::size_t>& indexForLocation)
        {
            for (auto i : indexForLocation)
            {
                target.push_back(source[i]);
            }
        }

        template<typename T, typename Compare, typename Other>
        static void coordinated_sort(std::vector<T>& toSort, Compare&& compareFunc, std::vector<Other>& coordinatedOther)
        {
            // Fill indices
            std::vector<std::size_t> is(toSort.size(), 0);
            std::iota(is.begin(), is.end(), 0);
            std::sort(is.begin(), is.end(), [&toSort, &compareFunc](auto i0, auto i1)
            {
                return compareFunc(toSort[i0], toSort[i1]);
            });
            std::vector<T> temp;
            // Apply permutation
            apply_permutation(toSort, temp, is);
            std::vector<Other> temp2;
            apply_permutation(coordinatedOther, temp2, is);
            toSort = std::move(temp);
            coordinatedOther = std::move(temp2);
        }

        class DidInsertVertex : public CGAL::Arr_observer<Traits::Arr>
        {
            bool m_didInsert = false;
            std::optional< Vertex_handle> m_created;
        public:
            void clear()
            {
                m_didInsert = false;
                m_created = {};
            }
            bool didInsert() const
            {
                return m_didInsert;
            }
            std::optional< Vertex_handle> created()const
            {
                return m_created;
            }
            void after_create_vertex(Vertex_handle vh) override
            {
                m_didInsert = true;
                m_created = vh;
            }
        };

        class SplitterObserver : public CGAL::Arr_observer<Traits::Arr>
        {
            Halfedge_handle startEdge;
            Halfedge_handle endEdge;
            Halfedge_handle cwDeleteEdge;
            Halfedge_handle ccwDeleteEdge;

            std::optional<Halfedge_handle> createdHalfEdge;
            std::optional<Vertex_handle> splitV;
            std::optional<std::size_t>splitId;
            Vertex_handle cwDeleteTarget;
            Vertex_handle ccwDeleteSrc;
            bool splittingCw = false;
            bool splittingCcw = false;
        public:
            SplitterObserver(Halfedge_handle startEdge,
                Halfedge_handle endEdge,
                Halfedge_handle cwDeleteEdge,
                Halfedge_handle ccwDeleteEdge) : startEdge(startEdge), endEdge(endEdge), cwDeleteEdge(cwDeleteEdge), ccwDeleteEdge(ccwDeleteEdge)
            {
                cwDeleteTarget = cwDeleteEdge->target();
                ccwDeleteSrc = ccwDeleteEdge->source();
            }

            std::map<std::size_t, std::pair<Halfedge_handle, Halfedge_handle>> splits;

            void after_create_edge(Halfedge_handle he) override
            {
                createdHalfEdge = he;
            }
            Halfedge_handle deleteEdgeCw() const
            {
                return cwDeleteEdge;
            }
            Halfedge_handle deleteEdgeCcw() const
            {
                return ccwDeleteEdge;
            }
            std::optional<Halfedge_handle> halfEdgeCreated()const
            {
                return createdHalfEdge;
            }

            void after_create_vertex(Vertex_handle) override {}
            void after_split_edge(Halfedge_handle he0, Halfedge_handle he1) override
            {
                splits[splitId.value()] = std::make_pair(he0, he1);
                if (splittingCw)
                {
                    std::cout << "[BboxMerge] CW edge was split\n";
                    splittingCw = false;
                    if (he0->target() == cwDeleteTarget || he0->source() == cwDeleteTarget)
                    {
                        cwDeleteEdge = he0->target() == cwDeleteTarget ? he0 : he0->twin();
                    }
                    else if (he1->target() == cwDeleteTarget || he1->source() == cwDeleteTarget)
                    {
                        cwDeleteEdge = he1->target() == cwDeleteTarget ? he1 : he1->twin();
                    }
                    else
                    {
                        throw std::runtime_error("DERP");
                    }
                }
                if (splittingCcw)
                {
                    std::cout << "[BboxMerge] CCW edge was split\n";
                    splittingCcw = false;
                    if (he0->source() == ccwDeleteSrc || he0->target() == ccwDeleteSrc)
                    {
                        ccwDeleteEdge = he0->source() == ccwDeleteSrc ? he0 : he0->twin();
                    }
                    else if (he1->source() == ccwDeleteSrc || he1->target() == ccwDeleteSrc)
                    {
                        ccwDeleteEdge = he1->source() == ccwDeleteSrc ? he1 : he1->twin();
                    }
                    else
                    {
                        throw std::runtime_error("DERP");
                    }
                }
            }
            void after_split_face(Face_handle, Face_handle, bool) override {}
            void before_split_edge(Halfedge_handle he, Vertex_handle v, const X_monotone_curve_2&,
                const X_monotone_curve_2&) override
            {
                splitId = he->data().id.value();
                splitV = v;
                if (he == cwDeleteEdge || he->twin() == cwDeleteEdge)
                {
                    splittingCw = true;
                }
                else if (he == ccwDeleteEdge || he->twin() == ccwDeleteEdge)
                {
                    splittingCcw = true;
                }
            }
            void before_split_face(Face_handle, Halfedge_handle) override {}
        };
        bool m_verbose = false;

        void compute_splitter(const std::vector<Traits::Point_2>& bbox, Traits::Segment_2& axisSegment, std::size_t& majorAxis, std::size_t& minorAxis)
        {
            auto l0 = Traits::approximate_length(bbox[0], bbox[1]);
            auto l1 = Traits::approximate_length(bbox[1], bbox[2]);
            majorAxis = l0 > l1 ? 0 : 1;
            minorAxis = 1 - majorAxis;
            auto dir = CGAL::midpoint(bbox[minorAxis + 2], bbox[minorAxis + 3]) - CGAL::midpoint(bbox[minorAxis], bbox[minorAxis + 1]);
            // Segment along major axis through middle of box.
            axisSegment = Traits::Make_segment_2()(
                CGAL::midpoint(bbox[minorAxis], bbox[minorAxis + 1]) - dir,
                CGAL::midpoint(bbox[minorAxis + 2], bbox[minorAxis + 3]) + dir
                );
            if (m_verbose) { std::cout << "Seg: " << axisSegment.left() << " " << axisSegment.right() << '\n'; }
        }

        template<typename Observer>
        static void handleNewEdge(Traits::Halfedge_handle he, Observer& observer)
        {
            observer.handleNewEdge(he->data().id.value(), he->source()->data().id, he->target()->data().id);
        }
        template<typename Observer>
        static void handleEdgeReplace(std::size_t e, Traits::Halfedge_handle he0, Traits::Halfedge_handle he1, Observer& observer)
        {
            observer.handleEdgeReplace(e, std::vector<std::size_t>{he0->data().id.value(), he1->data().id.value()});
        }
        static void printVert(Traits::Vertex_handle vert)
        {
            std::cout << "{" << vert->point().x() << "," << vert->point().y() << "}";
        }
        static void printEdge(Traits::Halfedge_handle edge)
        {
            printVert(edge->source());
            std::cout << ',';
            printVert(edge->target());
        }

        template<typename Observer>
        Traits::Halfedge_handle insert_bbox_segment(Traits::Arr& arrangement, Traits::Halfedge_handle startHe, Traits::Halfedge_handle endHe,
            const Traits::Point_2& intersectionStart, const Traits::Point_2& intersectionEnd, Observer& observer) const
        { //Construct
            auto startVs = std::make_pair(startHe->source(), startHe->target());
            auto endVs = std::make_pair(endHe->source(), endHe->target());

            auto e0 = startHe->data().id.value();
            auto e1 = endHe->data().id.value();

            DidInsertVertex tracker;
            tracker.attach(arrangement);
            Traits::Vertex_handle v0, v1;
            bool noRebuildE0 = false;
            bool noRebuildE1 = false;
            bool v0AlreadyExists = true;
            {
                auto pnt = Traits::reducedPrecision(intersectionStart);
                if (Traits::approximate_length(pnt, startHe->source()->point()) < 1.0)
                {
                    v0 = startHe->source();
                    noRebuildE0 = true;
                }
                else if (Traits::approximate_length(pnt, startHe->target()->point()) < 1.0)
                {
                    v0 = startHe->target();
                    noRebuildE0 = true;
                }
                else
                {
                    arrangement.remove_edge(startHe, false, false);
                    v0 = CGAL::insert_point(arrangement, Traits::reducedPrecision(intersectionStart));
                    v0AlreadyExists = !tracker.didInsert();
                }
            }


            tracker.clear();
            bool v1AlreadyExists = true;
            {
                auto pnt = Traits::reducedPrecision(intersectionEnd);
                if (Traits::approximate_length(pnt, endHe->source()->point()) < 1.0)
                {
                    v1 = endHe->source();
                    noRebuildE1 = true;
                }
                else if (Traits::approximate_length(pnt, endHe->target()->point()) < 1.0)
                {
                    v1 = endHe->target();
                    noRebuildE1 = true;
                }
                else
                {
                    arrangement.remove_edge(endHe, false, false);
                    v1 = CGAL::insert_point(arrangement, Traits::reducedPrecision(intersectionEnd));
                    v1AlreadyExists = !tracker.didInsert();
                }
            }
            if (m_verbose) { std::cout << "[BboxMerge] V0: " << v0->point() << " exists:" << v0AlreadyExists << '\n'; }
            if (m_verbose) { std::cout << "[BboxMerge] V1: " << v1->point() << " exists:" << v1AlreadyExists << '\n'; }
            if (!v0AlreadyExists)observer.handleNewVertex(v0->data().id);
            if (!v1AlreadyExists) observer.handleNewVertex(v1->data().id);
            // reconstrut
            if (v0AlreadyExists && !noRebuildE0)
            {
                auto newHe = arrangement.insert_at_vertices(Traits::Make_segment_2()(startVs.first->point(), startVs.second->point()), startVs.first, startVs.second);
                handleNewEdge(newHe, observer);
                observer.handleEdgeReplace(e0, std::vector<std::size_t>{newHe->data().id.value()});
            }
            else if (!v0AlreadyExists)
            {
                auto startHe0 = arrangement.insert_at_vertices(Traits::Make_segment_2()(startVs.first->point(), v0->point()), startVs.first, v0);
                auto startHe1 = arrangement.insert_at_vertices(Traits::Make_segment_2()(v0->point(), startVs.second->point()), v0, startVs.second);
                handleNewEdge(startHe0, observer);
                handleNewEdge(startHe1, observer);
                if (m_verbose)
                {
                    std::cout << "\tStart he0:";
                    printEdge(startHe0);
                    std::cout << "\n";
                    std::cout << "\tStart he1:";
                    printEdge(startHe1);
                    std::cout << "\n";
                }
                handleEdgeReplace(e0, startHe0, startHe1, observer);
            }

            if (v1AlreadyExists && !noRebuildE1)
            {
                auto newHe = arrangement.insert_at_vertices(Traits::Make_segment_2()(endVs.first->point(), endVs.second->point()), endVs.first, endVs.second);
                handleNewEdge(newHe, observer);
                observer.handleEdgeReplace(e1, std::vector<std::size_t>{newHe->data().id.value()});
            }
            else if (!v1AlreadyExists)
            {
                auto endHe0 = arrangement.insert_at_vertices(Traits::Make_segment_2()(endVs.first->point(), v1->point()), endVs.first, v1);
                auto endHe1 = arrangement.insert_at_vertices(Traits::Make_segment_2()(v1->point(), endVs.second->point()), v1, endVs.second);
                handleNewEdge(endHe0, observer);
                handleNewEdge(endHe1, observer);
                if (m_verbose)
                {
                    std::cout << "End he0:";
                    printEdge(endHe0);
                    std::cout << "\n";
                    std::cout << "End he1:";
                    printEdge(endHe1);
                    std::cout << "\n";
                }
                handleEdgeReplace(e1, endHe0, endHe1, observer);
            }

            // Add middle
            auto middleEdge = arrangement.insert_at_vertices(Traits::Make_segment_2()(v0->point(), v1->point()), v0, v1);
            if (middleEdge->source() != v0) throw std::runtime_error("TWISTED!");
            handleNewEdge(middleEdge, observer);
            if (m_verbose) 
            {
                std::cout << "[BboxMerge] Middle edge: ";
                printEdge(middleEdge);
                std::cout << '\n';
            }

            if(m_verbose)
            {
                std::size_t leftComplexity = 0, rightComplexity = 0;
                Traits::for_each_ccb_edge(middleEdge->face(), [&leftComplexity](auto face, auto&cmd)
                {
                    ++leftComplexity;
                });
                Traits::for_each_ccb_edge(middleEdge->twin()->face(), [&rightComplexity](auto face, auto&cmd)
                {
                    ++rightComplexity;
                });
                std::cout << "[bboxMerge] Left right face complexities: " << leftComplexity << " " << rightComplexity << '\n';
            }
            return middleEdge;
        }

        void find_edges_to_delete(Traits::Halfedge_handle middleEdge, std::vector<Traits::Halfedge_handle>& cwDelete,
            std::vector<Traits::Halfedge_handle>& ccwDelete)
        {
            // CCW direction
            {
                std::vector<Traits::NT> cumulativeApproxLength;
                {
                    auto curr = middleEdge->twin()->next();
                    while (curr != middleEdge->twin())
                    {
                        auto len = Traits::approximate_length(curr->source()->point(), curr->target()->point());
                        cumulativeApproxLength.push_back(len + (cumulativeApproxLength.empty() ? 0 : cumulativeApproxLength.back()));
                        curr = curr->next();
                    }
                }
                const auto totalLength = cumulativeApproxLength.back();

                // Keep track of chain of degree-2 vertices, delete them all when applicable.

                Traits::Halfedge_handle curr = middleEdge->twin()->next();
                Traits::Halfedge_handle chainStart = curr;
                std::size_t index = 0;
                while (curr != middleEdge->twin())
                {
                    if (cumulativeApproxLength[index] > totalLength*0.5)break;
                    curr = curr->next();
                    if (curr->source()->degree() > 2) chainStart = curr;
                    ++index;
                }
                while( curr->target()->degree() <= 2)
                {
                    curr = curr->next();
                    if (curr == middleEdge->twin()) throw std::runtime_error("Chain is invalid, central edge should have deg>2 vertices");
                }
                // Add edges before curr
                for (; chainStart != curr; chainStart = chainStart->next())
                {
                    ccwDelete.push_back(chainStart);
                }
                ccwDelete.push_back(curr);
            }
            // CW direction
            {
                std::vector<Traits::NT> cumulativeApproxLength;
                {
                    auto curr = middleEdge->next();
                    while (curr != middleEdge)
                    {
                        auto len = Traits::approximate_length(curr->source()->point(), curr->target()->point());
                        cumulativeApproxLength.push_back(len + (cumulativeApproxLength.empty() ? 0 : cumulativeApproxLength.back()));
                        curr = curr->next();
                    }
                }
                const auto totalLength = cumulativeApproxLength.back();

                Traits::Halfedge_handle curr = middleEdge->next();
                Traits::Halfedge_handle chainStart = curr;
                std::size_t index = 0;
                while (curr != middleEdge)
                {
                    if (cumulativeApproxLength[index] > totalLength*0.5)break;
                    curr = curr->next();
                    if (curr->source()->degree() > 2) chainStart = curr;
                    ++index;
                }
                while (curr->target()->degree() <= 2)
                {
                    curr = curr->next();
                    if (curr == middleEdge) throw std::runtime_error("Chain is invalid, central edge should have deg>2 vertices");
                }
                // Add edges before curr
                for (; chainStart != curr; chainStart = chainStart->next())
                {
                    cwDelete.push_back(chainStart);
                }
                cwDelete.push_back(curr);
            }
        }
        template<typename Observer = EdgeTrackingObserver>
        void notify_complement_reroute(const  std::vector<Traits::Halfedge_handle>& chain, Observer& observer)
        {
            {
                std::vector<std::size_t> rerouteEdges;
                auto start = chain.back()->next();
                if (chain.back()->face()->is_unbounded()) throw std::runtime_error("Face fucked up");
                while (start != chain.front())
                {
                    rerouteEdges.push_back(start->data().id.value());
                    start = start->next();
                }

                std::vector<std::size_t> deleteVertices;
                std::vector<std::size_t> deletedIds;
                for (auto el : chain)
                {
                    if (el->source() != chain.front()->source())
                    {
                        deleteVertices.push_back(el->source()->data().id);
                    }
                    deletedIds.push_back(el->data().id.value());
                }
                // We move in the other direction actually, so revert the ids (with appropriate ID's).
                SpecialEdgeIdFunctions::reverseCanonicalComplement(deletedIds);
                observer.handleEdgeReroute(deletedIds, rerouteEdges);
                // Delete intermediate vertices.
                for (auto v : deleteVertices)observer.handleVertexDelete(v);
            }
        }
    public:
        // Flag that we are going to try to remap faces. 
        bool recomputeAfterMerge() const
        {
            return true;
        }
        void set_verbose(bool val)
        {
            m_verbose = val;
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
            if (face->is_unbounded() || face->is_fictitious()) return;

            // Expect clean interior!
            using ListPolygon = std::vector<Traits::Point_2>;
            // Get face points
            std::vector<Traits::Point_2> points;
            std::vector<Traits::Halfedge_handle> hes;
            std::vector<std::size_t> degrees;
            std::size_t maxDegree = 0;
            Traits::for_each_ccb_edge(face, [&points, &hes, &degrees, &maxDegree](auto he, auto& cmd)
            {
                points.push_back(he->source()->point());
                hes.push_back(he);
                degrees.push_back(he->source()->degree());
                maxDegree = std::max<std::size_t>(maxDegree, he->source()->degree());
            });
            if (maxDegree == 2)
            {
                throw std::runtime_error("Invalid face, disconnected!");
            }
            { // Check if we only have 1 >degree 2  vertex
                std::size_t above2Count = 0;
                std::size_t index = 0;
                std::size_t highDegEl = 0;
                for (auto deg : degrees)
                {
                    if (deg > 2) {
                        above2Count++;
                        highDegEl = index;
                    }
                    ++index;
                }
                if (above2Count == 1)
                {
                    std::cout << "[BboxMerge] Collapsing face with single high degree vert \n";
                    // collapse to element
                    auto targetVert = hes[highDegEl]->source();
                    for (auto e : hes)
                    {
                        observer.handleEdgeCollapse(e->data().id.value(), targetVert->data().id);
                        if (e->source() != targetVert) observer.handleVertexMerge(e->source()->data().id, targetVert->data().id);
                        arrangement.remove_edge(e);
                    }
                    return;
                }
            }

            if (m_verbose)
            {
                std::cout << "[BboxMerge] Point count: " << points.size() << '\n';
                for (auto el : points)
                {
                    std::cout << el;
                    std::cout << '\n';
                }
            }
            //Compue CH
            ListPolygon CH;
            using QuadTraits = CGAL::Min_quadrilateral_default_traits_2<Traits::Kernel>;

            auto endIt = CGAL::ch_graham_andrew(points.begin(), points.end(), std::back_inserter(CH), CGAL::Convex_hull_traits_2<Traits::ArrTraits>{});
            if(m_verbose)
            {
                std::cout << "[BboxMerge] CH count: " << CH.size() << '\n';
                for (auto el : CH)
                {
                    std::cout << el;
                    std::cout << '\n';
                }
            }
            // Compute BBOX
            ListPolygon bbox;
            QuadTraits traits;
            auto bboxEnd = CGAL::min_rectangle_2(CH.begin(), CH.end(), std::back_inserter(bbox), traits);
            if(m_verbose)
            {
                std::cout << "[BboxMerge] bbox count: " << bbox.size() << '\n';
                for (auto el : bbox)
                {
                    std::cout << el;
                    std::cout << '\n';
                }
            }
            if (bbox.size() == 4)
            {
                bbox.push_back(bbox[0]);
            }
            for (auto& el : bbox)
            {
                el = Traits::reducedPrecision(el);
            }

            // Get major axis.
            Traits::Segment_2 axisSegment;
            std::size_t majorAxis, minorAxis;
            compute_splitter(bbox, axisSegment, majorAxis, minorAxis);

            std::vector<Traits::Point_2> intersections;
            std::vector<std::size_t> intersectionEdgeIndices;
            // Determine cuts through polygon //Note: special case if intersectin full segments!
            for (std::size_t i = 0; i < points.size(); ++i)
            {
                auto next = (i + 1) % points.size();
                auto polySeg = Traits::Make_segment_2()(points[i], points[next]);
                using Pnt = std::pair<Traits::ArrTraits::Point_2, Traits::ArrTraits::Multiplicity>;
                using Seg = Traits::ArrTraits::X_monotone_curve_2;
                typedef boost::variant < Pnt, Seg> Intersection_result;

                std::vector<Intersection_result> localIntersections;
                Traits::ArrTraits::Intersect_2 intersector = Traits::ArrTraits{}.intersect_2_object();
                intersector(axisSegment, polySeg, std::back_inserter(localIntersections));
                if (localIntersections.empty()) continue;
                for (const auto& inters : localIntersections)
                {
                    if (const auto* pnt = boost::get<Pnt>(&inters))
                    {
                        intersections.push_back(pnt->first);
                        intersectionEdgeIndices.push_back(i);
                    }
                    else if (const auto*seg = boost::get<Seg>(&inters))
                    {
                        intersections.push_back(seg->left());
                        intersections.push_back(seg->right());
                        intersectionEdgeIndices.push_back(i);
                        intersectionEdgeIndices.push_back(i);
                    }
                    else
                    {
                        throw std::runtime_error("No clue what this is");
                    }
                }
            }

            // Find edges to delete
            std::size_t edgeStart, edgeEnd;
            Traits::Point_2 intersectionStart, intersectionEnd;

            if (intersections.size() <= 1)
            {
                //err
                throw std::runtime_error("Expected at least 2 intersections");
            }
            else if (intersections.size() == 2)
            {
                //"Easy"
                edgeStart = intersectionEdgeIndices[0];
                intersectionStart = intersections[0];
                edgeEnd = intersectionEdgeIndices[1];
                intersectionEnd = intersections[1];
            }
            else
            {
                std::cout << "Selecting largest cut\n";
                Traits::ArrTraits::Line_2 line(bbox[minorAxis], bbox[minorAxis + 1]);
                Traits::ArrTraits::Compare_signed_distance_to_line_2 comparer;

                coordinated_sort(intersections, [&line, &comparer](const auto& pnt0, const auto& pnt1)
                {
                    return comparer(line, pnt0, pnt1);
                }, intersectionEdgeIndices);

                // Determine inside/outside.
                std::vector<CGAL::Bounded_side> midpointSides;
                // Pick largest
                Traits::NT maxLen = 0;
                std::optional<std::size_t> maxLenIndex;
                for (std::size_t i = 0; i < intersections.size() - 1; ++i)
                {
                    auto insideness = CGAL::bounded_side_2(points.begin(), points.end(), CGAL::midpoint(intersections[i], intersections[i + 1]), Traits::ArrTraits{});;
                    // If point is on boundary (i.e. segment intersection) or outside, this is not an interior segment
                    if (insideness == CGAL::ON_BOUNDARY || insideness == CGAL::ON_UNBOUNDED_SIDE) continue;
                    auto approxLen = Traits::approximate_length(intersections[i], intersections[i + 1]);
                    if (maxLen < approxLen)
                    {
                        maxLen = approxLen;
                        maxLenIndex = i;
                    }
                }
                if (!maxLenIndex.has_value()) throw std::runtime_error("All length zero segs?!?!");
                intersectionStart = intersections[maxLenIndex.value()];
                intersectionEnd = intersections[maxLenIndex.value() + 1];
                edgeStart = intersectionEdgeIndices[maxLenIndex.value()];
                edgeEnd = intersectionEdgeIndices[maxLenIndex.value() + 1];
            }

            if(m_verbose)
            {
                std::cout << "[BboxMerge] Start and end edges: " << edgeStart << ',' << edgeEnd << '\n';
                std::cout << "[BboxMerge] Total size: " << points.size() << '\n';
            }

            Traits::Halfedge_handle middleEdge = insert_bbox_segment(arrangement, hes[edgeStart], hes[edgeEnd],
                intersectionStart, intersectionEnd, observer);
            std::vector<Traits::Halfedge_handle> ccwDeleteHes, cwDeleteHes;
            find_edges_to_delete(middleEdge, cwDeleteHes, ccwDeleteHes);

            // Reroute by the 'complement' of the edges in the face
            notify_complement_reroute(cwDeleteHes, observer);
            notify_complement_reroute(ccwDeleteHes, observer);

            // Delete the edges in the arrangement
            for (auto el : ccwDeleteHes) arrangement.remove_edge(el);
            for (auto el : cwDeleteHes) arrangement.remove_edge(el);
        }
    };
}
#endif