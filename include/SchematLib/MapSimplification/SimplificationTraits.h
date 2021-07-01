#ifndef SCHEMATLIB_MAPSIMPLIFICATION_SIMPLIFICATIONTRAITS_H
#define SCHEMATLIB_MAPSIMPLIFICATION_SIMPLIFICATIONTRAITS_H
#include <queue>
#include "FaceMergeTraits.h"
#include <vector>
#include <GCpp/Math/Vector.h>

#include "SchematLib/MapSimplification/EdgeTrackingObserver.h"
namespace SchematLib::MapSimplification
{
    struct SimplificationTraits
    {
        // Setup CGAL types.
        using CGAL_NT = CGAL::Quotient<CGAL::MP_Float>;
        using Kernel = CGAL::Cartesian<CGAL_NT>;
        //using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
        using ArrTraits = CGAL::Arr_segment_traits_2<Kernel>;
        using Point_2 = typename ArrTraits::Point_2;
        using Segment_2 = typename ArrTraits::X_monotone_curve_2;

        // Custom data objects for the DCEL elements in the CGAL arrangement
        struct VertexData
        {
            // Assigned ID.
            std::size_t id = std::numeric_limits<std::size_t>::max();
            // Empty index: no vertex for it in the given Boost graph -> non-planarity vertex.
            std::optional<std::size_t> index{};
        };
        struct FaceData
        {
            std::optional<std::size_t> id{};
            // Marks that this faces is being deleted.
            bool beingDeleted = false;
            bool isInner = false;
            FaceData() {}
            FaceData(std::size_t fId) : id(fId) {}
        };
        struct EdgeData
        {
            std::optional<std::size_t> id{};

            friend std::ostream& operator<<(std::ostream& stream, const EdgeData& data)
            {
                stream << "EdgeData: id=";
                if (data.id.has_value()) stream << data.id.value();
                else stream << "unassigned";
                return stream;
            }
        };



        // Associate with the arrangement vertex and edges an ID, and with the faces the FaceData struct.
        using Arr = CGAL::Arrangement_with_history_2<ArrTraits, CGAL::Arr_extended_dcel<ArrTraits, VertexData, EdgeData, FaceData>>;
        using CGAL_Vertex_handle = typename Arr::Vertex_handle;
        using Face_iterator = typename Arr::Face_iterator;
        using Face_handle = Arr::Face_handle;
        using HE_handle = Arr::Halfedge_handle;
        using Vertex_handle = Arr::Vertex_handle;

        using PolygonTraits = Kernel;

        class IdGenerator
        {
            std::size_t m_currId;
        public:
            IdGenerator(std::size_t currId = 0) :m_currId(currId) {}
            std::size_t generateNext()
            {
                auto val = m_currId;
                ++m_currId;
                return val;
            }
        };
        struct GraphIdGenerator
        {
            IdGenerator faceIdGen;
            IdGenerator edgeIdGen;
            IdGenerator vertexIdGen;
        };

        class DuringRunChangeTracker : public CGAL::Arr_observer<Arr>
        {
            EdgeTrackingObserver* m_observer;
            GraphIdGenerator* m_idGenerators;
            std::optional<EdgeData> m_splitData;
        public:
            DuringRunChangeTracker(EdgeTrackingObserver* observer, GraphIdGenerator* idGenerators) :
                m_observer(observer),
                m_idGenerators(idGenerators)
            {}


            void after_create_vertex(Vertex_handle v) override
            {
                v->data().id = m_idGenerators->vertexIdGen.generateNext();
                m_observer->handleNewVertex(v->data().id);
            }

            void before_split_edge(Halfedge_handle he, Vertex_handle v, const X_monotone_curve_2& c0, const X_monotone_curve_2& c1) override
            {
                if (m_splitData.has_value()) throw std::runtime_error("Overwriting split data!");
                m_splitData = he->data();
            }

            void after_split_edge(Halfedge_handle he0, Halfedge_handle he1) override
            {
                if (!m_splitData.has_value()) throw std::runtime_error("No split data available!");
                // Generate ID's.
                he0->data().id = m_idGenerators->edgeIdGen.generateNext();
                he0->twin()->data().id = he0->data().id;
                he1->data().id = m_idGenerators->edgeIdGen.generateNext();
                he1->twin()->data().id = he1->data().id;
                // Notify observer
                m_observer->handleEdgeReplace(m_splitData.value().id.value(), { he0->data().id.value() ,he1->data().id.value() });
                // Unset split data
                m_splitData = {};
            }
        };

        static Point_2 gcppToCgalPosition(const GCpp::Math::Vec2<long double>& pos)
        {
            //return Point_2(pos.x(), pos.y());
            return Point_2(CGAL_NT(pos.x()), CGAL_NT(pos.y()));
        };

        static inline std::size_t getEdgeId(HE_handle handle)
        {
            return handle->data().id.value();
        }

        /**
         * \brief Reduce the precision of the point in the precise CGAL Kernel. Needed when computing new points, since otherwise
         * we can get runaway precision, bogging down all computations. We do this by intermediately casting the result to double
         * and reconstructing a point with these double values as input.
         * \param input The input point
         * \param output The output point with reduced precision
         */
        static void reducePrecision(const Point_2& input, Point_2& output) {
            output = Point_2(CGAL_NT(CGAL::to_double(input[0])), CGAL_NT(CGAL::to_double(input[1])));
        }

        template<typename EmbeddedGraph, typename SimplificationObserver>
        class SetupArrangement
        {
            using Vertex = typename EmbeddedGraph::vertex_descriptor;
            using Edge = typename EmbeddedGraph::edge_descriptor;
            /**
             * \brief Observer to assign edge IDs of the original graph to created (half)edges, without
             * interfering with the insertion process. Needs to take into account non-planarity that will split edges.
             * Notification of that is centralized in a later step.
             */
            class AssignIdsAtConstruction : public CGAL::Arr_observer<Arr>
            {
                const EmbeddedGraph& m_graph;
                std::size_t m_nextVertId;
                std::size_t m_nextEdgeId;
                // Original next-ids. Tells use which vertices/edges were already present.
                const std::size_t m_originalNextEdgeId;
                const std::size_t m_originalNextVertexId;
                std::size_t m_numFaces;
                SimplificationObserver* m_observer;

                std::size_t generateEdgeId()
                {
                    const auto val = m_nextEdgeId;
                    ++m_nextEdgeId;
                    return val;
                }
                std::size_t generateVertexId()
                {
                    const auto val = m_nextVertId;
                    ++m_nextVertId;
                    return val;
                }
                template<typename HeHandle>
                std::size_t setNewEdgeId(HeHandle handle)
                {
                    const auto eId = generateEdgeId();
                    handle->data().id = eId;
                    handle->twin()->data().id = eId;
                    return eId;
                }
            public:

                AssignIdsAtConstruction(const EmbeddedGraph& graph, SimplificationObserver* observer) :
                    m_graph(graph),
                    m_nextVertId(boost::get_property(m_graph, GCpp::DS::next_vertex_id_t{})),
                    m_nextEdgeId(boost::get_property(m_graph, GCpp::DS::next_edge_id_t{})),
                    m_originalNextVertexId(boost::get_property(m_graph, GCpp::DS::next_vertex_id_t{})),
                    m_originalNextEdgeId(boost::get_property(m_graph, GCpp::DS::next_edge_id_t{})),
                    m_numFaces(0),
                    m_observer(observer)
                {}

                std::size_t nextVertexId() const
                {
                    return m_nextVertId;
                }
                std::size_t nextEdgeId() const
                {
                    return m_nextEdgeId;
                }

                void after_create_vertex(typename ::CGAL::Arr_observer<Arr>::Vertex_handle vert) override
                {
                    if (vert->is_at_open_boundary()) return;
                    const auto vid = generateVertexId();
                    vert->set_data(VertexData{ vid,{} });
                    std::cout << "Vert at construction: " << vid << "\n";
                    m_observer->handleNewVertex(vid);
                }
                static auto lookup_edge(Vertex v0, Vertex v1, const EmbeddedGraph& graph)
                {
                    auto result = boost::lookup_edge(v0, v1, graph);
                    if (!result.second)
                    {
                        return boost::lookup_edge(v1, v0, graph);
                    }
                    return result;
                }
                void before_split_edge(typename ::CGAL::Arr_observer<Arr>::Halfedge_handle e, typename ::CGAL::Arr_observer<Arr>::Vertex_handle v,
                    const typename ::CGAL::Arr_observer<Arr>::X_monotone_curve_2 &c1, const typename ::CGAL::Arr_observer<Arr>::X_monotone_curve_2 &c2) override
                {
                    std::cout << "[Facemerger] Splitting edge " << e->data() << " with vert " << v->data().id << '\n';
                }
                void after_split_edge(typename ::CGAL::Arr_observer<Arr>::Halfedge_handle he0, typename ::CGAL::Arr_observer<Arr>::Halfedge_handle he1) override
                {


                    // Assign IDs if not present, but also assign a new ID for an original edge (otherwise, the mapping gets messed up).
                    if (!he0->data().id.has_value() || he0->data().id.value() < m_originalNextEdgeId)
                    {
                        const auto eid = setNewEdgeId(he0);
                        m_observer->handleNewEdge(eid, he0->source()->data().id, he0->target()->data().id);
                    }
                    if (!he1->data().id.has_value() || he1->data().id.value() < m_originalNextEdgeId)
                    {
                        const auto eid = setNewEdgeId(he1);
                        m_observer->handleNewEdge(eid, he1->source()->data().id, he1->target()->data().id);
                    }
                    std::cout << "[Facemerger] Split edge to  " << he0->data() << ", " << he1->data() << ", vs= "
                        << he0->source()->data().id << ',' << he0->target()->data().id << " | "
                        << he1->source()->data().id << ',' << he1->target()->data().id
                        << '\n';
                }
                void after_create_edge(typename ::CGAL::Arr_observer<Arr>::Halfedge_handle he) override
                {
                    if (he->is_fictitious()) return;
                    if (he->source()->is_at_open_boundary() || he->target()->is_at_open_boundary()) return;
                    const FaceMergeTraits::VertexData srcVertData = he->source()->data();
                    const FaceMergeTraits::VertexData targetVertData = he->target()->data();

                    // One of the two is newly added (i.e. vertex due to non-planarity) Happensw hen it has no
                    // index in the original graph
                    if (!srcVertData.index.has_value() || !targetVertData.index.has_value())
                    {
                        const auto eid = setNewEdgeId(he);
                        std::cout << "[FaceMerger] New edge added! " << eid << " with vs="
                            << he->source()->data().id << ',' << he->target()->data().id <<
                            std::endl;
                        m_observer->handleNewEdge(eid, srcVertData.id, targetVertData.id);
                        return;
                    }

                    auto[edge, exists] = lookup_edge(srcVertData.index.value(), targetVertData.index.value(), m_graph);
                    if (exists)
                    {
                        const auto eid = boost::get(GCpp::DS::edge_id_t{}, m_graph, edge);
                        he->data().id = eid;
                        he->twin()->data().id = eid;
                    }
                    else //Problem: this effectively splits an original edge into two pieces...
                    {
                        const auto eid = setNewEdgeId(he);
                        std::cout << "[FaceMerger] Edge got split!" << std::endl;
                        m_observer->handleNewEdge(eid, srcVertData.id, targetVertData.id);

                    }
                }

                void after_split_face(typename ::CGAL::Arr_observer<Arr>::Face_handle origAndNewFace,
                    typename ::CGAL::Arr_observer<Arr>::Face_handle newFace, bool is_hole) override
                {
                    newFace->data().id = m_numFaces;
                    ++m_numFaces;
                }
            };
            void notifyNonPlanaritySplits(const EmbeddedGraph& inputGraph, const Arr& arr, SimplificationObserver& observer) const
            {
                // Find full paths.
                for (auto it = arr.curves_begin(); it != arr.curves_end(); ++it)
                {
                    if (arr.number_of_induced_edges(it) == 1) continue;
                    std::cout << "Number of non-planar edges for curve: " << arr.number_of_induced_edges(it) << '\n';

                    std::set<std::size_t> vertexIndices;
                    {
                        for (auto testIt = arr.induced_edges_begin(it); testIt != arr.induced_edges_end(it); ++testIt)
                        {
                            auto heHandle = *testIt;
                            if (heHandle->source()->data().index.has_value()) vertexIndices.insert(heHandle->source()->data().index.value());
                            if (heHandle->target()->data().index.has_value()) vertexIndices.insert(heHandle->target()->data().index.value());
                        }
                    }
                    if (vertexIndices.size() != 2)
                    {
                        throw std::runtime_error("Malformed path");
                    }

                    std::size_t v0Index = *vertexIndices.begin();
                    std::size_t v1Index = *std::next(vertexIndices.begin());
                    auto[edge, exists] = boost::lookup_edge(v0Index, v1Index, inputGraph);
                    if (!exists)
                    {
                        std::tie(edge, exists) = boost::lookup_edge(v1Index, v0Index, inputGraph);
                        if (!exists) throw std::runtime_error("Could not deduce original edge from split edges");
                    }
                    std::vector<std::size_t> replacementEdges;
                    GCpp::Helpers::Iterators::for_each_it(arr.induced_edges_begin(it), arr.induced_edges_end(it), [&replacementEdges](auto eIt)
                    {
                        replacementEdges.push_back((*eIt)->data().id.value());
                    });
                    const auto origId = boost::get(GCpp::DS::edge_id_t{}, inputGraph, edge);
                    std::cout << "Replacing: " << origId << " with ";
                    for (auto edgeId : replacementEdges) std::cout << edgeId << " ";
                    std::cout << '\n';
                    // Assign the path to the edge
                    observer.handleEdgeReplace(origId, replacementEdges);
                }
            }
        public:

            void operator()(const EmbeddedGraph& inputGraph, Arr& arrangement, SimplificationObserver& observer, std::size_t& nextEdgeId, std::size_t& nextVertId, std::size_t& nextFaceId) const
            {
                namespace it = GCpp::Helpers::Iterators;
                auto cgalVertexPos = [&inputGraph](auto graphVertex)
                {
                    const auto pos = GCpp::DS::get_vertex_location(graphVertex, inputGraph);
                    //return Point_2(pos.x(), pos.y());
                    return Point_2(CGAL_NT(pos.x()), CGAL_NT(pos.y()));
                };
                // Add vertices. This way, we can easily set additional data.
                std::vector<CGAL_Vertex_handle> vertices;
                for (auto vertex : it::range(boost::vertices(inputGraph)))
                {
                    auto cgalVert = CGAL::insert_point(arrangement, cgalVertexPos(vertex));
                    // 
                    cgalVert->set_data(VertexData{ boost::get(GCpp::DS::vertex_id_t{}, inputGraph, vertex), vertices.size() });
                    vertices.push_back(cgalVert);
                }

                // Assigns edge and vertex ids for new elements in the arrangment
                AssignIdsAtConstruction assigner(inputGraph, &observer);
                assigner.attach(arrangement);

                //Add edges
                for (auto edge : it::range(boost::edges(inputGraph)))
                {
                    const auto src = boost::source(edge, inputGraph);
                    const auto target = boost::target(edge, inputGraph);
                    if (src >= vertices.size() || target >= vertices.size()) throw std::runtime_error("Vertex out of bounds");
                    auto srcPoint = cgalVertexPos(src);
                    auto targetPoint = cgalVertexPos(target);
                    Segment_2 seg(srcPoint, targetPoint);
                    // Observer AssignIdsAtConstruction will handle ID assignment.
                    CGAL::insert(arrangement, seg);
                }
                assigner.detach();

                nextEdgeId = assigner.nextEdgeId();
                nextVertId = assigner.nextVertexId();

                // Assign face IDs and meta data.
                std::size_t currentFace = 0;
                auto assignFaceIfUnassigned = [&currentFace](auto faceHandle)
                {
                    if (!faceHandle->has_outer_ccb()) return;
                    if (!faceHandle->data().id.has_value())
                    {
                        faceHandle->data().id = currentFace;
                        ++currentFace;
                    }
                };
                // Assign face IDs
                for (auto eit = arrangement.edges_begin(); eit != arrangement.edges_end(); ++eit)
                {
                    assignFaceIfUnassigned(eit->face());
                    assignFaceIfUnassigned(eit->twin()->face());
                }
                for (auto fit = arrangement.faces_begin(); fit != arrangement.faces_end(); ++fit)
                {
                    if (!fit->has_outer_ccb()) continue;

                    typename Arr::Hole_iterator holeIt;

                    for (holeIt = fit->holes_begin(); holeIt != fit->holes_end(); ++holeIt)
                    {
                        std::queue<typename Arr::Halfedge_iterator> toProcess;
                        std::set<std::size_t> seenIds;
                        seenIds.insert(fit->data().id.value());

                        // Runs over outer boundary of hole
                        typename Arr::Ccb_halfedge_circulator eIt = *holeIt;
                        do
                        {
                            auto twinFaceHandle = eIt->twin()->face();
                            if (twinFaceHandle == fit || seenIds.find(twinFaceHandle->data().id.value()) != seenIds.end())
                            {
                                ++eIt;
                                continue;
                            }
                            // Not yet seen face inside another face.
                            seenIds.insert(twinFaceHandle->data().id.value());
                            toProcess.push(eIt->twin());

                            ++eIt;
                        } while (eIt != *holeIt);
                        while (!toProcess.empty())
                        {
                            auto el = toProcess.front();
                            el->face()->data().isInner = true;
                            toProcess.pop();
                            auto ccb_it = el->ccb();
                            do
                            {
                                auto twinFaceHandle = eIt->twin()->face();
                                if (twinFaceHandle == fit || seenIds.find(twinFaceHandle->data().id.value()) != seenIds.end())
                                {
                                    ++eIt;
                                    continue;
                                }
                                // Not yet seen face inside another face.
                                seenIds.insert(twinFaceHandle->data().id.value());
                                toProcess.push(eIt->twin());
                                ++ccb_it;
                            } while (ccb_it != el->ccb());
                        }
                    }
                }

                // Notify of splits due to non-planarity.
                notifyNonPlanaritySplits(inputGraph, arrangement, observer);
            }
        };

    };
}
#endif
