#ifndef SCHEMATLIB_MAPSIMPLIFICATION_FACEMERGETRAITS_H
#define SCHEMATLIB_MAPSIMPLIFICATION_FACEMERGETRAITS_H
#ifdef foreach
#undef foreach
#pragma message("Undefining foreach")
#endif

#ifdef foreach
#pragma error("foreach should not be here!")
#endif
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h> // Simple kernel, replace when needed.
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <vector>
#include <GCpp/Math/Vector.h>
#include "SchematLib/MapSimplification/EdgeTrackingObserver.h"
namespace SchematLib::MapSimplification
{
#define SCHEMAT_MODIFIED_ARR
    // TODO: template this at some point
    struct FaceMergeTraits
    {
        using NT = long double;
        // Setup CGAL types.
        using CGAL_NT = CGAL::Quotient<CGAL::MP_Float>;
        using Kernel = CGAL::Cartesian<CGAL_NT>;
        


        // Custom data objects for the DCEL elements in the CGAL arrangement
        struct VertexData
        {
            // Assigned ID.
            std::size_t id = std::numeric_limits<std::size_t>::max();
            // Empty index: no vertex for it in the given Boost graph -> non-planarity vertex.
            std::optional<std::size_t> index;
        };
        struct FaceData
        {
            std::optional<std::size_t> id;
            // TODO: delete these flags.
            bool beingDeleted = false;
            bool isInner = false;
            FaceData() {}
            FaceData(std::size_t fId) : id(fId) {}
        };
        struct EdgeData
        {
            std::optional<std::size_t> id;

            friend std::ostream& operator<<(std::ostream& stream, const EdgeData& data)
            {
                stream << "EdgeData: id=";
                if (data.id.has_value()) stream << data.id.value();
                else stream << "unassigned";
                return stream;
            }
        };

#ifdef SCHEMAT_MODIFIED_ARR
        using ArrTraits = CGAL::Arr_segment_traits_2<Kernel>;
        // Add ID to curve
        using AddedCurveTraits = CGAL::Arr_curve_data_traits_2<ArrTraits, std::size_t>;
        using Point_2 = typename ArrTraits::Point_2;
        using Segment_2 = typename AddedCurveTraits::Curve_2;

        using Arr = CGAL::Arrangement_2<AddedCurveTraits, CGAL::Arr_extended_dcel<AddedCurveTraits, VertexData, EdgeData, FaceData>>;
#else
        //using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
        using ArrTraits = CGAL::Arr_segment_traits_2<Kernel>;

        // Add extra data to curve: the ID of the edge associated with it, if applicable.
        //using AddedCurveTraits = CGAL::Arr_curve_data_traits_2<ArrTraits, std::size_t>;
        //using Point_2 = typename AddedCurveTraits::Point_2;
        //using Segment_2 = typename AddedCurveTraits::Curve_2;

        using Point_2 = typename ArrTraits::Point_2;
        using Segment_2 = typename ArrTraits::X_monotone_curve_2;
        using Arr = CGAL::Arrangement_2<ArrTraits, CGAL::Arr_extended_dcel<ArrTraits, VertexData, EdgeData, FaceData>>;
#endif

        using CGAL_Vertex_handle = typename Arr::Vertex_handle;
        using Face_iterator = typename Arr::Face_iterator;
        using Face_handle = Arr::Face_handle;
        using HE_handle = Arr::Halfedge_handle;
        using Halfedge_handle = Arr::Halfedge_handle;
        using Vertex_handle = Arr::Vertex_handle;

        using PolygonTraits = Kernel;

        class IdGenerator
        {
            std::size_t m_currId;
            std::size_t m_startId;
        public:
            class Mark
            {
                friend class IdGenerator;
                std::size_t m_nextId;
                Mark(std::size_t nextId):m_nextId(nextId){}
            public:
                bool isNewIdComparedToThis(std::size_t id)const
                {
                    return id >= m_nextId;
                }
            };
            Mark createMark() const
            {
                return Mark(m_currId);
            }
            IdGenerator(std::size_t currId = 0) :m_currId(currId),m_startId(currId) {}
            /**
             * \brief Returns whether this generator is at 0 for the current id.
             * \return 
             */
            bool isAtZero() const
            {
                return m_currId == 0;
            }
            bool isAtInitial() const
            {
                return m_currId == m_startId;
            }
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

        using GcppVec2 = GCpp::Math::Vec2<NT>;

        static Point_2 gcppToCgalPosition(const GCpp::Math::Vec2<long double>& pos)
        {
            //return Point_2(pos.x(), pos.y());
            return Point_2(CGAL_NT(pos.x()), CGAL_NT(pos.y()));
        }
        static GCpp::Math::Vec2<NT> cgalPositionToGcpp(const Point_2& pos)
        {
            //return Point_2(pos.x(), pos.y());
            return GCpp::Math::Vec2<NT>(CGAL::to_double(pos.x()), CGAL::to_double(pos.y()));
        };
        static GCpp::Math::Vec2<NT> cgalPositionToGcpp(const Kernel::Vector_2& pos)
        {
            //return Point_2(pos.x(), pos.y());
            return GCpp::Math::Vec2<NT>(CGAL::to_double(pos.x()), CGAL::to_double(pos.y()));
        };

        static NT approximate_length(const Point_2& pos, const Point_2& pos2)
        {
            return (cgalPositionToGcpp(pos) - cgalPositionToGcpp(pos2)).length();
        }
        static NT approximate_length(const Kernel::Vector_2& vec)
        {
            return (cgalPositionToGcpp(vec)).length();
        }

        template<typename HEHandle>
        static inline bool isCanonical(HEHandle he)
        {
            return he->source()->data().id < he->target()->data().id;
        }
        template<typename HEHandle>
        static inline HEHandle getCanonical(HEHandle he)
        {
            if (isCanonical(he)) return he;
            return he->twin();
        }

        static inline std::size_t getEdgeId(HE_handle handle)
        {
            return handle->data().id.value();
        }

        static NT approximate_angle(const Kernel::Vector_2& p0, const Kernel::Vector_2& p1)
        {
            auto vec0 = cgalPositionToGcpp(p0);
            auto vec1 = cgalPositionToGcpp(p1);
            return std::acos(std::clamp<NT>(vec0 * vec1 / (vec0.length() * vec1.length()),0,1));
        }
        static NT approximate_cos_angle(const Kernel::Vector_2& p0, const Kernel::Vector_2& p1)
        {
            auto vec0 = cgalPositionToGcpp(p0);
            auto vec1 = cgalPositionToGcpp(p1);
            return vec0 * vec1 / (vec0.length() * vec1.length());
        }

        template<typename Arrangement>
        static void computeEdgeIdSet(const Arrangement& arrangement, std::set<std::size_t>& edgeIds, bool includeCanonical=true)
        {
            for (const auto& el : arrangement.edge_handles())
            {
                edgeIds.insert(el->data().id.value());
                if(includeCanonical)edgeIds.insert(el->twin()->data().id.value());
            }
        }
        template<typename Arrangement>
        static void computeVertexIdSet(const Arrangement& arrangement, std::set<std::size_t>& vertexIds)
        {
            for (const auto& el : arrangement.vertex_handles())
            {
                vertexIds.insert(el->data().id);
            }
        }
        template<typename Arrangement>
        static void computeFaceIdSet(const Arrangement& arrangement, std::set<std::size_t>& faceIds)
        {
            for (const auto& el : arrangement.face_handles())
            {
                faceIds.insert(el->data().id.value());
            }
            faceIds.insert(arrangement.unbounded_face()->data().id.value());
        }

        template<typename HalfedgeHandle>
        static void setEdgeId(HalfedgeHandle he, std::size_t id)
        {
            he->data().id = id;
            he->twin()->data().id = id;
        }

        static void computeBoundingBox(const Arr& arrangement, std::array<Point_2, 2>& bboxPoints)
        {
            Kernel::Compare_x_2 compareX;
            Kernel::Compare_y_2 compareY;
            CGAL_NT minX, maxX, minY, maxY;
            bool first = true;
            for (auto v : arrangement.vertex_handles())
            {
                if (first)
                {
                    minX = v->point().x();
                    maxX = v->point().x();
                    minY = v->point().y();
                    maxY = v->point().y();
                    first = false;
                }
                else
                {
                    if (minX > v->point().x()) minX = v->point().x();
                    if (minY > v->point().y()) minY = v->point().y();
                    if (maxX < v->point().x()) maxX = v->point().x();
                    if (maxY < v->point().y()) maxY = v->point().y();
                }
            }
            bboxPoints = { Point_2(minX,minY), Point_2(maxX, maxY) };
        }

        static void getEdgeIdToVerticesMapping(const Arr& inputGraph, std::map<std::size_t, std::pair<std::size_t, std::size_t>>& mapping)
        {
            auto getVId = [&inputGraph](auto v)
            {
                return v->data().id;
            };
            auto getEId = [&inputGraph](auto e)
            {
                return e->data().id.value();
            };
            for (auto e : inputGraph.edge_handles())
            {
                const auto srcId = getVId(e->source());
                const auto targetId = getVId(e->target());
                mapping[getEId(e)] = std::make_pair(std::min(srcId, targetId), std::max(srcId, targetId));
            }
        }

        template<typename Arrangement>
        static inline void verifyIds(const Arrangement& arrangement, bool canonicalEdgeIds=true)
        {
            // Vertices
            {
                std::set<std::size_t> vIds;
                for (const auto& v : arrangement.vertex_handles())
                {
                    const auto vId = v->data().id;
                    if (vIds.find(vId) != vIds.end())
                    {
                        std::cout << "Duplicate vertex ID " << vId << '\n';
                        throw std::runtime_error("Duplicate vert ID");
                    }
                    vIds.insert(v->data().id);
                }
            }
            // Edges
            {
                std::set<std::size_t> eIds;
                std::size_t eCount = 0;
                for (const auto& e : arrangement.edge_handles())
                {
                    if (!e->data().id.has_value()) throw std::runtime_error("Edge without ID");
                    if (!e->twin()->data().id.has_value())throw std::runtime_error("Edge without ID");
                    if(canonicalEdgeIds)
                    {
                        if (SpecialEdgeIdFunctions::flipCanonical(e->data().id.value()) != e->twin()->data().id.value()) throw std::runtime_error("Incongruency edge IDs");
                    }
                    else
                    {
                        if (e->data().id.value() != e->twin()->data().id.value()) throw std::runtime_error("Incongruency edge IDs");
                    }
                    
                    const auto eId = e->data().id.value();
                    if (eIds.find(eId) != eIds.end()) throw std::runtime_error("Duplicate edge");
                    eIds.insert(eId);
                    ++eCount;
                }
            }
            // Faces
            {
                if (!arrangement.unbounded_face()->data().id.has_value()) throw std::runtime_error("Unbounded face has no ID");
                if (arrangement.unbounded_face()->data().id.value() != 0) throw std::runtime_error("Invalid unbounded face ID");
                std::set<std::size_t> fIds;
                for (const auto& f : arrangement.face_handles())
                {
                    if (f->is_unbounded()) continue;
                    if (!f->data().id.has_value()) throw std::runtime_error("Face without id");
                    if (f->data().id.value() == 0) throw std::runtime_error("Invalid face ID, should not be zero");
                    auto[it, wasInserted] = fIds.insert(f->data().id.value());
                    if (!wasInserted) throw std::runtime_error("Duplicate ID");
                }
            }
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
        static Point_2 reducedPrecision(const Point_2& input) {
            return Point_2(CGAL_NT(CGAL::to_double(input[0])), CGAL_NT(CGAL::to_double(input[1])));
        }

        enum class IterationCommand
        {
            Break,
            Continue,
            Noop
        };

        struct ItCommand
        {
            bool m_break = false;
            bool m_continue = false;

            void reset() { m_break = false; m_continue = false; }
            void doBreak() { m_break = true; }
            void doContinue() { m_continue = true; }
            bool shouldBreak()const { return m_break; }
            bool shouldContinue()const { return m_continue; }
        };

        template<typename VertexHandle, typename Func>
        static void for_each_incident_edge(VertexHandle handle, Func&& func)
        {
            auto circ = handle->incident_halfedges();
            ItCommand cmd;
            do
            {
                func(circ, cmd);
                if (cmd.shouldBreak())break;

                ++circ;
            } while (circ != handle->incident_halfedges());
        }

        template<typename FaceHandle, typename Func>
        static void for_each_ccb_edge(FaceHandle handle, Func&& func)
        {
            auto circ = handle->outer_ccb();
            ItCommand cmd;
            do
            {
                func(circ, cmd);
                if (cmd.shouldBreak())break;

                ++circ;
            } while (circ != handle->outer_ccb());
        }

        struct ComputeComplexity
        {
            template<typename FaceHandle>
            std::size_t operator()(FaceHandle handle) const
            {
                std::size_t counter = 0;
                for_each_ccb_edge(handle,[&counter](auto handle, auto& cmd)
                {
                    ++counter;
                });
                return counter;
            }
        };

        template<typename Arrangement>
        static void verifyIntegrity(const Arrangement& arrangement, const EdgeTrackingObserver& obs)
        {
            std::set<std::size_t> edgeIdSet, vertexIdSet;
            computeEdgeIdSet(arrangement, edgeIdSet);
            computeVertexIdSet(arrangement, vertexIdSet);
            std::cout << "[VerifyIntegrity] Arrangement: #V=" << arrangement.number_of_vertices() << ",#E="
                << arrangement.number_of_edges() << ", #Vset=" << vertexIdSet.size() << ",#Eset=" << edgeIdSet.size() << '\n';
            obs.verifyIntegrity(edgeIdSet, vertexIdSet);
            // TODO: faces
        }

        //template<typename EmbeddedGraph, typename SimplificationObserver>
        //class SetupArrangement
        //{
        //    using Vertex = typename EmbeddedGraph::vertex_descriptor;
        //    using Edge = typename EmbeddedGraph::edge_descriptor;
        //    /**
        //     * \brief Observer to assign edge IDs of the original graph to created (half)edges, without
        //     * interfering with the insertion process. Needs to take into account non-planarity that will split edges.
        //     * Notification of that is centralized in a later step.
        //     */
        //    class AssignIdsAtConstruction : public CGAL::Arr_observer<Arr>
        //    {
        //        const EmbeddedGraph& m_graph;
        //        // Original next-ids. Tells use which vertices/edges were already present.
        //        SimplificationObserver* m_observer;

        //        template<typename HeHandle>
        //        std::size_t setNewEdgeId(HeHandle handle)
        //        {
        //            const auto eId = m_idGenerator->edgeIdGen.generateNext();
        //            handle->data().id = eId;
        //            handle->twin()->data().id = eId;
        //            return eId;
        //        }

        //        GraphIdGenerator* m_idGenerator;

        //        std::optional<std::size_t> insertingEdge;

        //        std::optional<std::size_t> splittingEdgeId;
        //        bool wasSplit = false;

        //        std::vector<std::size_t> insertedSplitIds;
        //    public:

        //        AssignIdsAtConstruction(const EmbeddedGraph& graph, SimplificationObserver* observer, GraphIdGenerator* idGenerator) :
        //            m_graph(graph),
        //            m_idGenerator(idGenerator),
        //            m_observer(observer)
        //        {}

        //        void setCurrentlyInsertingEdge(std::size_t edgeId)
        //        {
        //            if (!insertedSplitIds.empty())
        //            {
        //                m_observer->handleEdgeReroute(std::vector<std::size_t>{ insertingEdge.value() }, insertedSplitIds);
        //            }
        //            insertingEdge = edgeId;
        //            wasSplit = false;
        //            insertedSplitIds.clear();
        //        }
        //        void clearInsertingEdge()
        //        {
        //            insertingEdge = {};
        //        }


        //        void after_create_vertex(typename ::CGAL::Arr_observer<Arr>::Vertex_handle vert) override
        //        {
        //            // This should always be a non-planarity vertex.
        //            if (vert->is_at_open_boundary()) return;
        //            //std::cout << "[AfterCreateVertex]\n";
        //            const auto vid = m_idGenerator->vertexIdGen.generateNext();
        //            vert->set_data(VertexData{ vid,{} });
        //            //std::cout << "Vert at construction: " << vid << "\n";
        //            m_observer->handleNewVertex(vid);
        //        }
        //        static auto lookup_edge(Vertex v0, Vertex v1, const EmbeddedGraph& graph)
        //        {
        //            auto result = boost::lookup_edge(v0, v1, graph);
        //            if (!result.second)
        //            {
        //                return boost::lookup_edge(v1, v0, graph);
        //            }
        //            return result;
        //        }
        //        // When inserting an edge, splits should only be called on already present edges.
        //        void before_split_edge(typename ::CGAL::Arr_observer<Arr>::Halfedge_handle e, typename ::CGAL::Arr_observer<Arr>::Vertex_handle v,
        //            const typename ::CGAL::Arr_observer<Arr>::X_monotone_curve_2 &c1, const typename ::CGAL::Arr_observer<Arr>::X_monotone_curve_2 &c2) override
        //        {
        //            splittingEdgeId = e->data().id.value();
        //            wasSplit = true;
        //            //std::cout << "[BeforeSplitEdge]\n";
        //            //std::cout << "[Facemerger] Splitting edge " << e->data() << " with vert " << v->data().id << '\n';
        //        }
        //        void after_split_edge(typename ::CGAL::Arr_observer<Arr>::Halfedge_handle he0, typename ::CGAL::Arr_observer<Arr>::Halfedge_handle he1) override
        //        {
        //            if (!splittingEdgeId.has_value()) throw std::runtime_error("Splitting edge data was not set");
        //            //std::cout << "[AfterSplitEdge]\n";
        //            // Assign IDs if not present, but also assign a new ID for an original edge (otherwise, the mapping gets messed up).
        //            std::size_t eId0 = setNewEdgeId(he0);
        //            std::size_t eId1 = setNewEdgeId(he1);
        //            m_observer->handleNewEdge(eId0, he0->source()->data().id, he0->target()->data().id);
        //            m_observer->handleNewEdge(eId1, he1->source()->data().id, he1->target()->data().id);
        //            /*std::cout << "[Facemerger] Split edge to  " << he0->data() << ", " << he1->data() << ", vs= "
        //                << he0->source()->data().id << ',' << he0->target()->data().id << " | "
        //                << he1->source()->data().id << ',' << he1->target()->data().id
        //                << '\n';*/
        //            m_observer->handleEdgeReroute({ splittingEdgeId.value() }, { eId0, eId1 });
        //            splittingEdgeId = {};
        //        }
        //        void after_create_edge(typename ::CGAL::Arr_observer<Arr>::Halfedge_handle he) override
        //        {
        //            //std::cout << "[AfterCreateEdge]\n";
        //            if (he->is_fictitious()) return;
        //            if (he->source()->is_at_open_boundary() || he->target()->is_at_open_boundary()) return;
        //            const FaceMergeTraits::VertexData srcVertData = he->source()->data();
        //            const FaceMergeTraits::VertexData targetVertData = he->target()->data();

        //            // One of the two is newly added (i.e. vertex due to non-planarity) Happensw hen it has no
        //            // index in the original graph
        //            //if (!srcVertData.index.has_value() || !targetVertData.index.has_value())
        //            if (wasSplit)
        //            {
        //                const auto eid = setNewEdgeId(he);
        //                /*std::cout << "[FaceMerger] New edge added! " << eid << " with vs="
        //                    << he->source()->data().id << ',' << he->target()->data().id <<
        //                    std::endl;*/
        //                    //if (!wasSplit) throw std::runtime_error("Expecting split");
        //                insertedSplitIds.push_back(eid);
        //                m_observer->handleNewEdge(eid, srcVertData.id, targetVertData.id);
        //                return;
        //            }
        //            auto[edge, exists] = lookup_edge(srcVertData.index.value(), targetVertData.index.value(), m_graph);
        //            if (exists)
        //            {
        //                const auto eid = boost::get(GCpp::DS::edge_id_t{}, m_graph, edge);
        //                he->data().id = eid;
        //                he->twin()->data().id = eid;
        //            }
        //            else //Problem: this effectively splits an original edge into two pieces...
        //            {
        //                const auto eid = setNewEdgeId(he);
        //            }
        //            m_observer->handleNewEdge(he->data().id.value(), srcVertData.id, targetVertData.id);
        //        }

        //        void after_split_face(typename ::CGAL::Arr_observer<Arr>::Face_handle origAndNewFace,
        //            typename ::CGAL::Arr_observer<Arr>::Face_handle newFace, bool is_hole) override
        //        {
        //            newFace->data().id = m_idGenerator->faceIdGen.generateNext();
        //        }
        //    };

        //public:

        //    void operator()(const EmbeddedGraph& inputGraph, Arr& arrangement, SimplificationObserver& observer, GraphIdGenerator& idGenerator) const
        //    {
        //        // Unbounded face gets the zero index.
        //        arrangement.unbounded_face()->data().id = 0;
        //        idGenerator.faceIdGen = IdGenerator(1);
        //        idGenerator.vertexIdGen = IdGenerator(boost::get_property(inputGraph, GCpp::DS::next_vertex_id_t{}));
        //        idGenerator.edgeIdGen = IdGenerator(boost::get_property(inputGraph, GCpp::DS::next_edge_id_t{}));

        //        namespace it = GCpp::Helpers::Iterators;
        //        auto cgalVertexPos = [&inputGraph](auto graphVertex)
        //        {
        //            const auto pos = GCpp::DS::get_vertex_location(graphVertex, inputGraph);
        //            //return Point_2(pos.x(), pos.y());
        //            return Point_2(CGAL_NT(pos.x()), CGAL_NT(pos.y()));
        //        };
        //        // Add vertices. This way, we can easily set additional data.
        //        std::vector<CGAL_Vertex_handle> vertices;
        //        for (auto vertex : it::range(boost::vertices(inputGraph)))
        //        {
        //            auto cgalVert = CGAL::insert_point(arrangement, cgalVertexPos(vertex));
        //            //
        //            const auto vId = GCpp::DS::getVertexId(inputGraph, vertex);
        //            cgalVert->set_data(FaceMergeTraits::VertexData{ vId, vertices.size() });
        //            observer.handleNewVertex(vId);
        //            vertices.push_back(cgalVert);
        //        }

        //        // Assigns edge and vertex ids for new elements in the arrangment
        //        AssignIdsAtConstruction assigner(inputGraph, &observer, &idGenerator);
        //        assigner.attach(arrangement);

        //        //Add edges
        //        std::size_t edgeCount = 0;
        //        std::size_t edgeTotal = boost::num_edges(inputGraph);
        //        const std::size_t percInc = 5;
        //        std::size_t nextPerc = percInc;
        //        for (auto edge : it::range(boost::edges(inputGraph)))
        //        {
        //            const auto eId = GCpp::DS::getEdgeId(inputGraph, edge);
        //            const auto src = boost::source(edge, inputGraph);
        //            const auto target = boost::target(edge, inputGraph);
        //            if (src >= vertices.size() || target >= vertices.size()) throw std::runtime_error("Vertex out of bounds");
        //            auto srcPoint = cgalVertexPos(src);
        //            auto targetPoint = cgalVertexPos(target);
        //            Segment_2 seg(srcPoint, targetPoint);
        //            assigner.setCurrentlyInsertingEdge(eId);
        //            // Observer AssignIdsAtConstruction will handle ID assignment.
        //            //std::cout << "[SetupArr] Creating edge " << edgeCount << '\n';
        //            CGAL::insert(arrangement, seg);
        //            ++edgeCount;
        //            if (edgeCount * 100 > nextPerc * edgeTotal)
        //            {
        //                std::cout << "At " << nextPerc << "% of edges\n";
        //                nextPerc += percInc;
        //            }
        //        }
        //        assigner.detach();

        //        // Assign face IDs and meta data.
        //        auto assignFaceIfUnassigned = [&idGenerator, &observer](auto faceHandle)
        //        {
        //            if (!faceHandle->has_outer_ccb()) return;
        //            if (!faceHandle->data().id.has_value())
        //            {
        //                faceHandle->data().id = idGenerator.faceIdGen.generateNext();
        //                observer.handleNewFace(faceHandle->data().id.value());
        //            }
        //        };
        //        for (auto eit = arrangement.edges_begin(); eit != arrangement.edges_end(); ++eit)
        //        {
        //            assignFaceIfUnassigned(eit->face());
        //            assignFaceIfUnassigned(eit->twin()->face());
        //        }
        //        FaceMergeTraits::verifyIds(arrangement);
        //        verifyIntegrity(arrangement, observer);
        //    }
        //};
        struct Make_segment_2
        {
            Segment_2 operator()(const Point_2& p0, const Point_2& p1) const
            {
                return Segment_2(CGAL::Arr_segment_2<Kernel>(p0, p1));
            }
        };

        struct SortVertexByPoint
        {
            bool operator()(Vertex_handle v0, Vertex_handle v1) const
            {
                return v0->point() < v1->point();
            }
        };

        template<typename EmbeddedGraph>
        class SetupModifiedArrangement
        {
        public:
            using SimplificationObserver = EdgeTrackingObserver;
            using Traits = FaceMergeTraits;
        private:
            using CGAL_NT = Traits::CGAL_NT;
            using Vertex = typename EmbeddedGraph::vertex_descriptor;
            using Edge = typename EmbeddedGraph::edge_descriptor;
            /**
             * \brief Observer to assign edge IDs of the original graph to created (half)edges, without
             * interfering with the insertion process. Needs to take into account non-planarity that will split edges.
             * Notification of that is centralized in a later step.
             */
            class AssignIdsAtConstruction : public CGAL::Arr_observer<Traits::Arr>
            {
            public:
                using Arr = Traits::Arr;

            private:
                const EmbeddedGraph& m_graph;
                // Original next-ids. Tells use which vertices/edges were already present.
                SimplificationObserver* m_observer;

                template<typename HeHandle>
                std::size_t setNewEdgeId(HeHandle handle)
                {
                    const auto eId = m_idGenerator->edgeIdGen.generateNext();
                    handle->data().id = eId;
                    handle->twin()->data().id = eId;
                    return eId;
                }

                Traits::GraphIdGenerator* m_idGenerator;

                std::optional<std::size_t> insertingEdge;

                std::optional<std::size_t> splittingEdgeId;
                bool wasSplit = false;

                std::vector<std::size_t> insertedSplitIds;

                std::map<std::size_t, std::set<Arr::Halfedge_handle>> splitEdges;

                std::size_t edgesCreated = 0;
                std::size_t edgesSplit = 0;
                // MAps edge ID to its vert IDs, ordered by geometry, lexicographically.
                std::map<std::size_t, std::pair<std::size_t, std::size_t>> m_edgeIdToVertIdGeomLexicographically;
                // Track which vertices got their ID.
                std::set<Vertex_handle> idWasSet;

                bool vertexHasId(Vertex_handle vertex)
                {
                    return idWasSet.find(vertex) != idWasSet.end();
                }
                std::size_t setIdIfNotSetAndNotify(Vertex_handle vertex, std::size_t id, bool verifySame = true)
                {
                    if (!vertexHasId(vertex))
                    {
                        vertex->data().id = id;
                        idWasSet.insert(vertex);
                        if (!m_notifyNewOnly || inputVertexIds.find(id) == inputVertexIds.end())
                        {
                            m_observer->handleNewVertex(vertex->data().id);
                        }
                    }
                    else if (verifySame && id != vertex->data().id)
                    {
                        throw std::runtime_error("Incompatible ids");
                    }
                    return vertex->data().id = id;;
                }
                std::pair<std::size_t, bool> generateIdIfNotSetAndNotify(Vertex_handle vertex)
                {
                    bool wasGenerated = false;
                    if (!vertexHasId(vertex))
                    {
                        vertex->data().id = m_idGenerator->vertexIdGen.generateNext();
                        idWasSet.insert(vertex);
                        wasGenerated = true;
                        m_observer->handleNewVertex(vertex->data().id);
                    }
                    return std::make_pair(vertex->data().id, wasGenerated);
                }
                void notifyEdgeCreate(Halfedge_handle he, std::size_t id)
                {
                    he->data().id = id;
                    he->twin()->data().id = SpecialEdgeIdFunctions::flipCanonical(id);
                    if(!m_notifyNewOnly || inputEdgeIds.find(id) == inputEdgeIds.end())
                    {
                        m_observer->handleNewEdge(id, he->source()->data().id, he->target()->data().id);
                    }
                }
                bool m_notifyNewOnly;
                std::set<std::size_t> inputVertexIds;
                std::set<std::size_t> inputEdgeIds;
            public:

                /**
                 * \brief 
                 * \param graph The input graph, not planar perse.
                 * \param observer 
                 * \param idGenerator 
                 * \param edgeToVertGeomLexico 
                 * \param notifyNewOnly Only notify on new objects that were not present in graph.
                 */
                AssignIdsAtConstruction(const EmbeddedGraph& graph, SimplificationObserver* observer, Traits::GraphIdGenerator* idGenerator,
                    const std::map<std::size_t, std::pair<std::size_t, std::size_t>>& edgeToVertGeomLexico, bool notifyNewOnly=false) :
                    m_graph(graph),
                    m_idGenerator(idGenerator),
                    m_observer(observer),
                    m_edgeIdToVertIdGeomLexicographically(edgeToVertGeomLexico),
                    m_notifyNewOnly(notifyNewOnly)
                {
                    if(notifyNewOnly)
                    {
                        GCpp::DS::computeVertexIdSet(graph, inputVertexIds);
                        for (const auto& kv : edgeToVertGeomLexico)
                        {
                            inputEdgeIds.insert(kv.first);
                        }
                    }
                    
                }

                void setCurrentlyInsertingEdge(std::size_t edgeId)
                {
                    if (!insertedSplitIds.empty())
                    {
                        m_observer->handleEdgeReroute(std::vector<std::size_t>{ insertingEdge.value() }, insertedSplitIds);
                    }
                    insertingEdge = edgeId;
                    wasSplit = false;
                    insertedSplitIds.clear();
                }
                void clearInsertingEdge()
                {
                    insertingEdge = {};
                }

                /**
                 * \brief Notify the observer of all the edge changes after planarization.
                 */
                void notifyEdgeChanges()
                {
                    std::cout << "Created edges: " << edgesCreated << ", split edges: " << edgesSplit << '\n';
                    // Reconstruct which edge gets mapped where.
                    for (auto kv : splitEdges)
                    {
                        // Edge ID of edge in input!
                        const auto eId = kv.first;

                        if (kv.second.size() == 1) //Maps to self
                        {
                            auto el = *kv.second.begin();
                            auto startVert = el->source();
                            auto endVert = el->target();
                            if (el->source()->point() > el->target()->point())
                            {
                                // Swap handles
                                std::swap(startVert, endVert);
                                el = el->twin();
                            }

                            auto verts = m_edgeIdToVertIdGeomLexicographically[eId];
                            setIdIfNotSetAndNotify(startVert, verts.first);
                            setIdIfNotSetAndNotify(endVert, verts.second);
                            notifyEdgeCreate(el, eId);
                        }
                        else
                        {
                            std::cout << "Edge count: " << kv.second.size() << '\n';
                            // Map vertices to connected edges, sort the vertices by point lexicographically.
                            std::map<Vertex_handle, std::vector<Halfedge_handle>, SortVertexByPoint> vertToEdge;
                            for (auto e : kv.second)
                            {
                                if (vertToEdge.find(e->source()) == vertToEdge.end()) vertToEdge[e->source()] = {};
                                if (vertToEdge.find(e->target()) == vertToEdge.end()) vertToEdge[e->target()] = {};
                                vertToEdge[e->source()].push_back(e);
                                vertToEdge[e->target()].push_back(e);
                            }
                            // Check our mapping makes sense.
                            if (vertToEdge.begin()->second.size() != 1) throw std::runtime_error("Invalid edge");
                            if (std::prev(vertToEdge.end())->second.size() != 1) throw std::runtime_error("Invalid edge");

                            // Assign vertices
                            auto verts = m_edgeIdToVertIdGeomLexicographically[kv.first];
                            setIdIfNotSetAndNotify(vertToEdge.begin()->first, verts.first);
                            setIdIfNotSetAndNotify(std::prev(vertToEdge.end())->first, verts.second);

                            // Notify observer: cannot go via CGAL, since this edge actually never existed in the arrangement.
                            if(!m_notifyNewOnly) m_observer->handleNewEdge(eId, verts.first, verts.second);

                            for (auto it = std::next(vertToEdge.begin()); it != std::prev(vertToEdge.end()); ++it)
                            {
                                // Result unused....
                                auto[id, wasGenerated] = generateIdIfNotSetAndNotify(it->first);

                            }

                            std::vector<std::size_t> subsegEdgeIDs;
                            Halfedge_handle prevEdge;
                            for (auto it = vertToEdge.begin(); it != std::prev(vertToEdge.end()); ++it)
                            {
                                if (it == vertToEdge.begin())
                                {
                                    prevEdge = it->second[0];
                                    if(it->first->data().id != prevEdge->source()->data().id)
                                    {
                                        prevEdge = prevEdge->twin();
                                    }
                                    const auto newEId = m_idGenerator->edgeIdGen.generateNext();
                                    notifyEdgeCreate(prevEdge, newEId);
                                    subsegEdgeIDs.push_back(newEId);
                                }
                                else
                                {
                                    if (it->second.size() != 2) throw std::runtime_error("invalid path");

                                    Halfedge_handle edge;
                                    if (prevEdge == it->second[0] || prevEdge == it->second[0]->twin()) edge = it->second[1];
                                    else if (prevEdge == it->second[1] || prevEdge == it->second[1]->twin()) edge = it->second[0];
                                    else throw std::runtime_error("Disjointness in path!");
                                    if(edge->source() != prevEdge->target())
                                    {
                                        edge = edge->twin();
                                    }
                                    const auto newEId = m_idGenerator->edgeIdGen.generateNext();
                                    notifyEdgeCreate(edge, newEId);
                                    subsegEdgeIDs.push_back(newEId);
                                    prevEdge = edge;
                                }
                                std::cout << "Edge between " << prevEdge->source()->data().id << " and " << prevEdge->target()->data().id << "\n";
                            }
                            // Notify of reroute
                            m_observer->handleEdgeReroute(std::vector<std::size_t>{eId}, subsegEdgeIDs);
                        }
                    }
                    for (auto v : this->arrangement()->vertex_handles())
                    {
                        if (!vertexHasId(v))
                        {
                            generateIdIfNotSetAndNotify(v);
                        }
                    }
                }


                void after_create_vertex(Vertex_handle vert) override
                {
                    // This should always be a non-planarity vertex.
                    //if (vert->is_at_open_boundary()) return;
                    //// Get established ID or generate a new one
                    //std::size_t vid = vert->point().data().has_value() ? vert->point().data().value() : m_idGenerator->vertexIdGen.generateNext();
                    //vert->set_data(Traits::VertexData{ vid,{} });

                    //if(vert->point().data().has_value())
                    //{
                    //    std::cout << "preassigned vertex\n";
                    //}
                    ////std::cout << "Vert at construction: " << vid << "\n";
                    //m_observer->handleNewVertex(vid);

                    // Will be assigned at notifyEdgeChanges()
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
                // When inserting an edge, splits should only be called on already present edges.
                void before_split_edge(typename ::CGAL::Arr_observer<Arr>::Halfedge_handle e, typename ::CGAL::Arr_observer<Arr>::Vertex_handle v,
                    const typename ::CGAL::Arr_observer<Arr>::X_monotone_curve_2 &c1, const typename ::CGAL::Arr_observer<Arr>::X_monotone_curve_2 &c2) override
                {
                    const auto eId = e->curve().data();
                    // Remove the halfedge as a split, since we are adding the splitted edges later.
                    if (splitEdges.find(eId) != splitEdges.end())
                    {
                        splitEdges[eId].erase(FaceMergeTraits::getCanonical(e));
                    }
                    //std::cout << "[BeforeSplitEdge]\n";
                    //std::cout << "[Facemerger] Splitting edge " << e->data() << " with vert " << v->data().id << '\n';
                }
                void after_split_edge(typename ::CGAL::Arr_observer<Arr>::Halfedge_handle he0, typename ::CGAL::Arr_observer<Arr>::Halfedge_handle he1) override
                {
                    if (!splittingEdgeId.has_value()) throw std::runtime_error("Splitting edge data was not set");

                    const auto eId = he0->curve().data();
                    if (splitEdges.find(eId) == splitEdges.end()) splitEdges[eId] = {};
                    // These edges we will process later.
                    splitEdges[eId].insert(FaceMergeTraits::getCanonical(he0));
                    splitEdges[eId].insert(FaceMergeTraits::getCanonical(he1));
                    ++edgesSplit;
                }
                void after_create_edge(typename ::CGAL::Arr_observer<Arr>::Halfedge_handle he) override
                {
                    //std::cout << "[AfterCreateEdge]\n";
                    if (he->is_fictitious()) return;
                    if (he->source()->is_at_open_boundary() || he->target()->is_at_open_boundary()) return;
                    const auto eId = he->curve().data();
                    if (splitEdges.find(eId) == splitEdges.end()) splitEdges[eId] = {};
                    splitEdges[eId].insert(he);

                    he->data().id = eId;
                    he->twin()->data().id = SpecialEdgeIdFunctions::flipCanonical(eId);
                    ++edgesCreated;
                }
            };
            void verifyVerticesConnectivity(const EmbeddedGraph& graph, const Traits::Arr& arrangement) const
            {
                /*using Vert = typename boost::graph_traits<EmbeddedGraph>::vertex_descriptor;
                std::map<std::size_t, Vert> mapping;
                GCpp::DS::computeVertexIdToVertexMap(graph, mapping);
                for(auto v: arrangement.vertex_handles())
                {
                    auto vid = v->data().id;
                    std::set<Vert> connectedVerts;
                    for(auto e : GCpp::Helpers::Iterators::range(boost::out_edges(mapping[vid], graph)))
                    {
                        auto target = boost::target(e, graph);
                        connectedVerts.insert(GCpp::DS::getVertexId(graph, target));
                    }
                    if (connectedVerts.find(vid) != connectedVerts.end()) throw std::runtime_error("Invalid vertex connection");
                    
                    std::set<Vert> verified;
                    std::size_t count = 0;
                    FaceMergeTraits::for_each_incident_edge(v,[&graph,&connectedVerts,&verified,&count](auto eh, auto& cmd)
                    {
                        auto other = eh->source()->data().id;
                        if (connectedVerts.find(other) == connectedVerts.end())throw std::runtime_error("illegal connection");
                        verified.insert(other);
                        ++count;
                    });
                    if(count < connectedVerts.size())
                    {
                        throw std::runtime_error("Missing connections");
                    }
                    if(verified.size() != connectedVerts.size())
                    {
                        throw std::runtime_error("Vertex duplicate ID!");
                    }
                }*/
            }
        public:

            void operator()(const EmbeddedGraph& inputGraph, Traits::Arr& arrangement, SimplificationObserver& observer, Traits::GraphIdGenerator& idGenerator,bool notifyNewOnly=false) const
            {
                // Unbounded face gets the zero index.
                arrangement.unbounded_face()->data().id = 0;
                observer.handleNewFace(0);
                // Modify to accomodate unbounded face.
                if(idGenerator.faceIdGen.isAtZero())
                {
                    idGenerator.faceIdGen = Traits::IdGenerator(1);
                }
                idGenerator.vertexIdGen = Traits::IdGenerator(boost::get_property(inputGraph, GCpp::DS::next_vertex_id_t{}));
                idGenerator.edgeIdGen = Traits::IdGenerator(boost::get_property(inputGraph, GCpp::DS::next_edge_id_t{}));

                namespace it = GCpp::Helpers::Iterators;
                auto cgalVertexPos = [&inputGraph](auto graphVertex)
                {
                    const auto pos = GCpp::DS::get_vertex_location(graphVertex, inputGraph);
                    //return Point_2(pos.x(), pos.y());
                    return Traits::Kernel::Point_2(CGAL_NT(pos.x()), CGAL_NT(pos.y()));
                };
                // Add vertices. This way, we can easily set additional data.
                std::vector<Traits::Point_2> vertices;
                for (auto vertex : it::range(boost::vertices(inputGraph)))
                {
                    vertices.push_back(Traits::Point_2(cgalVertexPos(vertex)));
                    // Don't need to notify, will happen during construction.
                }
                std::vector<Traits::Segment_2> segments;

                std::map<std::size_t, std::pair<std::size_t, std::size_t>> edgetToVertLexico;

                auto& inputGraphMod = const_cast<EmbeddedGraph&>(inputGraph);
                auto prevUseCanonical = inputGraphMod.m_use_canonical_edges_ids;
                inputGraphMod.m_use_canonical_edges_ids = true;
                std::cout << "Edges in the graph: " << boost::num_edges(inputGraph) << "\n";
                for (auto edge : it::range(boost::edges(inputGraph)))
                {
                    auto eId = GCpp::DS::getEdgeId(inputGraph, edge);
                    const auto src = boost::source(edge, inputGraph);
                    const auto target = boost::target(edge, inputGraph);
                    if (src >= vertices.size() || target >= vertices.size()) throw std::runtime_error("Vertex out of bounds");
                    auto srcPoint = vertices[src];
                    auto targetPoint = vertices[target];

                    // Get vertex ids.
                    auto minVId = GCpp::DS::getVertexId(inputGraph, src);
                    auto maxVId = GCpp::DS::getVertexId(inputGraph, target);
                    if (targetPoint < srcPoint) //Sort IDs by their points, lexicographically.
                    {
                        std::swap(minVId, maxVId);
                        eId = SpecialEdgeIdFunctions::flipCanonical(eId);
                    }
                    edgetToVertLexico[eId] = std::make_pair(minVId, maxVId);
                    Traits::Segment_2 seg(Traits::Kernel::Segment_2(srcPoint, targetPoint), eId);
                    segments.push_back(seg);
                }
                inputGraphMod.m_use_canonical_edges_ids = prevUseCanonical;
                // Assigns edge and vertex ids for new elements in the arrangment
                AssignIdsAtConstruction assigner(inputGraph, &observer, &idGenerator, edgetToVertLexico, notifyNewOnly);
                assigner.attach(arrangement);
                CGAL::insert(arrangement, segments.begin(), segments.end());
                assigner.notifyEdgeChanges();
                assigner.detach();

                // Assign face IDs and meta data.
                auto assignFaceIfUnassigned = [&idGenerator, &observer](auto faceHandle)
                {
                    if (!faceHandle->has_outer_ccb()) return;
                    if (!faceHandle->data().id.has_value())
                    {
                        faceHandle->data().id = idGenerator.faceIdGen.generateNext();
                        observer.handleNewFace(faceHandle->data().id.value());
                    }
                };
                for (auto eit = arrangement.edges_begin(); eit != arrangement.edges_end(); ++eit)
                {
                    assignFaceIfUnassigned(eit->face());
                    assignFaceIfUnassigned(eit->twin()->face());
                }
                FaceMergeTraits::verifyIds(arrangement);
                FaceMergeTraits::verifyIntegrity(arrangement, observer);
                verifyVerticesConnectivity(inputGraph, arrangement);
            }
        };

        template<typename EmbeddedGraph>
        using SetupArrangement = SetupModifiedArrangement<EmbeddedGraph>;
    };
}
#endif