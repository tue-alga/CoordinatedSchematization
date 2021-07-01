#ifndef SCHEMATLIB_MAPSIMPLIFICATION_INEFFICIENTARRANGEMENTSETUP_H
#define SCHEMATLIB_MAPSIMPLIFICATION_INEFFICIENTARRANGEMENTSETUP_H
#include "FaceMergeTraits.h"
namespace SchematLib::MapSimplification
{
    template<typename EmbeddedGraph, typename SimplificationObserver>
    class InefficientSetupArrangement
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

            GraphIdGenerator* m_idGenerator;

            std::optional<std::size_t> insertingEdge;

            std::optional<std::size_t> splittingEdgeId;
            bool wasSplit = false;

            std::vector<std::size_t> insertedSplitIds;
        public:

            AssignIdsAtConstruction(const EmbeddedGraph& graph, SimplificationObserver* observer, GraphIdGenerator* idGenerator) :
                m_graph(graph),
                m_idGenerator(idGenerator),
                m_observer(observer)
            {}

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


            void after_create_vertex(typename ::CGAL::Arr_observer<Arr>::Vertex_handle vert) override
            {
                // This should always be a non-planarity vertex.
                if (vert->is_at_open_boundary()) return;
                //std::cout << "[AfterCreateVertex]\n";
                const auto vid = m_idGenerator->vertexIdGen.generateNext();
                vert->set_data(VertexData{ vid,{} });
                //std::cout << "Vert at construction: " << vid << "\n";
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
            // When inserting an edge, splits should only be called on already present edges.
            void before_split_edge(typename ::CGAL::Arr_observer<Arr>::Halfedge_handle e, typename ::CGAL::Arr_observer<Arr>::Vertex_handle v,
                const typename ::CGAL::Arr_observer<Arr>::X_monotone_curve_2 &c1, const typename ::CGAL::Arr_observer<Arr>::X_monotone_curve_2 &c2) override
            {
                splittingEdgeId = e->data().id.value();
                wasSplit = true;
                //std::cout << "[BeforeSplitEdge]\n";
                //std::cout << "[Facemerger] Splitting edge " << e->data() << " with vert " << v->data().id << '\n';
            }
            void after_split_edge(typename ::CGAL::Arr_observer<Arr>::Halfedge_handle he0, typename ::CGAL::Arr_observer<Arr>::Halfedge_handle he1) override
            {
                if (!splittingEdgeId.has_value()) throw std::runtime_error("Splitting edge data was not set");
                //std::cout << "[AfterSplitEdge]\n";
                // Assign IDs if not present, but also assign a new ID for an original edge (otherwise, the mapping gets messed up).
                std::size_t eId0 = setNewEdgeId(he0);
                std::size_t eId1 = setNewEdgeId(he1);
                m_observer->handleNewEdge(eId0, he0->source()->data().id, he0->target()->data().id);
                m_observer->handleNewEdge(eId1, he1->source()->data().id, he1->target()->data().id);
                /*std::cout << "[Facemerger] Split edge to  " << he0->data() << ", " << he1->data() << ", vs= "
                    << he0->source()->data().id << ',' << he0->target()->data().id << " | "
                    << he1->source()->data().id << ',' << he1->target()->data().id
                    << '\n';*/
                m_observer->handleEdgeReroute({ splittingEdgeId.value() }, { eId0, eId1 });
                splittingEdgeId = {};
            }
            void after_create_edge(typename ::CGAL::Arr_observer<Arr>::Halfedge_handle he) override
            {
                //std::cout << "[AfterCreateEdge]\n";
                if (he->is_fictitious()) return;
                if (he->source()->is_at_open_boundary() || he->target()->is_at_open_boundary()) return;
                const FaceMergeTraits::VertexData srcVertData = he->source()->data();
                const FaceMergeTraits::VertexData targetVertData = he->target()->data();

                // One of the two is newly added (i.e. vertex due to non-planarity) Happensw hen it has no
                // index in the original graph
                //if (!srcVertData.index.has_value() || !targetVertData.index.has_value())
                if (wasSplit)
                {
                    const auto eid = setNewEdgeId(he);
                    /*std::cout << "[FaceMerger] New edge added! " << eid << " with vs="
                        << he->source()->data().id << ',' << he->target()->data().id <<
                        std::endl;*/
                        //if (!wasSplit) throw std::runtime_error("Expecting split");
                    insertedSplitIds.push_back(eid);
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
                }
                m_observer->handleNewEdge(he->data().id.value(), srcVertData.id, targetVertData.id);
            }

            void after_split_face(typename ::CGAL::Arr_observer<Arr>::Face_handle origAndNewFace,
                typename ::CGAL::Arr_observer<Arr>::Face_handle newFace, bool is_hole) override
            {
                newFace->data().id = m_idGenerator->faceIdGen.generateNext();
            }
        };

    public:

        void operator()(const EmbeddedGraph& inputGraph, Arr& arrangement, SimplificationObserver& observer, GraphIdGenerator& idGenerator) const
        {
            // Unbounded face gets the zero index.
            arrangement.unbounded_face()->data().id = 0;
            idGenerator.faceIdGen = IdGenerator(1);
            idGenerator.vertexIdGen = IdGenerator(boost::get_property(inputGraph, GCpp::DS::next_vertex_id_t{}));
            idGenerator.edgeIdGen = IdGenerator(boost::get_property(inputGraph, GCpp::DS::next_edge_id_t{}));

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
                const auto vId = GCpp::DS::getVertexId(inputGraph, vertex);
                cgalVert->set_data(FaceMergeTraits::VertexData{ vId, vertices.size() });
                observer.handleNewVertex(vId);
                vertices.push_back(cgalVert);
            }

            // Assigns edge and vertex ids for new elements in the arrangment
            AssignIdsAtConstruction assigner(inputGraph, &observer, &idGenerator);
            assigner.attach(arrangement);

            //Add edges
            std::size_t edgeCount = 0;
            std::size_t edgeTotal = boost::num_edges(inputGraph);
            const std::size_t percInc = 5;
            std::size_t nextPerc = percInc;
            for (auto edge : it::range(boost::edges(inputGraph)))
            {
                const auto eId = GCpp::DS::getEdgeId(inputGraph, edge);
                const auto src = boost::source(edge, inputGraph);
                const auto target = boost::target(edge, inputGraph);
                if (src >= vertices.size() || target >= vertices.size()) throw std::runtime_error("Vertex out of bounds");
                auto srcPoint = cgalVertexPos(src);
                auto targetPoint = cgalVertexPos(target);
                Segment_2 seg(srcPoint, targetPoint);
                assigner.setCurrentlyInsertingEdge(eId);
                // Observer AssignIdsAtConstruction will handle ID assignment.
                //std::cout << "[SetupArr] Creating edge " << edgeCount << '\n';
                CGAL::insert(arrangement, seg);
                ++edgeCount;
                if (edgeCount * 100 > nextPerc * edgeTotal)
                {
                    std::cout << "At " << nextPerc << "% of edges\n";
                    nextPerc += percInc;
                }
            }
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
            verifyIntegrity(arrangement, observer);
        }
    };

}
#endif
