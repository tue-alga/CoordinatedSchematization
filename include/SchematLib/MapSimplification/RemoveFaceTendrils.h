#ifndef SCHEMATLIB_MAPSIMPLIFICATION_REMOVEFACETENDRILS_H
#define SCHEMATLIB_MAPSIMPLIFICATION_REMOVEFACETENDRILS_H
#include "FaceMergeTraits.h"
#include <boost/graph/lookup_edge.hpp>
#include <GCpp/DS/BoostEmbeddedGraph.h>

#include "SimplficiationObserver.h"
#include "FaceMergeStrategies.h"

namespace SchematLib::MapSimplification
{
    /**
     * \brief Merges faces of an embeddedgraph by applying the given MergeStrategy
     * \tparam EmbeddedGraph Embedded graph type. Should be boost::graph compatible and GCpp Graph compatible.
     * \tparam SelectPredicate Predicate for selecting faces that need merging
     * \tparam OrderStrategy The order of merging marked faces
     * \tparam MergeStrategy Strategy to use for merging faces that were marked as needing merging
     */
    template<
        typename EmbeddedGraph,
        typename SimplificationObserver = TrivialSimplificationObserver
    >
        class RemoveFaceTendrils
    {
    public:
        using Vertex = typename boost::graph_traits<EmbeddedGraph>::vertex_descriptor;
        using Edge = typename boost::graph_traits<EmbeddedGraph>::edge_descriptor;
        using Point = decltype(GCpp::DS::get_vertex_location(std::declval<Vertex>(), std::declval<const EmbeddedGraph&>()));
        using NT = decltype(std::declval<const Point&>().x());
    private:
        // Setup CGAL types.
        using CGAL_NT = FaceMergeTraits::CGAL_NT;
        using Kernel = FaceMergeTraits::Kernel;
        //using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
        using Point_2 = FaceMergeTraits::Point_2;
        using Segment_2 = FaceMergeTraits::Segment_2;

        // Associate with the arrangement vertex and edges an ID, and with the faces the FaceData struct.
        using Arr = FaceMergeTraits::Arr;
        using CGAL_Vertex_handle = typename Arr::Vertex_handle;
        using Face_iterator = typename Arr::Face_iterator;

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
                vert->set_data(FaceMergeTraits::VertexData{ vid,{} });
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

        /**
         * \brief Assigns IDs to vertices and edges during the merging steps
         */
        class AssignVertexId : public CGAL::Arr_observer<Arr>
        {
            std::size_t m_nextVertexId = 0;
            std::size_t m_nextEdgeId = 0;
            std::size_t generateEdgeId()
            {
                const auto val = m_nextEdgeId;
                ++m_nextEdgeId;
                return val;
            }
            std::size_t generateVertexId()
            {
                const auto val = m_nextVertexId;
                ++m_nextVertexId;
                return val;
            }
            SimplificationObserver* m_observer;
        public:
            using Obs = CGAL::Arr_observer<Arr>;
            using HEHandle = typename ::CGAL::Arr_observer<Arr>::Halfedge_handle;
            using VHandle = typename CGAL::Arr_observer<Arr>::Vertex_handle;
            using XmonCurve = typename Obs::X_monotone_curve_2;

            AssignVertexId(std::size_t nextVertexId, std::size_t nextEdgeId, SimplificationObserver& observer) :
                m_nextVertexId(nextVertexId),
                m_nextEdgeId(nextEdgeId),
                m_observer(&observer)
            {}
            // Only called on newly added vertices!
            void after_create_vertex(typename CGAL::Arr_observer<Arr>::Vertex_handle vertexHandle) override
            {
                if (vertexHandle->is_at_open_boundary()) return;
                const auto vid = generateVertexId();
                vertexHandle->set_data(FaceMergeTraits::VertexData{ vid, {} });
                std::cout << "[FM-merge-step] Constructed new vertex: " << vid << '\n';
            }
            std::size_t nextVertId() const { return m_nextVertexId; }
            std::size_t nextEdgeId() const { return m_nextEdgeId; }

            void after_split_edge(typename ::CGAL::Arr_observer<Arr>::Halfedge_handle he0, typename ::CGAL::Arr_observer<Arr>::Halfedge_handle he1) override
            {
                //TODO: reinstate later if needed!!!!
                //throw std::runtime_error("Not expecting splits!");
                // For now, regenerate ids
                after_create_edge(he0);
                after_create_edge(he1);
            }

            void before_remove_vertex(VHandle vertexHandle) override
            {
                //m_observer->handleVertexDelete(vertexHandle->data().id);
            }

            void after_create_edge(typename CGAL::Arr_observer<Arr>::Halfedge_handle heHandle) override
            {
                const auto eid = generateEdgeId();
                heHandle->data().id = eid;
                heHandle->twin()->data().id = SpecialEdgeIdFunctions::flipCanonical(eid);
                std::cout << "[FM-merge-step] Constructed new edge: " << heHandle->data() << ", vs=" << heHandle->source()->data().id << ',' << heHandle->target()->data().id << '\n';
            }
            void before_split_face(typename Obs::Face_handle fHandle, typename Obs::Halfedge_handle heHandle)override
            {
                //std::cout << "[FM-merge-step] Splitting face?\n";
            }

            void after_merge_edge(typename CGAL::Arr_observer<Arr>::Halfedge_handle) override
            {
                //std::cout << "[FM-merge-step] Merging edge?\n";
            }
            void after_merge_face(
                typename CGAL::Arr_observer<Arr>::Face_handle) override
            {
                //std::cout << "[FM-merge-step] Merging face?\n";
            }
            void after_merge_inner_ccb(
                typename CGAL::Arr_observer<Arr>::Face_handle,
                typename CGAL::Arr_observer<Arr>::Ccb_halfedge_circulator) override
            {
                //std::cout << "[FM-merge-step] Merging inner CBB?\n";
            }
            void after_merge_outer_ccb(
                typename CGAL::Arr_observer<Arr>::Face_handle,
                typename CGAL::Arr_observer<Arr>::Ccb_halfedge_circulator) override
            {
                //std::cout << "[FM-merge-step] Merging outer CBB?\n";
            }
            void after_remove_inner_ccb(
                typename CGAL::Arr_observer<Arr>::Face_handle) override
            {
                //std::cout << "[FM-merge-step] Removing inner CBB?\n";
            }
            void after_remove_outer_ccb(
                typename CGAL::Arr_observer<Arr>::Face_handle) override
            {
                //std::cout << "[FM-merge-step] Removing inner CBB?\n";
            }
        };

        /**
         * \brief Tracks the changes when creating and merging faces.
         * Ensures that any face handle we originally provided as associated with an ID still points
         * to the ''original'' face, albeit with potentially slight modifications.
         */
        class FaceChangeTracker : public CGAL::Arr_observer<Arr>
        {
            std::map<std::size_t, typename Arr::Face_handle>& m_handles;
        public:
            FaceChangeTracker(std::map<std::size_t, typename Arr::Face_handle>& handles) : m_handles(handles) {}

            void before_merge_face(typename ::CGAL::Arr_observer<Arr>::Face_handle f1,
                typename ::CGAL::Arr_observer<Arr>::Face_handle f2,
                typename ::CGAL::Arr_observer<Arr>::Halfedge_handle he) override
            {
                // When merging, make sure we preserve the face ID of the original not-to-be-merged face.
                if (!f1->data().id.has_value() || !f2->data().id.has_value()) return;
                if (f1->data().beingDeleted)
                {
                    f1->data().id = f2->data().id;
                    f1->data().beingDeleted = false;
                }
                else
                {
                    f2->data().id = f1->data().id;
                    f2->data().beingDeleted = false;
                }
            }
            void after_merge_face(typename ::CGAL::Arr_observer<Arr>::Face_handle face) override
            {
                if (!face->data().id.has_value()) return;
                auto val = face->data().id.value();
                // Upddate the handle if it was referenced before.
                if (m_handles.find(val) == m_handles.end()) return;
                m_handles[val] = face;
            }

            void after_split_face(typename ::CGAL::Arr_observer<Arr>::Face_handle f1, typename ::CGAL::Arr_observer<Arr>::Face_handle f2, bool isHole) override
            {
                // We only split faces to delete them afterwards.
                f1->data().beingDeleted = true;
                f1->data().id = 0;
                f2->data().beingDeleted = true;
                f2->data().id = 0;
            }
        };
        using GCppVector2 = decltype(GCpp::DS::get_vertex_location(std::declval<Vertex>(), std::declval<const EmbeddedGraph&>()));

        void setupArrangement(const EmbeddedGraph& inputGraph, Arr& arrangement, SimplificationObserver& observer, std::size_t& nextEdgeId, std::size_t& nextVertId, std::size_t& nextFaceId) const
        {
            namespace it = GCpp::Helpers::Iterators;
            FaceMergeTraits::SetupArrangement<EmbeddedGraph> setup;
            FaceMergeTraits::GraphIdGenerator gen;
            gen.faceIdGen = FaceMergeTraits::IdGenerator(observer.maxFId + 1);
            setup(inputGraph, arrangement, observer, gen, true);


            nextEdgeId = gen.edgeIdGen.generateNext();
            nextVertId = gen.vertexIdGen.generateNext();
            nextFaceId = gen.faceIdGen.generateNext();

            EmbeddedGraph test;
            convertArrangementToGraph(arrangement, test, nextVertId, nextEdgeId);

            observer.verifyIntegrity(test);
        }

        // TODO: mapping of edges/vertices relative to original graph
        void convertArrangementToGraph(const Arr& arrangement, EmbeddedGraph& output, std::size_t nextVertId, std::size_t nextEdgeId, SimplificationObserver* obs = nullptr) const
        {
            namespace it = GCpp::Helpers::Iterators;

            // Setup initial size of graph
            output = EmbeddedGraph(arrangement.number_of_vertices());
            output.m_use_canonical_edges_ids = true;
            boost::get_property(output, GCpp::DS::next_edge_id_t{}) = nextEdgeId;
            boost::get_property(output, GCpp::DS::next_vertex_id_t{}) = nextVertId;

            // Function for creating an embeddedgraph position from a CGAL Point.
            auto cgalVertexToOwn = [](const auto& vertex)
            {
                //return Models::Point(CGAL::to_double(vertex[0]), CGAL::to_double(vertex[1])); //Hopefully?
                return Models::Point(CGAL::to_double(vertex[0]), CGAL::to_double(vertex[1])); //Hopefully?
                //return Models::Point(vertex[0], vertex[1]);
            };
            // Verify unique vert ID's on vertices.
            verifyUniqueVertIds(arrangement);
            verifyUniqueEdgeIds(arrangement);

            std::size_t vertexNum = 0;
            std::map<std::size_t, std::size_t> vertexToIndex;

            auto hasVertex = [&vertexToIndex](std::size_t val)
            {
                return vertexToIndex.find(val) != vertexToIndex.end();
            };
            std::cout << "[FM] Adding verts: ";
            for (auto vHandle = arrangement.vertices_begin(); vHandle != arrangement.vertices_end(); ++vHandle)
            {
                if (vHandle->is_at_open_boundary()) continue; //Ignore vertices at infinity
                if (vHandle->is_isolated())
                {
                    std::cout << "[FM] Isolated vertex removed" << std::endl;
                    if (obs)
                    {
                        obs->handleVertexDelete(vHandle->data().id);
                        continue;
                    }
                    else
                    {
                        throw std::runtime_error("Isolated vertex in arrangement?");
                    }
                }
                boost::put(GCpp::DS::vertex_location_t{}, output, vertexNum, cgalVertexToOwn(vHandle->point()));
                boost::put(GCpp::DS::vertex_id_t{}, output, vertexNum, vHandle->data().id);

                if (hasVertex(vHandle->data().id))
                {
                    throw std::runtime_error("Duplicate vertex ID detected!");
                }
                vertexToIndex[vHandle->data().id] = vertexNum;
                std::cout << vHandle->data().id << " ";
                ++vertexNum;
            }
            std::cout << '\n';
            std::cout << "Assigned " << vertexNum << " vertices, in arrangement we have " << arrangement.number_of_vertices() << '\n';

            for (std::size_t i = arrangement.number_of_vertices() - 1; i >= vertexNum; --i)
            {
                boost::remove_vertex(i, output);
            }

            if (!GCpp::DS::hasUniqueVertexIds(output)) throw std::runtime_error("Non-unique vertex IDs before adding edges");
            // Create edges
            for (auto edgeHandle = arrangement.edges_begin(); edgeHandle != arrangement.edges_end(); ++edgeHandle)
            {
                if (edgeHandle->is_fictitious()) continue;

                const auto srcId = edgeHandle->source()->data().id;
                const auto targetId = edgeHandle->target()->data().id;
                if (!hasVertex(srcId) || !hasVertex(targetId)) throw std::runtime_error("Could not find vertex for edge");
                const auto edgeId = edgeHandle->data().id.value();
                const auto src = vertexToIndex.at(srcId);
                const auto target = vertexToIndex.at(targetId);
                // TODO: assumption vertex index and handle in boost is the same
                auto[edge, wasAdded] = boost::add_edge(src, target, output);
                assert(wasAdded);
                boost::put(GCpp::DS::edge_id_t{}, output, edge, edgeId);
            }
            if (!GCpp::DS::hasUniqueVertexIds(output)) throw std::runtime_error("Non-unique vertex IDs");
            if (!GCpp::DS::hasUniqueEdgeIds(output)) throw std::runtime_error("Non-unique edge IDs");

        }

        void verifyUniqueVertIds(const Arr& arr) const
        {
            std::unordered_set<std::size_t> ids;
            for (auto v : GCpp::Helpers::Iterators::range(arr.vertices_begin(), arr.vertices_end()))
            {
                const auto vId = v.data().id;
                if (ids.find(vId) != ids.end()) throw std::runtime_error("Duplicate vert ID!");
                ids.insert(vId);
            }
        }
        void verifyUniqueEdgeIds(const Arr& arr) const
        {
            std::unordered_set<std::size_t> ids;
            for (auto e : GCpp::Helpers::Iterators::range(arr.edges_begin(), arr.edges_end()))
            {
                if (!e.data().id.has_value()) throw std::runtime_error("Edge without ID!");
                const auto eId = e.data().id.value();
                if (ids.find(eId) != ids.end()) throw std::runtime_error("Duplicate edge ID!");
                ids.insert(eId);
            }
        }

    public:
        RemoveFaceTendrils(){}

        void operator()(EmbeddedGraph& inputGraph, SimplificationObserver& observer)//, std::map<Vertex, Vertex>& mergeGroups, std::map<Vertex,Point>& newVertexLocations) const
        {
            runWithRecompute(inputGraph, observer);
        }
        void runWithRecompute(EmbeddedGraph& inputGraph, SimplificationObserver& observer)
        {
            std::cout << "#############################################\n";
            std::cout << "##################  RemoveFaceTendrils #############\n";
            std::cout << "#############################################\n";
            namespace it = GCpp::Helpers::Iterators;

            // Construct arrangement for input graph
            Arr arrangement;
            std::size_t nextEdgeId, nextVertId, nextFaceId;
            // Setup the arrangement
            setupArrangement(inputGraph, arrangement, observer, nextEdgeId, nextVertId, nextFaceId);

            verifyUniqueVertIds(arrangement);
            verifyUniqueEdgeIds(arrangement);

            // Run over the faces
            for (auto iter = arrangement.faces_begin(); iter != arrangement.faces_end(); ++iter)
            {
                // If the faces doesn't have a ccb, it is not a 'real' face.
                if (!iter->has_outer_ccb()) continue;
                if (iter->is_unbounded()) continue;

                typename Arr::Face_handle handle = iter;
                Operators::RemoveInteriorFaceTendrils tendrilRemover;
                std::set<std::size_t> affectedVertices, affectedEdges;
                tendrilRemover(arrangement, handle, affectedVertices, affectedEdges);
                // Notify observer that all these are mapped to the face.
                for (auto v : affectedVertices) observer.handleVertexToFaceMerge(v, handle->data().id.value());
                for (auto e : affectedEdges) observer.handleEdgeToFaceMerge(e, handle->data().id.value());
            }
            convertArrangementToGraph(arrangement, inputGraph, nextVertId, nextEdgeId, &observer);
        }
    };
}
#endif