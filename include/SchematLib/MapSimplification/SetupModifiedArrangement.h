#ifndef SCHEMATLIB_MAPSIMPLIFICATION_SETUPMODIFIEDARRANGEMENT_H
#define SCHEMATLIB_MAPSIMPLIFICATION_SETUPMODIFIEDARRANGEMENT_H
#include "FaceMergeTraits.h"
namespace SchematLib::MapSimplification
{
    // Try out here first
    struct ModifiedFaceMergeTraits : public FaceMergeTraits
    {
        using ArrTraits = CGAL::Arr_segment_traits_2<Kernel>;
        // Add ID to curve
        using AddedCurveTraits = CGAL::Arr_curve_data_traits_2<ArrTraits, std::size_t>;
        using Point_2 = typename ArrTraits::Point_2;
        using Segment_2 = typename AddedCurveTraits::Curve_2;

        using Arr = CGAL::Arrangement_2<AddedCurveTraits, CGAL::Arr_extended_dcel<AddedCurveTraits, VertexData, EdgeData, FaceData>>;
    };

    template<typename EmbeddedGraph>
    class SetupModifiedArrangement
    {
    public:
        using SimplificationObserver = EdgeTrackingObserver;
        using Traits = ModifiedFaceMergeTraits;
    private:
        using CGAL_NT = Traits::CGAL_NT;
        using Vertex = typename EmbeddedGraph::vertex_descriptor;
        using Edge = typename EmbeddedGraph::edge_descriptor;
        /**
         * \brief Observer to assign edge IDs of the original graph to created (half)edges, without
         * interfering with the insertion process. Needs to take into account non-planarity that will split edges.
         * Notification of that is centralized in a later step.
         */
        class AssignIdsAtConstruction : public CGAL::Arr_observer<ModifiedFaceMergeTraits::Arr>
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
            std::size_t setIdIfNotSetAndNotify(Vertex_handle vertex, std::size_t id, bool verifySame=true)
            {
                if(!vertexHasId(vertex))
                {
                    vertex->data().id = id;
                    idWasSet.insert(vertex);
                    m_observer->handleNewVertex(vertex->data().id);
                }
                else if (verifySame && id != vertex->data().id) 
                {
                    throw std::runtime_error("Incompatible ids");
                }
                return vertex->data().id = id;;
            }
            std::pair<std::size_t,bool> generateIdIfNotSetAndNotify(Vertex_handle vertex)
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
                he->twin()->data().id = id;
                m_observer->handleNewEdge(id, he->source()->data().id, he->target()->data().id);
            }
        public:
            
            AssignIdsAtConstruction(const EmbeddedGraph& graph, SimplificationObserver* observer, Traits::GraphIdGenerator* idGenerator,
                const std::map<std::size_t, std::pair<std::size_t, std::size_t>>& edgeToVertGeomLexico) :
                m_graph(graph),
                m_idGenerator(idGenerator),
                m_observer(observer),
                m_edgeIdToVertIdGeomLexicographically(edgeToVertGeomLexico)
            {
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

            struct SortVertexByPoint
            {
                bool operator()(Vertex_handle v0, Vertex_handle v1) const
                {
                    return v0->point() < v1->point();
                }
            };

            /**
             * \brief Notify the observer of all the edge changes after planarization.
             */
            void notifyEdgeChanges()
            {
                std::cout << "Created edges: " << edgesCreated << ", split edges: " << edgesSplit << '\n';
                // Reconstruct which edge gets mapped where.
                for(auto kv : splitEdges)
                {
                    const auto eId = kv.first;
                    if(kv.second.size() == 1) //Maps to self
                    {
                        auto el = *kv.second.begin();

                        auto startVert = el->source();
                        auto endVert = el->target();
                        if(el->source()->point() > el->target()->point())
                        {
                            // Swap handles
                            std::swap(startVert, endVert);
                        }

                        auto verts = m_edgeIdToVertIdGeomLexicographically[kv.first];
                        setIdIfNotSetAndNotify(startVert, verts.first);
                        setIdIfNotSetAndNotify(endVert, verts.second);
                        notifyEdgeCreate(el, kv.first);
                    }
                    else
                    {
                        std::cout << "Edge count: " << kv.second.size() << '\n';
                        // Map vertices to connected edges, sort the vertices by point lexicographically.
                        std::map<Vertex_handle, std::vector<Halfedge_handle>, SortVertexByPoint> vertToEdge;
                        for(auto e : kv.second)
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

                        // Notify observer
                        m_observer->handleNewEdge(eId, verts.first, verts.second);

                        for (auto it = std::next(vertToEdge.begin()); it != std::prev(vertToEdge.end()); ++it)
                        {
                            auto [id,wasGenerated] = generateIdIfNotSetAndNotify(it->first);
                        }

                        std::vector<std::size_t> subsegEdgeIDs;
                        Halfedge_handle prevEdge;
                        for(auto it = vertToEdge.begin(); it != std::prev(vertToEdge.end()); ++it)
                        {
                            if(it ==vertToEdge.begin())
                            {
                                prevEdge = it->second[0];
                                const auto newEId = m_idGenerator->edgeIdGen.generateNext();
                                notifyEdgeCreate(prevEdge, newEId);
                                subsegEdgeIDs.push_back(newEId);
                            }
                            else
                            {
                                if (it->second.size() != 2) throw std::runtime_error("invalid path");

                                Halfedge_handle edge;
                                if (prevEdge == it->second[0]) edge = it->second[1];
                                else if (prevEdge == it->second[1]) edge = it->second[0];
                                else throw std::runtime_error("Disjointness in path!");
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
                for(auto v : this->arrangement()->vertex_handles())
                {
                    if(!vertexHasId(v))
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
                ++edgesCreated;
            }
        };

    public:

        void operator()(const EmbeddedGraph& inputGraph, Traits::Arr& arrangement, SimplificationObserver& observer, Traits::GraphIdGenerator& idGenerator) const
        {
            //TODO Cleanup later
            bool doCanonical = true;

            // Unbounded face gets the zero index.
            arrangement.unbounded_face()->data().id = 0;
            idGenerator.faceIdGen = Traits::IdGenerator(1);
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
                if(targetPoint < srcPoint) //Sort IDs by their points, lexicographically.
                {
                    std::swap(minVId, maxVId);
                    eId = SpecialEdgeIdFunctions::flipCanonical(eId);
                }
                edgetToVertLexico[eId] = std::make_pair(minVId, maxVId);
                Traits::Segment_2 seg(Traits::Kernel::Segment_2(srcPoint, targetPoint), eId);
                segments.push_back(seg);
            }
            // Assigns edge and vertex ids for new elements in the arrangment
            AssignIdsAtConstruction assigner(inputGraph, &observer, &idGenerator, edgetToVertLexico);
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
        }
    };

}
#endif