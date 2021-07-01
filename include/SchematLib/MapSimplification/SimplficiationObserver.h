#ifndef SCHEMATLIB_MAPSIMPLIFICATION_SIMPLIFICATIONOBSERVER_H
#define SCHEMATLIB_MAPSIMPLIFICATION_SIMPLIFICATIONOBSERVER_H
#include <vector>
#include <iostream>
namespace SchematLib::MapSimplification
{
    /**
     * \brief Observer for structural changes of a map simplification.
     * Key assumption: initially, edge IDs and vertex IDs are in the range of the number of edges
     * resp. vertices.
     * Stub implementation that does nothing.
     */
    struct TrivialSimplificationObserver
    {
        TrivialSimplificationObserver(){}
        TrivialSimplificationObserver(const std::set<std::size_t>& originalVertexIds, const std::map<std::size_t, std::pair<std::size_t, std::size_t>> & originalEdgeIds){}

        void handleInitialEdgeCount(std::size_t edgeCount){}

        void handleInitialVertexCount(std::size_t vertexCount) {}

        void handleNewEdge(std::size_t edgeId, std::size_t srcVert, std::size_t targetVert){}

        void handleEdgeReplace(std::size_t edgeId, const std::vector<std::size_t>& replacementEdges) {}

        void handleNewVertex(std::size_t vertexId)
        {
        }
        void handleNewFace(std::size_t faceId)
        {
            
        }
        void handleEdgeToFaceConnection(std::size_t edgeId, std::size_t leftF, std::size_t rightF)
        {

        }
        void handleFaceMerge(std::size_t face, std::size_t targetFace)
        {
            
        }
        void handleFaceCollapseToVertex(std::size_t face, std::size_t vertex)
        {
            
        }
        void handleEdgeToFaceMerge(std::size_t edge, std::size_t face){}
        void handleVertexToFaceMerge(std::size_t vertex, std::size_t face) {}

        void handleVertexDelete(std::size_t vertex)
        {
        }
        /**
         * \brief Handles an edge collapse, where the given edge (ID) is completely merge to a vertex
         * \param edge 
         * \param newVertex The target merge vertex
         */
        void handleEdgeCollapse(std::size_t edge, std::size_t newVertex){}

        /**
         * \brief Handles an edge merge, where the edge is represented by the new edge.
         * \param edge 
         * \param newEdge 
         */
        void handleEdgeEdgeMerge(std::size_t edge, std::size_t newEdge) {}

        /**
         * \brief Handles the case where a vertex is inserted in an edge, splitting the edge in two sub edges
         * \param edge The edge to be split
         * \param newVertex 
         * \param newStartEdge 
         * \param newEndEdge 
         */
        void handleVertexInsert(std::size_t edge, std::size_t newVertex, std::size_t newStartEdge, std::size_t newEndEdge) {}

        /**
         * \brief Handles the case where a vertex is merged to another. Note that edges that will potentially be merged to the vertex emit
         * their own events.
         * \param vertex The original vertex
         * \param mergeVertex The merge vertex
         */
        void handleVertexMerge(std::size_t vertex, std::size_t mergeVertex)
        {
            
        }

        /**
         * \brief Handles the case when an edge is deleted and has no representation anymore
         * \param edge 
         */
        void handleEdgeDelete(std::size_t edge){}
        /**
         * \brief Handles the event when we reindex an edge.
         * \param edgeId The original edge ID
         * \param newEdgeId The new edge ID.
         */
        void handleEdgeReindex(std::size_t edgeId, std::size_t newEdgeId) {}

        void handleEdgeReroute(const std::vector<std::size_t>& originalRout, const std::vector<std::size_t>& newRoute) {}
        void handleEdgeRerouteVertexRoute(const std::vector<std::size_t>& originalRout, const std::vector<std::size_t>& newRoute) {}
    };
    struct LoggingSimplificationObserver
    {
        static inline const char* prefix = "[LoggingSimplificationObserver]";

        struct LogLine
        {
            LogLine()
            {
                std::cout << LoggingSimplificationObserver::prefix;
            }
            template<typename T>
            LogLine& operator<<(const T& t)
            {
                std::cout << t;
                return *this;
            }
            template<typename T>
            LogLine& operator<<(const std::vector<T>& t)
            {
                for(auto el : t)
                {
                    std::cout << ' ';
                    std::cout << el;
                }
                return *this;
            }
            ~LogLine() {
                std::cout << '\n';
            }
        };
        LoggingSimplificationObserver() {}
        LoggingSimplificationObserver(const std::set<std::size_t>& originalVertexIds, const std::map<std::size_t, std::pair<std::size_t, std::size_t>> & originalEdgeIds) {}

        void handleInitialEdgeCount(std::size_t edgeCount)
        {
            LogLine() << " InitialEdgeCount: " << edgeCount;
        }

        void handleInitialVertexCount(std::size_t vertexCount)
        {
            LogLine() << " InitialVertexCount: " << vertexCount;
        }

        void handleNewEdge(std::size_t edgeId, std::size_t srcVert, std::size_t targetVert)
        {
            LogLine() << " NewEdge: " << edgeId << ", src&target:" << srcVert << ',' << targetVert;
        }

        void handleNewVertex(std::size_t vertexId)
        {
            LogLine() << " NewVertex: " << vertexId ;
        }

        void handleEdgeReplace(std::size_t edgeId, const std::vector<std::size_t>& replacementEdges)
        {
            LogLine() << " EdgeReplace: " << edgeId << replacementEdges;
        }

        /**
         * \brief Handles an edge collapse, where the given edge (ID) is completely merge to a vertex
         * \param edge
         * \param newVertex The target merge vertex
         */
        void handleEdgeCollapse(std::size_t edge, std::size_t newVertex)
        {
            LogLine() << "EdgeCollapse: edge " << edge << " to vertex " << newVertex;
        }
        void handleVertexDelete(std::size_t vertex)
        {
            LogLine() << "Deleting vertex: " << vertex ;
        }

        /**
         * \brief Handles an edge merge, where the edge is represented by the new edge.
         * \param edge
         * \param newEdge
         */
        void handleEdgeEdgeMerge(std::size_t edge, std::size_t newEdge)
        {
            LogLine() << "EdgeEdgeMerge: " << edge << " to " << newEdge;
        }

        /**
         * \brief Handles the case where a vertex is inserted in an edge, splitting the edge in two sub edges
         * \param edge The edge to be split
         * \param newVertex
         * \param newStartEdge
         * \param newEndEdge
         */
        void handleVertexInsert(std::size_t edge, std::size_t newVertex, std::size_t newStartEdge, std::size_t newEndEdge)
        {
            LogLine() << "VertexInsert: " << edge << ": " << newVertex << " splitting edge into " << newStartEdge << " and " << newEndEdge;
        }

        /**
         * \brief Handles the case when an edge is deleted and has no representation anymore
         * \param edge
         */
        void handleEdgeDelete(std::size_t edge)
        {
            LogLine() << "EdgeDelete: " << edge;
        }
        /**
         * \brief Handles the event when we reindex an edge.
         * \param edgeId The original edge ID
         * \param newEdgeId The new edge ID.
         */
        void handleEdgeReindex(std::size_t edgeId, std::size_t newEdgeId)
        {
            LogLine() << "EdgeReindex: " << edgeId << " to " << newEdgeId;
        }
        void handleEdgeReroute(const std::vector<std::size_t>& originalRout, const std::vector<std::size_t>& newRoute)
        {
            LogLine() << "Rerouting " << originalRout << " to " << newRoute;
        }
        void handleEdgeRerouteVertexRoute(const std::vector<std::size_t>& originalRout, const std::vector<std::size_t>& newRoute)
        {
            LogLine() << "Rerouting " << originalRout << " to vertex route " << newRoute;
        }
        /**
         * \brief Handles the case where a vertex is merged to another. Note that edges that will potentially be merged to the vertex emit
         * their own events.
         * \param vertex The original vertex
         * \param mergeVertex The merge vertex
         */
        void handleVertexMerge(std::size_t vertex, std::size_t mergeVertex)
        {
            LogLine() << "Merging vertex  " << vertex << " to vertex  " << mergeVertex;
        }
    };
}
#endif