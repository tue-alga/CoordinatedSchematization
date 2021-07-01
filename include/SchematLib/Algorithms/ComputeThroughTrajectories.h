#ifndef SCHEMATLIB_ALGORITHMS_COMPUTETHROUGHTRAJECTORIES_H
#define SCHEMATLIB_ALGORITHMS_COMPUTETHROUGHTRAJECTORIES_H
#include <filesystem>
#include <map>
#include <GCpp/DS/Trajectory.h>

#include "SchematLib/IO/GraphTrajectoryCsvReader.h"
#include "SchematLib/IO/ReadGeneralization.h"
#include "SchematLib/IO/ReadWriteOsmMap.h"
#include "SchematLib/Models/EmbeddedGraph.h"
#include <iostream>

namespace SchematLib::Algorithms
{
    struct Subtrajectory
    {
        std::string sourceTrajectory;
        // First index of this subtrajectory
        std::size_t startEdgeIndex;
        // Exclusive end index of this subtrajectory.
        std::size_t endEdgeIndex;
        // In base graph!
        Models::UndirectedEmbeddedGraph::vertex_descriptor startV;
        Models::UndirectedEmbeddedGraph::vertex_descriptor endV;
    };
    struct OffroadElement
    {
        using Vertex = Models::UndirectedEmbeddedGraph::vertex_descriptor;
        std::string sourceTrajectory;
        // First index of this subtrajectory
        std::size_t startEdgeIndex;
        // Exclusive end index of this subtrajectory.
        std::size_t endEdgeIndex;
        // In base graph!
        Models::UndirectedEmbeddedGraph::vertex_descriptor startV;
        Models::UndirectedEmbeddedGraph::vertex_descriptor endV;
        // Start and end connected to subgraph.
        bool startConnected = false;
        bool endConnected = false;

        void reset()
        {
            startConnected = false;
            endConnected = false;
            startEdgeIndex = 0;
            endEdgeIndex = 0;
        }

        void setConnectivity(bool start, bool end)
        {
            startConnected = start;
            endConnected = end;
        }
        void setVertices(Vertex start, Vertex end)
        {
            startV = start;
            endV = end;
        }
        void setEdgeIndices(std::size_t start, std::size_t end)
        {
            startEdgeIndex=start;
            endEdgeIndex = end;
        }

        std::size_t size() const
        {
            return endEdgeIndex - startEdgeIndex; //+1 ? 
        }
        bool isThrough() const { return startConnected && endConnected; }
        bool isFullyInside() const { return !startConnected && !endConnected; }
    };

    class ComputeThroughTrajectories
    {
    public:
        using Trajectory = GCpp::DS::Trajectory<std::size_t, std::string, std::vector>;
        using SubtrajectoryType = std::variant<Subtrajectory, OffroadElement>;
        
    private:
        void readTrajectories(const std::filesystem::path& trajectoriesPath, std::map<std::string, Trajectory>& output)
        {
            struct TrajectoryLoader
            {
                using Trajectory = GCpp::DS::Trajectory<std::size_t, std::string, std::vector>;
                std::map<std::string, Trajectory>* m_target = nullptr;
                TrajectoryLoader() {}
                TrajectoryLoader(std::map<std::string, GCpp::DS::Trajectory<std::size_t, std::string, std::vector>>* target) :m_target(target) {}
                TrajectoryLoader& operator=(Trajectory&& traj)
                {
                    (*m_target)[traj.trajectoryData()] = std::move(traj);
                    return *this;
                }
                TrajectoryLoader& operator*()
                {
                    return *this;
                }
            };
            TrajectoryLoader loader(&output);
            SchematLib::IO::GraphTrajectoryCsvReader reader;
            reader.read(trajectoriesPath, loader);
        }

        struct GraphData
        {
            using Vertex = Models::UndirectedEmbeddedGraph::vertex_descriptor;
            using Edge = Models::UndirectedEmbeddedGraph::edge_descriptor;
            // Generalization data
            std::set<std::size_t> edgeIdSet;
            std::set<std::size_t> vertexIdSet;
            // Base map for comparison, with id to edge mapping
            const Models::UndirectedEmbeddedGraph* m_baseGraph = nullptr;
            std::map<std::size_t, Models::UndirectedEmbeddedGraph::edge_descriptor> edgeIdToEdgeMap;

            void compute(const Models::UndirectedEmbeddedGraph& generalization, const Models::UndirectedEmbeddedGraph& base)
            {
                GCpp::DS::computeEdgeIdSet(generalization, edgeIdSet);
                GCpp::DS::computeVertexIdSet(generalization, vertexIdSet);
                GCpp::DS::computeEdgeIdToEdgeMap(base, edgeIdToEdgeMap);
                m_baseGraph = &base;
            }
            // Vertex on base map
            Vertex source(Edge e)const
            {
                return boost::source(e, *m_baseGraph);
            }
            // Vertex on base map
            Vertex target(Edge e)const
            {
                return boost::target(e, *m_baseGraph);
            }
            // Resulting edge in basemap
            Models::UndirectedEmbeddedGraph::edge_descriptor edgeForIdInBase(std::size_t edgeId) const
            {
                return edgeIdToEdgeMap.at(edgeId);
            }
            // On selection
            bool hasEdge(std::size_t edgeId) const
            {
                return edgeIdSet.find(edgeId) != edgeIdSet.end();
            }
            // On selection
            bool hasVertex(std::size_t vertexId) const
            {
                return vertexIdSet.find(vertexId) != vertexIdSet.end();
            }
            // On basemap
            std::optional<Vertex> findCommonVertex(Edge e0, Edge e1) const
            {
                if (source(e0) == source(e1) || source(e0) == target(e1)) return { source(e0) };
                if (target(e0) == source(e1) || target(e0) == target(e1)) return { target(e0) };
                return {};
            }
            // On basemap
            std::optional<Vertex> otherVert(Edge e, Vertex oneOfEdge) const
            {
                if (source(e) == oneOfEdge) return { target(e) };
                if (target(e) == oneOfEdge) return { source(e) };
                return {};
            }


            std::pair<bool, bool> edgeIdVertexPresence(std::size_t  edge) const
            {
                auto e = edgeForIdInBase(edge);
                auto srcId = GCpp::DS::getVertexId(*m_baseGraph, boost::source(e, *m_baseGraph));
                auto targetId = GCpp::DS::getVertexId(*m_baseGraph, boost::target(e, *m_baseGraph));
                return std::make_pair(hasVertex(srcId), hasVertex(targetId));
            }
            /**
             * \brief Both vertices are present in graph or not.
             * \param edge 
             * \return 
             */
            bool edgeIdHasPresentVertices(std::size_t edge) const
            {
                auto pair = edgeIdVertexPresence(edge);
                return pair.first && pair.second;
            }
        };

        void computeTrajectoryVertices(const Trajectory& trajectory, const GraphData& data, std::vector<Models::UndirectedEmbeddedGraph::vertex_descriptor>& vertices)
        {
            if(trajectory.size() == 1)
            {
                auto edge= data.edgeIdToEdgeMap.at(*trajectory.begin());
                vertices = { boost::source(edge,*data.m_baseGraph), boost::target(edge,*data.m_baseGraph) };
            }
            using It = std::decay_t<decltype(trajectory.end())>;
            It curr = trajectory.begin();
            It next = std::next(curr);
            It end = trajectory.end();
            for(; next != end; ++curr, ++next)
            {
                if (*curr == *next) continue;
                // Should be a common vertex between them.
                auto commonVert = data.findCommonVertex(data.edgeIdToEdgeMap.at(*curr), data.edgeIdToEdgeMap.at(*next));
                if (!commonVert.has_value()) throw std::runtime_error("Disconnected trajectory");
                auto vert = commonVert.value();

                vertices.push_back(vert);
                // Reconstruct backward
                while(true)
                {
                    auto otherVert = data.otherVert(data.edgeIdToEdgeMap.at(*curr), vertices.back());
                    if (!otherVert.has_value()) throw std::runtime_error("Incorrect path chain");
                    vertices.push_back(otherVert.value());
                    if (curr == trajectory.begin()) break;
                    --curr;
                }
                std::reverse(vertices.begin(), vertices.end());
                curr = next;
                // Reconstruct forward
                while(curr != trajectory.end())
                {
                    auto otherVert = data.otherVert(data.edgeIdToEdgeMap.at(*curr), vertices.back());
                    if (!otherVert.has_value()) throw std::runtime_error("Incorrect path chain");
                    vertices.push_back(otherVert.value());
                    ++curr;
                }

                break;
            }
        }

        void segmentTrajectory(const Trajectory& trajectory, const Models::UndirectedEmbeddedGraph& graph,
            const GraphData& data,
            std::vector<SubtrajectoryType>& segments)
        {
            const auto trajectoryId = trajectory.trajectoryData();
            using It = std::decay_t<decltype(trajectory.end())>;
            It curr = trajectory.begin();
            It end = trajectory.end();
            std::size_t index = 0;
            std::vector<Models::UndirectedEmbeddedGraph::vertex_descriptor> trajectoryVertices;
            computeTrajectoryVertices(trajectory, data, trajectoryVertices);

            std::vector<OffroadElement> offroadElements;

            OffroadElement element{ trajectoryId, 0,0 };

            for (; curr != end; ++curr, ++index)
            {
                const auto eId = *curr;
                // Edge is offroad.
                if (!data.hasEdge(eId))
                {
                    // Existing subtrajectory
                    if(element.size() > 0)
                    {
                        // Check if vertex cuts through present graph.
                        if(data.hasVertex(GCpp::DS::getVertexId(*data.m_baseGraph,trajectoryVertices[index+1])))
                        {
                            element.endConnected = true;
                            element.endEdgeIndex = index + 1;
                            element.endV = trajectoryVertices[index + 1];
                            offroadElements.push_back(element);
                            element.reset();
                            std::cout << "[ComputeThroughTrajectories]\tAdding element" << '\n';
                        }
                        else
                        {
                            element.endEdgeIndex = index + 1;
                            element.endV = trajectoryVertices[index+1];
                        }
                    }
                    // Non-existing subtrajectory
                    else
                    {
                        // Check if is singleton
                        if (data.edgeIdHasPresentVertices(eId))
                        {
                            element.setConnectivity(true, true);
                            element.setEdgeIndices(index, index + 1);
                            element.setVertices(trajectoryVertices[index], trajectoryVertices[index + 1]);
                            offroadElements.push_back(element);
                            // Reset
                            element.reset();
                        }
                        // Start new segment
                        else
                        {
                            element.setConnectivity(data.hasVertex(GCpp::DS::getVertexId(*data.m_baseGraph, trajectoryVertices[index])), false);
                            element.startV = trajectoryVertices[index];
                            element.endV = trajectoryVertices[index+1];
                            element.setEdgeIndices(index, index + 1);
                        }
                    }
                }
                else
                {
                    if(element.size() > 0)
                    {
                        offroadElements.push_back(element);
                        std::cout << "[ComputeThroughTrajectories]\tAdding element" << '\n';
                        element.reset();
                    }
                }
            }
            if (element.size() > 0)
            {
                offroadElements.push_back(element);
                std::cout << "[ComputeThroughTrajectories]\tAdding element" << '\n';
            }
            for(const auto& el: offroadElements)
            {
                std::cout << "Conn:" << el.startConnected << ',' << el.endConnected << '\n';
                if(el.startConnected)
                {
                    const auto vId = GCpp::DS::getVertexId(*data.m_baseGraph, el.startV);
                    if (!data.hasVertex(vId)) throw std::runtime_error("Stopped on vertex not in selection");
                }
                if(el.endConnected)
                {
                    const auto vId = GCpp::DS::getVertexId(*data.m_baseGraph, el.endV);
                    if (!data.hasVertex(vId)) throw std::runtime_error("Stopped on vertex not in selection");
                }
            }

            if(offroadElements.size() > 0)
            {
                // Interleave subtrajectories with offroad elements
                if(offroadElements.front().startEdgeIndex > 0)
                {
                    // Insert initial element
                    segments.push_back(
                        Subtrajectory{
                            trajectoryId, 0, offroadElements.front().startEdgeIndex,
                            trajectoryVertices.front(), trajectoryVertices[offroadElements.front().startEdgeIndex]
                        }
                    );
                }
                auto it = offroadElements.begin();
                for(;it != offroadElements.end(); ++it)
                {
                    if(segments.size() > 0 && std::holds_alternative<OffroadElement>(segments.back()))
                    {
                        // Determine between subtrajectory.
                        const OffroadElement& prev = std::get<OffroadElement>(segments.back());
                        if(prev.endEdgeIndex != it->startEdgeIndex)
                        {
                            // Add subtrajectory.
                            Subtrajectory subtrajectory{ it->sourceTrajectory, prev.endEdgeIndex, it->startEdgeIndex,
                                trajectoryVertices[prev.endEdgeIndex + 1], trajectoryVertices[it->startEdgeIndex]
                            };
                            segments.push_back(subtrajectory);
                        }
                        segments.push_back(*it);
                    }
                    else
                    {
                        segments.push_back(*it);
                    }
                }
            }
            else
            {
                segments.push_back(
                    Subtrajectory{
                        trajectoryId, 0, trajectory.size(),
                        trajectoryVertices.front(), trajectoryVertices.back()
                    }
                );
            }
        }

        void processTrajectories(const std::map<std::string, Trajectory>& trajectories, const Models::UndirectedEmbeddedGraph& graph,
            const GraphData& data,
            std::map<std::string,std::vector<SubtrajectoryType>>& segmentation)
        {
            for(const auto& kv: trajectories)
            {
                if (kv.second.size() == 0) continue;
                segmentation[kv.first] = {};
                segmentTrajectory(kv.second, graph, data, segmentation[kv.first]);
            }
        }
    public:

        void operator()(const std::filesystem::path& mapFile, const std::filesystem::path& generalizationFile, const std::filesystem::path& mmTrajectories,
            std::map<std::string, std::vector<SubtrajectoryType>>& segmentation, std::map<std::string, Trajectory>& trajectories)
        {
            // Read trajectories;
            readTrajectories(mmTrajectories, trajectories);
            // Read graph
            Models::UndirectedEmbeddedGraph generalization, map;
            IO::ReadGeneralization reader;
            reader(generalizationFile, generalization);

            IO::ReadOsmMap undirectedReader;
            int epsg;
            undirectedReader.readCustomFormat(mapFile, map, epsg);

            GraphData data;
            data.compute(generalization, map);
            processTrajectories(trajectories, generalization, data, segmentation);
        }

        void operator()(const std::filesystem::path& mapFile, const std::filesystem::path& generalizationFile, const std::filesystem::path& mmTrajectories,
            std::map<std::string, std::vector<SubtrajectoryType>>& segmentation)
        {
            // Read trajectories;
            std::map<std::string, Trajectory> trajectories;
            this->operator()(mapFile, generalizationFile, mmTrajectories, segmentation, trajectories);
        }

        void operator()(const std::filesystem::path& mapFile, const std::filesystem::path& generalizationFile, const std::filesystem::path& mmTrajectories,
            const std::filesystem::path& outputFile)
        {
            std::map<std::string, std::vector<SubtrajectoryType>> segmentation;
            // Read trajectories;
            std::map<std::string, Trajectory> trajectories;
            this->operator()(mapFile, generalizationFile, mmTrajectories, segmentation, trajectories);
            //readTrajectories(mmTrajectories, trajectories);
            //// Read graph
            //Models::UndirectedEmbeddedGraph generalization, map;
            //IO::ReadGeneralization reader;
            //reader(generalizationFile, generalization);

            //IO::ReadOsmMap undirectedReader;
            //int epsg;
            //undirectedReader.readCustomFormat(mapFile, map, epsg);

            //GraphData data;
            //data.compute(generalization, map);
            //processTrajectories(trajectories, generalization, data, offroadElements);

            std::size_t totalTrajectoryEdges = 0;
            for(const auto& kv: trajectories)
            {
                totalTrajectoryEdges += kv.second.size();
            }

            std::size_t offroadEdgeCount = 0;
            std::ofstream output(outputFile.string());
            for(const auto& kv: segmentation)
            {
                const auto& trajectoryId = kv.first;
                for(const auto& segmentationElement: kv.second)
                {
                    if (std::holds_alternative<OffroadElement>(segmentationElement))
                    {
                        const auto& el = std::get<OffroadElement>(segmentationElement);
                        output << el.sourceTrajectory << ' ' << trajectories[el.sourceTrajectory].size() << ' ' << el.startEdgeIndex << ' ' << el.startConnected << ' ' << el.endConnected << ' ' << (el.endEdgeIndex - el.startEdgeIndex);
                        const auto& traj = trajectories[el.sourceTrajectory];
                        for (auto i = el.startEdgeIndex; i < el.endEdgeIndex; ++i) {
                            output << ' ' << traj.at(i);
                        }
                        output << '\n';
                        offroadEdgeCount += (el.endEdgeIndex - el.startEdgeIndex);
                    }
                    else
                    {
                        // TODO: maybe write.
                    }
                }
                
            }
            std::cout << "Offroad els: " << offroadEdgeCount << "/" << totalTrajectoryEdges << '\n';
        }
    };
}
#endif