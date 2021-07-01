#ifndef SCHEMATLIB_TRAJECTORYSCHEMATIZATION_OFFROADTRAITS_H
#define SCHEMATLIB_TRAJECTORYSCHEMATIZATION_OFFROADTRAITS_H
#include <filesystem>
#include <fstream>
#include <GCpp/DS/Trajectory.h>
#include <GCpp/DS/BoostEmbeddedGraph.h>
namespace SchematLib::TrajectorySchematization
{
    template<typename EmbeddedGraph>
    struct OffroadTraits
    {
        using Vertex = typename EmbeddedGraph::vertex_descriptor;
        using Point = std::decay_t<decltype(GCpp::DS::get_vertex_location(std::declval<Vertex>(), std::declval<EmbeddedGraph>()))>;
        struct OffroadData
        {
            std::string sourceTrajectoryId;
            std::size_t segmentIndex = std::numeric_limits<std::size_t>::max();
            bool startOnBoundary = false;
            bool endOnBoundary = false;
            std::size_t startVId = std::numeric_limits<std::size_t>::max();
            std::size_t endVId = std::numeric_limits<std::size_t>::max();
        };
        using OffroadTrajectory = GCpp::DS::Trajectory<std::size_t, OffroadData>;

        void readOffroadTrajectories(const std::filesystem::path& path, std::vector<OffroadTrajectory>& trajectories)
        {
            std::ifstream stream(path.string());
            std::string line;
            while(std::getline(stream,line))
            {
                std::stringstream lineStream(line);
                trajectories.push_back({});
                auto& traj = trajectories.back();
                auto& trajectoryData = traj.trajectoryData();
                std::size_t size = 0;
                stream >> trajectoryData.sourceTrajectoryId >> trajectoryData.segmentIndex
                    >> trajectoryData.startOnBoundary >> trajectoryData.endOnBoundary
                    >> trajectoryData.startVId >> trajectoryData.endVId
                    >> size;
                for(std::size_t i = 0; i < size; ++i)
                {
                    std::size_t eId;
                    stream >> eId;
                    traj.push_back(eId);
                }
            }
        }
        void writeOffroadTrajectories(const std::filesystem::path& path, const std::vector<OffroadTrajectory>& trajectories)
        {
            std::ofstream stream(path.string());
            const char sp = ' ';
            for(const auto& traj: trajectories)
            {
                const auto& trajectoryData = traj.trajectoryData();
                stream << trajectoryData.sourceTrajectoryId << sp << trajectoryData.segmentIndex
                    << sp << trajectoryData.startOnBoundary << sp << trajectoryData.endOnBoundary
                    << sp << trajectoryData.startVId << sp << trajectoryData.endVId
                    << sp << traj.size();
                for(auto el : traj)
                {
                    stream << sp << el;
                }
                stream << '\n';
            }
        }
    };
}
#endif
