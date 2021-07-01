#ifndef SCHEMATLIB_IO_GRAPHTRAJECTORYCSVREADER_H
#define SCHEMATLIB_IO_GRAPHTRAJECTORYCSVREADER_H
#include <filesystem>
#include <fstream>
#include <sstream>
#include <GCpp/DS/Trajectory.h>
namespace SchematLib::IO
{
    class GraphTrajectoryCsvReader
    {
        std::size_t m_trajectoriesToRead = std::numeric_limits<std::size_t>::max();
    public:
        void setTrajectoriesToRead(std::size_t trajectoriesToRead)
        {
            m_trajectoriesToRead = trajectoriesToRead;
        }
        std::size_t trajectoriesToRead() const
        {
            return m_trajectoriesToRead;
        }
        using Trajectory = GCpp::DS::Trajectory<std::size_t, std::string,std::vector>;
        template<typename OutputIt>
        void read(const std::filesystem::path& path, OutputIt trajectories)
        {
            std::ifstream stream(path.c_str());
            if (!stream.is_open()) throw std::runtime_error("Could not open" + path.string());

            std::string line;
            // Read header line
            std::getline(stream, line);
            std::size_t readCount = 0;
            while(std::getline(stream,line))
            {
                Trajectory traj;
                std::stringstream lineStream(line);
                lineStream >> traj.trajectoryData();
                std::size_t elementCount = 0;
                lineStream >> elementCount;
                std::size_t id;
                for(auto i = 0; i < elementCount; ++i)
                {
                    lineStream >> id;
                    traj.push_back(id);
                }
                *trajectories = std::move(traj);
                ++readCount;
                if (readCount == m_trajectoriesToRead) break;
            }
        }
    };
}
#endif
