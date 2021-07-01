#ifndef SCHEMATLIB_IO_READRAWTRAJECTORIES_H
#define SCHEMATLIB_IO_READRAWTRAJECTORIES_H
#include <filesystem>
#include <fstream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <GCpp/DS/Trajectory.h>


#include "SchematLib/Models/EmbeddedGraph.h"
namespace SchematLib::IO
{
    struct ReadMMTrajectories
    {
        using MMTrajectory = GCpp::DS::Trajectory<std::size_t, std::string, std::vector>;
        void operator()(const std::filesystem::path& trajectoryInputFilePath, std::vector<MMTrajectory>& output)
        {
            auto ext = trajectoryInputFilePath.extension();
            if(ext.string() == ".boosttxt")
            {
                readBoostTxt(trajectoryInputFilePath, output);
            }
            else if(ext.string() == ".intermediate")
            {
                readIntermediate(trajectoryInputFilePath, output);
            }
            else
            {
                throw std::runtime_error("Unsupported extension " + ext.string() + " for ReadMMTrajectories");
            }
        }
        void readIntermediate(const std::filesystem::path& trajectoryInputFilePath, std::vector<MMTrajectory>& output)
        {
            std::ifstream inputStream(trajectoryInputFilePath.string());
            if (!inputStream.is_open()) throw std::runtime_error("Could not open input file: " + trajectoryInputFilePath.string());
            auto clean_id = [](std::string& id)
            {
                for (std::size_t i = 0; i < id.size(); ++i)
                {
                    if (std::isspace(id[i]))
                    {
                        id[i] = '_';
                    }
                    if (id[i] == ';')
                    {
                        id[i] = '_';
                    }
                }
            };
            std::string line;
            while(true)
            {
                if (!std::getline(inputStream, line)) break;

                std::stringstream lineStream(line);

                /*if (i * 100 > currentTargetPerc * size)
                {
                    std::cout << "[MapGeneralizer]\tLoading at " << currentTargetPerc << "%, " << (i + 1) << "/" << size << "\n";
                    currentTargetPerc += percInc;
                }*/

                output.emplace_back();
                auto& trajectory = output.back();
                // Read from line
                std::getline(lineStream, trajectory.trajectoryData(), ';');
                clean_id(trajectory.trajectoryData());

                /*if (currentTargetPerc == 100)
                {
                    std::cout << "[MapGeneralizer]\t" << (i + 1) << "/" << size << ":" << trajectory.trajectoryData() << "\n";
                }*/


                std::size_t trajectorySize = 0;
                using ElementData = std::remove_cv_t<std::remove_reference_t<decltype(trajectory.at(0))>>;
                {
                    std::size_t data;
                    while (lineStream >> data)
                    {
                        trajectory.push_back(data);
                        ++trajectorySize;
                    }
                }
            }
        }
        void readBoostTxt(const std::filesystem::path& trajectoryInputFilePath, std::vector<MMTrajectory>& output)
        {
            std::ifstream inputStream(trajectoryInputFilePath.string());
            if (!inputStream.is_open()) throw std::runtime_error("Could not open input file: " + trajectoryInputFilePath.string());

            // Input stream
            boost::archive::text_iarchive in(inputStream);

            // Execute mapmatching
            std::string mapFileString;
            in >> boost::serialization::make_nvp("map", mapFileString);
            std::size_t size = 0;
            in >> boost::serialization::make_nvp("count", size);

            // Output trajectory edge indices in generalization graph (embedded output.
            std::size_t uniqueEdgesPerTrajectory = 0;
            std::size_t duplicateEdgesPerTrajectory = 0;
            std::size_t maxDegeneracy = 0;
            // Stream from input file
            for (auto i = 0; i < size; ++i)
            {
                output.emplace_back();
                auto& trajectory = output.back();
                in >> boost::serialization::make_nvp("trajectoryName", trajectory.trajectoryData());
                std::size_t trajectorySize;
                in >> boost::serialization::make_nvp("trajectorySize", trajectorySize);
                using ElementData = std::remove_cv_t<std::remove_reference_t<decltype(trajectory.at(0))>>;
                for (auto j = 0; j < trajectorySize; ++j)
                {
                    std::size_t data = 0;
                    in >> data;
                    trajectory.push_back(data);
                }
            }
        }
    };
}
#endif