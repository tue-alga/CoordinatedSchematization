#ifndef SCHEMATLIB_TRAJECTORYSCHEMATIZATION_TRAJECTORYBUNDLETRAITS_H
#define SCHEMATLIB_TRAJECTORYSCHEMATIZATION_TRAJECTORYBUNDLETRAITS_H
#include <vector>
namespace SchematLib::TrajectorySchematization
{

    struct TrajectoryBundle
    {
        std::vector<std::size_t> edges;
        std::size_t matchedTrajectoriesCount = 0;
        std::size_t supportClass = 0;
        long double bundleLength = 0;
    };

    struct TrajectoryBundleTraits
    {
        using Bundle = TrajectoryBundle;
    };
}

#endif
