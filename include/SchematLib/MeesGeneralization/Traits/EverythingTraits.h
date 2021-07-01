#ifndef EVERYTHINGTRAITS_H
#define EVERYTHINGTRAITS_H
#include "tulib/HereTrajectoryTraits.h"
//#include "LAMetroTrajectoryTraits.h"
#include "GeometryBackendTraits.h"

//I don't think this class is used anymore.

//Qt classes cannot be templated. Therefore, all qt objects inherit some things like the Number Type and other traits from here.
struct EverythingTraits
{
	// Specializations for the Commit2Data raw probe format
	using TrajectoryTraits = here::c2d::raw::TabularTrajectoryTraits;
	using ProbeTraits = typename TrajectoryTraits::ProbeTraits;
	using GeometryTraits = GeometryKernel::TulibGeometryKernel;
	typedef GeometryTraits::NT NT;
};

#endif //EVERYTHINGTRAITS_H
