#ifndef TULIB_GRAPHHOPPER_H
#define TULIB_GRAPHHOPPER_H

#include "tulib/TrajectoryTraits.h"
#include "GraphHopperProbeTraits.h"
#include "tulib/ColumnarTrajectory.h"
#include "tulib/TabularTrajectory.h"

namespace here {
	namespace c2d {
		namespace raw {

			typedef here::c2d::raw::ProbeTraits ProbeTraits;
			constexpr static int SplitByFieldIdx = ProbeTraits::ProbeColumns::PROBE_ID;
			constexpr static int SortByFieldIdx = ProbeTraits::ProbeColumns::PROBE_ID;
			using columnar_trajectory_type = ColumnarTrajectory<string, long double, long double>;
			using tabular_trajectory_type = TabularTrajectory<string, long double, long double>;

			using ColumnarTrajectoryTraits = _TrajectoryTraits<ProbeTraits, SplitByFieldIdx, SortByFieldIdx, columnar_trajectory_type>;
			using TabularTrajectoryTraits = _TrajectoryTraits<ProbeTraits, SplitByFieldIdx, SortByFieldIdx, tabular_trajectory_type>;
		}  // namespace raw
	}  // namespace c2d
}  // namespace here

#endif //TULIB_LAMETRO_H