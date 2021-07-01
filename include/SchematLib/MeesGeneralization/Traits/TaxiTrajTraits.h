/*
Copyright (C) 2018-2019 HERE Global B.V. and its affiliate(s).
All rights reserved.

This software and other materials contain proprietary information
controlled by HERE and are protected by applicable copyright legislation.
Any use and utilization of this software and other materials and
disclosure to any third parties is conditional upon having a separate
agreement with HERE for the access, use, utilization or disclosure of this
software. In the absence of such agreement, the use of the software is not
allowed.
 */

 //
 // Created by onur on 10/12/18.
 //

#ifndef TULIB_TAXITRAJECTORYTRAITSNOHEADER_H
#define TULIB_TAXITRAJECTORYTRAITSNOHEADER_H

#include "tulib/TrajectoryTraits.h"
#include "TaxiProbeTraits.h"
#include "tulib/ColumnarTrajectory.h"
#include "tulib/TabularTrajectory.h"

namespace here {
	namespace c2d {
		namespace raw {

			typedef here::c2d::raw::ProbeTraits ProbeTraits;
			constexpr static int SplitByFieldIdx = ProbeTraits::ProbeColumns::PROBE_ID;
			constexpr static int SortByFieldIdx = ProbeTraits::ProbeColumns::SAMPLE_DATE;
			using columnar_trajectory_type = ColumnarTrajectory<string, ProbeParseDate, double, double>;
			using tabular_trajectory_type = TabularTrajectory<string, ProbeParseDate, double, double>;

			using ColumnarTrajectoryTraits = _TrajectoryTraits<ProbeTraits, SplitByFieldIdx, SortByFieldIdx, columnar_trajectory_type>;
			using TabularTrajectoryTraits = _TrajectoryTraits<ProbeTraits, SplitByFieldIdx, SortByFieldIdx, tabular_trajectory_type>;
		}  // namespace raw
	}  // namespace c2d
}  // namespace here

#endif //TULIB_HERETRAJECTORYTRAITS_H
