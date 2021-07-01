/*! @file OutlierDetectionTraits.h
 *  @brief  Defines the consistency check used by the output sensitive outlier detection
 *  @authors Mees van de Kerkhof (m.a.vandekerkhof@uu.nl)
 */

#ifndef TULIB_OUTLIERDETECTIONTRAITS_H
#define TULIB_OUTLIERDETECTIONTRAITS_H

//#include "tulib/geom/Interface.h"
#include "tulib/geo/geo.h"
#include <cmath>
#include <vector>

namespace tulib_algorithms
{

/*!
 *
 * @tparam GeometryTraits
 */
	template <class GeometryTraits>
	class DistanceOverTimeCheck {
	private:
		typedef typename GeometryTraits::NT NT;
		typedef std::tuple<NT,NT, time_t> InputTuple;
		NT threshold;

	public:
		DistanceOverTimeCheck() = default;
		/*!
		*@param InThreshold
		*/
		DistanceOverTimeCheck(NT InThreshold) : threshold(InThreshold) {};

		/*!
		*@param i1
		*@param i2
		*/
		bool operator()(InputTuple i1, InputTuple i2) {
			auto lat0 = std::get<0>(i1);
			auto lat1 = std::get<0>(i2);
			auto lon0 = std::get<1>(i1);
			auto lon1 = std::get<1>(i2);
			auto t0 = std::get<2>(i1);
			auto t1 = std::get<2>(i2);
			auto len = distance_exact(lat0, lon0, lat1, lon1);
			/*cout << "--Consistency Check--" << endl;
			cout << "Lat/Lon1: (" << lat0 << "," << lon0 << ");" << endl;
			cout << "Lat/Lon2: (" << lat1 << "," << lon1 << ");" << endl;
			cout << "Spatial Distance: " << len << "m;" << endl;
			cout << "Time 1: " << t0 << endl;
			cout << "Time 2: " << t1 << endl;
			cout << "Temporal Distance: " << abs(t1 - t0) << "s;" << endl;
			cout << "Average speed: " << len / abs(t1 - t0) << "m/s;" << endl << endl;*/
			return len / abs(t1 - t0) <= threshold;
		}
	};
}

#endif //TULIB_OUTLIERDETECTIONTRAITS_H
