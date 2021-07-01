#ifndef TULIB_GRAPHHOPPERPROBETRAITS_H
#define TULIB_GRAPHHOPPERPROBETRAITS_H

#include "tulib/io/csv/ParseDate.h"
#include "tulib/io/csv/CategoricalField.h"
#include "tulib/io/csv/csv.h"
#include "tulib/ProbeTraits.h"

namespace here {
	namespace c2d {
		namespace raw {

			enum InputColumns {
				_PROBE_ID, _LAT, _LON,
			};

			// Fields of interest: all
			enum ProbeColumns {
				PROBE_ID, LAT, LON,
			};

			class ProbeParseDate : public ParseDate {
			public:
				explicit ProbeParseDate(std::time_t ts = 0, string date_format = "%Y-%m-%d %H:%M:%S")
					:ParseDate(ts, std::move(date_format)) { }
			};

			class ProviderCategoricalField : public CategoricalField<std::string, ProviderCategoricalField> {};

			typedef csv<std::tuple<string, long double, long double>,
				_PROBE_ID, _LAT, _LON> ProbeCsv;

			typedef typename ProbeCsv::value_type ProbePoint;

			typedef _ProbeTraits<ProbeColumns, ProbeParseDate, ProbeCsv, ProbePoint> ProbeTraits;

		}  // namespace raw
	}  // namespace c2d
}  // namespace here

#endif //TULIB_GRAPHHOPPERPROBETRAITS_H