#pragma once
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
 // Created by onur on 9/23/18.
 //

#ifndef TULIB_TAXITRAITSNOHEADER_H
#define TULIB_TAXITRAITSNOHEADER_H

#include "tulib/io/csv/ParseDate.h"
#include "tulib/io/csv/CategoricalField.h"
#include "tulib/io/csv/csv.h"
#include "tulib/ProbeTraits.h"

namespace here {
	namespace c2d {
		namespace raw {

			enum InputColumns {
				_PROBE_ID, _SAMPLE_DATE, _LON, _LAT
			};

			// Fields of interest: all
			enum ProbeColumns {
				PROBE_ID, SAMPLE_DATE, LON, LAT
			};

			class ProbeParseDate : public ParseDate {
			public:
				explicit ProbeParseDate(std::time_t ts = 0, string date_format = "%Y-%m-%d %H:%M:%S")
					:ParseDate(ts, std::move(date_format)) { }
			};

			class ProviderCategoricalField : public CategoricalField<std::string, ProviderCategoricalField> {};

			typedef csv<std::tuple<string, ProbeParseDate, double, double>,
				_PROBE_ID, _SAMPLE_DATE, _LON, _LAT> ProbeCsv;

			typedef typename ProbeCsv::value_type ProbePoint;

			typedef _ProbeTraits<ProbeColumns, ProbeParseDate, ProbeCsv, ProbePoint, false> ProbeTraits;

		}  // namespace raw
	}  // namespace c2d
}  // namespace here

#endif //TULIB_HEREPROBETRAITSNOHEADER_H
