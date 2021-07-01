#pragma once

#include <iostream>
#include "tulib/utils/Asserts.h"

template <class Kernel, class Norm >
class Discrete_Frechet {
private:
	typedef typename Kernel::NT NT;
	template <class InputIterator,
		typename = tulib_core::requires_random_access_iterator<InputIterator>,
		typename = tulib_core::requires_tulib_point<Kernel,
		typename InputIterator::value_type>>
		NT distance(InputIterator iter_a, InputIterator iter_b) {
		Norm norm;
		typename Kernel::Tulib_Vector v = *iter_b - *iter_a;
		return norm(v);
	}
public:
	template <class InputIterator,
		typename = tulib_core::requires_random_access_iterator<InputIterator>,
		typename = tulib_core::requires_tulib_point<Kernel,
		typename InputIterator::value_type>>
		NT operator()(InputIterator polyline_a_first, InputIterator polyline_a_beyond,
			InputIterator polyline_b_first, InputIterator polyline_b_beyond) {
		std::size_t size_polyline_b = std::distance(polyline_b_first, polyline_b_beyond);
		std::vector<NT> dp_row(size_polyline_b);
		std::fill(std::begin(dp_row), std::begin(dp_row) + size_polyline_b, -1);
		InputIterator it_a = polyline_a_first;
		while (it_a != polyline_a_beyond) {
			InputIterator it_b = polyline_b_first;
			std::size_t j = 0;
			typename Kernel::NT previous = -1, current = -1;
			while (it_b != polyline_b_beyond) {
				if ((it_a == polyline_a_first) && (it_b == polyline_b_first)) {
					current = distance(it_a, it_b);
					dp_row[0] = current;
				}
				else if ((it_a != polyline_a_first) && (it_b == polyline_b_first)) {
					current = std::max(dp_row[0], distance(it_a, it_b));
					dp_row[0] = current;
				}
				else if ((it_a == polyline_a_first) && (it_b != polyline_b_first)) {
					current = std::max(previous, distance(it_a, it_b));
					dp_row[j - 1] = previous;
				}
				else {
					current = std::max(std::min(
						std::min(dp_row[j], previous),
						dp_row[j - 1]),
						distance(it_a, it_b));
					dp_row[j - 1] = previous;
				}
				previous = current;
				j++;
				it_b++;
			}
			dp_row[j - 1] = previous;
			it_a++;
		}
		NT dfd = dp_row.back();
		NT n = 1 / static_cast<NT>(Norm::P);
		return  std::pow(dfd, n);
	}
};