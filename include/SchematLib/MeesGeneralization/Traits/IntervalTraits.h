/*! @file IntervalTraits.h
 *  @brief A collection of traits for Nontransitive Physics Outlier Detection
 *  @authors Mees van de Kerkhof (m.a.vandekerkhof@uu.nl)
 */

#ifndef TULIB_INTERVALTRAITS_H
#define TULIB_INTERVALTRAITS_H

#include <vector>
#include <stack>
#include "tulib/geo/geo.h"

namespace tulib_algorithms {

	template <class NT>
	struct SpeedInterval
	{
		NT min, max;
		__int64 prev; //Maintain the index of the previous succesful node for fast back-propagation. Is signed so -1 can imply no predecessor
	public:
		SpeedInterval(NT InMin, NT InMax, __int64 InPrev) { min = InMin; max = InMax; prev = InPrev; }
	};
	template <class Interval>
	struct IntervalComparer
	{
		bool operator()	(Interval i1, Interval i2) {
			if (i1.min == i2.min)
				return i1.max > i2.max;
			else return i1.min < i2.min;
		}
	};

	template <class Interval>
	class SpeedIntervalMerger 
	{
	public:
		////Better version of merge function which I can't get to work unfortunately
		//template <class InputIterator, class OutputIterator, typename = tulib_core::requires_random_access_iterator <InputIterator>,
		//	typename = tulib_core::requires_output_iterator <OutputIterator>,
		//	typename = tulib_core::requires_equality<typename InputIterator::value_type,
		//	typename OutputIterator::value_type::value_type >, typename = tulib_core::requires_equality<typename InputIterator::value_type, Interval>>
		//void operator() (InputIterator first, InputIterator beyond, OutputIterator result, bool skip_sorting = false)
		//{
		//	cout << "merging " << beyond - first << " intervals" << endl;
		//	if (first == beyond)
		//		return;
		//	//Merge algorithm adapted from https://www.geeksforgeeks.org/merging-intervals/
		//	if (!skip_sorting) //The merging algorithm assumes the intervals have been sorted on starting times. If this is not yet the case, they can be sorted here. This increases running time of this function from O(n) to O(n log n)
		//	{
		//		IntervalComparer<Interval> compare_interval;
		//		std::sort(first, beyond, compare_interval);
		//	}
		//	std::stack<Interval> s;
		//	s.push(*first);
		//	for (auto it = std::next(first, 1); it != beyond; it++)
		//	{
		//		auto top_interval = s.top();
		//		if (top_interval.max < it->min)
		//		{
		//			//Intervals do not overlap. New interval becomes the new top
		//			s.push(*it);
		//		}
		//		else if (top_interval.max < it->max) 
		//		{
		//			s.pop();
		//			//Intervals overlap and need to be merged; to help with backtracing we keep both, but remove the overlapping part. If they have the same prev reference, they are merged.
		//			if (top_interval.prev != it->prev)
		//			{
		//				auto new_interval = Interval(top_interval.min, it->min, top_interval.prev);
		//				s.push(new_interval);
		//				s.push(*it);
		//			}
		//			else
		//			{
		//				auto new_interval = Interval(top_interval.min, it->max, top_interval.prev);
		//				s.push(new_interval);
		//			}
		//		}
		//		//Otherwise this interval is redundant and should not be added
		//	}
		//	//Put the intervals in a vector in the correct order.
		//	while (!s.empty())
		//	{
		//		*result = s.top();
		//		s.pop();
		//	}
		//}

		
			std::vector<Interval> operator() (std::vector<Interval> input, bool skip_sorting = false)
			{
				std::vector<Interval> result;
				if (input.size() == 0)
					return result;

				//Merge algorithm adapted from https://www.geeksforgeeks.org/merging-intervals/
				if (!skip_sorting) //The merging algorithm assumes the intervals have been sorted on starting times. If this is not yet the case, they can be sorted here. This increases running time of this function from O(n) to O(n log n)
				{
				IntervalComparer<Interval> compare_interval;
				std::sort(input.begin(), input.end(), compare_interval);
				}

				std::stack<Interval> s;
				s.push(input.front());
				for (auto it = std::next(input.begin(), 1); it != input.end(); it++)
				{
					auto top_interval = s.top();
					if (top_interval.max < it->min)
					{
						//Intervals do not overlap. New interval becomes the new top
						s.push(*it);
					}
					else if (top_interval.max < it->max)
					{
						s.pop();
						//Intervals overlap and need to be merged; to help with backtracing we keep both, but remove the overlapping part. If they have the same prev reference, they are merged.
						if (top_interval.prev != it->prev)
						{
							auto new_interval = Interval(top_interval.min, it->min, top_interval.prev);
							s.push(new_interval);
							s.push(*it);
						}
						else
						{
							auto new_interval = Interval(top_interval.min, it->max, top_interval.prev);
							s.push(new_interval);
						}
					}
				//Otherwise this interval is redundant and should not be added
				}

				//Put the intervals in a vector in the correct order.
				while (!s.empty())
				{
					result.insert(result.begin(), s.top());
					s.pop();
				}
				return result;
			}
	};

	template <class Interval>
	class SpeedIntervalBacktracker
	{
	public:
		size_t operator()(std::vector<Interval> input)
		{
			return input.front().prev;
		}
	};

	template <class NT, class Interval>
	class BSpeedBAccPropagator 
	{
	private:
		NT min_speed, max_speed, min_acc, abs_min_acc, max_acc, detourfactor;
		NT sqrt_min_acc, sqrt_max_acc, sqrt_abs_min_acc, sqrt_2;
	public:
		BSpeedBAccPropagator<NT, Interval>() {};
		BSpeedBAccPropagator<NT, Interval>(NT InMin_speed, NT InMax_speed, NT InMin_acc, NT InMax_acc, NT InDetourFactor = 1)
		{
			min_speed = InMin_speed;
			max_speed = InMax_speed;
			min_acc = InMin_acc;
			abs_min_acc = abs(min_acc);
			max_acc = InMax_acc;
			detourfactor = InDetourFactor;
			sqrt_min_acc = sqrt(min_acc);
			sqrt_abs_min_acc = sqrt(abs_min_acc);
			sqrt_max_acc = sqrt(max_acc);
			sqrt_2 = sqrt(2);
		}

		//[Deprecated] Propagator with a distance map as argument, so only new distances get computed. 
		std::vector<Interval> operator()(std::vector<Interval> input, tuple<NT,NT,time_t> u, tuple<NT,NT, time_t> v, __int64 prev_index, __int64 next_index, map<tuple<__int64, __int64>, NT> &distancemap)
		{
			std::vector<Interval> result;
			NT min_distance;
			tuple<__int64, __int64> key(prev_index, next_index);
			auto dist = distancemap.find(key);
			if (dist != distancemap.end())
			{
				min_distance = get<1>(*dist);
			}
			else {
				NT u_x = get<0>(u);
				NT u_y = get<1>(u);

				NT v_x = get<0>(v);
				NT v_y = get<1>(v);
				min_distance = distance_exact(u_x, u_y, v_x, v_y);
				distancemap[key] = min_distance;
			}
			time_t u_t = get<2>(u);
			time_t v_t = get<2>(v);
			//Determine the trunc velocities, outside of which the next probe is definitely unreachable.
			NT max_distance = min_distance * detourfactor;
			auto time_lapsed = v_t - u_t;
			NT acc_delta = max_acc - min_acc;
			NT min_avg_speed = min_distance / time_lapsed;
			NT max_avg_speed = max_distance / time_lapsed;
			//Test if driving full duration at max speed is enough (otherwise the points are incompatible no matter the speed)
			if (time_lapsed * max_speed < min_distance)
				return result;
			NT minimum_consistent_speed = max_speed - sqrt(2)*sqrt(max_acc)*sqrt(time_lapsed * max_speed - min_distance);
			//Test if this case actually applies. 
			if ((max_speed - minimum_consistent_speed) / max_acc > time_lapsed)
				minimum_consistent_speed = (2 * min_distance - max_acc * time_lapsed*time_lapsed) / (2 * time_lapsed);
			//Test if driving full duration at min speed does not exceed (otherwise the points are incompatible no matter the speed)
			if (min_speed * time_lapsed > max_distance)
				return result;
			NT abs_min_acc = abs(min_acc);
			NT maximum_consistent_speed = min_speed + sqrt(2)*sqrt(abs_min_acc)*sqrt(max_distance - time_lapsed * min_speed);
			//Test if this case actually applies. 
			if ((min_speed - maximum_consistent_speed) / min_acc > time_lapsed)
				maximum_consistent_speed = (2 * max_distance - abs_min_acc * time_lapsed * time_lapsed - 2 * min_acc * time_lapsed * time_lapsed) / (2 * time_lapsed);
			if (minimum_consistent_speed > maximum_consistent_speed)
				return result;
			for (auto it = input.begin(); it != input.end(); it++)
			{
				//Truncate the interval to allowable values
				NT i_min = max(it->min, minimum_consistent_speed);
				NT i_max = min(it->max, maximum_consistent_speed);
				if (i_min > i_max) //If the truncated interval is empty there is no need to propagate
					continue; 

				//Either one of the 2 special cases applies, or the min/max of the interval are the min/max possible speed;
				NT new_min_speed = min_speed;
				NT new_max_speed = max_speed;

				bool done = false;
				//Find the minimum reachable speed at v_t if we leave u with the highest speed in this interval.
				//Case 1: Accelerate until some time t_star, then decellerate until v_t
				if (min_avg_speed - 0.5*max_acc*time_lapsed <= i_max && i_max <= min_avg_speed - 0.5*min_acc*time_lapsed)
				{
					NT t_star = v_t - sqrt((max_acc*time_lapsed*time_lapsed + 2 * time_lapsed*i_max - 2 * min_distance) / acc_delta);
					if (u_t <= t_star && t_star <= v_t && i_max + (t_star - u_t)*max_acc <= max_speed)
					{
						new_min_speed = i_max + time_lapsed * min_acc + acc_delta * (t_star - u_t);
						if (new_min_speed < min_speed)
						{
							new_min_speed = min_speed;
						}
						else
							done = true;
					}
				}

				//Case 2: Accelerate until top speed, keep going at max_speed until t_star, then decelerate until v_t
				if (!done && i_max >= max_speed - sqrt(2 * max_acc*(max_speed*time_lapsed - min_distance)) && i_max <= max_speed && i_max >= max_speed - max_acc * time_lapsed)
				{
					NT t_star = v_t - sqrt(((max_speed - i_max)*(max_speed - i_max) + 2 * max_acc*(min_distance - max_speed * time_lapsed)) / (max_acc*min_acc));
					if (u_t <= t_star && t_star <= v_t)
					{
						new_min_speed = max_speed + min_acc * (v_t - t_star);
						if (new_min_speed < min_speed)
						{
							new_min_speed = min_speed;
						}
					}
				}

				done = false;
				//Find the maximum reachable speed at v_t if we leave u with the lowest speed in this interval.
				//Case 1: Decelerate until some time t_star, then acccellerate until v_t
				if (max_avg_speed - 0.5*max_acc*time_lapsed <= i_min && i_min <= max_avg_speed - 0.5*min_acc*time_lapsed)
				{
					NT t_star = v_t - sqrt((min_acc*time_lapsed*time_lapsed + 2 * time_lapsed*i_min - 2 * max_distance) / (-1* acc_delta));
					if (u_t <= t_star && t_star <= v_t && i_min + (t_star - u_t) * min_acc >= min_speed)
					{
						new_max_speed = i_min + time_lapsed * max_acc + -1 * acc_delta * (t_star - u_t);
						if (new_max_speed > max_speed)
						{
							new_max_speed = max_speed;
						}
						else
							done = true;
					}
				}

				//Case 2: Decelerate until min speed, keep going at min_speed until t_star, then accelerate until v_t
				if (!done && min_speed <= i_min && i_min <= min_speed - min_acc * time_lapsed)
				{
					NT t_star = v_t - sqrt(((min_speed - i_min)*(min_speed - i_min) + 2 * min_acc*(max_distance - min_speed * time_lapsed)) / (max_acc*min_acc));
					if (u_t <= t_star && t_star <= v_t)
					{
						new_max_speed = min_speed + max_acc * (v_t - t_star);
						if (new_max_speed > max_speed)
						{
							new_max_speed = max_speed;
						}
					}
				}
				Interval new_i(new_min_speed, new_max_speed, prev_index);
				result.push_back(new_i);
			}
			return result;
		}

		//Version with precomputed distance; this is the version used by the greedy heuristics.
		std::vector<Interval> operator()(std::vector<Interval> input, tuple<NT, NT, time_t> u, tuple<NT, NT, time_t> v, __int64 prev_index, NT min_distance)
		{
			std::vector<Interval> result;
			
			time_t u_t = get<2>(u);
			time_t v_t = get<2>(v);
			//Determine the trunc velocities, outside of which the next probe is definitely unreachable.
			NT max_distance = min_distance * detourfactor;
			auto time_lapsed = v_t - u_t;
			NT acc_delta = max_acc - min_acc;
			NT min_avg_speed = min_distance / time_lapsed;
			NT max_avg_speed = max_distance / time_lapsed;
			//Test if driving full duration at max speed is enough (otherwise the points are incompatible no matter the speed)
			if (time_lapsed * max_speed < min_distance)
				return result;
			NT minimum_consistent_speed = max_speed - sqrt_2*sqrt_max_acc*sqrt(time_lapsed * max_speed - min_distance);
			//Test if this case actually applies. 
			if ((max_speed - minimum_consistent_speed) / max_acc > time_lapsed)
				minimum_consistent_speed = (2 * min_distance - max_acc * time_lapsed*time_lapsed) / (2 * time_lapsed);
			//Test if driving full duration at min speed does not exceed (otherwise the points are incompatible no matter the speed)
			if (min_speed * time_lapsed > max_distance)
				return result;
			NT maximum_consistent_speed = min_speed + sqrt_2*sqrt_abs_min_acc*sqrt(max_distance - time_lapsed * min_speed);
			//Test if this case actually applies. 
			if ((min_speed - maximum_consistent_speed) / min_acc > time_lapsed)
				maximum_consistent_speed = (2 * max_distance - abs_min_acc * time_lapsed * time_lapsed - 2 * min_acc * time_lapsed * time_lapsed) / (2 * time_lapsed);
			
			if (minimum_consistent_speed > maximum_consistent_speed)
			{
				throw std::runtime_error("minimum_consistent_speed is bigger than maximum_consistent_speed");
				return result;
			}
			for (auto it = input.begin(); it != input.end(); it++)
			{
				//Truncate the interval to allowable values
				NT i_min = max(it->min, minimum_consistent_speed);
				NT i_max = min(it->max, maximum_consistent_speed);
				if (i_min > i_max) //If the truncated interval is empty there is no need to propagate
					continue;

				//Either one of the 2 special cases applies, or the min/max of the interval are the min/max possible speed;
				NT new_min_speed = min_speed;
				NT new_max_speed = max_speed;

				bool done = false;
				//Find the minimum reachable speed at v_t if we leave u with the highest speed in this interval.
				//Case 1: Accelerate until some time t_star, then decellerate until v_t
				if (min_avg_speed - 0.5*max_acc*time_lapsed <= i_max && i_max <= min_avg_speed - 0.5*min_acc*time_lapsed)
				{
					NT t_star = v_t - sqrt((max_acc*time_lapsed*time_lapsed + 2 * time_lapsed*i_max - 2 * min_distance) / acc_delta);
					if (u_t <= t_star && t_star <= v_t && i_max + (t_star - u_t)*max_acc <= max_speed)
					{
						
						new_min_speed = i_max + time_lapsed * min_acc + acc_delta * (t_star - u_t);

						if (new_min_speed < min_speed)
						{
							new_min_speed = min_speed;
						}
						else
							done = true;
					}
				}

				//Case 2: Accelerate until top speed, keep going at max_speed until t_star, then decelerate until v_t
				if (!done && i_max >= max_speed - sqrt(2 * max_acc*(max_speed*time_lapsed - min_distance)) && i_max <= max_speed && i_max >= max_speed - max_acc * time_lapsed)
				{
					NT t_star = v_t - sqrt(((max_speed - i_max)*(max_speed - i_max) + 2 * max_acc*(min_distance - max_speed * time_lapsed)) / (max_acc*min_acc));
					if (u_t <= t_star && t_star <= v_t)
					{
						new_min_speed = max_speed + min_acc * (v_t - t_star);

						if (new_min_speed < min_speed)
						{
							new_min_speed = min_speed;
						}
					}
				}

				done = false;
				//Find the maximum reachable speed at v_t if we leave u with the lowest speed in this interval.
				//Case 1: Decelerate until some time t_star, then acccellerate until v_t
				if (max_avg_speed - 0.5*max_acc*time_lapsed <= i_min && i_min <= max_avg_speed - 0.5*min_acc*time_lapsed)
				{
					NT t_star = v_t - sqrt((min_acc*time_lapsed*time_lapsed + 2 * time_lapsed*i_min - 2 * max_distance) / (-1 * acc_delta));
					if (u_t <= t_star && t_star <= v_t && i_min + (t_star - u_t) * min_acc >= min_speed)
					{
						
						new_max_speed = i_min + time_lapsed * max_acc + -1 * acc_delta * (t_star - u_t);

						if (new_max_speed > max_speed)
						{
							new_max_speed = max_speed;
						}
						else
							done = true;
					}
				}

				//Case 2: Decelerate until min speed, keep going at min_speed until t_star, then accelerate until v_t
				if (!done && min_speed <= i_min && i_min <= min_speed - min_acc * time_lapsed)
				{
					NT t_star = v_t - sqrt(((min_speed - i_min)*(min_speed - i_min) + 2 * min_acc*(max_distance - min_speed * time_lapsed)) / (max_acc*min_acc));
					if (u_t <= t_star && t_star <= v_t)
					{
						new_max_speed = min_speed + max_acc * (v_t - t_star);

						if (new_max_speed > max_speed)
						{
							new_max_speed = max_speed;
						}
					}
				}
				Interval new_i(new_min_speed, new_max_speed, prev_index);
				result.push_back(new_i);
			}
			return result;
		}
		
		//Version with precomputed distance and square roots, this is the version used by the optimal algorithm
		std::vector<Interval> operator()(std::vector<Interval> input, tuple<NT, NT, time_t> u, tuple<NT, NT, time_t> v, __int64 prev_index, NT min_distance, NT min_const_sqrt, NT max_const_sqrt)
		{
			std::vector<Interval> result;

			time_t u_t = get<2>(u);
			time_t v_t = get<2>(v);
			//Determine the trunc velocities, outside of which the next probe is definitely unreachable.
			NT max_distance = min_distance * detourfactor;
			auto time_lapsed = v_t - u_t;
			NT acc_delta = max_acc - min_acc;
			NT min_avg_speed = min_distance / time_lapsed;
			NT max_avg_speed = max_distance / time_lapsed;
			//Test if driving full duration at max speed is enough (otherwise the points are incompatible no matter the speed)
			if (time_lapsed * max_speed < min_distance)
				return result;
			NT minimum_consistent_speed = max_speed - sqrt_2 * sqrt_max_acc*min_const_sqrt;
			//Test if this case actually applies. 
			if ((max_speed - minimum_consistent_speed) / max_acc > time_lapsed)
				minimum_consistent_speed = (2 * min_distance - max_acc * time_lapsed*time_lapsed) / (2 * time_lapsed);
			//Test if driving full duration at min speed does not exceed (otherwise the points are incompatible no matter the speed)
			if (min_speed * time_lapsed > max_distance)
				return result;
			NT maximum_consistent_speed = min_speed + sqrt_2 * sqrt_abs_min_acc*max_const_sqrt;
			//Test if this case actually applies. 
			if ((min_speed - maximum_consistent_speed) / min_acc > time_lapsed)
				maximum_consistent_speed = (2 * max_distance - abs_min_acc * time_lapsed * time_lapsed - 2 * min_acc * time_lapsed * time_lapsed) / (2 * time_lapsed);
			if (minimum_consistent_speed > maximum_consistent_speed)
				return result;
			for (auto it = input.begin(); it != input.end(); it++)
			{
				//Truncate the interval to allowable values
				NT i_min = max(it->min, minimum_consistent_speed);
				NT i_max = min(it->max, maximum_consistent_speed);
				if (i_min > i_max) //If the truncated interval is empty there is no need to propagate
					continue;

				//Either one of the 2 special cases applies, or the min/max of the interval are the min/max possible speed;
				NT new_min_speed = min_speed;
				NT new_max_speed = max_speed;

				bool done = false;
				//Find the minimum reachable speed at v_t if we leave u with the highest speed in this interval.
				//Case 1: Accelerate until some time t_star, then decellerate until v_t
				if (min_avg_speed - 0.5*max_acc*time_lapsed <= i_max && i_max <= min_avg_speed - 0.5*min_acc*time_lapsed)
				{
					NT t_star = v_t - sqrt((max_acc*time_lapsed*time_lapsed + 2 * time_lapsed*i_max - 2 * min_distance) / acc_delta);
					if (u_t <= t_star && t_star <= v_t && i_max + (t_star - u_t)*max_acc <= max_speed)
					{
						new_min_speed = i_max + time_lapsed * min_acc + acc_delta * (t_star - u_t);
						if (new_min_speed < min_speed)
						{
							new_min_speed = min_speed;
						}
						else
							done = true;
					}
				}

				//Case 2: Accelerate until top speed, keep going at max_speed until t_star, then decelerate until v_t
				if (!done && i_max >= max_speed - sqrt(2 * max_acc*(max_speed*time_lapsed - min_distance)) && i_max <= max_speed && i_max >= max_speed - max_acc * time_lapsed)
				{
					NT t_star = v_t - sqrt(((max_speed - i_max)*(max_speed - i_max) + 2 * max_acc*(min_distance - max_speed * time_lapsed)) / (max_acc*min_acc));
					if (u_t <= t_star && t_star <= v_t)
					{
						new_min_speed = max_speed + min_acc * (v_t - t_star);
						if (new_min_speed < min_speed)
						{
							new_min_speed = min_speed;
						}
					}
				}

				done = false;
				//Find the maximum reachable speed at v_t if we leave u with the lowest speed in this interval.
				//Case 1: Decelerate until some time t_star, then acccellerate until v_t
				if (max_avg_speed - 0.5*max_acc*time_lapsed <= i_min && i_min <= max_avg_speed - 0.5*min_acc*time_lapsed)
				{
					NT t_star = v_t - sqrt((min_acc*time_lapsed*time_lapsed + 2 * time_lapsed*i_min - 2 * max_distance) / (-1 * acc_delta));
					if (u_t <= t_star && t_star <= v_t && i_min + (t_star - u_t) * min_acc >= min_speed)
					{
						new_max_speed = i_min + time_lapsed * max_acc + -1 * acc_delta * (t_star - u_t);
						if (new_max_speed > max_speed)
						{
							new_max_speed = max_speed;
						}
						else
							done = true;
					}
				}

				//Case 2: Decelerate until min speed, keep going at min_speed until t_star, then accelerate until v_t
				if (!done && min_speed <= i_min && i_min <= min_speed - min_acc * time_lapsed)
				{
					NT t_star = v_t - sqrt(((min_speed - i_min)*(min_speed - i_min) + 2 * min_acc*(max_distance - min_speed * time_lapsed)) / (max_acc*min_acc));
					if (u_t <= t_star && t_star <= v_t)
					{
						new_max_speed = min_speed + max_acc * (v_t - t_star);
						if (new_max_speed > max_speed)
						{
							new_max_speed = max_speed;
						}
					}
				}
				Interval new_i(new_min_speed, new_max_speed, prev_index);
				result.push_back(new_i);
			}
			return result;
		}
	};

	template <class NT, class Interval>
	class BSpeedBAccBackPropagator
	{
		NT min_speed, max_speed, min_acc, abs_min_acc, max_acc, detourfactor, sqrt_min_acc, sqrt_abs_min_acc, sqrt_max_acc, sqrt_2;
	public:
		BSpeedBAccBackPropagator<NT, Interval>() {};
		BSpeedBAccBackPropagator<NT, Interval>(NT InMin_speed, NT InMax_speed, NT InMin_acc, NT InMax_acc, NT InDetourFactor = 1)
		{
			min_speed = InMin_speed;
			max_speed = InMax_speed;
			min_acc = InMin_acc;
			abs_min_acc = abs(min_acc);
			max_acc = InMax_acc;
			detourfactor = InDetourFactor;
			sqrt_min_acc = sqrt(min_acc);
			sqrt_abs_min_acc = sqrt(abs_min_acc);
			sqrt_max_acc = sqrt(max_acc);
			sqrt_2 = sqrt(2);
		}
		Interval operator()(Interval input, tuple<NT, NT, time_t> u, tuple<NT, NT, time_t> v, __int64 prev_index, NT min_distance,NT slack)
		{
			time_t u_t = get<2>(u);
			time_t v_t = get<2>(v);
			auto time_lapsed = v_t - u_t;
			NT max_distance = min_distance * detourfactor;
			auto outmax = input.max;
			auto outmin = input.min;
			auto acc_delta = max_acc - min_acc;
			
			auto vOut = outmin;
			auto distance = max_distance;

			NT orig_max = min_speed;
			//get the upper bound on the interval
			bool good_t_star = false;

			if (vOut == min_speed)
			{
			bool check = false;
			if (2 * distance*abs_min_acc - 2 * time_lapsed*min_speed*abs_min_acc >= 0)
			{
				orig_max = min_speed + sqrt(2 * abs_min_acc*distance - 2 * abs_min_acc*time_lapsed*min_speed);
				if ((min_speed - orig_max) / min_acc > time_lapsed)
					orig_max = (2 * distance - abs_min_acc * time_lapsed * time_lapsed - 2 * min_acc * time_lapsed * time_lapsed) / (2 * time_lapsed);
				check = true;
			}
			if (0.5*abs_min_acc*time_lapsed*time_lapsed + min_speed * time_lapsed == distance)
			{
				orig_max = min_speed + abs_min_acc * time_lapsed;

				check = true;
			}
			if (!check)
				throw std::runtime_error("No cases apply during backpropagation");
			good_t_star = check;
			}
			
			if (!good_t_star && ((0.5 * (2 * vOut*time_lapsed - max_acc * time_lapsed*time_lapsed) - distance <= slack && distance - 0.5 * (2 * vOut * time_lapsed - min_acc * time_lapsed * time_lapsed <= slack)) || distance - 0.5 * (2 * vOut * time_lapsed - max_acc * time_lapsed * time_lapsed) <= slack))
			{
				if ((-2 * min_acc * vOut * time_lapsed + 2 * max_acc * vOut * time_lapsed + min_acc * min_acc * time_lapsed * time_lapsed - min_acc * max_acc*time_lapsed*time_lapsed + 2 * min_acc * distance - 2 * max_acc*distance) <= slack)
				{
					orig_max = vOut - min_acc * time_lapsed;
					NT t_star = time_lapsed - sqrt((max_acc*time_lapsed*time_lapsed + 2 * time_lapsed*orig_max - 2 * min_distance) / acc_delta);

					if (t_star >= 0 && orig_max + t_star * max_acc <= max_speed)
						good_t_star = true;
				}
				
				if (!good_t_star)
				{
					orig_max = vOut - min_acc * time_lapsed + sqrt(-2 * min_acc * vOut * time_lapsed + 2 * max_acc * vOut * time_lapsed + min_acc * min_acc * time_lapsed * time_lapsed - min_acc * max_acc*time_lapsed*time_lapsed + 2 * min_acc * distance - 2 * max_acc*distance);
					NT t_star = time_lapsed - sqrt((max_acc*time_lapsed*time_lapsed + 2 * time_lapsed*orig_max - 2 * min_distance) / acc_delta);
					if (t_star >= 0 && orig_max + t_star * max_acc <= max_speed)
						good_t_star = true;
				}
				if (!good_t_star || orig_max < min_speed || orig_max > max_speed)
				{
					NT other_max = vOut - min_acc * time_lapsed - sqrt(-2 * min_acc * vOut * time_lapsed + 2 * max_acc * vOut * time_lapsed + min_acc * min_acc * time_lapsed * time_lapsed - min_acc * max_acc*time_lapsed*time_lapsed + 2 * min_acc * distance - 2 * max_acc*distance);
					if (other_max > min_speed && other_max < max_speed)
					{
						orig_max = other_max;
						NT t_star = time_lapsed - sqrt((max_acc*time_lapsed*time_lapsed + 2 * time_lapsed*orig_max - 2 * min_distance) / acc_delta);
						if (t_star <= 0 && orig_max + t_star * max_acc <= max_speed)
							good_t_star = false;
					}

				}
			}
			if (
				!good_t_star &&
				(
				(max_acc >= ((min_speed * min_speed - 2 * min_speed * max_speed + max_speed * max_speed) / (2 * max_speed * time_lapsed - 2 * distance)) && vOut <= max_speed && min_speed <= vOut && distance <= max_speed * time_lapsed && ((-max_acc * vOut*vOut + 2 * max_acc*vOut*max_speed - max_acc * max_speed * max_speed) / (-min_speed * min_speed + 2 * min_speed * max_speed - max_speed * max_speed + 2 * max_acc * max_speed * time_lapsed - 2 * max_acc * distance)) <= min_acc && min_acc <= ((-vOut * vOut + 2 * vOut * max_speed - max_speed * max_speed) / (2 * max_speed*time_lapsed - 2 * distance)))
					||
					(vOut <= max_speed && min_speed <= vOut && distance <= max_speed * time_lapsed && max_acc > 0 && max_acc <= ((min_speed * min_speed - 2 * min_speed * max_speed + max_speed * max_speed) / (2 * max_speed * time_lapsed - 2 * distance)) && min_acc <= ((-vOut * vOut + 2 * vOut * max_speed - max_speed * max_speed) / (2 * max_speed * time_lapsed - 2 * distance)))
					||
					(vOut <= max_speed && min_speed <= vOut && distance <= max_speed * time_lapsed && ((-vOut * vOut + 2 * vOut * max_speed - max_speed * max_speed) / (2 * max_speed * time_lapsed - 2 * distance)) <= max_acc && max_acc < 0 && ((-vOut * vOut + 2 * vOut * max_speed - max_speed * max_speed) / (2 * max_speed * time_lapsed - 2 * distance)) <= min_acc && min_acc < max_acc)
					)
				)
			{
				NT squirt = sqrt((max_acc * vOut * vOut - 2 * max_acc * vOut * max_speed + max_acc * max_speed * max_speed + 2 * min_acc * max_acc * max_speed * time_lapsed - 2 * min_acc * max_acc * distance) / min_acc);
				if (squirt < slack)
					squirt = 0;

				orig_max = max_speed - squirt;
				NT t_star = time_lapsed - sqrt(((max_speed - orig_max)*(max_speed - orig_max) + 2 * max_acc*(min_distance - max_speed * time_lapsed)) / (max_acc*min_acc));
				if (t_star >= 0 && abs(min_speed + (time_lapsed - t_star) * max_acc - vOut) <= slack)
					good_t_star = true;
			}
			
			if (!good_t_star)
			{
				bool check = false;
				if (2 * distance*abs_min_acc - 2 * time_lapsed*min_speed*abs_min_acc >= 0)
				{
					orig_max = min_speed + sqrt(2 * abs_min_acc*distance - 2 * abs_min_acc*time_lapsed*min_speed);
					if ((min_speed - orig_max) / min_acc > time_lapsed)
					{
						orig_max = (2 * distance - abs_min_acc * time_lapsed * time_lapsed - 2 * min_acc * time_lapsed * time_lapsed) / (2 * time_lapsed);
					}

					check = true;
				}
				if (0.5*abs_min_acc*time_lapsed*time_lapsed + min_speed * time_lapsed == distance)
				{
					orig_max = min_speed + abs_min_acc * time_lapsed;
					check = true;
				}
				if (!check)
					throw std::runtime_error("No cases apply during backpropagation");
				good_t_star = check;
			}
			good_t_star = false;
			vOut = outmax;
			distance = min_distance;

			NT orig_min = max_speed + 1;
			if (vOut == max_speed)
			{
			orig_min = (2 * distance - time_lapsed * time_lapsed * max_acc) / (2 * time_lapsed);
			if ((max_speed - orig_min) / max_acc < time_lapsed)
			{
				orig_min = vOut - sqrt(-2 * max_acc * distance + 2 * time_lapsed *  vOut *  max_acc);
			}
			good_t_star = true;
			}

			if (!good_t_star && ((((1 / 2) * (2 * vOut * time_lapsed - max_acc * time_lapsed * time_lapsed) - distance <= slack && distance - (1 / 2) * (2 * vOut * time_lapsed - min_acc * time_lapsed * time_lapsed)) <= slack && min_acc < max_acc && time_lapsed > 0) || (min_acc < max_acc && time_lapsed > 0 && (1 / 2) * (2 * vOut * time_lapsed - min_acc * time_lapsed * time_lapsed)) - distance <= slack))
			{
				if (2 * min_acc * vOut * time_lapsed - 2 * max_acc * vOut * time_lapsed - min_acc * max_acc * time_lapsed * time_lapsed + max_acc * max_acc * time_lapsed * time_lapsed - 2 * min_acc * distance + 2 * max_acc * distance <= slack)
					orig_min = vOut - max_acc * time_lapsed;
				else
				orig_min = vOut - max_acc * time_lapsed - sqrt(2 * min_acc * vOut * time_lapsed - 2 * max_acc * vOut * time_lapsed - min_acc * max_acc * time_lapsed * time_lapsed + max_acc * max_acc * time_lapsed * time_lapsed - 2 * min_acc * distance + 2 * max_acc * distance);
				
				if (orig_min < min_speed || orig_min > max_speed)
				{
					NT other_min = vOut - max_acc * time_lapsed + sqrt(2 * min_acc * vOut * time_lapsed - 2 * max_acc * vOut * time_lapsed - min_acc * max_acc * time_lapsed * time_lapsed + max_acc * max_acc * time_lapsed * time_lapsed - 2 * min_acc * distance + 2 * max_acc * distance);
					if (other_min > min_speed && other_min < max_speed)
						orig_min = other_min;
				}
				NT t_star = time_lapsed - sqrt((min_acc*time_lapsed*time_lapsed + 2 * time_lapsed*orig_min - 2 * max_distance) / (-1 * acc_delta));
				if (t_star >= 0 && orig_min + t_star * min_acc >= min_speed)
					good_t_star = true;
			}
			if (!good_t_star && ((time_lapsed > 0 && max_acc >= (2 * min_speed * vOut - vOut * vOut - 2 * min_speed * max_speed + max_speed * max_speed) / (2 * min_speed * time_lapsed - 2 * distance) && vOut <= max_speed && min_speed <= vOut && distance <= min_speed * time_lapsed && 0 <= min_acc && min_acc <= (max_acc * min_speed * min_speed - 2 * max_acc * min_speed * max_speed + max_acc * max_speed * max_speed) / (min_speed * min_speed - 2 * min_speed * vOut + vOut * vOut + 2 * max_acc * min_speed * time_lapsed - 2 * max_acc * distance))
				||
				(time_lapsed > 0 && distance >= min_speed * time_lapsed && max_acc >= (-1 * min_speed * min_speed + 2 * min_speed * vOut - vOut * vOut) / (2 * min_speed * time_lapsed - 2 * distance) && vOut <= max_speed && min_speed <= vOut && (max_acc * min_speed * min_speed - 2 * max_acc * min_speed * max_speed + max_acc * max_speed * max_speed) / (min_speed * min_speed - 2 * min_speed * vOut + vOut * vOut + 2 * max_acc * min_speed * time_lapsed - 2 * max_acc * distance) <= min_acc && min_acc <= 0)
				||
				(time_lapsed > 0 && distance >= min_speed * time_lapsed && vOut <= max_speed && min_speed <= vOut && 0 <= min_acc && min_acc <= max_acc && 0 <= max_acc && max_acc <= (-1 * min_speed * min_speed + 2 * min_speed * vOut - vOut * vOut) / (2 * min_speed * time_lapsed - 2 * distance))
				||
				(time_lapsed >= 0 && vOut <= max_speed && min_speed <= vOut && distance <= min_speed * time_lapsed && 0 <= max_acc && max_acc <= (2 * min_speed * vOut - vOut * vOut - 2 * min_speed * max_speed + max_speed * max_speed) / (2 * min_speed * time_lapsed - 2 * distance) && 0 <= min_acc && min_acc <= max_acc)))
			{
				NT squirt = sqrt((min_acc * min_speed * min_speed - 2 * min_acc * min_speed * vOut + min_acc * vOut * vOut + 2 * min_acc * max_acc * min_speed * time_lapsed - 2 * min_acc * max_acc * distance) / max_acc);
				if (squirt < slack)
					squirt = 0;
				orig_min = min_speed + squirt;
				NT t_star = time_lapsed - sqrt(((min_speed - orig_min)*(min_speed - orig_min) + 2 * min_acc*(max_distance - min_speed * time_lapsed)) / (max_acc*min_acc));
				if (t_star >= 0 && abs(max_speed + min_acc * (time_lapsed - t_star) - vOut) <= slack)
					good_t_star = true;
			}

			if (!good_t_star)
			{
				orig_min = (2 * distance - time_lapsed * time_lapsed * max_acc) / (2 * time_lapsed);

				if ((max_speed - orig_min) / max_acc < time_lapsed)
				{
					orig_min = vOut - sqrt(-2 * max_acc * distance + 2 * time_lapsed *  vOut *  max_acc);
				}
				good_t_star = true;
			}
			if (orig_max > max_speed)
				orig_max = max_speed;
			if (orig_min < min_speed)
				orig_min = min_speed;
			return Interval(orig_min, orig_max, prev_index);
		}
	};
	template <class NT>
	class DistanceFunction
	{
	public:
		NT operator()(tuple<NT, NT, time_t> u, tuple<NT, NT, time_t> v)
		{
			auto ux = get<0>(u);
			auto uy = get<1>(u);
			auto vx = get<0>(v);
			auto vy = get<1>(v);
			return distance_exact(ux, uy, vx, vy);
		}
	};

	template <class NT>
	class MinSqrtPrecomputer
	{
	private:
		NT max_speed;
	public:
		MinSqrtPrecomputer<NT>() {};
		MinSqrtPrecomputer(NT InMaxSpeed)
		{
			max_speed = InMaxSpeed;
		}
	
		NT operator()(tuple<NT,NT,time_t> u, tuple<NT,NT,time_t> v, NT distance)
		{
			auto tu = get<2>(u);
			auto tv = get<2>(v);
			return sqrt((tv-tu)*max_speed - distance);
		}
	};

	template <class NT>
	class MaxSqrtPrecomputer
	{
	private:
		NT min_speed, detourfactor;
	public:
		MaxSqrtPrecomputer<NT>() {};
		MaxSqrtPrecomputer<NT>(NT InMinSpeed, NT InDetourFactor)
		{
			min_speed = InMinSpeed;
			detourfactor = InDetourFactor;
		}

		NT operator()(tuple<NT, NT, time_t> u, tuple<NT, NT, time_t> v, NT distance)
		{
			auto tu = get<2>(u);
			auto tv = get<2>(v);
			return sqrt(detourfactor * distance - (tv-tu) * min_speed);
		}
	};
}

template <class NT>
struct SpeedIntervalTraits {
	typedef tulib_algorithms::SpeedInterval<NT> Interval;
	tulib_algorithms::BSpeedBAccPropagator<NT, Interval> propagate;
	tulib_algorithms::SpeedIntervalMerger<Interval> merge;
	tulib_algorithms::SpeedIntervalBacktracker<Interval> backtrack;
	tulib_algorithms::DistanceFunction<NT> distance;
	tulib_algorithms::MinSqrtPrecomputer<NT> get_min_sqrt;
	tulib_algorithms::MaxSqrtPrecomputer<NT> get_max_sqrt;
	tulib_algorithms::BSpeedBAccBackPropagator<NT, Interval> back_propagate;
	NT min_speed, max_speed, min_acc, max_acc, detourfactor;

	SpeedIntervalTraits() {};
	SpeedIntervalTraits(NT InMin_speed, NT InMax_speed, NT InMin_acc, NT InMax_acc, NT InDetourfactor = 1)
	{
		min_speed = InMin_speed;
		max_speed = InMax_speed;
		min_acc = InMin_acc;
		max_acc = InMax_acc;
		detourfactor = InDetourfactor;
		propagate = tulib_algorithms::BSpeedBAccPropagator<NT,Interval>(min_speed, max_speed, min_acc, max_acc, detourfactor);
		back_propagate = tulib_algorithms::BSpeedBAccBackPropagator<NT, Interval>(min_speed, max_speed, min_acc, max_acc, detourfactor);
		merge = tulib_algorithms::SpeedIntervalMerger<Interval>();
		backtrack = tulib_algorithms::SpeedIntervalBacktracker<Interval>();
		distance = tulib_algorithms::DistanceFunction<NT>();
		get_min_sqrt = tulib_algorithms::MinSqrtPrecomputer<NT>(max_speed);
		get_max_sqrt = tulib_algorithms::MaxSqrtPrecomputer<NT>(min_speed, detourfactor);
	}

	Interval base_interval() { return Interval(min_speed, max_speed, -1); }
};

#endif //TULIB_INTERVALTRAITS_H