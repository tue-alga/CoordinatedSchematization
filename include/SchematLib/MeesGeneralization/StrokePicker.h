#pragma once
#include <algorithm>
#include <map>
#include <vector>
namespace SchematLib::MeesGeneralization {
template <class NT>
class StrokePicker {
public:
	void Solution1(std::vector<std::map<size_t, size_t>> const &trafficflows,
	               std::vector<NT> const &base_scores,
	               std::vector<NT> const &updated_scores,
	               std::vector<NT> const &stroke_lengths,
	               NT budget,
	               std::vector<bool> &stroke_taken) {
		// Tuple of score, length and index in list
		std::vector<std::tuple<NT, NT, size_t>> sortable_base_scores;
		auto len_it = stroke_lengths.begin();
		NT remaining_budget = budget;
		stroke_taken = std::vector<bool>(base_scores.size(), false);
		for (auto score_it = base_scores.begin(); score_it != base_scores.end(); score_it++) {
			sortable_base_scores.push_back(std::make_tuple(*score_it, *len_it, distance(base_scores.begin(), score_it)));
			++len_it;
		}
		std::sort(sortable_base_scores.begin(), sortable_base_scores.end(), std::greater<std::tuple<NT, NT, size_t>>());

		std::vector<std::tuple<NT, NT, size_t>> sortable_neighborhood;
		// Get the first stroke by solely looking at GP.
		bool selected_top = false;
		auto ss_it = sortable_base_scores.begin();  // Create the base score iterator
		for (; ss_it != sortable_base_scores.end(); ss_it++) {
			size_t curr_id = distance(sortable_base_scores.begin(), ss_it);

			const auto [score, length, index] = *ss_it;
			if (length <= remaining_budget) {
				if (selected_top) {
					sortable_neighborhood.push_back(make_tuple(base_scores[curr_id], stroke_lengths[curr_id], curr_id));
					++ss_it;  // Make sure that the iterator is at the next unpicked stroke
					break;
				}
				remaining_budget -= length;
				stroke_taken[index] = true;
				selected_top = true;
				for (auto tf : trafficflows[curr_id]) {
					size_t outstroke = tf.first;
					sortable_neighborhood.push_back(make_tuple(updated_scores[outstroke], stroke_lengths[outstroke], outstroke));
				}
			}
		}

		while (sortable_neighborhood.size() > 0) {
			std::sort(sortable_neighborhood.begin(), sortable_neighborhood.end(), std::greater<std::tuple<NT, NT, size_t>>());
			std::vector<std::tuple<NT, NT, size_t>> new_candidates;
			// Find the best stroke in the candidate set
			for (auto sn_it = sortable_neighborhood.begin(); sn_it != sortable_neighborhood.end(); sn_it++) {
				const auto [score, length, curr_id] = *sn_it;
				// size_t curr_id = get<2>(*sn_it);
				if (length <= remaining_budget && !stroke_taken[curr_id]) {
					remaining_budget -= length;
					stroke_taken[curr_id] = true;

					for (auto downstroke : trafficflows[curr_id]) {
						auto down_id = downstroke.first;
						new_candidates.push_back(make_tuple(updated_scores[down_id], stroke_lengths[down_id], down_id));
					}
					break;
				}
			}
			// Add the unpicked stroke with the highest GP
			for (ss_it; ss_it != sortable_base_scores.end(); ss_it++) {
				size_t curr_id = get<2>(*ss_it);
				if (((get<1>(*ss_it)) <= remaining_budget) && (!stroke_taken[curr_id])) {
					new_candidates.push_back(make_tuple(base_scores[curr_id], stroke_lengths[curr_id], curr_id));
					ss_it++;
					break;
				}
			}
			sortable_neighborhood = new_candidates;
		}
	}

	void Solution1fixed(std::vector<std::map<size_t, size_t>> const &trafficflows,
	                    std::vector<NT> const &base_scores,
	                    std::vector<NT> const &stroke_lengths,
	                    NT budget,
	                    std::vector<NT> const &trafficsums,
	                    std::vector<bool> &stroke_taken) {
		std::vector<std::tuple<NT, NT, size_t>> sortable_scores;  // Tuple contains score, length, and the index
		auto len_it = stroke_lengths.begin();
		NT remaining_budget = budget;
		stroke_taken = std::vector<bool>(base_scores.size(), false);
		std::vector<NT> influence =
		    std::vector<NT>(base_scores.size(), 0);  // Downstroke influence will be collected in this vector

		// Initially there is no influence, so we only consider the base scores.
		for (auto score_it = base_scores.begin(); score_it != base_scores.end(); score_it++) {
			sortable_scores.push_back(std::make_tuple(*score_it, *len_it, distance(base_scores.begin(), score_it)));
			len_it++;
		}
		std::sort(sortable_scores.begin(),
		          sortable_scores.end(),
		          std::greater<std::tuple<NT, NT, size_t>>());  // The strokes are now sorted on GP

		// Pick the first
		std::vector<std::tuple<NT, NT, size_t>> sortable_neighborhood;
		// Get the first stroke by solely looking at GP.
		bool selected_top = false;
		auto ss_it = sortable_scores.begin();  // Create the base score iterator
		for (ss_it; ss_it != sortable_scores.end(); ss_it++) {
			// size_t curr_id = distance(sortable_scores.begin(), ss_it);
			size_t curr_id = get<2>(*ss_it);
			if (get<1>(*ss_it) <= remaining_budget) {
				if (selected_top) {
					// Besides the downstrokes of the selected stroke, we also add the second-best stroke on GP.
					sortable_neighborhood.push_back(
					    make_tuple(base_scores[curr_id] + influence[curr_id], stroke_lengths[curr_id], curr_id));
					ss_it++;  // Make sure that the iterator is at the next unpicked stroke
					break;
				} else {
					remaining_budget -= get<1>(*ss_it);
					stroke_taken[curr_id] = true;
					selected_top = true;
					// When we select a stroke, we add its influence to its outstrokes, and then add them to the neighborhood

					// for (int i = 0; i < trafficflows.size(); i++)
					//{
					for (auto outflow_it = trafficflows[curr_id].begin(); outflow_it != trafficflows[curr_id].end();
					     outflow_it++) {
						auto outstroke_id = get<0>(*outflow_it);
						auto n = get<1>(*outflow_it);
						influence[outstroke_id] += (n / trafficsums[curr_id]) * base_scores[curr_id];
						sortable_neighborhood.push_back(make_tuple(
						    base_scores[outstroke_id] + influence[outstroke_id], stroke_lengths[outstroke_id], outstroke_id));
					}
					//}
					// for (auto tf : trafficflows[curr_id])
					//{
					//	size_t outstroke = get<0>(tf);
					//	sortable_neighborhood.push_back(make_tuple(updated_scores[outstroke], stroke_lengths[outstroke],
					// outstroke));
					//}
				}
			}
		}

		while (sortable_neighborhood.size() > 0) {
			std::sort(sortable_neighborhood.begin(), sortable_neighborhood.end(), std::greater<std::tuple<NT, NT, size_t>>());
			std::vector<std::tuple<NT, NT, size_t>> new_candidates;
			// Find the best stroke in the candidate set
			for (auto sn_it = sortable_neighborhood.begin(); sn_it != sortable_neighborhood.end(); sn_it++) {
				size_t curr_id = get<2>(*sn_it);
				if (get<1>(*sn_it) <= remaining_budget && !stroke_taken[curr_id]) {
					remaining_budget -= get<1>(*sn_it);
					stroke_taken[curr_id] = true;

					for (auto downstroke : trafficflows[curr_id]) {
						// Add the influence of the stroke
						auto outstroke_id = get<0>(downstroke);
						auto n = get<1>(downstroke);
						influence[outstroke_id] += (n / trafficsums[curr_id]) * base_scores[curr_id];
						new_candidates.push_back(make_tuple(
						    base_scores[outstroke_id] + influence[outstroke_id], stroke_lengths[outstroke_id], outstroke_id));

						// auto down_id = get<0>(downstroke);

						// new_candidates.push_back(make_tuple(updated_scores[down_id], stroke_lengths[down_id], down_id));
					}
					break;
				}
			}

			// Remake the sortable scores set
			sortable_scores.clear();
			len_it = stroke_lengths.begin();
			auto inf_it = influence.begin();
			for (auto score_it = base_scores.begin(); score_it != base_scores.end(); score_it++) {
				sortable_scores.push_back(
				    std::make_tuple(*score_it + *inf_it, *len_it, distance(base_scores.begin(), score_it)));
				len_it++;
				inf_it++;
			}
			std::sort(sortable_scores.begin(),
			          sortable_scores.end(),
			          std::greater<std::tuple<NT, NT, size_t>>());  // The strokes are now sorted on (partial) IP

			// Add the unpicked stroke with the highest GP
			for (ss_it = sortable_scores.begin(); ss_it != sortable_scores.end(); ss_it++) {
				size_t curr_id = get<2>(*ss_it);
				if (((get<1>(*ss_it)) <= remaining_budget) && (!stroke_taken[curr_id])) {
					new_candidates.push_back(
					    make_tuple(base_scores[curr_id] + influence[curr_id], stroke_lengths[curr_id], curr_id));
					ss_it++;
					break;
				}
			}
			sortable_neighborhood = new_candidates;
		}
	}

	void Solution2(const std::vector<NT> &stroke_scores,
	               const std::vector<NT> &stroke_lengths,
	               NT budget,
	               std::vector<bool> &stroke_taken) {
		const auto stroke_count = stroke_scores.size();
		std::vector<std::tuple<NT, NT, size_t>> sortable_scores;
		NT remaining_budget = budget;
		stroke_taken = std::vector<bool>(stroke_count, false);
		for (std::size_t i = 0; i < stroke_count; ++i) {
			sortable_scores.push_back(std::make_tuple(stroke_scores[i], stroke_lengths[i], i));
		}
		// Greedily select highest score strokes first.
		std::sort(sortable_scores.begin(), sortable_scores.end(), std::greater<std::tuple<NT, NT, size_t>>());
		for (const auto &[score, length, id] : sortable_scores) {
			if (length > remaining_budget)
				continue;
			remaining_budget -= length;
			stroke_taken[id] = true;
		}
	}
};
}  // namespace SchematLib::MeesGeneralization