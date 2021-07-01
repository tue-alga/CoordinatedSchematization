#ifndef SCHEMATLIB_TRAJECTORYSCHEMATIZATION_BIDIRECTIONALBRUTEFORCETRAJECTORYBUNDLER_H
#define SCHEMATLIB_TRAJECTORYSCHEMATIZATION_BIDIRECTIONALBRUTEFORCETRAJECTORYBUNDLER_H
#include <GCpp/DS/BoostEmbeddedGraph.h>
#include <GCpp/DS/Heap.h>
#include <SchematLib/IO/settings_io.h>

#include <optional>
#include <queue>
#include <set>
#include <tuple>
#include <type_traits>
#include <vector>

#include "TrajectoryBundleTraits.h"

namespace SchematLib::TrajectorySchematization {
struct IsSuperSequence {
	template <typename It, typename It2>
	bool operator()(It seq0Start, It seq0End, It2 seq1Start, It2 seq1End) const {
		auto size0 = std::distance(seq0Start, seq0End);
		auto size1 = std::distance(seq1Start, seq1End);
		// First sequence is smaller, so cannot be super sequence
		if (size0 < size1)
			return false;
		// Find start in first sequence of second
		auto start0 = std::find(seq0Start, seq0End, *seq1Start);
		// If not found, or subsequence cannot fully match seq1, seq0 cannot be a supersequence.
		if (start0 == seq0End || std::distance(start0, seq0End) < size1)
			return false;
		// Run over the sequences, comparing elements
		auto curr1 = seq1Start;
		auto curr0 = start0;
		for (; curr0 != seq0End && curr1 != seq1End; ++curr0, ++curr1) {
			if (*curr0 != *curr1)
				return false;
		}
		return curr1 == seq1End;
	}
};

template <typename EmbeddedGraph>
class BidirectionalBruteforceTrajectoryBundler {
public:
	enum class WeighingScheme { Length, Support, LengthAndSupport, Custom };
	struct Settings {
		// std::size_t m_minimum_bundle_element_complexity = 3;

		std::size_t m_max_bundle_set_complexity = 10;

		Models::NT m_minimum_bundle_length = 1000;  // meters assumed

		// Minimum match count for first class.
		std::size_t m_min_matched_trajectories = 10;

		bool m_priotizie_by_length = false;

		bool m_verbose = false;

		bool m_weigh_edges_by_bundles = false;

		bool m_superbundles_only = false;

		// Capture overlapping trajectories that differ by at most the given amount
		bool m_match_intermediate_to_bundle = false;
		Models::NT m_max_intermediate_bundle_difference = 500;

		Models::NT m_weight_exponent = 1.0;

		WeighingScheme m_weighingScheme = WeighingScheme::Length;
	};

private:
	using SEF = MapSimplification::SpecialEdgeIdFunctions;

	struct TrajectoryEdge {
		std::size_t trajectory;
		std::size_t edgeNum;
		bool operator<(const TrajectoryEdge& other) const {
			return trajectory == other.trajectory ? edgeNum < other.edgeNum : trajectory < other.trajectory;
		}
		bool operator==(const TrajectoryEdge& other) const {
			return trajectory == other.trajectory && edgeNum == other.edgeNum;
		}
	};

	using Edge = typename EmbeddedGraph::edge_descriptor;
	using GT = boost::graph_traits<EmbeddedGraph>;
	using Vertex = typename EmbeddedGraph::vertex_descriptor;
	struct PathEdge {
		enum class Explore { Source, Target, Both };
		Edge edge;
		Explore toExplore;
		Vertex vertexToExplore(const EmbeddedGraph& graph) const {
			if (toExplore == Explore::Source)
				return boost::source(edge, graph);
			if (toExplore == Explore::Target)
				return boost::target(edge, graph);
			throw std::runtime_error("Vertex to explore only works for Source/Target");
		}
	};
	// Semantically tag the weight to avoid confusion with IDs of std::size_t
	struct EdgeWeight {
		std::size_t value;
		bool operator<(const EdgeWeight& other) const { return value < other.value; }
		bool operator>(const EdgeWeight& other) const { return value > other.value; }
		bool operator==(const EdgeWeight& other) const { return value == other.value; }
	};
	enum class Direction { Forward = 0, Backward = 1 };
	/**
	 * \brief Matcher object for a trajectory, in a particular direction.
	 */
	struct TrajectoryMatch {
		std::size_t trajectoryIndex;
		std::size_t startEdgeIndex;
		std::size_t endEdgeIndex;
		Models::NT pathLength = 0;

		bool matches(Direction direction,
		             std::size_t edgeId,
		             const std::vector<std::vector<std::size_t>>& trajectories) const {
			if (direction == Direction::Backward) {
				if (startEdgeIndex == 0)
					return false;
				return trajectories[trajectoryIndex][startEdgeIndex - 1] == edgeId;
			} else {
				if (endEdgeIndex == trajectories[trajectoryIndex].size() - 1)
					return false;
				return trajectories[trajectoryIndex][endEdgeIndex + 1] == edgeId;
			}
		}
		void extend(Direction direction) {
			if (direction == Direction::Backward)
				--startEdgeIndex;
			else
				++endEdgeIndex;
		}
	};

	using MatchingMap = std::map<std::size_t, TrajectoryMatch>;

	struct ComputationData {
		// Edge IDs of the given graph
		std::unordered_set<std::size_t> edgeIds;

		std::map<std::size_t, EdgeWeight> edgeWeights;
		// Indexed by edgeID, then trajectory index.
		std::map<std::size_t, std::map<std::size_t, TrajectoryEdge>> edgeIdToTrajectoryEdges;
		std::map<std::size_t, std::set<std::size_t>> edgeIdToTrajectoryNum;
		std::map<std::size_t, typename GT::edge_descriptor> edgeIdToEdgeMap;
		std::map<std::size_t, std::size_t> edgeIdToBundleCount;
		std::vector<Models::NT> m_trajectoryLengths;
		const EmbeddedGraph& graph;
		ComputationData(const EmbeddedGraph& inputGraph) : graph(inputGraph) { GCpp::DS::computeEdgeIdSet(graph, edgeIds); }

		Models::NT get_trajectory_length(std::size_t trajectoryNum) const { return m_trajectoryLengths[trajectoryNum]; }

		void precomputeComputationData(const std::vector<std::vector<std::size_t>>& trajectories) {
			GCpp::DS::computeEdgeIdToEdgeMap(graph, edgeIdToEdgeMap);
			std::size_t trajectoryNum = 0;
			for (const auto& traj : trajectories) {
				std::size_t currentEdge = 0;
				Models::NT totalLength = 0;
				for (auto e : traj) {
					if (edgeIdToTrajectoryNum.find(e) == edgeIdToTrajectoryNum.end()) {
						edgeIdToTrajectoryNum[e] = {};
					}
					edgeIdToTrajectoryNum[e].insert(trajectoryNum);
					if (edgeIds.find(e) == edgeIds.end())
						throw std::runtime_error("Trajectory off-road. Clean your data first!");

					if (edgeWeights.find(e) == edgeWeights.end())
						edgeWeights[e] = EdgeWeight{0};
					++edgeWeights[e].value;
					if (edgeIdToTrajectoryEdges.find(e) == edgeIdToTrajectoryEdges.end())
						edgeIdToTrajectoryEdges[e] = {};
					edgeIdToTrajectoryEdges[e][trajectoryNum] = TrajectoryEdge{trajectoryNum, currentEdge};
					totalLength += getWeightedEdgeLength(e, false);
					++currentEdge;
				}
				m_trajectoryLengths.push_back(totalLength);
				++trajectoryNum;
			}
		}

		std::size_t bundleCount(std::size_t edge_id) const {
			return edgeIdToBundleCount.find(edge_id) == edgeIdToBundleCount.end() ? 0 : edgeIdToBundleCount.at(edge_id);
		}

		Models::NT getWeightedEdgeLength(std::size_t edgeId, bool weight_by_bundles = false) const {
			auto weight = GCpp::DS::get_linear_edge_length(edgeIdToEdgeMap.at(edgeId), graph);
			if (weight_by_bundles) {
				return weight / (Models::NT{1} + Models::NT(bundleCount(edgeId)));
			}
			return weight;
		}

		void updatePriority(const std::vector<std::size_t>& bundleEdges, const MatchingMap& trajectories) {
			for (auto e : bundleEdges) {
				edgeWeights[e] = EdgeWeight{edgeWeights[e].value - trajectories.size()};
				if (edgeIdToBundleCount.find(e) == edgeIdToBundleCount.end())
					edgeIdToBundleCount[e] = 0;
				edgeIdToBundleCount[e] += 1;
			}
		}

		/**
		 * \brief Update edge ID to trajectories mapping by removing the given trajectories for the given edges.
		 * \param edges The dges
		 * \param trajectories The trajectories.
		 */
		void updateEdgeIdToTrajectoryMapping(const std::vector<std::size_t>& edges, const MatchingMap& trajectories) {
			for (auto edge : edges) {
				auto& trajectorySet = edgeIdToTrajectoryNum[edge];
				auto& edgeToTrajectoryEdges = edgeIdToTrajectoryEdges[edge];
				for (const auto& [trajId, trajMatching] : trajectories) {
					trajectorySet.erase(trajId);
					edgeToTrajectoryEdges.erase(trajId);
				}
			}
		}
	};

	static bool isConnected(typename EmbeddedGraph::vertex_descriptor v,
	                        typename EmbeddedGraph::edge_descriptor e,
	                        const EmbeddedGraph& g) {
		return boost::source(e, g) == v || boost::target(e, g) == v;
	}

	struct ProcessState {
		std::vector<PathEdge> currentPath;
		MatchingMap activeTrajectories;
		typename EmbeddedGraph::edge_descriptor edge;
		typename EmbeddedGraph::vertex_descriptor endVertex;
		std::size_t edgeId;
		Models::NT pathLength;
		std::size_t depth;
	};

	bool isBetter(const ProcessState& s0, const ProcessState& s1) const {
		return isBetter(std::make_pair(s0.pathLength, s0.activeTrajectories.size()),
		                std::make_pair(s1.pathLength, s1.activeTrajectories.size()));
	}
	bool isBetterBundle(const TrajectoryBundle& b0, const TrajectoryBundle& b1) const {
		return isBetter(std::make_pair(b0.bundleLength, b0.matchedTrajectoriesCount),
		                std::make_pair(b1.bundleLength, b1.matchedTrajectoriesCount));
	}
	bool isBetter(const std::pair<Models::NT, std::size_t>& lenAndSize0,
	              const std::pair<Models::NT, std::size_t>& lenAndSize1) const {
		auto weight = [this](Models::NT val, bool oneMinus = false) {
			return std::pow<Models::NT>(
			    val, oneMinus ? 1.0 - m_currentSettings.m_weight_exponent : m_currentSettings.m_weight_exponent);
		};
		switch (m_currentSettings.m_weighingScheme) {
			case WeighingScheme::Length: return lenAndSize0.first > lenAndSize1.first;
			case WeighingScheme::Support: return lenAndSize0.second > lenAndSize1.second;
			case WeighingScheme::LengthAndSupport:
				return lenAndSize0.second * lenAndSize0.first > lenAndSize1.second * lenAndSize1.first;
			case WeighingScheme::Custom:
				return weight(lenAndSize0.first) * weight(lenAndSize0.second, true) >
				       weight(lenAndSize1.first) * weight(lenAndSize1.second, true);
			default: throw std::runtime_error("Invalid scheme!");
		}
	}

	/**
	 * \brief Search for a bundle starting with the edge, given as an ID. Outputs bundle, its score and matched
	 * trajectories. Returns whether a bundle was found \param edgeID ID of the edge to start at \param data Computation
	 * data that is needed \param graph The graph \param outputBundle The output bundle \param outputScore The output
	 * score \param matchedTrajectories The output matched trajectories. \return Whether a bundle was found satisfying the
	 * criteria.
	 */
	bool searchForBundle(std::size_t edgeID,
	                     const ComputationData& data,
	                     const std::vector<std::vector<std::size_t>>& mmTrajectories,
	                     std::size_t min_bundle_count,
	                     TrajectoryBundle& bundle,
	                     MatchingMap& matchedTrajectories) const {
		if (data.edgeIdToTrajectoryNum.find(edgeID) == data.edgeIdToTrajectoryNum.end()) {
			if (m_currentSettings.m_verbose) {
				std::cout << "\tNo trajectories\n";
			}
			return false;
		}
		if (m_currentSettings.m_verbose) {
			std::cout << "###############################################\n";
			std::cout << "###############  New edge " << edgeID << "  ##########################\n";
			std::cout << "###############################################\n";
		}

		namespace it = GCpp::Helpers::Iterators;

		auto searchV = boost::target(data.edgeIdToEdgeMap.at(edgeID), data.graph);

		ProcessState initState{std::vector<PathEdge>{PathEdge{data.edgeIdToEdgeMap.at(edgeID), PathEdge::Explore::Both}},
		                       {},
		                       data.edgeIdToEdgeMap.at(edgeID),
		                       searchV,
		                       edgeID,
		                       data.getWeightedEdgeLength(edgeID, m_currentSettings.m_weigh_edges_by_bundles),
		                       0};

		// Initial matched trajectories.
		for (auto el : data.edgeIdToTrajectoryNum.at(edgeID)) {
			const auto edgeIndex = data.edgeIdToTrajectoryEdges.at(edgeID).at(el).edgeNum;
			initState.activeTrajectories[el] = TrajectoryMatch{
			    el, edgeIndex, edgeIndex, data.getWeightedEdgeLength(edgeID, m_currentSettings.m_weigh_edges_by_bundles)};
		}
		if (m_currentSettings.m_verbose) {
			for (auto d = 0; d < initState.depth; ++d) std::cout << "##|";
			std::cout << "Init bundle of length " << initState.pathLength
			          << ", matches: " << initState.activeTrajectories.size() << "/" << initState.activeTrajectories.size()
			          << '\n';
		}

		std::queue<ProcessState> elementsToProcess;
		std::optional<ProcessState> bestState;
		if (initState.activeTrajectories.size() >= min_bundle_count &&
		    initState.pathLength >= m_currentSettings.m_minimum_bundle_length) {
			bestState = initState;
		}
		elementsToProcess.push(initState);
		// Keep going while our score improves.
		while (!elementsToProcess.empty()) {
			auto state = elementsToProcess.front();
			elementsToProcess.pop();

			for (auto e : GCpp::Helpers::Iterators::range(boost::out_edges(state.endVertex, data.graph))) {
				auto eID = GCpp::DS::getEdgeId(data.graph, e);
				// Don't go to self
				if (eID == state.edgeId || eID == SEF::flipCanonical(state.edgeId))
					continue;
				// Skip edges without trajectories.
				if (data.edgeIdToTrajectoryNum.find(eID) == data.edgeIdToTrajectoryNum.end())
					continue;

				// Maybe weigh edges.
				const auto edgeLength = data.getWeightedEdgeLength(eID, m_currentSettings.m_weigh_edges_by_bundles);
				const auto newPathLength =
				    state.pathLength + data.getWeightedEdgeLength(eID, m_currentSettings.m_weigh_edges_by_bundles);

				ProcessState newState = state;
				newState.currentPath.push_back(PathEdge{e, PathEdge::Explore::Both});
				newState.activeTrajectories = {};
				for (const auto& kv : state.activeTrajectories) {
					if (!kv.second.matches(Direction::Forward, eID, mmTrajectories)) {
						continue;
					}
					auto copyVal = kv.second;
					copyVal.extend(Direction::Forward);
					copyVal.pathLength += edgeLength;
					newState.activeTrajectories[kv.first] = copyVal;
				}
				newState.edge = e;
				newState.edgeId = eID;
				++newState.depth;
				newState.endVertex = boost::source(e, data.graph) != state.endVertex ? boost::source(e, data.graph)
				                                                                     : boost::target(e, data.graph);
				newState.pathLength += data.getWeightedEdgeLength(eID, m_currentSettings.m_weigh_edges_by_bundles);

				if (newState.activeTrajectories.size() >= min_bundle_count) {
					if (m_currentSettings.m_verbose) {
						for (auto d = 0; d < newState.depth; ++d) std::cout << "##|";
						std::cout << "Checking bundle of length " << newState.pathLength
						          << ", matches: " << newState.activeTrajectories.size() << "/" << state.activeTrajectories.size()
						          << '\n';
					}
					if (newState.pathLength >= m_currentSettings.m_minimum_bundle_length) {
						if (!bestState.has_value() || this->isBetter(newState, bestState.value())) {
							bestState = newState;
						}
					}

					elementsToProcess.push(newState);
				} else {
					if (m_currentSettings.m_verbose) {
						for (auto d = 0; d < newState.depth; ++d) std::cout << "##|";
						std::cout << "Too litle matches: bundle of length " << newState.pathLength
						          << ", matches: " << newState.activeTrajectories.size() << "/" << state.activeTrajectories.size()
						          << '\n';
					}
				}
			}
		}
		if (!bestState.has_value()) {
			return false;
		}
		matchedTrajectories = bestState.value().activeTrajectories;
		bundle.matchedTrajectoriesCount = matchedTrajectories.size();
		bundle.edges.clear();
		bundle.bundleLength = bestState.value().pathLength;
		for (auto el : bestState.value().currentPath) {
			bundle.edges.push_back(GCpp::DS::getEdgeId(data.graph, el.edge));
		}
		// std::cout << "\tBundle with " << bundle.matchedTrajectoriesCount << " matched, length " << bundle.bundleLength <<
		// '\n';
		// Succesful path of required minimum length found.
		return bestState.value().pathLength >= m_currentSettings.m_minimum_bundle_length &&
		       matchedTrajectories.size() >= min_bundle_count;
	}

	/**
	 * \brief Check if b0 is a superbundle of b1, that is, b0 contains all edges of b1 (in order).
	 * \param b0 The first bundle
	 * \param b1 The second bundle
	 * \return b0 is a superbundle of b1 (or not)
	 */
	static bool isSuperBundle(const TrajectoryBundle& b0, const TrajectoryBundle& b1) {
		if (b0.edges.size() < b1.edges.size())
			return false;
		auto start = std::find(b0.edges.begin(), b0.edges.end(), b1.edges.front());
		if (start == b0.edges.end() || std::distance(start, b0.edges.end()) < b1.edges.size())
			return false;
		auto b1It = b1.edges.begin();
		for (; start != b0.edges.end() && b1It != b1.edges.end(); ++start, ++b1It) {
			if (*b1It != *start)
				return false;
		}
		return b1It == b1.edges.end();
	}

	// Settings object
	Settings m_currentSettings;

public:
	void setSettings(const Settings& settings) { m_currentSettings = settings; }
	const Settings& settings() const { return m_currentSettings; }
#pragma region Parameters
	void set_weigh_edges_by_bundles(bool val) { m_currentSettings.m_weigh_edges_by_bundles = val; }
	bool get_weigh_edges_by_bundles() const { return m_currentSettings.m_weigh_edges_by_bundles; }
	void set_superbundles_only(bool val) { m_currentSettings.m_superbundles_only = val; }
	void set_weighing_scheme(WeighingScheme scheme) { m_currentSettings.m_weighingScheme = scheme; }
	WeighingScheme weighing_scheme() const { return m_currentSettings.m_weighingScheme; }
	bool get_superbundles_only() const { return m_currentSettings.m_superbundles_only; }

	void set_verbose(bool value) { m_currentSettings.m_verbose = value; }
	bool get_verbose() const { return m_currentSettings.m_verbose; }
	void set_priotizie_by_length(bool val) { m_currentSettings.m_priotizie_by_length = val; }
	bool get_priotizie_by_length() const { return m_currentSettings.m_priotizie_by_length; }
	void set_minimum_bundle_length(const Models::NT& minimum_bundle_length) {
		m_currentSettings.m_minimum_bundle_length = minimum_bundle_length;
	}
	Models::NT get_minimum_bundle_length() const { return m_currentSettings.m_minimum_bundle_length; }
	void set_max_bundle_set_complexity(const std::size_t& max_bundle_complexity) {
		m_currentSettings.m_max_bundle_set_complexity = max_bundle_complexity;
	}
	std::size_t get_max_bundle_set_complexity() const { return m_currentSettings.m_max_bundle_set_complexity; }
	void set_weight_exponent(const Models::NT& exponent) { m_currentSettings.m_weight_exponent = exponent; }
	Models::NT weight_exponent() const { return m_currentSettings.m_weight_exponent; }
	void set_min_matched_trajectories(const std::size_t& min_match) {
		m_currentSettings.m_min_matched_trajectories = min_match;
	}
	std::size_t get_min_matched_trajectories() const { return m_currentSettings.m_min_matched_trajectories; }
#pragma endregion

	/**
	 * \brief
	 * \param graph Input simplified graph
	 * \param trajectories Given as list of edge IDs.
	 * \param output Vector of trajectory bundles, that contain score and bundle as list of edge IDs
	 */
	void operator()(const EmbeddedGraph& graph,
	                const std::vector<std::vector<std::size_t>>& trajectories,
	                std::vector<TrajectoryBundle>& output) const {
		// Compute data
		ComputationData data(graph);
		data.precomputeComputationData(trajectories);
		std::cout << "[BruteforceBundling] Starting bundling\n";

		print_max_trajectory(data, trajectories);

		if (!graph.m_use_canonical_edges_ids)
			throw std::runtime_error("Not using canonical edge IDS!");

		verify_connectivity(data, trajectories);

		// Allow getting smaller, increasing the class
		std::size_t current_bundle_count = m_currentSettings.m_min_matched_trajectories;
		std::size_t currentClass = 0;

		while (output.size() < m_currentSettings.m_max_bundle_set_complexity) {
			TrajectoryBundle bundle, bestBundle;
			MatchingMap matchedTrajectories, bestMatchedTrajectories;
			const std::size_t edgeNum = boost::num_edges(data.graph);
			std::size_t currEdge = 0;
			std::size_t perc = 5;
			const std::size_t percInc = 5;
			// Go over all directed edges (2 per undirected)
			for (auto e : GCpp::Helpers::Iterators::range(boost::edges(data.graph))) {
				bundle = {};
				matchedTrajectories = {};
				if (searchForBundle(GCpp::DS::getEdgeId(data.graph, e),
				                    data,
				                    trajectories,
				                    current_bundle_count,
				                    bundle,
				                    matchedTrajectories)) {
					if (isBetterBundle(bundle, bestBundle)) {
						bestBundle = bundle;
						bestMatchedTrajectories = matchedTrajectories;
					}
				}
				++currEdge;
				if (currEdge * 100 > edgeNum * perc) {
					std::cout << "[BruteforceBundling] At " << output.size() << ", " << perc << "%\n";
					perc += percInc;
				}
			}
			if (bestBundle.bundleLength > m_currentSettings.m_minimum_bundle_length) {
				bestBundle.supportClass = currentClass;
				data.updateEdgeIdToTrajectoryMapping(bestBundle.edges, bestMatchedTrajectories);
				data.updatePriority(bestBundle.edges, bestMatchedTrajectories);
				// TODO: cross class superbundles?
				if (m_currentSettings.m_superbundles_only) {
					bool shouldBeAdded = true;
					// Check if current is a subbundle or superbundle of someone else.
					for (std::size_t i = 0; i < output.size(); ++i) {
						if (isSuperBundle(output[i], bestBundle)) {
							shouldBeAdded = false;
							break;
						} else if (isSuperBundle(bestBundle, output[i])) {
							// Replace in output
							output[i] = bestBundle;
							shouldBeAdded = false;
							break;
						}
						// TODO: track where they go?
					}
					if (shouldBeAdded) {
						output.push_back(bestBundle);
					}
				} else {
					output.push_back(bestBundle);
				}
			} else {
				++currentClass;
				current_bundle_count /= 2;
				if (current_bundle_count == 1)
					break;
			}
		}
		std::cout << "[BruteforceBundling] Found " << output.size() << "/" << m_currentSettings.m_max_bundle_set_complexity
		          << " bundles \n";
		for (const auto& bundle : output) {
			std::cout << "[BruteforceBundling] Bundle with length " << bundle.bundleLength << " and support "
			          << bundle.matchedTrajectoriesCount << '\n';
		}
	}

private:
	void print_max_trajectory(const ComputationData& data,
	                          const std::vector<std::vector<std::size_t>>& trajectories) const {
		Models::NT maxLength = 0;
		for (const auto& traj : trajectories) {
			Models::NT localLength = 0;
			for (auto e : traj) {
				localLength += GCpp::DS::get_linear_edge_length(data.edgeIdToEdgeMap.at(e), data.graph);
			}
			maxLength = std::max(maxLength, localLength);
		}
		std::cout << "[BruteforceBundling] Max length trajectory" << maxLength << "\n";
	}

	void verify_connectivity(const ComputationData& data,
	                         const std::vector<std::vector<std::size_t>>& trajectories) const {  // Verify connectivity
		for (const auto& t : trajectories) {
			for (std::size_t i = 1; i < t.size(); ++i) {
				auto e0 = data.edgeIdToEdgeMap.at(t[i - 1]);
				auto e1 = data.edgeIdToEdgeMap.at(t[i]);
				if (boost::source(e1, data.graph) != boost::target(e0, data.graph)) {
					std::cout << GCpp::DS::getEdgeId(data.graph, e0) << "|" << GCpp::DS::getEdgeId(data.graph, e1) << '\n';
					std::cout << "Inv:" << SEF::flipCanonical(GCpp::DS::getEdgeId(data.graph, e0)) << "|"
					          << SEF::flipCanonical(GCpp::DS::getEdgeId(data.graph, e1));
					throw std::runtime_error("Disconnected path");
				}
			}
		}
	}
};
}  // namespace SchematLib::TrajectorySchematization

namespace SchematLib::IO {
template <typename EmbeddedGraph>
void read_settings(
    boost::property_tree::ptree& data,
    typename TrajectorySchematization::BidirectionalBruteforceTrajectoryBundler<EmbeddedGraph>::Settings& settings) {
#define READ_MEM(member) settings.m_member = data.get<decltype(settings.member)>(#member)
	READ_MEM(min_matched_trajectories);
	READ_MEM(minimum_bundle_length);
	READ_MEM(verbose);
}
}  // namespace SchematLib::IO
#endif
