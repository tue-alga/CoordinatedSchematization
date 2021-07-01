#ifndef MEESGENERALIZATION_RPRPG_H    
#define MEESGENERALIZATION_RPRPG_H
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <boost/property_tree/ptree.hpp>

#include <boost/property_tree/json_parser.hpp>


#include "GreedyTrajectoryPicker.h"
#include "MostCoveredHeuristic.h"
#include "SecondaryAlgs/ComponentCounter.h"

namespace SchematLib::MeesGeneralization{

    //TODO move
    

    template<typename NT>
    struct RoutePreservingGeneralizationSettings
    {
        NT bbox[4]{ 0 };
        NT min_length;
        NT min_unique_length;
        NT max_point_error;
        NT max_frechet_error;
        NT budget_percentage;
        int heuristic=-1;
        bool count_bad_leafs = false;
        bool debug_messages = false;
    };

    template<typename NT>
    struct RoutePreservingRoadGeneralizationStats
    {
        // Total length (units of spatial reference) of graph.
        NT graph_length = 0;
        // Budget for retaining
        NT budget = 0;
        // 
        NT subgraphlength = 0;
        NT avgplyedge = 0;
        NT avgplymeter = 0;
        std::size_t num_picked = 0;
        NT length_picked = 0;
        std::size_t num_components = 0;
        std::size_t num_leaf_nodes = 0;
        std::optional<std::size_t> num_bad_leafs;
        NT biggest_component_edge_length = 0;
        std::size_t debug_bigcomp_idx = 0;
    };

    template<typename EmbeddedGraph, typename NT, typename EdgeTrajectory>
    class RoutePreservingGen
    {
        using Edge = typename EmbeddedGraph::Edge;

        // Define type preconditions
        static_assert(std::is_convertible_v<Edge, decltype(*std::declval<const EdgeTrajectory&>().begin())>, "Trajectory should be iterable with result Edge of the graph");

        bool m_debugMessages = false;
        bool m_countBadLeafs = false;
        bool m_rescore = false;
        bool m_removeWorst = false;
        NT m_budgetPercentage = 0.05;
    public:
        enum class Heuristic
        {
            MostCovered,
            GreedyTrajectoryPicker
        };
    private:
        Heuristic m_heuristic = Heuristic::MostCovered;

        void writeResult(const std::string& resultsFile, const RoutePreservingRoadGeneralizationStats<NT>& stats, std::vector<std::size_t>& pickedEdges)
        {
            boost::property_tree::ptree output;

            constexpr const char* GRAPH_LENGTH = "graph_length";

            output.put(GRAPH_LENGTH, stats.graph_length);
            output.put("budget", stats.budget);
            output.put("subgraph_length", stats.subgraphlength);
            output.put("edge_avg_subgraph_ply", stats.avgplyedge);
            output.put("length_avg_subgraph_ply ", stats.avgplymeter);
            output.put("number_of_trajectories_picked", stats.num_picked);
            output.put("length_of_trajectories_picked", stats.length_picked);
            output.put("number_of_connected_components", stats.num_components);
            output.put("number_of_leaf_nodes", stats.num_leaf_nodes);
            if (stats.num_bad_leafs.has_value())
                output.put("number_of_bad_leaf_nodes", stats.num_bad_leafs.value());
            output.put("length_of_biggest_component ", stats.biggest_component_edge_length);
            output.put("debug_bigcomp_idx",stats.debug_bigcomp_idx);

            // 

            std::ofstream outputStream(resultsFile);
            if (!outputStream.is_open()) throw std::runtime_error("Could not open output stream");
            boost::property_tree::write_json(outputStream, output);

        }

        /**
         * \brief 
         * \return Whether the current trajectory should be kept 
         */
        bool countBadLeafsForTrajectory(const std::vector<std::size_t>& fulleulerian, const EmbeddedGraph& graph, std::vector<size_t>& endpointverts)
        {
            auto fue_it = fulleulerian.begin();
            auto e1 = *std::next(graph.edges().begin(), *fue_it);
            size_t e2_index = *fue_it;
            while (e2_index == *fue_it)
                fue_it++;
            auto e2 = *std::next(graph.edges().begin(), *fue_it);
            auto e2v1 = graph.source(e2);
            auto e2v2 = graph.target(e2);
            if (graph.source(e1) == e2v1 || graph.source(e1) == e2v2)
                endpointverts.push_back(graph.target(e1));
            else if (graph.target(e1) == e2v1 || graph.target(e1) == e2v2)
                endpointverts.push_back(graph.source(e1));
            else {
                if (m_debugMessages) std::cout << "[RPRPG] Could not determine endpoints of trajectory " << "\n";// << the_name << "." << endl;
                return false;
            }
            auto fue_back_it = fulleulerian.end();
            fue_back_it--;
            auto back_e1 = *(std::next(graph.edges().begin(), *fue_back_it));
            size_t back_e2_index = *fue_back_it;
            while (back_e2_index == *fue_back_it)
                fue_back_it--;
            auto back_e2 = *std::next(graph.edges().begin(), *fue_back_it);
            auto back_e2v1 = graph.source(back_e2);
            auto back_e2v2 = graph.target(back_e2);

            if (graph.source(back_e1) == back_e2v1 || graph.source(back_e1) == back_e2v2)
                endpointverts.push_back(graph.target(back_e1));
            else if (graph.target(back_e1) == back_e2v1 || graph.target(back_e1) == back_e2v2)
                endpointverts.push_back(graph.source(back_e1));
            else {
                endpointverts.pop_back();
                if (m_debugMessages) std::cout << "[RPRPG] Could not determine endpoints of trajectory \n";// << the_name << "." << endl;
                return false;
            }

            //Find the endpoints created by a trajectory doubling back
            if (fulleulerian.size() > 2)
            {
                auto fulleuler_it = fulleulerian.begin();
                size_t prev_pref_ref = *fulleuler_it;
                ++fulleuler_it;
                size_t prev_ref = *fulleuler_it;
                ++fulleuler_it;
                for (; fulleuler_it != fulleulerian.end(); ++fulleuler_it)
                {
                    if (prev_ref == *fulleuler_it && prev_pref_ref != prev_ref)
                    {
                        auto prev_edge = graph.edge(prev_ref);
                        auto prev_prev_edge = graph.edge(prev_pref_ref);
                        size_t prev_edge_v1 = graph.source(prev_edge);
                        size_t prev_edge_v2 = graph.target(prev_edge);
                        if (prev_edge_v1 == graph.source(prev_prev_edge) || prev_edge_v1 == graph.target(prev_prev_edge))
                            endpointverts.push_back(prev_edge_v2);
                        else if (prev_edge_v2 == graph.source(prev_prev_edge) || prev_edge_v2 == graph.target(prev_prev_edge))
                            endpointverts.push_back(prev_edge_v1);
                        else
                        {
                            throw std::runtime_error("Impossible trajectory");
                        }
                    }
                    prev_pref_ref = prev_ref;
                    prev_ref = *fulleuler_it;
                }
            }
            return true;
        }
    public:
        void setRemoveWorst(const bool& removeWorst)
        {
            m_removeWorst = removeWorst;
        }
        bool removeWorst() const
        {
            return m_removeWorst;
        }

        void setHeuristic(const Heuristic& heuristic)
        {
            m_heuristic = heuristic;
        }
        Heuristic heuristic() const
        {
            return m_heuristic;
        }

        NT  budgetPercentage() const { return m_budgetPercentage; }

        void setBudgetPercentage(NT budgetPercentage)
        {
            m_budgetPercentage = budgetPercentage;
        }
        void setRescore(const bool& rescore)
        {
            m_rescore = rescore;
        }
        bool rescore() const
        {
            return m_rescore;
        }
        void setDebugMessage(bool value)
        {
            m_debugMessages = value;
        }
        bool debugMessages() const
        {
            return m_debugMessages;
        }
        void setCountBadLeafs(const bool& countBadLeafs)
        {
            m_countBadLeafs = countBadLeafs;
        }
        bool countBadLeafs() const
        {
            return m_countBadLeafs;
        }



        //void operator()(const EmbeddedGraph& graph, const std::vector<EdgeTrajectory>& trajectories, const RoutePreservingGeneralizationSettings<NT>& settings, const std::string& resultFile)
        
        /**
         * \brief 
         * \param graph 
         * \param trajectories 
         * \param stats 
         * \param pickedEdges List of edge indices of the selected edges.
         */
        void operator()(const EmbeddedGraph& graph, const std::vector<EdgeTrajectory>& trajectories, RoutePreservingRoadGeneralizationStats<NT>& stats, std::vector<std::size_t>& pickedEdges)
        {
            using namespace std;

            const auto numberOfEdges = graph.num_edges();

            std::vector<std::vector<size_t>> euleriantrajectories;
            std::vector<size_t> covered_counts(numberOfEdges, 0);

            std::vector<size_t> endpointverts;

            if (m_debugMessages) std::cout << "[RPRPG] Loading trajectories." << std::endl;
            std::size_t cnt = 0;
            std::size_t printEvery = 50;
            std::size_t totalCover = 0; //For DEBUG!
            for (const auto& trajectory: trajectories)
            //for (auto& p : std::filesystem::directory_iterator(infolder))
            {
                ++cnt;
                if (m_debugMessages && cnt % printEvery) std::cout << "[RPRPG] At " << static_cast<double>( cnt )*100. / static_cast<double>(trajectories.size()) << std::endl;
                // READ TRAJECTORY AND FILTER ON SOME TERMS

                /*std::ifstream eulerfile(p.path());
                std::string line;
                std::getline(eulerfile, line);
                std::istringstream iss(line);
                string nametag, the_name, lengthtag, pointerrortag, frechettag, uniquelengthtag;
                NT the_length, point_error, frechet_error, unique_length;
                if (!(iss >> nametag >> the_name >> lengthtag >> the_length >> pointerrortag >> point_error >> frechettag >> frechet_error)) { cout << "Error parsing file" << endl; continue; }
                if ((iss >> uniquelengthtag >> unique_length)) { if (min_unique_length != 0 && unique_length < min_unique_length) continue; }
                if ((min_length != 0 && the_length < min_length) || (max_point_error != 0 && point_error > max_point_error) || (max_frechet_error != 0 && frechet_error > max_frechet_error))
                    continue;*/

                std::vector<size_t> eulerian; // Sorted unique edges
                std::set<size_t> set_eulerian; //(Sorted) edges, same as previous but the set.
                std::vector<size_t> fulleulerian; // Full path
                for(const auto& e: trajectory)
                {
                    auto index = graph.getProperty(edge_index_t{}, e);
                    set_eulerian.insert(index);
                    fulleulerian.push_back(index);
                }
                // Copy ordered unique edges to eulerian
                std::copy(set_eulerian.begin(), set_eulerian.end(), std::back_inserter(eulerian));
                
                // Covered length, not "actual" length.
                const NT trajlength = std::accumulate(set_eulerian.begin(), set_eulerian.end(), static_cast<NT>(0), [&graph](NT accum, const std::size_t& e)
                {
                    return accum + graph.getProperty(edge_length_t{}, graph.edge(e));
                });

                if (m_countBadLeafs)
                {
                    countBadLeafsForTrajectory(fulleulerian, graph, endpointverts);
                }
                for (const auto& e : eulerian)
                {
                    covered_counts[e] += 1;
                    ++totalCover;
                }
                euleriantrajectories.push_back(eulerian);
                //std::cout << "[RPRPG] Trajectory size: set:" << eulerian.size() << "=?" << set_eulerian.size() << ", full:" << fulleulerian.size() << std::endl;
            }
            if (m_debugMessages) std::cout << "[RPRPG] Processed trajectory files. Total edges covering sth: " << totalCover << std::endl;

            const NT graphLength = std::accumulate(graph.edges().begin(), graph.edges().end(), static_cast<NT>(0), [&graph](NT accum, const Edge& e)
            {
                return accum + graph.getProperty(edge_length_t{}, e);
            });
            /*NT graph_length = 0;
            for (auto edg : edges)
                graph_length += edg.length;*/
            if (m_debugMessages) std::cout << "[RPRPG] graph length: " << graphLength << std::endl;
            stats.graph_length = graphLength;

            //return 0;
            NT one_percent = stats.graph_length * 0.01;
            NT budget = one_percent * m_budgetPercentage;
            if (m_debugMessages) std::cout << "[RPRPG] Budget: " << budget << std::endl;
            
            std::vector<bool> pickeroutput;
            switch(m_heuristic)
            {
            case Heuristic::GreedyTrajectoryPicker:
                {
                GreedyTrajectoryPicker<EmbeddedGraph, NT, std::vector<size_t>> picker;
                picker.setRescore(m_rescore);
                picker.setRemoveWorst(m_removeWorst);
                picker.setDebugMessage(debugMessages());
                pickeroutput = picker(euleriantrajectories, graph, budget);
                }
                break;
            case Heuristic::MostCovered:
                {
                MostCoveredHeuristic<EmbeddedGraph,NT> mostcovered;
                //mostcovered.setDebugMessage(debugMessages());
                auto edgeSet = mostcovered(euleriantrajectories.begin(), euleriantrajectories.end(), graph, budget);
                pickeroutput.resize(numberOfEdges, false);
                if (m_debugMessages) std::cout << "[RPRPG] Mostcovered edgeset size:" << edgeSet.size() << std::endl;
                for (const auto& el : edgeSet) pickeroutput[el] = true;
                }
                break;
            default:break;
            }
            if (m_debugMessages) std::cout << "[RPRPG] Heuristic applied\n[RPRPG] Computing subgraph length and ply " << std::endl;
            //
            NT subgraphply = 0;
            NT subgraphplylength = 0;
            auto edges_it = graph.edges().begin();
            size_t index = 0;
            std::vector<size_t> pickedcoveredcounts(graph.num_edges(),0);
            size_t highest_coverage = 0;
            size_t picked_edge_count = 0;
            for (auto pit = pickeroutput.begin(); pit != pickeroutput.end(); ++pit)
            {
                if (*pit)
                {
                    picked_edge_count++;
                    stats.subgraphlength += graph.getProperty(edge_length_t{},*edges_it);
                    subgraphply += covered_counts[index];
                    subgraphplylength += covered_counts[index] * graph.getProperty(edge_length_t{}, *edges_it);
                }
                edges_it++;
                index++;
            }


            if (m_debugMessages) std::cout << "[RPRPG] Computing picked lengths" << std::endl;
            // Compute picking stats
            for (const auto& eulerTraj : euleriantrajectories)
            {
                bool picked = true;
                NT runninglength = 0;
                for (const auto& edgeIndex: eulerTraj)
                {
                    if (!pickeroutput[edgeIndex])
                    {
                        picked = false;
                        break;
                    }
                    runninglength += graph.getProperty(edge_length_t{}, graph.edge(edgeIndex));
                }
                // Not picked, so look at next
                if (!picked) continue;
                
                ++stats.num_picked;
                stats.length_picked += runninglength;
                for (const auto& edge: eulerTraj)
                {
                    pickedcoveredcounts[edge] += 1;
                    if (pickedcoveredcounts[edge] > highest_coverage)
                        highest_coverage = pickedcoveredcounts[edge];
                }
            }

            if (m_debugMessages) std::cout << "[RPRPG] Computing components" << std::endl;
            SecondaryAlgs::ComponentCounter<NT> componentcounter;
            auto count_results = componentcounter(graph, pickeroutput.begin(), pickeroutput.end());
            stats.num_components = count_results.num_connected_components;
            stats.num_leaf_nodes = count_results.num_leaf_nodes;

            size_t num_bad_leafs = 0;
            if (m_debugMessages) std::cout << "[RPRPG] \tNumber of leaf nodes: " << count_results.leafnodes.size() << std::endl;
            if (m_countBadLeafs)
            {
                for (const auto& leafNode: count_results.leafnodes)
                {
                    if (std::find(endpointverts.begin(), endpointverts.end(), leafNode) == endpointverts.end())
                        num_bad_leafs++;
                }
            }
            if (m_debugMessages) std::cout << "[RPRPG] Computing biggest component" << std::endl;
            std::vector<bool> biggest_component;
            size_t biggest_component_size = 0;
            size_t indx = 0;
            if (m_debugMessages) std::cout << "[RPRPG] \tComponent ids count: " << count_results.componentids.size() << std::endl;
            stats.debug_bigcomp_idx = count_results.largest_component_idx;
            if (m_debugMessages) std::cout << "[RPRPG] \tLargest component id: " << count_results.largest_component_idx << std::endl;
            for (auto count_it = count_results.componentids.begin(); count_it != count_results.componentids.end(); ++count_it)
            {
                if (*count_it == stats.debug_bigcomp_idx)
                {
                    biggest_component.push_back(true);
                    biggest_component_size++;
                    stats.biggest_component_edge_length += graph.getProperty(edge_length_t{}, graph.edge(indx));
                }
                else biggest_component.push_back(false);

                indx++;
            }

            stats.graph_length = graphLength;
            stats.budget = budget;
            stats.avgplyedge = subgraphply / static_cast<NT>(picked_edge_count);
            stats.avgplymeter = subgraphplylength / stats.subgraphlength;
            if (m_countBadLeafs)
                stats.num_bad_leafs = num_bad_leafs;

            for(auto i = 0; i < pickeroutput.size(); ++i)
            {
                if (pickeroutput[i]) pickedEdges.push_back(i);
            }
        }

        void operator()(const EmbeddedGraph& graph, const std::vector<EdgeTrajectory>& trajectories, const std::string& resultFile)
        {
            OutputStats stats;
            std::vector<std::size_t> edges;
            this->operator()(graph, trajectories, stats, edges);
            if (debugMessages()) std::cout << "[RPRPG] Writing result ot file " << resultFile << std::endl;
            writeResult(resultFile, stats, edges);
        }
    };
}
#endif