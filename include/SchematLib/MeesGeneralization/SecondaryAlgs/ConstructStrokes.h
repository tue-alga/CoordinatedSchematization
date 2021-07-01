#pragma once
#include <vector>
#include <tulib/geom/GeometryInterface.h>
#include <tulib/geo/geo.h>
//#include "DummyGraph.h"

template <class NT>
class ConstructStrokes
{
public: 
	template<class VertexIterator, class EdgeIterator>
	void operator()(VertexIterator v_begin, VertexIterator v_beyond, EdgeIterator e_begin, EdgeIterator e_beyond, NT angle_diff_threshold, std::vector<int> &stroke_ids, std::vector<std::vector<size_t>> &strokes)
	{
		//Start by initializing all stroke ids to -1;
		for (auto e_it = e_begin; e_it != e_beyond; e_it++)
			stroke_ids.push_back(-1);
		//size_t stroke_count = 0;
		size_t vertex_iteration = 0;
		//Go over all vertices and divide their incident edges into strokes
		for (auto v_it = v_begin; v_it != v_beyond; v_it++)
		{

			vertex_iteration++;

			auto incident_edges = v_it->edges;
			//if (incident_edges.size() > 4)
			//	throw "degree > 4";
			std::vector<size_t> opposite_verts;
			for (auto ie_it = incident_edges.begin(); ie_it != incident_edges.end(); ie_it++)
			{
				auto edg = e_begin + *ie_it;
				if (get<0>(edg->vertices) == distance(v_begin, v_it))
					opposite_verts.push_back(get<1>(edg->vertices));
				else if (get<1>(edg->vertices) == distance(v_begin, v_it))
					opposite_verts.push_back(get<0>(edg->vertices));
				else throw "bad vertex";
			}
			std::vector<std::tuple<size_t, size_t, NT>> angles;
			//Now opposite_verts is indexed like incident_edges.
			//Compute the turning angles between each pair of edges
			//std::vector<NT> angle;
			for (auto op_it = opposite_verts.begin(); op_it != opposite_verts.end(); op_it++)
				for (auto op_sub_it = op_it + 1; op_sub_it != opposite_verts.end(); op_sub_it++)
				{
					auto vert_1 = v_begin + *op_it;
					auto vert_2 = v_begin + *op_sub_it;
					auto bearing_1 = bearing_exact(vert_1->lat, vert_1->lon, v_it->lat, v_it->lon);
					auto bearing_2 = bearing_exact(v_it->lat, v_it->lon, vert_2->lat, vert_2->lon);
					//Get the difference in angle between the bearings. Because the bearings are mod 360, we must check both ways of going from one angle to the other
					//Start by ensuring bearing_1 <= bearing_2
					if (bearing_2 < bearing_1)
					{
						auto tmp = bearing_1;
						bearing_1 = bearing_2;
						bearing_2 = tmp;
					}
					auto bearing_diff = bearing_2 - bearing_1;
					if (bearing_diff > 180)
						bearing_diff = 360 - bearing_diff;
					////Now, take the two versions of the distance
					//auto bearing_diff_inner = bearing_2 - bearing_1;
					//auto bearing_diff_outer = bearing_1 + (360 - bearing_2);
					//NT bearing_diff = min(bearing_diff_inner, bearing_diff_outer);
					if (bearing_diff <= angle_diff_threshold)
					{
						auto edg1 = *(incident_edges.begin() + std::distance(opposite_verts.begin(), op_it));
						auto edg2 = *(incident_edges.begin() + std::distance(opposite_verts.begin(), op_sub_it));
						angles.push_back(make_tuple(edg1, edg2, bearing_diff));
					}
				}
			std::sort(angles.begin(), angles.end(), [](std::tuple<size_t, size_t, NT> a, std::tuple<size_t, size_t, NT> b) {
				return get<2>(a) < get<2>(b); 
			});
			std::vector<size_t> strokified_edges;
			for (auto angle_it = angles.begin(); angle_it != angles.end(); angle_it++)
			{
				auto e1 = get<0>(*angle_it);
				auto e2 = get<1>(*angle_it);
				//Make sure e1 <= e2 so we only have 3 cases to consider later
				if (stroke_ids[e1] > stroke_ids[e2])
				{
					auto tmp = e1;
					e1 = e2;
					e2 = tmp;
				}
				//If one of the edges was already covered, skip the angle
				if (std::find(strokified_edges.begin(), strokified_edges.end(), e1) != strokified_edges.end() || std::find(strokified_edges.begin(), strokified_edges.end(), e2) != strokified_edges.end())
					continue;
				
				//See if a new stroke should be constructed
				if (stroke_ids[e1] == -1)
				{
					if (stroke_ids[e2] == -1)
					{
						//Create a new stroke
						stroke_ids[e1] = strokes.size();
						stroke_ids[e2] = strokes.size();
						std::vector<size_t> stroke;
						stroke.push_back(e1);
						stroke.push_back(e2);
						strokes.push_back(stroke);
						//stroke_count++;
					}
					else
					{
						//Add e1 to e2's stroke
						stroke_ids[e1] = stroke_ids[e2];
						strokes[stroke_ids[e2]].push_back(e1);
					}
				}
				else
				{
					//e1 != -1 && e2 != -1
					if (stroke_ids[e1] == stroke_ids[e2])
					{
						strokified_edges.push_back(e1);
						strokified_edges.push_back(e2);
						continue;
					}
					
					//Delete the later stroke and merge it into the earlier stroke

					//Save the index from where we need to start updating strokes
					size_t old_e2_key = stroke_ids[e2];

					//Add stroke 2 to stroke 1
					strokes[stroke_ids[e1]].insert(strokes[stroke_ids[e1]].end(), strokes[old_e2_key].begin(), strokes[old_e2_key].end());

					//Update the refs so all edges can find their new stroke
					for (auto s_it = strokes[old_e2_key].begin(); s_it != strokes[old_e2_key].end(); s_it++)
						stroke_ids[*s_it] = stroke_ids[e1];

					//Delete the old stroke

					//cout << "vertex " <<vertex_iteration << " old_e2_key : " << old_e2_key << " strokes_size "<< strokes.size() << " " << endl;
					strokes.erase(strokes.begin() + old_e2_key);

					//Because a stroke has been erased, all of the stroke ids that come later have been invalidated;
					for (auto s_it = (strokes.begin() + old_e2_key); s_it != strokes.end(); s_it++)
					{
						for (auto s_sub_it = s_it->begin(); s_sub_it != s_it->end(); s_sub_it++)
							stroke_ids[*s_sub_it] = std::distance(strokes.begin(), s_it);
					}
					//Get the last stroke
					auto debugstroke = strokes[strokes.size() - 1];
					if (vertex_iteration == 111157)
						for (auto s : debugstroke)
							cout << s << endl;
				}
				strokified_edges.push_back(e1);
				strokified_edges.push_back(e2);
			}
			//If any edges are not yet part of a stroke, they will be made into new strokes
			for (auto e_it = incident_edges.begin(); e_it != incident_edges.end(); e_it++)
			{
				if (std::find(strokified_edges.begin(), strokified_edges.end(), *e_it) == strokified_edges.end() && stroke_ids[*e_it] == -1)
				{
					//Create a new stroke
					stroke_ids[*e_it] = strokes.size();
					std::vector<size_t> stroke;
					stroke.push_back(*e_it);
					strokes.push_back(stroke);
					//stroke_count++;
				}
			}

			

			//if (true)//if (debug_vert)
			//{
			//	//Debug step: After every step, validate that every stroke id is valid
			//	for (auto id : stroke_ids)
			//		if (id != -1 && id >= strokes.size())
			//			cout << "Stroke id " << id << " is a bad stroke id!" << " vertex iteration = " << vertex_iteration << endl;
			//	//Verify that the stroke ids and strokes match up
			//	for (auto strok_it = strokes.begin(); strok_it != strokes.end(); strok_it++)
			//	{
			//		//cout << "distance " << distance(strokes.begin(), strok_it) << endl;
			//		for (auto id_it = strok_it->begin(); id_it != strok_it->end(); id_it++)
			//			if (stroke_ids[*id_it] != distance(strokes.begin(), strok_it))
			//				cout << "Bad stroke / id matchup at vertex iteration " << vertex_iteration << ": Stroke " << distance(strokes.begin(), strok_it) << " has edge " << *id_it << " that has stroke_id " << stroke_ids[*id_it] << endl;
			//	}
			//	//Verify that all stroke ref counts are correct;
			//	std::vector<size_t> strokeref_counts = std::vector<size_t>(stroke_ids.size(), 0);
			//	for (auto id : stroke_ids)
			//	{
			//		if (id != -1)
			//			strokeref_counts[id] = strokeref_counts[id] + 1;
			//	}
			//	for (int i = 0; i < strokes.size(); i++)
			//		if (strokes[i].size() != strokeref_counts[i])
			//			cout << "Bad strokeref counts at iteration " << vertex_iteration << endl;
			//}
		}
	}
	
	template<class VertexIterator, class EdgeIterator, class SuperVertexIterator, class SuperEdgeIterator>
	void supergraph_strokes(VertexIterator v_begin, EdgeIterator e_begin, SuperVertexIterator s_v_begin, SuperVertexIterator s_v_beyond, SuperEdgeIterator s_e_begin, SuperEdgeIterator s_e_beyond, NT angle_diff_threshold, std::vector<int> &stroke_ids, std::vector<std::vector<size_t>> &strokes)
	{
		//Start by initializing all stroke ids to -1;
		for (auto e_it = s_e_begin; e_it != s_e_beyond; e_it++)
			stroke_ids.push_back(-1);

		size_t vertex_iteration = 0;
		//Go over all supervertices and divide their incident edges into strokes
		for (auto v_it = s_v_begin; v_it != s_v_beyond; v_it++)
		{
			bool debug_vert = false;
			vertex_iteration++;
			//Get the base vertex
			auto base_vert_id = v_it->base_vertex;
			if (base_vert_id == 53608)
				debug_vert = true;
			auto base_vert = v_begin + base_vert_id;
			//auto incident_edges = base_vert.edges;

			auto incident_edges = v_it->incident_superedges;
			std::vector<size_t> opposite_verts;
			for (auto ie_it = incident_edges.begin(); ie_it != incident_edges.end(); ie_it++)
			{
				auto edg = s_e_begin + *ie_it;
				size_t relevant_edge_id;
				if (get<0>(edg->endpoint_super_vertices) == distance(s_v_begin, v_it))
					relevant_edge_id = edg->component_edges.front();
				else if (get<1>(edg->endpoint_super_vertices) == distance(s_v_begin, v_it))
					relevant_edge_id = edg->component_edges.back();
				else throw "bad supervertex";

				auto relevant_edge = *(e_begin + relevant_edge_id);

					if (get<0>(relevant_edge.vertices) == base_vert_id)
						opposite_verts.push_back(get<1>(relevant_edge.vertices));
					else if (get<1>(relevant_edge.vertices) == base_vert_id)
						opposite_verts.push_back(get<0>(relevant_edge.vertices));
					else throw "bad superedge";
			}
			if (debug_vert)
			{
				cout << "opposite verts: " << endl;
				for (auto opver : opposite_verts)
					cout << opver << endl;
			}
			std::vector<std::tuple<size_t, size_t, NT>> angles;
			//Now opposite_verts is indexed like incident_edges.
			//Compute the turning angles between each pair of edges
			//std::vector<NT> angle;
			for (auto op_it = opposite_verts.begin(); op_it != opposite_verts.end(); op_it++)
				for (auto op_sub_it = op_it + 1; op_sub_it != opposite_verts.end(); op_sub_it++)
				{
					//To get the accurate bearings, we have to find the relevant component edges out of the superedges

					auto vert_1 = v_begin + *op_it;
					auto vert_2 = v_begin + *op_sub_it;
					if (debug_vert)
					{
						cout << "getting angle for vertices " << *op_it << " and " << *op_sub_it << endl;
						cout << "vert1 lat lon : " << vert_1->lat << " " << vert_1->lon << " " << endl;
						cout << "vert2 lat lon : " << vert_2->lat << " " << vert_2->lon << " " << endl;
						cout << "base lat lon : " << base_vert->lat << " " << base_vert->lon << endl;
						cout << "Projected:" << endl;
						LocalCoordinateReference<NT> ref(52.1389, 4.27673);
						auto pt = ref.project(52.331, 4.91239);

						auto pt1 = ref.project(vert_1->lat, vert_1->lon);
						cout << "vert1: " << get<0>(pt1) << " " << (get<1>(pt) - get<1>(pt1)) << " (" << get<1>(pt1) << ")" << endl;
						auto pt2 = ref.project(vert_2->lat, vert_2->lon);
						cout << "vert2: " << get<0>(pt2) << " " << (get<1>(pt) - get<1>(pt2)) << " (" << get<1>(pt2) << ")" << endl;
						auto ptbase = ref.project(base_vert->lat, base_vert->lon);
						cout << "base: " << get<0>(ptbase) << " " << (get<1>(pt) - get<1>(ptbase)) << " (" << get<1>(ptbase) << ")" << endl;
					}
					auto bearing_1 = bearing_exact(vert_1->lat, vert_1->lon, base_vert->lat, base_vert->lon);
					auto bearing_2 = bearing_exact(base_vert->lat, base_vert->lon, vert_2->lat, vert_2->lon);
					if (debug_vert)
					{
						cout << "bearing_1 (vert1->base): " << bearing_1 << endl;
						cout << "bearing_2 (base->vert2): " << bearing_2 << endl;
					}
					//Get the difference in angle between the bearings. Because the bearings are mod 360, we must check both ways of going from one angle to the other
					//Start by ensuring bearing_1 <= bearing_2
					if (bearing_2 < bearing_1)
					{
						auto tmp = bearing_1;
						bearing_1 = bearing_2;
						bearing_2 = tmp;
					}
					auto bearing_diff = bearing_2 - bearing_1;
					if (bearing_diff > 180)
						bearing_diff = 360 - bearing_diff;
					if (debug_vert)
					{
						cout << "bearing diff: " << bearing_diff << endl;
					}
					////Now, take the two versions of the distance
					//auto bearing_diff_inner = bearing_2 - bearing_1;
					//auto bearing_diff_outer = bearing_1 + (360 - bearing_2);
					//NT bearing_diff = min(bearing_diff_inner, bearing_diff_outer);
					if (bearing_diff <= angle_diff_threshold)
					{
						auto edg1 = *(incident_edges.begin() + std::distance(opposite_verts.begin(), op_it));
						auto edg2 = *(incident_edges.begin() + std::distance(opposite_verts.begin(), op_sub_it));
						angles.push_back(make_tuple(edg1, edg2, bearing_diff));
					}
				}
			std::sort(angles.begin(), angles.end(), [](std::tuple<size_t, size_t, NT> a, std::tuple<size_t, size_t, NT> b) {
				return get<2>(a) < get<2>(b);
			});
			if (debug_vert)
			{
				cout << "sorted angles: " << endl;
				for (auto ang : angles)
					cout << get<0>(ang) << " " << get<1>(ang) << " " << get<2>(ang) << endl;
				cout << "pause for breakpoint" << endl;
			}
			std::vector<size_t> strokified_edges;
			for (auto angle_it = angles.begin(); angle_it != angles.end(); angle_it++)
			{
				auto e1 = get<0>(*angle_it);
				auto e2 = get<1>(*angle_it);
				//Make sure e1 <= e2 so we only have 3 cases to consider later
				if (stroke_ids[e1] > stroke_ids[e2])
				{
					auto tmp = e1;
					e1 = e2;
					e2 = tmp;
				}
				//If one of the edges was already covered, skip the angle
				if (std::find(strokified_edges.begin(), strokified_edges.end(), e1) != strokified_edges.end() || std::find(strokified_edges.begin(), strokified_edges.end(), e2) != strokified_edges.end())
					continue;

				//See if a new stroke should be constructed
				if (stroke_ids[e1] == -1)
				{
					if (stroke_ids[e2] == -1)
					{
						//Create a new stroke
						stroke_ids[e1] = strokes.size();
						stroke_ids[e2] = strokes.size();
						std::vector<size_t> stroke;
						stroke.push_back(e1);
						stroke.push_back(e2);
						strokes.push_back(stroke);
						//stroke_count++;
					}
					else
					{
						//Add e1 to e2's stroke
						stroke_ids[e1] = stroke_ids[e2];
						strokes[stroke_ids[e2]].push_back(e1);
					}
				}
				else
				{
					//e1 != -1 && e2 != -1
					if (stroke_ids[e1] == stroke_ids[e2])
					{
						strokified_edges.push_back(e1);
						strokified_edges.push_back(e2);
						continue;
					}

					//Delete the later stroke and merge it into the earlier stroke

					//Save the index from where we need to start updating strokes
					size_t old_e2_key = stroke_ids[e2];

					//Add stroke 2 to stroke 1
					strokes[stroke_ids[e1]].insert(strokes[stroke_ids[e1]].end(), strokes[old_e2_key].begin(), strokes[old_e2_key].end());

					//Update the refs so all edges can find their new stroke
					for (auto s_it = strokes[old_e2_key].begin(); s_it != strokes[old_e2_key].end(); s_it++)
						stroke_ids[*s_it] = stroke_ids[e1];

					//Delete the old stroke

					//cout << "vertex " <<vertex_iteration << " old_e2_key : " << old_e2_key << " strokes_size "<< strokes.size() << " " << endl;
					strokes.erase(strokes.begin() + old_e2_key);

					//Because a stroke has been erased, all of the stroke ids that come later have been invalidated;
					for (auto s_it = (strokes.begin() + old_e2_key); s_it != strokes.end(); s_it++)
					{
						for (auto s_sub_it = s_it->begin(); s_sub_it != s_it->end(); s_sub_it++)
							stroke_ids[*s_sub_it] = std::distance(strokes.begin(), s_it);
					}
				}
				strokified_edges.push_back(e1);
				strokified_edges.push_back(e2);
			}
			//If any edges are not yet part of a stroke, they will be made into new strokes
			for (auto e_it = incident_edges.begin(); e_it != incident_edges.end(); e_it++)
			{
				if (std::find(strokified_edges.begin(), strokified_edges.end(), *e_it) == strokified_edges.end() && stroke_ids[*e_it] == -1)
				{
					//Create a new stroke
					stroke_ids[*e_it] = strokes.size();
					std::vector<size_t> stroke;
					stroke.push_back(*e_it);
					strokes.push_back(stroke);
					//stroke_count++;
				}
			}



			//if (true)//if (debug_vert)
			//{
			//	//Debug step: After every step, validate that every stroke id is valid
			//	for (auto id : stroke_ids)
			//		if (id != -1 && id >= strokes.size())
			//			cout << "Stroke id " << id << " is a bad stroke id!" << " vertex iteration = " << vertex_iteration << endl;
			//	//Verify that the stroke ids and strokes match up
			//	for (auto strok_it = strokes.begin(); strok_it != strokes.end(); strok_it++)
			//	{
			//		//cout << "distance " << distance(strokes.begin(), strok_it) << endl;
			//		for (auto id_it = strok_it->begin(); id_it != strok_it->end(); id_it++)
			//			if (stroke_ids[*id_it] != distance(strokes.begin(), strok_it))
			//				cout << "Bad stroke / id matchup at vertex iteration " << vertex_iteration << ": Stroke " << distance(strokes.begin(), strok_it) << " has edge " << *id_it << " that has stroke_id " << stroke_ids[*id_it] << endl;
			//	}
			//	//Verify that all stroke ref counts are correct;
			//	std::vector<size_t> strokeref_counts = std::vector<size_t>(stroke_ids.size(), 0);
			//	for (auto id : stroke_ids)
			//	{
			//		if (id != -1)
			//			strokeref_counts[id] = strokeref_counts[id] + 1;
			//	}
			//	for (int i = 0; i < strokes.size(); i++)
			//		if (strokes[i].size() != strokeref_counts[i])
			//			cout << "Bad strokeref counts at iteration " << vertex_iteration << endl;
			//}
		}
	}
};