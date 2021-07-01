/*! @file BronKerbosch.h
 *  @brief  The Bron-Kerbosch algorithm for finding maximal cliques
 *  @authors Mees van de Kerkhof (m.a.vandekerkhof@uu.nl)
 */

 /*We use the Bron-Kerbosch algorithm to find maximal cliques in the consistency graph we create when doing Physics-Outlier Detection*/

#ifndef TULIB_BRONKERBOSCH_H
#define TULIB_BRONKERBOSCH_H

#include "tulib/utils/Requirements.h"

using namespace std;
/*!
 * @brief a collection of algorithms in tulib
 */
namespace tulib_algorithms
{
	/*!
	* @tparam GeometryTraits
	* @tparam TrajectoryTraits
	* @tparam TestType
	*/
	template <class GeometryTraits, class TrajectoryTraits, class TestType>
	class BronKerbosch {


	private:
		typedef typename GeometryTraits::NT NT;
		typedef TestType TEST;
		typedef typename GeometryTraits::Tulib_Point Point;
		typedef typename TrajectoryTraits::ProbeTraits::ProbeColumns Columns;
		TEST test;

		NT threshold;

		/*!
		* @tparam InputIterator
		*/
		template <class InputIterator>
		struct GraphNode
		{
		public:
			InputIterator probe; //Pointer to the underlying probe. 
			 //Length of the largest consistent trajectory ending in this probe
			vector<GraphNode> neighborhood;

			/*!
			*@param probeptr
			*/
			GraphNode(InputIterator InProbe) {
				probe = InProbe;
			}
		};
	public:

		/*!
		*@param InThreshold
		*/
		BronKerbosch(NT InThreshold) { threshold = InThreshold; test = TEST(InThreshold); };

		/*!
 *
 * @tparam InputIterator
 * @tparam OutputIterator
 * @param first
 * @param beyond
 * @param result
 */
		template <class InputIterator, class OutputIterator,
			typename = tulib_core::requires_random_access_iterator <InputIterator>,
			typename = tulib_core::requires_output_iterator <OutputIterator>>
			void operator()(InputIterator first, InputIterator beyond, OutputIterator result)
		{
			vector<GraphNode> vertices;
			//We start by creating the graph
			for (auto it = first; it != beyond; it++)
				vertices.push_back(GraphNode(it));
			auto verts_it = next(vertices.begin(),1);
			for (auto it = next(first, 1); it != beyond; it++)
			{
				auto verts_sub_it = vertices.begin();
				for (auto sub_it = first; sub_it != it; sub_it++)
				{
					if (test(*sub_it, *it))
					{
						verts_it->neighborhood.push_back(*verts_sub_it);
						verts_sub_it->neighborhood.push_back(*verts_it);
					}
					verts_sub_it++;
				}
				verts_it++;
			}


		} //operator()

		void _bronkerbosch(vector<GraphNode> P, vector<GraphNode> X, vector<GraphNode> R, vector<vector<GraphNode>> &maximalcliques)
		{
			if (P.size() == 0 && X.size() == 0)
			{
				maximalcliques.push_back(R);
				return;
			}
			GraphNode pivot;
			if (P.size() != null)
			{
				pivot = P.front();
			}
			else pivot = X.front();
			std::vector<GraphNode> Pprime(P.size());
			auto it = set_difference(P.begin(), P.end(), pivot.neighborhood.begin(), pivot.neighborhood.end(), Pprime.begin());
			Pprime.resize(it - Pprime.begin());
			for (auto it = Pprime.begin(); it != Pprime.end(); i++)
			{
				vector<GraphNode> recursiveR = copy(R);
				recursiveR.insert
			}
		}
	}; //class BronKerbosch
}

#endif //TULIB_BRONKERBOSCH_H
