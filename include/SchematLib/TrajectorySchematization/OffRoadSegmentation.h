#ifndef SCHEMATLIB_TRAJECTORYSCHEMATIZATION_OFFROADSEGMENTATION_H
#define SCHEMATLIB_TRAJECTORYSCHEMATIZATION_OFFROADSEGMENTATION_H
#include <optional>
#include <set>
#include <type_traits>
#include <tuple>

namespace SchematLib::TrajectorySchematization
{
    struct EdgeBased{};
    struct IndexBased{};
    template<typename Graph, typename Approach>
    class OffRoadSegmentation;
    /**
     * \brief Retrieves segments of a given (mapmatched) trajectory that go off-road (i.e. over edges not in the network).
     * \tparam Graph The graph type of the network
     */
    template<typename Graph>
    class OffRoadSegmentation<Graph, EdgeBased>
    {
        std::set<typename Graph::edge_descriptor> m_edges;

        bool hasEdge(typename Graph::edge_descriptor edge) const
        {
            return m_edges.find(edge) != m_edges.end();
        }
    public:
        OffRoadSegmentation(const std::set<typename Graph::edge_descriptor>& edges): m_edges(edges){}

        /**
         * \brief Segments the trajectory, returning the off-road sections as pairs of iterators in the trajectory.
         * Note that the end iterator is exclusive!
         * \tparam Trajectory The trajectory type. Should be iterable and return values that are convertible to graph edges
         * \tparam OutputIterator The output iterator type. Should be assignable with pairs of iterators of the trajectory
         * \param trajectory The trajectory
         * \param output The output iterator
         */
        template<typename Trajectory, typename OutputIterator,
            typename = std::enable_if_t<std::is_assignable_v<decltype(*std::declval<OutputIterator&>()),
                std::pair<decltype(std::declval<const Trajectory&>().begin()), decltype(std::declval<const Trajectory&>().begin())>>
                >
            >
        void operator()(const Trajectory& trajectory, OutputIterator output)
        {
            using It = decltype(trajectory.begin());
            // Start of segment that goes off-road.
            std::optional<It> segmentStart;
            for(It it = trajectory.begin(); it != trajectory.end(); ++it)
            {
                if(segmentStart.has_value() && hasEdge(*it))
                {
                    *output = std::make_pair(segmentStart.value(), it);
                    // Reset  the start
                    segmentStart = {};
                }
                else if(!segmentStart.has_value() && !hasEdge(*it))
                {
                    segmentStart = it;
                }
            }
            if(segmentStart.has_value())
            {
                *output = std::make_pair(segmentStart.value(), trajectory.end());
            }
        }
    };
    /**
     * \brief Retrieves segments of a given (mapmatched) trajectory that go off-road (i.e. over edges not in the network).
     * \tparam Graph The graph type of the network
     */
    template<typename Graph>
    class OffRoadSegmentation<Graph, IndexBased>
    {
        std::set<std::size_t> m_edges;

        bool hasEdgeIndex(std::size_t edgeIndex) const
        {
            return m_edges.find(edgeIndex) != m_edges.end();
        }
    public:
        OffRoadSegmentation(const std::set<std::size_t>& edges) : m_edges(edges) {}

        /**
         * \brief Segments the trajectory, returning the off-road sections as pairs of iterators in the trajectory.
         * Note that the end iterator is exclusive!
         * \tparam Trajectory The trajectory type. Should be iterable and return values that are convertible to graph edges
         * \tparam OutputIterator The output iterator type. Should be assignable with pairs of iterators of the trajectory
         * \param trajectory The trajectory
         * \param output The output iterator
         */
        template<typename Trajectory, typename OutputIterator,
            typename = std::enable_if_t<std::is_assignable_v<decltype(*std::declval<OutputIterator&>()),
            std::pair<decltype(std::declval<const Trajectory&>().begin()), decltype(std::declval<const Trajectory&>().begin())>>
            >,
            typename = std::enable_if_t<std::is_convertible_v<decltype(*std::declval<const Trajectory&>().begin()), std::size_t>>
            >
            void operator()(const Trajectory& trajectory, OutputIterator output)
        {
            using It = decltype(trajectory.begin());
            // Start of segment that goes off-road.
            std::optional<It> segmentStart;
            for (It it = trajectory.begin(); it != trajectory.end(); ++it)
            {
                if (segmentStart.has_value() && hasEdgeIndex(*it))
                {
                    *output = std::make_pair(segmentStart.value(), it);
                    // Reset  the start
                    segmentStart = {};
                }
                else if (!segmentStart.has_value() && !hasEdgeIndex(*it))
                {
                    segmentStart = it;
                }
            }
            if (segmentStart.has_value())
            {
                *output = std::make_pair(segmentStart.value(), trajectory.end());
            }
        }
    };

}
#endif