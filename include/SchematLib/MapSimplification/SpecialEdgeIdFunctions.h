#ifndef SCHEMATLIB_MAPSIMPLIFICATION_SPECIALEDGDEIDFUNCTIONS_H
#define SCHEMATLIB_MAPSIMPLIFICATION_SPECIALEDGDEIDFUNCTIONS_H
#include <vector>
namespace SchematLib::MapSimplification {
    struct SpecialEdgeIdFunctions
    {
        static constexpr std::size_t NONCANONICAL_BIT = (std::size_t{ 1 } << (sizeof(std::size_t) * 8 - 1));
        static inline std::size_t getEdgeId(std::size_t id)
        {
            return id & (~(NONCANONICAL_BIT));
        }
        static inline bool isNonCanonical(std::size_t id)
        {
            return (id & NONCANONICAL_BIT ) > 0;
        }
        static inline bool isCanonical(std::size_t id)
        {
            return !isNonCanonical(id);
        }
        static inline std::size_t makeNonCanonical(std::size_t id)
        {
            return getEdgeId(id) | NONCANONICAL_BIT;
        }
        static inline std::size_t makeCanonical(std::size_t id)
        {
            return getEdgeId(id); //Same, since non-canonical bit is filtered out
        }
        static inline std::size_t flipCanonical(std::size_t id)
        {
            return isCanonical(id) ? makeNonCanonical(id) : makeCanonical(id);
        }
        static inline void reverseCanonicalComplement(std::vector<std::size_t>& ids)
        {
            std::reverse(ids.begin(), ids.end());
            for (auto& id : ids)
            {
                id = isCanonical(id) ? makeNonCanonical(id) : makeCanonical(id);
            }
        }
    };
}
#endif