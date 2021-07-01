#ifndef SCHEMATLIB_MEESGENERALIZATION_DS_MULTIPROPERTYMAP_H
#define SCHEMATLIB_MEESGENERALIZATION_DS_MULTIPROPERTYMAP_H
#include <tuple>

namespace SchematLib::MeesGeneralization::DS
{
    namespace detail
    {
        template <class T, class Tuple>
        struct Index;

        template <class T, class... Types>
        struct Index<T, std::tuple<T, Types...>> {
            static const std::size_t value = 0;
        };

        template <class T, class U, class... Types>
        struct Index<T, std::tuple<U, Types...>> {
            static const std::size_t value = 1 + Index<T, std::tuple<Types...>>::value;
        };
    }

    template<typename Tags, typename Types>
    struct MultipropertyMap;

    template<typename...Tags, typename...Ts>
    struct MultipropertyMap<std::tuple<Tags...>, std::tuple<Ts...>>
    {
        std::tuple<std::unordered_map<std::size_t, Ts>...> m_maps;
    public:
        template<typename Tag>
        using TagDataType = std::tuple_element_t <
            detail::Index<Tag, std::tuple<Tags...>>::value,
            std::tuple<Ts...>
        >;

        template<typename Tag>
        const TagDataType<Tag>& get(Tag, std::size_t id) const
        {
            return std::get<detail::Index<Tag, std::tuple<Tags...>>::value>(m_maps).at(id);
        }
        template<typename Tag>
        TagDataType<Tag>& get(Tag, std::size_t id)
        {
            return std::get<detail::Index<Tag, std::tuple<Tags...>>::value>(m_maps).at(id);
        }
        template<typename Tag>
        void put(Tag, std::size_t id, const std::tuple_element_t<detail::Index<Tag, std::tuple<Tags...>>::value, std::tuple<Ts...>>& value)
        {
            std::get<detail::Index<Tag, std::tuple<Tags...>>::value>(m_maps).at(id) = value;
        }
    };
}
#endif