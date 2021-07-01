#ifndef SCHEMATLIB_MEESGENERALIZATION_SECONDARYALGS_RANGEITERATOR_H
#define SCHEMATLIB_MEESGENERALIZATION_SECONDARYALGS_RANGEITERATOR_H
namespace SchematLib::MeesGeneralization::SecondaryAlgs
{
    template<typename Iterator>
    class Range
    {
        Iterator m_begin, m_end;
    public:
        Range(Iterator begin, Iterator end):m_begin(begin), m_end(end){}
        Range(std::pair<Iterator, Iterator> rangePair) :m_begin(rangePair.first), m_end(rangePair.second) {}

        Iterator begin() const
        {
            return m_begin;
        }
        Iterator end() const
        {
            return m_end;
        }
    };
    template<typename Iterator>
    Range<Iterator> range(Iterator begin, Iterator end)
    {
        return Range(begin, end);
    }
}

#endif