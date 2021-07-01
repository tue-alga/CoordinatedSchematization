#ifndef SCHEMATLIB_IO_SETTINGS_IO_H
#define SCHEMATLIB_IO_SETTINGS_IO_H
#include <boost/property_tree/json_parser.hpp>
namespace SchematLib::IO {
template <typename Settings>
void read_settings(boost::property_tree::ptree& data, Settings& output) {}

template <typename Settings>
void write_settings(boost::property_tree::ptree& data, Settings& output) {}
}  // namespace SchematLib::IO
#endif
