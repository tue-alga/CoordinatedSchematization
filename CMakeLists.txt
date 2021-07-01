cmake_minimum_required(VERSION 3.16)

# The project name
set(PROJ_NAME SchematLib)

# Setup the project
project(${PROJ_NAME} LANGUAGES C CXX)

option(SCHEMATLIB_CREATE_FORMAT_TARGET "Create format target" OFF)

# Sources
file(GLOB_RECURSE SRC CONFIGURE_DEPENDS src/**/*.cpp src/*.cpp)

file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS include/**/*.h include/**/*.cpp)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SRC} ${HEADERS})

message(STATUS "[SchematLib] Toolchain: ${CMAKE_TOOLCHAIN_FILE}")

include(FetchContent)
include(cmake/InstallDependency.cmake)

# Set downloaded dependency install directory
set(_DEPENDENCY_DIRECTORY ${CMAKE_BINARY_DIR}/dependencies)
# Add dependencies
find_package(Boost 1.75 REQUIRED COMPONENTS serialization date_time timer graph)

# Add the dependency directory to our prefix path for find_package in config mode.
list(APPEND CMAKE_PREFIX_PATH ${_DEPENDENCY_DIRECTORY})

## Geometry Core++ library
FetchContent_Declare(
  gcpp-lib
  GIT_REPOSITORY https://gitlab.tue.nl/20184217/gcpp-lib.git
  GIT_TAG        783a869a2bd34d2f45fee910803d81a41e948a41 # Master 
)
INSTALL_DEPENDENCY(
    FETCHCONTENT_NAME gcpp-lib
    INSTALL_DIRECTORY ${_DEPENDENCY_DIRECTORY}
)
find_package(GeometryCore++Lib CONFIG REQUIRED)

# Movetk
# todo(bram): currently pulling custom fork, change this when movetk becomes more stable.
FetchContent_Declare(
  movetk
  GIT_REPOSITORY https://github.com/bacusters/movetk.git
  GIT_TAG        2627444624d5197be8e0630c2eef4d7b28b0abb3
)
INSTALL_DEPENDENCY(
    CLEAN
    FETCHCONTENT_NAME movetk
    INSTALL_DIRECTORY ${_DEPENDENCY_DIRECTORY}
)
find_package(movetk CONFIG REQUIRED)

# TODO: fix a version
find_package(CGAL CONFIG REQUIRED)
find_package(GMP REQUIRED)

#Add the library
add_library(${PROJ_NAME} ${SRC} ${HEADERS})

# Enable C++17
target_compile_features(${PROJ_NAME} PUBLIC cxx_std_17)

target_compile_definitions(${PROJ_NAME} PUBLIC _USE_MATH_DEFINES)

#Include appropriate directories
target_include_directories(${PROJ_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link libraries
target_link_libraries(${PROJ_NAME} PUBLIC GeometryCore++Lib::GeometryCore++Lib movetk CGAL::CGAL ${GMP_LIBRARY} Boost::boost)

#add_subdirectory(SchematLibCommon)

if(SCHEMATLIB_CREATE_FORMAT_TARGET)
find_program(CLANG_FORMAT_EXEC clang-format)
add_custom_target(${PROJ_NAME}_format
	COMMAND ${CLANG_FORMAT_EXEC} --style=file -i ${SRC} ${HEADERS}
)
endif()

# option(SCHEMATLIB_BUILD_UI "Build UI" OFF)
# if(SCHEMATLIB_BUILD_UI)
#     add_subdirectory(SchematLibUI)
# endif()

# option(SCHEMATLIB_BUILD_CLI "Build CLI" OFF)
# if(SCHEMATLIB_BUILD_CLI)
#     add_subdirectory(SchematLibCLI)
# endif()