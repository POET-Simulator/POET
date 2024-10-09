cmake_minimum_required(VERSION 3.14)

project(POET
  LANGUAGES CXX C
  DESCRIPTION "A coupled reactive transport simulator")

set(CMAKE_CXX_FLAGS "-O3 -march=native -mtune=native")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include("CMake/POET_Scripts.cmake")
list(APPEND CMAKE_MODULE_PATH "${POET_SOURCE_DIR}/CMake")

get_poet_version()

find_package(MPI REQUIRED)

find_package(RRuntime REQUIRED)

add_subdirectory(src)

option(POET_PREPROCESS_BENCHS "Preprocess benchmarks" ON)

if (POET_PREPROCESS_BENCHS)
  add_subdirectory(bench)
endif()

# as tug will also pull in doctest as a dependency
set(TUG_ENABLE_TESTING OFF CACHE BOOL "" FORCE)

add_subdirectory(ext/tug EXCLUDE_FROM_ALL)
add_subdirectory(ext/iphreeqc EXCLUDE_FROM_ALL)

option(POET_ENABLE_TESTING "Build test suite for POET" OFF)

if (POET_ENABLE_TESTING)
  add_subdirectory(test)
endif()

option(BUILD_DOC "Build documentation with doxygen" OFF)

add_subdirectory(docs)
