cmake_minimum_required(VERSION 3.20)

project(POET
  LANGUAGES CXX C
  DESCRIPTION "A coupled reactive transport simulator")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}'.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
    STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

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
add_subdirectory(ext/litephreeqc EXCLUDE_FROM_ALL)

option(POET_ENABLE_TESTING "Build test suite for POET" OFF)

if (POET_ENABLE_TESTING)
  add_subdirectory(test)
endif()

option(BUILD_DOC "Build documentation with doxygen" OFF)

add_subdirectory(docs)
