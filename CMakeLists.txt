cmake_minimum_required(VERSION 3.14)
include(ExternalProject)
project(POET
  LANGUAGES CXX C
  DESCRIPTION "A coupled reactive transport simulator")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


include("CMake/POET_Scripts.cmake")
include(FetchContent)

list(APPEND CMAKE_MODULE_PATH "${POET_SOURCE_DIR}/CMake")

get_poet_version()

find_package(MPI REQUIRED)
find_package(RRuntime REQUIRED)

set(USE_NAA ON)
FetchContent_Declare(naa-communication-prototype
  GIT_REPOSITORY gfz:naaice/naa-communication-prototype
  GIT_TAG        cmake-version)

if(USE_NAA)
  FetchContent_MakeAvailable(naa-communication-prototype)
  set(SOURCE_IP "10.3.10.41")
  set(TARGET_IP "10.3.10.42")
endif()


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
