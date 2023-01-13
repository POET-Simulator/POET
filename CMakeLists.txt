# Version 3.9+ offers new MPI package variables
cmake_minimum_required(VERSION 3.9)

project(POET CXX C)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include("CMake/POET_Scripts.cmake")
list(APPEND CMAKE_MODULE_PATH "${POET_SOURCE_DIR}/CMake")

# set(GCC_CXX_FLAGS "-D STRICT_R_HEADERS") add_definitions(${GCC_CXX_FLAGS})

find_package(MPI REQUIRED)

find_package(RRuntime REQUIRED)

add_subdirectory(src)
add_subdirectory(R_lib)
add_subdirectory(data)
add_subdirectory(app)
add_subdirectory(bench/dolo_diffu_inner)

add_subdirectory(ext/tug EXCLUDE_FROM_ALL)
add_subdirectory(ext/phreeqcrm EXCLUDE_FROM_ALL)

option(BUILD_DOC "Build documentation with doxygen" OFF)

if(BUILD_DOC)
  add_subdirectory(docs)
endif(BUILD_DOC)
