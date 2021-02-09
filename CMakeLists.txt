# Version 3.9+ offers new MPI package variables
cmake_minimum_required(VERSION 3.9)

project(POET VERSION 0.1 LANGUAGES CXX C)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

list(APPEND CMAKE_MODULE_PATH "${POET_SOURCE_DIR}/CMake")

# set(GCC_CXX_FLAGS "-D STRICT_R_HEADERS")
# add_definitions(${GCC_CXX_FLAGS})

find_package(MPI REQUIRED)

find_package(R REQUIRED)
find_package(Rcpp REQUIRED)
find_package(RInside REQUIRED)

add_subdirectory(src)
add_subdirectory(R_lib)
add_subdirectory(data)
