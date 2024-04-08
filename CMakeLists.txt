cmake_minimum_required(VERSION 3.14)

project(POET
  LANGUAGES CXX C
  DESCRIPTION "A coupled reactive transport simulator")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include("CMake/POET_Scripts.cmake")
list(APPEND CMAKE_MODULE_PATH "${POET_SOURCE_DIR}/CMake")

get_poet_version()

# set(GCC_CXX_FLAGS "-D STRICT_R_HEADERS") add_definitions(${GCC_CXX_FLAGS})

find_package(MPI REQUIRED)

find_package(RRuntime REQUIRED)

# add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
# add_link_options(-fsanitize=address)

add_subdirectory(src)
add_subdirectory(bench)

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