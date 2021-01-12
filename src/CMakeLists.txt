# Version 3.9+ offers new MPI package variables
cmake_minimum_required(VERSION 3.9)

project(POET VERSION 0.1)

# Not needed until now
# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(MPI REQUIRED)

set(GCC_CXX_FLAGS "-D STRICT_R_HEADERS")
add_definitions(${GCC_CXX_FLAGS})

# prepare R environment (Rcpp + RInside)
find_program(R_EXE "R")

# search for R executable, R header file and library path
if(R_EXE)
  execute_process(COMMAND ${R_EXE} RHOME
    OUTPUT_VARIABLE R_ROOT_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  find_path(R_INCLUDE_DIR R.h
    HINTS ${R_ROOT_DIR}
    PATHS /usr/inlcude /usr/local/include /usr/share
    PATH_SUFFIXES include/R R/include
  )

  find_library(R_LIBRARY R
    HINTS ${R_ROOT_DIR}/lib
  )
else()
  message(FATAL_ERROR "No R runtime found!")
endif()

set(R_LIBRARIES ${R_LIBRARY})
set(R_INCLUDE_DIRS ${R_INCLUDE_DIR})


# find Rcpp include directory
execute_process(COMMAND echo "cat(find.package('Rcpp'))"
  COMMAND ${R_EXE} --vanilla --slave
  RESULT_VARIABLE RCPP_NOT_FOUND
  ERROR_QUIET
  OUTPUT_VARIABLE RCPP_PATH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(RCPP_NOT_FOUND)
  message(FATAL_ERROR "Rcpp not found!")
endif()

find_path(R_Rcpp_INCLUDE_DIR Rcpp.h
  HINTS ${RCPP_PATH}
  PATH_SUFFIXES include)

list(APPEND R_INCLUDE_DIRS ${R_Rcpp_INCLUDE_DIR})


# find RInside libraries and include path
execute_process(COMMAND echo "cat(find.package('RInside'))"
  COMMAND ${R_EXE} --vanilla --slave
  RESULT_VARIABLE RINSIDE_NOT_FOUND
  ERROR_QUIET
  OUTPUT_VARIABLE RINSIDE_PATH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(RInside_NOT_FOUND)
  message(FATAL_ERROR "RInside not found!")
endif()

find_library(R_RInside_LIBRARY libRInside.so
  HINTS ${RINSIDE_PATH}/lib)

list(APPEND R_LIBRARIES ${R_RInside_LIBRARY})

find_path(R_RInside_INCLUDE_DIR RInside.h
  HINTS ${RINSIDE_PATH}
  PATH_SUFFIXES include)

list(APPEND R_INCLUDE_DIRS ${R_RInside_INCLUDE_DIR})

#include found directories for the whole scope (will be changed with modularization)
include_directories(${MPI_CXX_INCLUDE_DIRS})

#define program libraries

add_library(POET_Libs STATIC util/RRuntime.cpp model/Grid.cpp model/ChemSim.cpp model/ChemMaster.cpp model/ChemWorker.cpp DHT/DHT_Wrapper.cpp DHT/DHT.cpp util/Parser.cpp model/TransportSim.cpp)
target_include_directories(POET_Libs PUBLIC ${R_INCLUDE_DIRS})
target_link_libraries(POET_Libs ${R_LIBRARIES} MPI::MPI_CXX crypto)

#add_library(DHT OBJECT DHT.cpp dht_wrapper.cpp)
#target_link_libraries(DHT crypto R_Wrapper)

#add_library(Worker OBJECT worker.cpp)
#target_link_libraries(Worker ${R_LIBRARIES} MPI::MPI_CXX R_Wrapper)

add_executable(poet kin.cpp)
target_link_libraries(poet PUBLIC MPI::MPI_CXX POET_Libs)
