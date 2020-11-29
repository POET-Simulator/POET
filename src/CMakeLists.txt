# this is the current cmake version distributed by debian stable (buster)
cmake_minimum_required(VERSION 3.13)

project(Poet VERSION 0.0)

add_executable(poet kin.cpp)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(MPI)

list(APPEND R_DEPS "Rcpp" "RInside")

# prepare R environment (Rcpp + RInside)
find_program(R_EXE "R")

if(R_EXE)
  execute_process(COMMAND ${R_EXE} RHOME
    OUTPUT_VARIABLE R_ROOT_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  find_path(R_INCLUDE_DIR R.h
    HINTS ${R_ROOT_DIR}
    PATHS /usr/inlcude /usr/local/include
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

#foreach(_r_deps ${R_DEPS})
#  execute_process(COMMAND echo "cat(find.package('${_r_deps}'))"
#    COMMAND ${R_EXE} --vanilla --slave
#    RESULT_VARIABLE _rd
#    ERROR_QUIET
#    OUTPUT_VARIABLE _r_deps_path
#    OUTPUT_STRIP_TRAILING_WHITESPACE
#  )
#
#  if(_rd)
#    message(FATAL_ERROR "${_r_deps} not found!")
#  endif()
#
#  find_library(R_${_r_deps}_LIBRARY
#    lib${_r_deps}.so
#    HINTS ${_r_deps_path}/lib
#  )
#
#  if(R_${_r_deps}_LIBRARY) 
#    list(APPEND R_LIBRARIES ${R_${_r_deps}_LIBRARY})
#  else()
#    message(STATUS "lol")
#    unset(${R_${_r_deps}_LIBRARY})
#  endif()
#  
#endforeach()


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

#execute_process(
#  COMMAND ${R_EXE} "--slave" "--vanilla" "-e" "stopifnot(require('Rcpp'));cat(Rcpp:::Rcpp.system.file('include'))"
#  OUTPUT_VARIABLE LIBRCPP_INCLUDE_DIRS
#)
#
#execute_process(
#  COMMAND ${R_EXE} "--slave" "--vanilla" "-e" "stopifnot(require('RInside'));cat(RInside:::RInsideCFlags())"
#  OUTPUT_VARIABLE LIBRINSIDE_INCLUDE_DIRS
#)
#
#message(STATUS ${LIBRCPP_INCLUDE_DIRS})
#message(STATUS ${LIBRINSIDE_INCLUDE_DIRS})

# may find a better way ...
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBRINSIDE_INCLUDE_DIRS}")

include_directories(${R_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})

add_library(Worker STATIC worker.cpp)

add_library(DHT STATIC DHT.cpp dht_wrapper.cpp)
target_link_libraries(DHT crypto)

add_library(R_UTILS STATIC r_utils.cpp)

target_link_libraries(poet PUBLIC Worker DHT R_UTILS ${MPI_C_LIBRARIES} ${R_LIBRARIES})
