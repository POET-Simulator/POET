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

# find_library(R_Rcpp_LIBRARY Rcpp.so
#   HINTS ${RCPP_PATH}/libs)

# list(APPEND R_LIBRARIES ${R_Rcpp_LIBRARY})

find_path(R_Rcpp_INCLUDE_DIR Rcpp.h
  HINTS ${RCPP_PATH}
  PATH_SUFFIXES include)

list(APPEND R_INCLUDE_DIRS ${R_Rcpp_INCLUDE_DIR})