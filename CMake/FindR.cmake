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