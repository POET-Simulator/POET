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