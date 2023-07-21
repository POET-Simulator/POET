file(GLOB_RECURSE poet_lib_SRC
  CONFIGURE_DEPENDS
  "*.cpp" "*.c")

find_library(MATH_LIBRARY m)
find_library(CRYPTO_LIBRARY crypto)

add_library(poet_lib ${poet_lib_SRC})
target_include_directories(poet_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(poet_lib PUBLIC
  MPI::MPI_CXX ${MATH_LIBRARY} RRuntime PhreeqcRM tug)

target_compile_definitions(poet_lib PUBLIC STRICT_R_HEADERS OMPI_SKIP_MPICXX)

mark_as_advanced(PHREEQCRM_BUILD_MPI PHREEQCRM_DISABLE_OPENMP)
set(PHREEQCRM_DISABLE_OPENMP ON CACHE BOOL "" FORCE)

option(POET_DHT_DEBUG "Build with DHT debug info" OFF)

if(POET_DHT_DEBUG)
  target_compile_definitions(poet_lib PRIVATE DHT_STATISTICS)
endif()
