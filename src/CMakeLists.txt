file(GLOB_RECURSE poet_lib_SRC
  CONFIGURE_DEPENDS
  "*.cpp" "*.c")

configure_file(poet.h.in poet.h)

add_executable(poet ${poet_lib_SRC})
target_link_libraries(poet 
  MPI::MPI_CXX 
  ${MATH_LIBRARY} 
  RRuntime 
  PhreeqcRM 
  tug
)
target_include_directories(poet PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")

target_compile_definitions(poet PRIVATE STRICT_R_HEADERS OMPI_SKIP_MPICXX)

mark_as_advanced(PHREEQCRM_BUILD_MPI PHREEQCRM_DISABLE_OPENMP)
set(PHREEQCRM_DISABLE_OPENMP ON CACHE BOOL "" FORCE)

option(POET_DHT_DEBUG "Build with DHT debug info" OFF)

if(POET_DHT_DEBUG)
  target_compile_definitions(poet PRIVATE DHT_STATISTICS)
endif()

option(POET_PHT_ADDITIONAL_INFO "Enables additional information in the PHT" OFF)

if (POET_PHT_ADDITIONAL_INFO)
  target_compile_definitions(poet PRIVATE POET_PHT_ADD)
endif()

install(TARGETS poet DESTINATION bin)