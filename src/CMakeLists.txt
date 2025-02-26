set(POET_SOURCES 
  Base/Grid.cpp
  Base/SimParams.cpp
  Chemistry/ChemistryModule.cpp
  Chemistry/MasterFunctions.cpp 
  Chemistry/WorkerFunctions.cpp 
  Chemistry/SurrogateModels/DHT_Wrapper.cpp 
  Chemistry/SurrogateModels/HashFunctions.cpp 
  DataStructures/Field.cpp
  Transport/DiffusionModule.cpp
  Transport/AdvectionModule.cpp
)

# option(POET_USE_DHT_MPI "Use MPI for DHT" OFF)
#
# if (NOT POET_USE_DHT_MPI)
#   target_compile_definitions(poetlib PUBLIC POET_DHT_UCX)
#   target_link_libraries(poetlib PUBLIC DHT_UCX)
# else()
#   target_link_libraries(poetlib PUBLIC DHT_MPI)
# endif()
#
# target_compile_definitions(poetlib PUBLIC STRICT_R_HEADERS OMPI_SKIP_MPICXX)
#
# mark_as_advanced(PHREEQCRM_BUILD_MPI PHREEQCRM_DISABLE_OPENMP)
# set(PHREEQCRM_DISABLE_OPENMP ON CACHE BOOL "" FORCE)
#
# option(POET_DHT_DEBUG "Build with DHT debug info" OFF)
#
# if(POET_DHT_DEBUG)
#   target_compile_definitions(poetlib PRIVATE DHT_STATISTICS)
# endif()
#
# option(POET_PHT_ADDITIONAL_INFO "Enables additional information in the PHT" OFF)
#
# if (POET_PHT_ADDITIONAL_INFO)
#   target_compile_definitions(poetlib PRIVATE POET_PHT_ADD)
# endif()
#
# if (OpenMP_FOUND)
#   target_link_libraries(poetlib PUBLIC OpenMP::OpenMP_CXX)
# endif()

file(READ "${PROJECT_SOURCE_DIR}/R_lib/kin_r_library.R" R_KIN_LIB )

configure_file(poet.hpp.in poet.hpp @ONLY)

add_executable(poet_coarse poet.cpp ${POET_SOURCES})
target_link_libraries(poet_coarse PRIVATE 
  MPI::MPI_C 
  ${MATH_LIBRARY} 
  RRuntime 
  PhreeqcRM 
  tug
  OpenSSL::Crypto
  LMPI_DHT_OLD
)
target_include_directories(poet_coarse PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
install(TARGETS poet_coarse DESTINATION bin)

add_executable(poet_fine poet.cpp ${POET_SOURCES})
target_link_libraries(poet_fine PRIVATE 
  MPI::MPI_C 
  ${MATH_LIBRARY} 
  RRuntime 
  PhreeqcRM 
  tug
  OpenSSL::Crypto
  LFINE_LOCK
)
target_include_directories(poet_fine PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
install(TARGETS poet_fine DESTINATION bin)

add_executable(poet_no poet.cpp ${POET_SOURCES})
target_link_libraries(poet_no PRIVATE 
  MPI::MPI_C 
  ${MATH_LIBRARY} 
  RRuntime 
  PhreeqcRM 
  tug
  OpenSSL::Crypto
  LMPI_DHT
)
target_include_directories(poet_no PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
install(TARGETS poet_no DESTINATION bin)
