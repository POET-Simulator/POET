add_library(poetlib
  Base/Grid.cpp
  Base/SimParams.cpp
  Chemistry/ChemistryModule.cpp
  Chemistry/MasterFunctions.cpp 
  Chemistry/WorkerFunctions.cpp 
  Chemistry/SurrogateModels/DHT_Wrapper.cpp 
  Chemistry/SurrogateModels/HashFunctions.cpp 
  Chemistry/SurrogateModels/InterpolationModule.cpp 
  Chemistry/SurrogateModels/ProximityHashTable.cpp
  DataStructures/Field.cpp
  Transport/DiffusionModule.cpp
  Transport/AdvectionModule.cpp
)

target_link_libraries(poetlib PUBLIC
  MPI::MPI_C 
  ${MATH_LIBRARY} 
  RRuntime 
  PhreeqcRM 
  tug
  DHT_UCX
)

target_compile_definitions(poetlib PUBLIC STRICT_R_HEADERS OMPI_SKIP_MPICXX)

mark_as_advanced(PHREEQCRM_BUILD_MPI PHREEQCRM_DISABLE_OPENMP)
set(PHREEQCRM_DISABLE_OPENMP ON CACHE BOOL "" FORCE)

option(POET_DHT_DEBUG "Build with DHT debug info" OFF)

if(POET_DHT_DEBUG)
  target_compile_definitions(poetlib PRIVATE DHT_STATISTICS)
endif()

option(POET_PHT_ADDITIONAL_INFO "Enables additional information in the PHT" OFF)

if (POET_PHT_ADDITIONAL_INFO)
  target_compile_definitions(poetlib PRIVATE POET_PHT_ADD)
endif()

if (OpenMP_FOUND)
  target_link_libraries(poetlib PUBLIC OpenMP::OpenMP_CXX)
endif()

file(READ "${PROJECT_SOURCE_DIR}/R_lib/kin_r_library.R" R_KIN_LIB )

configure_file(poet.hpp.in poet.hpp @ONLY)

add_executable(poet poet.cpp)
target_link_libraries(poet PRIVATE poetlib MPI::MPI_C RRuntime)
target_include_directories(poet PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")

install(TARGETS poet DESTINATION bin)
