option(POET_USE_PRM_BACKEND "Use PhreeqcRM parallelization instead of poet's." OFF)
option(POET_DHT_DEBUG "Build with DHT debug info" OFF)

list(APPEND CHEM_MODEL_SRC "ChemistryModule.cpp" )

mark_as_advanced(PHREEQCRM_BUILD_MPI PHREEQCRM_DISABLE_OPENMP)
set(PHREEQCRM_DISABLE_OPENMP ON CACHE BOOL "" FORCE)

if(POET_USE_PRM_BACKEND)
  set(PHREEQCRM_BUILD_MPI ON CACHE BOOL "" FORCE)
  target_compile_definitions(poet_lib PRIVATE POET_USE_PRM)
else()
  set(PHREEQCRM_BUILD_MPI OFF CACHE BOOL "" FORCE)
list(APPEND CHEM_MODEL_SRC "WorkerFunctions.cpp" "MasterFunctions.cpp" "DHT_Wrapper.cpp" "HashFunctions.cpp" "DaosKeyValue.c" "DHT.c")
endif()

# ML: Fuer Nico: Includes + linking fuer Daos bereitstellen
find_library(DAOS_LIB libdaos.so
  PATH_SUFFIXES lib lib64
)
find_path(DAOS_INCLUDE daos.h)

add_library(DAOS INTERFACE IMPORTED)
set_target_properties(DAOS PROPERTIES
  INTERFACE_LINK_LIBRARIES "${DAOS_LIB}"
  INTERFACE_INCLUDE_DIRECTORIES "${DAOS_INCLUDE}"
)

add_library(ChemistryModule ${CHEM_MODEL_SRC})
target_include_directories(ChemistryModule PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(ChemistryModule
  PUBLIC MPI::MPI_CXX PhreeqcRM
  PUBLIC MPI::MPI_CXX PhreeqcRM DAOS
  PRIVATE ${MATH_LIBRARY}
  )

target_compile_definitions(ChemistryModule PUBLIC OMPI_SKIP_MPICXX)

if(POET_DHT_DEBUG)
  target_compile_definitions(ChemistryModule PRIVATE DHT_STATISTICS)
endif()

if(POET_USE_PRM_BACKEND)
  target_compile_definitions(ChemistryModule PUBLIC POET_USE_PRM)
endif()
