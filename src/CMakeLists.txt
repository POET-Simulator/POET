add_library(POETLib)

target_sources(POETLib 
  PRIVATE
  Init/InitialList.cpp 
  Init/GridInit.cpp 
  Init/DiffusionInit.cpp
  Init/ChemistryInit.cpp
  DataStructures/Field.cpp
  Transport/DiffusionModule.cpp
  Chemistry/ChemistryModule.cpp
  Chemistry/MasterFunctions.cpp
  Chemistry/WorkerFunctions.cpp
  Chemistry/SurrogateModels/DHT_Wrapper.cpp
  Chemistry/SurrogateModels/DHT.c
  Chemistry/SurrogateModels/HashFunctions.cpp
  Chemistry/SurrogateModels/InterpolationModule.cpp
  Chemistry/SurrogateModels/ProximityHashTable.cpp
) 

target_include_directories(
 POETLib PUBLIC
 "${CMAKE_CURRENT_SOURCE_DIR}"
)

option(USE_AI_SURROGATE "Compiles the AI surrogate functions with Python.h integration" OFF)
if(USE_AI_SURROGATE)
  add_definitions(-DUSE_AI_SURROGATE)
  message("Building with AI surrogate functions")
  message(" - Needs Python with Numpy and Keras and Threads libraries")
  # make sure to use the python installation from the conda environment
  if(DEFINED ENV{CONDA_PREFIX})
    set(Python3_EXECUTABLE "$ENV{CONDA_PREFIX}/bin/python3")
  endif()
  find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
  
  find_package(Threads REQUIRED)
  
  set_source_files_properties(
    Chemistry/SurrogateModels/AI_functions.cpp
    PROPERTIES COMPILE_FLAGS "-O3 -march=native -mtune=native"
  )

  target_sources(POETLib
   PRIVATE
   Chemistry/SurrogateModels/AI_functions.cpp
  )
  
  target_include_directories(POETLib PUBLIC
   "${Python3_INCLUDE_DIRS}"
   "${Python3_NumPy_INCLUDE_DIRS}"
  )
  
  target_link_libraries(POETLib
   PUBLIC "${Python3_LIBRARIES}"
   PUBLIC Threads::Threads pthread
  )
endif() # USE_AI_SURROGATE


target_link_libraries(
 POETLib
 PUBLIC RRuntime
 PUBLIC IPhreeqcPOET
 PUBLIC tug
 PUBLIC MPI::MPI_C
)

# Define Python API version
target_compile_definitions(POETLib PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)


include(FetchContent)
FetchContent_Declare(
  cli11
  QUIET 
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG v2.4.2
)

FetchContent_MakeAvailable(cli11)

# add_library(poetlib
#   Base/Grid.cpp
#   Base/SimParams.cpp
#   Chemistry/ChemistryModule.cpp
#   Chemistry/MasterFunctions.cpp 
#   Chemistry/WorkerFunctions.cpp 
#   Chemistry/SurrogateModels/DHT_Wrapper.cpp 
#   Chemistry/SurrogateModels/DHT.c 
#   Chemistry/SurrogateModels/HashFunctions.cpp 
#   Chemistry/SurrogateModels/InterpolationModule.cpp 
#   Chemistry/SurrogateModels/ProximityHashTable.cpp
#   Transport/DiffusionModule.cpp
# )

# target_link_libraries(poetlib PUBLIC
#   DataStructures
#   Init
#   MPI::MPI_C 
#   ${MATH_LIBRARY} 
#   RRuntime 
#   PhreeqcRM
#   tug
# )

target_compile_definitions(POETLib PUBLIC STRICT_R_HEADERS OMPI_SKIP_MPICXX)

# mark_as_advanced(PHREEQCRM_BUILD_MPI PHREEQCRM_DISABLE_OPENMP)
# set(PHREEQCRM_DISABLE_OPENMP ON CACHE BOOL "" FORCE)

# option(POET_DHT_DEBUG "Build with DHT debug info" OFF)

# if(POET_DHT_DEBUG)
#   target_compile_definitions(poetlib PRIVATE DHT_STATISTICS)
# endif()

# option(POET_PHT_ADDITIONAL_INFO "Enables additional information in the PHT" OFF)

# if (POET_PHT_ADDITIONAL_INFO)
#   target_compile_definitions(poetlib PRIVATE POET_PHT_ADD)
# endif()

file(READ "${PROJECT_SOURCE_DIR}/R_lib/kin_r_library.R" R_KIN_LIB )
file(READ "${PROJECT_SOURCE_DIR}/R_lib/init_r_lib.R" R_INIT_LIB)
file(READ "${PROJECT_SOURCE_DIR}/R_lib/ai_surrogate_model_functions.R" R_AI_SURROGATE_LIB)

configure_file(poet.hpp.in poet.hpp @ONLY)


add_executable(poet poet.cpp)
target_link_libraries(poet PRIVATE POETLib MPI::MPI_C RRuntime CLI11::CLI11)
target_include_directories(poet PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")

add_executable(poet_init initializer.cpp)
target_link_libraries(poet_init PRIVATE POETLib RRuntime CLI11::CLI11)
target_include_directories(poet_init PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")

install(TARGETS poet poet_init DESTINATION bin)
