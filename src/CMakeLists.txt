option(POET_DAOSKV "Build with DAOS, using the KV-API" OFF)
option(POET_DAOSOB "Build with DAOS, using the Object-API" OFF)

file(GLOB_RECURSE poet_lib_SRC
  CONFIGURE_DEPENDS
  "*.cpp" "*.c")

if(POET_DAOSKV OR POET_DAOSOB)
  #list(REMOVE_ITEM poet_lib_SRC "${CMAKE_CURRENT_SOURCE_DIR}/ChemistryModule/SurrogateModels/DHT.c")
  list(REMOVE_ITEM poet_lib_SRC "${CMAKE_CURRENT_SOURCE_DIR}/ChemistryModule/SurrogateModels/DHT_Wrapper.cpp")

  if (POET_DAOSKV)
    message(NOTICE "Using DAOS KV")
    list(REMOVE_ITEM poet_lib_SRC "${CMAKE_CURRENT_SOURCE_DIR}/ChemistryModule/DAOS/OB/DaosKeyValue.c")
  else () 
    message(NOTICE "Using DAOS OB")
    list(REMOVE_ITEM poet_lib_SRC "${CMAKE_CURRENT_SOURCE_DIR}/ChemistryModule/DAOS/KV/DaosKeyValue.c")
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
    set(PHREEQCRM_BUILD_MPI OFF CACHE BOOL "" FORCE)

else()
  list(REMOVE_ITEM poet_lib_SRC "${CMAKE_CURRENT_SOURCE_DIR}/ChemistryModule/DAOS/DHT_Wrapper.cpp")
  list(REMOVE_ITEM poet_lib_SRC "${CMAKE_CURRENT_SOURCE_DIR}/ChemistryModule/DAOS/KV/DaosKeyValue.c")
  list(REMOVE_ITEM poet_lib_SRC "${CMAKE_CURRENT_SOURCE_DIR}/ChemistryModule/DAOS/OB/DaosKeyValue.c")
endif()

foreach(src IN LISTS poet_lib_SRC)
  message(STATUS "  ${src}")
endforeach()

add_library(poet_lib ${poet_lib_SRC})
target_include_directories(poet_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(poet_lib PUBLIC
  MPI::MPI_CXX ${MATH_LIBRARY} RRuntime PhreeqcRM tug)

  if(POET_DAOSKV OR POET_DAOSOB)
  target_link_libraries(poet_lib PUBLIC DAOS)
  endif()

target_compile_definitions(poet_lib PUBLIC STRICT_R_HEADERS OMPI_SKIP_MPICXX)

mark_as_advanced(PHREEQCRM_BUILD_MPI PHREEQCRM_DISABLE_OPENMP)
set(PHREEQCRM_DISABLE_OPENMP ON CACHE BOOL "" FORCE)

option(POET_DHT_DEBUG "Build with DHT debug info" OFF)

if(POET_DHT_DEBUG)
  target_compile_definitions(poet_lib PRIVATE DHT_STATISTICS)
endif()

option(POET_PHT_ADDITIONAL_INFO "Enables additional information in the PHT" OFF)

if (POET_PHT_ADDITIONAL_INFO)
  target_compile_definitions(poet_lib PRIVATE POET_PHT_ADD)
endif()
