file(GLOB surrogate_models_SRC
  CONFIGURE_DEPENDS
  "*.cpp" "*.c")

find_library(MATH_LIBRARY m)

add_library(surrogate_models ${surrogate_models_SRC})
target_include_directories(surrogate_models PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(surrogate_models PUBLIC
    MPI::MPI_CXX ${MATH_LIBRARY})
target_compile_definitions(surrogate_models PUBLIC OMPI_SKIP_MPICXX)

option(POET_PHT_ADDITIONAL_INFO "Enables additional information in the PHT" OFF)

if (POET_PHT_ADDITIONAL_INFO)
  target_compile_definitions(surrogate_models PRIVATE POET_PHT_ADD)
endif()

if(POET_DHT_DEBUG)
  target_compile_definitions(surrogate_models PRIVATE DHT_STATISTICS)
endif()
