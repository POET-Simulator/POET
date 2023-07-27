add_subdirectory(DataStructures)

file(GLOB poet_lib_SRC
  CONFIGURE_DEPENDS
  "*.cpp" "*.c")

find_library(MATH_LIBRARY m)
find_library(CRYPTO_LIBRARY crypto)

add_library(poet_lib ${poet_lib_SRC})
target_include_directories(poet_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(poet_lib PUBLIC
  MPI::MPI_CXX ${MATH_LIBRARY} RRuntime tug PhreeqcRM DataStructures ChemistryModule)
target_compile_definitions(poet_lib PUBLIC STRICT_R_HEADERS OMPI_SKIP_MPICXX)

add_subdirectory(ChemistryModule)
