add_executable(poet poet.cpp)
target_link_libraries(poet PRIVATE POET_Model POET_Util MPI::MPI_C)
target_compile_definitions(poet PRIVATE OMPI_SKIP_MPICXX)

add_subdirectory(DHT)
add_subdirectory(model)
add_subdirectory(util)
