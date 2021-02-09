configure_file(poet.h.in poet.h)

add_executable(poet poet.cpp)
target_include_directories(poet PUBLIC "${PROJECT_BINARY_DIR}/src")
target_link_libraries(poet PUBLIC POET_Model POET_Util MPI::MPI_C)
target_compile_definitions(poet PRIVATE OMPI_SKIP_MPICXX)

install(TARGETS poet DESTINATION bin)

add_subdirectory(DHT)
add_subdirectory(model)
add_subdirectory(util)
