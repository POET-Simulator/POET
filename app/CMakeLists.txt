get_poet_version()

configure_file(poet.h.in poet.h)

add_executable(poet poet.cpp)
target_include_directories(poet PUBLIC "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(poet PUBLIC poet_lib MPI::MPI_C)
#target_compile_definitions(poet PRIVATE OMPI_SKIP_MPICXX)

install(TARGETS poet DESTINATION bin)
