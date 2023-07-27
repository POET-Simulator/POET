configure_file(poet.h.in poet.h)

if(POET_USE_PRM_BACKEND)
  set(poet_SRC poet_prm.cpp)
else()
  set(poet_SRC poet.cpp)
endif()

add_executable(poet ${poet_SRC})
target_include_directories(poet PUBLIC "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(poet PUBLIC poet_lib MPI::MPI_CXX)
#target_compile_definitions(poet PRIVATE OMPI_SKIP_MPICXX)

install(TARGETS poet DESTINATION bin)
