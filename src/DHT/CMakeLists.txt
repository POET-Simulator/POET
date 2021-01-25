add_library(DHT SHARED DHT.c DHT.h)
target_include_directories(DHT PRIVATE ${MPI_C_INCLUDE_DIRS})
target_link_libraries(DHT MPI::MPI_C)
target_compile_definitions(DHT PUBLIC OMPI_SKIP_MPICXX)

add_library(DHT_Wrapper STATIC DHT_Wrapper.cpp DHT_Wrapper.h)
target_include_directories(DHT_Wrapper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(DHT_Wrapper DHT m crypto POET_Util)