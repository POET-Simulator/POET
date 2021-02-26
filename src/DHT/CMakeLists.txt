add_library(DHT DHT.c DHT.h)
target_include_directories(DHT PUBLIC ${MPI_C_INCLUDE_DIRS})
target_link_libraries(DHT PRIVATE MPI::MPI_C)
target_compile_definitions(DHT PUBLIC OMPI_SKIP_MPICXX)

option(DHT_Debug "Toggle output of stastic table for each iteration" OFF)

if (DHT_Debug)
  target_compile_definitions(DHT PUBLIC DHT_STATISTICS)
endif()

find_library(MATH_LIBRARY m)
find_library(CRYPTO_LIBRARY crypto)

add_library(DHT_Wrapper DHT_Wrapper.cpp DHT_Wrapper.h)
target_include_directories(DHT_Wrapper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(DHT_Wrapper PRIVATE DHT ${MATH_LIBRARY} ${CRYPTO_LIBRARY} POET_Util)
