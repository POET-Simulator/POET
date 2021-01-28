add_library(POET_Model ChemSim.cpp ChemSim.h ChemMaster.cpp ChemWorker.cpp TransportSim.cpp Grid.cpp)
target_include_directories(POET_Model PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MPI_C_INCLUDE_DIRS})
target_link_libraries(POET_Model PRIVATE POET_Util MPI::MPI_C PUBLIC DHT_Wrapper)
target_compile_definitions(POET_Model PRIVATE OMPI_SKIP_MPICXX)