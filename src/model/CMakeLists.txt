add_library(POET_Model STATIC ChemSim.cpp ChemSim.h ChemMaster.cpp ChemWorker.cpp TransportSim.cpp Grid.cpp)
target_include_directories(POET_Model PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(POET_Model POET_Util DHT_Wrapper)