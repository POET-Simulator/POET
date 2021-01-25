add_executable(poet kin.cpp)
target_link_libraries(poet PRIVATE POET_Model POET_Util)

add_subdirectory(DHT)
add_subdirectory(model)
add_subdirectory(util)
