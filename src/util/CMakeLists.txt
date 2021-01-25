add_library(POET_Util STATIC RRuntime.cpp Parser.cpp)
target_include_directories(POET_Util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${R_INCLUDE_DIRS})
target_link_libraries(POET_Util ${R_LIBRARIES})