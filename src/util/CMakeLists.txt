add_library(POET_Util RRuntime.cpp SimParams.cpp)
target_include_directories(POET_Util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${R_INCLUDE_DIRS})
target_link_libraries(POET_Util ${R_LIBRARIES})
target_compile_definitions(POET_Util PUBLIC STRICT_R_HEADERS)