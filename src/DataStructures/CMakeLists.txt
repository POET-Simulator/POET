add_library(DataStructures OBJECT 
    Field.cpp
)

target_link_libraries(DataStructures
    PUBLIC RRuntime
)

target_include_directories(DataStructures
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)