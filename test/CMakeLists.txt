file(GLOB test_SRC
  CONFIGURE_DEPENDS
  "*.cpp" "*.c")

add_executable(testPOET ${test_SRC})
target_link_libraries(testPOET doctest poet_lib)

get_filename_component(TEST_RInsideSourceFile "RInsidePOET_funcs.R" REALPATH)
configure_file(testDataStructures.hpp.in testDataStructures.hpp)
target_include_directories(testPOET PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")

add_custom_target(check
    COMMAND $<TARGET_FILE:testPOET>
    DEPENDS testPOET
)
