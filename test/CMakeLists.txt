file(GLOB test_SRC
  CONFIGURE_DEPENDS
  "*.cpp" "*.c")

find_package(doctest REQUIRED)

add_executable(testPOET ${test_SRC})
target_link_libraries(testPOET doctest::doctest POETLib)
target_include_directories(testPOET PRIVATE "${PROJECT_SOURCE_DIR}/src")

get_filename_component(TEST_RInsideSourceFile "RInsidePOET_funcs.R" REALPATH)
configure_file(testDataStructures.hpp.in testDataStructures.hpp)
target_include_directories(testPOET PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")

add_custom_target(check
    COMMAND $<TARGET_FILE:testPOET>
    DEPENDS testPOET
)
