file(GLOB test_SRC
  CONFIGURE_DEPENDS
  "*.cpp" "*.c")

add_executable(testPOET ${test_SRC})
target_link_libraries(testPOET doctest poetlib)
target_include_directories(testPOET PRIVATE "${PROJECT_SOURCE_DIR}/src")

get_filename_component(TEST_RInsideSourceFile "RInsidePOET_funcs.R" REALPATH)
get_filename_component(TEST_SampleInputScript "${PROJECT_SOURCE_DIR}/test/advection/dolo_adv.R" REALPATH)
configure_file(testInputFiles.hpp.in InputFiles.hpp)
target_include_directories(testPOET PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")

add_custom_target(check
    COMMAND $<TARGET_FILE:testPOET>
    DEPENDS testPOET
)

add_executable(advection advection/testAdvection.cpp)
target_include_directories(advection PRIVATE "${CMAKE_CURRENT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(advection poetlib)
