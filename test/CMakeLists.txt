file(GLOB test_SRC
  CONFIGURE_DEPENDS
  "*.cpp" "*.c")

add_executable(testPOET ${test_SRC})
target_link_libraries(testPOET doctest poet_lib)

add_custom_target(check
    COMMAND $<TARGET_FILE:testPOET>
    DEPENDS testPOET
)
