file(GLOB INIT_SRCS CONFIGURE_DEPENDS "initializer/*.cpp")

add_executable(poet_initializer ${INIT_SRCS})
target_link_libraries(poet_initializer RRuntime tug)