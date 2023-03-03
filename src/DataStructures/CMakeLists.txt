file(GLOB DataStructures_SRC
  CONFIGURE_DEPENDS
  "*.cpp" "*.c")

add_library(DataStructures ${DataStructures_SRC})
target_include_directories(DataStructures PUBLIC ${PROJECT_SOURCE_DIR}/include)
