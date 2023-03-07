find_package(Doxygen)

if(DOXYGEN_FOUND)
  set(DOXYGEN_GENERATE_HTML YES)
  set(DOXYGEN_EXCLUDE_PATTERNS ${PROJECT_SOURCE_DIR}/*/argh.hpp)
  set(DOXYGEN_USE_MDFILE_AS_MAINPAGE ${PROJECT_SOURCE_DIR}/README.md)

  # See https://www.doxygen.nl/manual/customize.html
  set(DOXYGEN_DISABLE_INDEX NO)
  set(DOXYGEN_GENERATE_TREEVIEW YES)
  set(DOXYGEN_FULL_SIDEBAR YES)

  doxygen_add_docs(doxygen
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/README.md
    ${PROJECT_SOURCE_DIR}/docs/Input_Scripts.md
    ${PROJECT_SOURCE_DIR}/docs/Output.md
    COMMENT "Generate html pages")
endif()
