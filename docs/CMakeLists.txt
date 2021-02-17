find_package(Doxygen)
if(DOXYGEN_FOUND)

  set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/./Doxyfile.in)
  set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/Doxyfile)

  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

  add_custom_target(doc_doxygen ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating documantation with doxygen"
    VERBATIM)

  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxygen/ DESTINATION docs)
else(DOXYGEN_FOUND)
  message("Doxygen not found. Please install.")
endif(DOXYGEN_FOUND)
