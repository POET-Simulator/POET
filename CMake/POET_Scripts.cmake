# Set or get version
macro(get_POET_version)
  if(EXISTS ${PROJECT_SOURCE_DIR}/.git)
    find_program(GIT_EXECUTABLE git DOC "git executable")
    mark_as_advanced(GIT_EXECUTABLE)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      OUTPUT_VARIABLE POET_GIT_BRANCH
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe --always
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      OUTPUT_VARIABLE POET_GIT_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT POET_GIT_BRANCH STREQUAL "master")
      set(POET_VERSION "${POET_GIT_BRANCH}/${POET_GIT_VERSION}")
    else()
      set(POET_VERSION "${POET_GIT_VERSION}")
    endif()
  elseif(EXISTS ${PROJECT_SOURCE_DIR}/.svn)
    file(STRINGS .gitversion POET_VERSION)
  else()
    set(POET_VERSION "0.1")
  endif()

  message(STATUS "Configuring POET version ${POET_VERSION}")
endmacro(get_POET_version)
