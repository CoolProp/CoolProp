
#if(NOT COOLPROP_DEBUG_CMAKE)
#  message(FATAL_ERROR "Debug option missing, aborting.")
#else
if(NOT COOLPROP_SOURCES_ROOT)
  message(FATAL_ERROR "Root directory missing, aborting.")
elseif(NOT COOLPROP_VERSION)
  message(FATAL_ERROR "Version information missing, aborting.")
elseif(NOT COOLPROP_INPUT_FILE)
  message(FATAL_ERROR "Input file missing, aborting.")
elseif(NOT COOLPROP_OUTPUT_FILE)
  message(FATAL_ERROR "Output file missing, aborting.")
endif()

if(NOT GIT_FOUND)
  if(COOLPROP_DEBUG_CMAKE)
    find_package(Git)
  else()
    find_package(Git QUIET)
  endif()
endif()
if(GIT_FOUND)
  # Get branch (refspec)
  execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY ${COOLPROP_SOURCES_ROOT}
      OUTPUT_VARIABLE COOLPROP_GIT_BRANCH
      OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  # Get hash 
  #execute_process(
  #    COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
  execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
      WORKING_DIRECTORY ${COOLPROP_SOURCES_ROOT}
      OUTPUT_VARIABLE COOLPROP_GIT_HASH
      OUTPUT_STRIP_TRAILING_WHITESPACE
  )
else()
  if(COOLPROP_DEBUG_CMAKE)
    message(STATUS "Could not find Git - using dummy variables.")
  endif()
  set(COOLPROP_GIT_BRANCH "COOLPROP_GIT_BRANCH-NOTFOUND")
  set(COOLPROP_GIT_HASH "COOLPROP_GIT_HASH-NOTFOUND")
endif()

if(COOLPROP_GIT_HASH)
  configure_file("${COOLPROP_INPUT_FILE}" "${COOLPROP_INPUT_FILE}.tmp" @ONLY)
  file(SHA1 "${COOLPROP_INPUT_FILE}.tmp" COOLPROP_NEW_FILE_HASH)
  if(EXISTS "${COOLPROP_OUTPUT_FILE}")
    file(SHA1 "${COOLPROP_OUTPUT_FILE}" COOLPROP_OLD_FILE_HASH)
  else()
    set(COOLPROP_OLD_FILE_HASH "COOLPROP_OLD_FILE_HASH-NOTFOUND")
  endif()
  if("${COOLPROP_NEW_FILE_HASH}" STREQUAL "${COOLPROP_OLD_FILE_HASH}")
    if(COOLPROP_DEBUG_CMAKE)
      message(STATUS "${COOLPROP_OUTPUT_FILE} is up to date")
    endif()
  else()
    #configure_file("${COOLPROP_INPUT_FILE}" "${COOLPROP_OUTPUT_FILE}" @ONLY)
    file(RENAME "${COOLPROP_INPUT_FILE}.tmp" "${COOLPROP_OUTPUT_FILE}")
    if(COOLPROP_DEBUG_CMAKE)
      message(STATUS "${COOLPROP_OUTPUT_FILE} written to file")
    endif()
  endif()
else()
  message(WARNING "Could not get git hash, leaving '${COOLPROP_OUTPUT_FILE}' untouched")
endif()
