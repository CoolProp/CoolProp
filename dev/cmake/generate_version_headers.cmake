

#option(COOLPROP_USE_RPAVLIK_CMAKE OFF) #use CMake modules from https://github.com/rpavlik/cmake-modules
#  # Generate the version information headers, forces CMake to rerun after git commit
#  include(GetGitRevisionDescription)
#  get_git_head_revision(COOLPROP_GIT_REF COOLPROP_GIT_HASH)
#  git_local_changes(COOLPROP_GIT_CHG)
#  #message(STATUS "COOLPROP_GIT_REF: ${COOLPROP_GIT_REF}")
#  #message(STATUS "COOLPROP_GIT_HASH: ${COOLPROP_GIT_HASH}")
#  #message(STATUS "COOLPROP_GIT_CHNG: ${COOLPROP_GIT_CHNG}")
#  if(COOLPROP_GIT_HASH)
#    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake/cpversion.h.in"
#        "${CMAKE_CURRENT_SOURCE_DIR}/include/cpversion.h"
#        @ONLY)
#  else()
#    message(STATUS "Could not get the git hash, leaving '${CMAKE_CURRENT_SOURCE_DIR}/include/cpversion.h' untouched.")
#  endif()

option(COOLPROP_DEBUG_CMAKE OFF)

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
  file(SHA1 "${COOLPROP_OUTPUT_FILE}" COOLPROP_OLD_FILE_HASH)
  if("${COOLPROP_NEW_FILE_HASH}" STREQUAL "${COOLPROP_NEW_FILE_HASH}")
    if(COOLPROP_DEBUG_CMAKE)
      message(STATUS "${COOLPROP_OUTPUT_FILE} is up to date")
    endif()
  else()
    #configure_file("${COOLPROP_INPUT_FILE}" "${COOLPROP_OUTPUT_FILE}" @ONLY)
    file(COPY "${COOLPROP_INPUT_FILE}.tmp" DESTINATION "${COOLPROP_OUTPUT_FILE}")
    if(COOLPROP_DEBUG_CMAKE)
      message(STATUS "${COOLPROP_OUTPUT_FILE} written to file")
    endif()
  endif()
else()
  if(COOLPROP_DEBUG_CMAKE)
    message(WARNING "Could not get git hash, leaving '${COOLPROP_OUTPUT_FILE}' untouched")
  endif()
endif()
