# Try to find the build flags to compile octave shared objects (oct and mex files)
# Once done this will define
#
# OCTAVE_VERSION - Version of Octave
# OCTAVE_FOUND - if Coin3d is found
# OCTAVE_CXXFLAGS - extra flags
# OCTAVE_INCLUDE_DIRS - include directories
# OCTAVE_LINK_DIRS - link directories
# OCTAVE_LIBRARY_RELEASE - the release version
# OCTAVE_LIBRARY_DEBUG - the debug version
# OCTAVE_LIBRARY - a default library, with priority debug.

cmake_minimum_required(VERSION 2.8)

IF (WIN32)
  IF( "$ENV{OCTAVE_ROOT}" STREQUAL "" )
    message(FATAL_ERROR "On windows, environmental variable OCTAVE_ROOT must be set to folder containing folders bin, include, etc. for octave")
  ENDIF()
ENDIF()

set(OCTAVE_BIN)
IF( "$ENV{OCTAVE_ROOT}" STREQUAL "" )
ELSE()
    set(OCTAVE_BIN $ENV{OCTAVE_ROOT}/bin)
ENDIF()

# use mkoctfile
set(MKOCTFILE_EXECUTABLE MKOCTFILE_EXECUTABLE-NOTFOUND)
find_program(MKOCTFILE_EXECUTABLE 
             NAME mkoctfile 
             PATHS ${OCTAVE_BIN})
mark_as_advanced(MKOCTFILE_EXECUTABLE)

if(MKOCTFILE_EXECUTABLE)
  set(OCTAVE_FOUND 1)
  message(STATUS "Found mkoctfile executable")

  execute_process(
    COMMAND ${MKOCTFILE_EXECUTABLE} -p ALL_CXXFLAGS
    OUTPUT_VARIABLE _mkoctfile_cppflags
    RESULT_VARIABLE _mkoctfile_failed)
  string(REGEX REPLACE "[\r\n]" " " _mkoctfile_cppflags "${_mkoctfile_cppflags}")
  execute_process(
    COMMAND ${MKOCTFILE_EXECUTABLE} -p INCFLAGS
    OUTPUT_VARIABLE _mkoctfile_includedir
    RESULT_VARIABLE _mkoctfile_failed)
  string(REGEX REPLACE "[\r\n]" " " _mkoctfile_includedir "${_mkoctfile_includedir}")
  string(REGEX REPLACE "-I" " " _mkoctfile_includedir "${_mkoctfile_includedir}")
  execute_process(
    COMMAND ${MKOCTFILE_EXECUTABLE} -p ALL_LDFLAGS
    OUTPUT_VARIABLE _mkoctfile_ldflags
    RESULT_VARIABLE _mkoctfile_failed)
  string(REGEX REPLACE "[\r\n]" " " _mkoctfile_ldflags "${_mkoctfile_ldflags}")
  execute_process(
    COMMAND ${MKOCTFILE_EXECUTABLE} -p LFLAGS
    OUTPUT_VARIABLE _mkoctfile_ldirs
    RESULT_VARIABLE _mkoctfile_failed)
    
  string(REGEX REPLACE "[\r\n]" " " _mkoctfile_ldirs "${_mkoctfile_ldirs}")
  string(REGEX REPLACE "-L" "" _mkoctfile_ldirs "${_mkoctfile_ldirs}")
  
  separate_arguments(_mkoctfile_ldirs)
    
  execute_process(
    COMMAND ${MKOCTFILE_EXECUTABLE} -p LIBS
    OUTPUT_VARIABLE _mkoctfile_libs
    RESULT_VARIABLE _mkoctfile_failed)
  string(REGEX REPLACE "[\r\n]" " " _mkoctfile_libs "${_mkoctfile_libs}")
  execute_process(
    COMMAND ${MKOCTFILE_EXECUTABLE} -p OCTAVE_LIBS
    OUTPUT_VARIABLE _mkoctfile_octlibs
    RESULT_VARIABLE _mkoctfile_failed)
  string(REGEX REPLACE "[\r\n]" " " _mkoctfile_octlibs "${_mkoctfile_octlibs}")
  set(_mkoctfile_libs "${_mkoctfile_libs} ${_mkoctfile_octlibs}")
    
  string(REGEX MATCHALL "(^| )-l([./+-_\\a-zA-Z]*)" _mkoctfile_libs "${_mkoctfile_libs}")
  string(REGEX REPLACE "(^| )-l" "" _mkoctfile_libs "${_mkoctfile_libs}")

#~   string(REGEX MATCHALL "(^| )-L([./+-_\\a-zA-Z]*)" _mkoctfile_ldirs "${_mkoctfile_ldirs}")
#~   string(REGEX REPLACE "(^| )-L" "" _mkoctfile_ldirs "${_mkoctfile_ldirs}")

  string(REGEX REPLACE "(^| )-l([./+-_\\a-zA-Z]*)" " " _mkoctfile_ldflags "${_mkoctfile_ldflags}")
  string(REGEX REPLACE "(^| )-L([./+-_\\a-zA-Z]*)" " " _mkoctfile_ldflags "${_mkoctfile_ldflags}")

  if (WIN32)
    
    string(REGEX REPLACE "Program Files " "Program~Files~" _mkoctfile_ldirs "${_mkoctfile_ldirs}")
    separate_arguments(_mkoctfile_ldirs)
    string(REGEX REPLACE "Program Files " "Program~Files~" _mkoctfile_includedir "${_mkoctfile_includedir}")
    separate_arguments(_mkoctfile_includedir)
    
    set(includes)
    foreach(ITR ${_mkoctfile_includedir})
      #string(REGEX REPLACE "~" " " ITR ${ITR})
      string(REGEX REPLACE "\"" "" ITR ${ITR})
      list(APPEND includes ${ITR})
    endforeach()
    set(_mkoctfile_includedir ${includes})
    
    set(libs)
    foreach(ITR ${_mkoctfile_ldirs})
      #string(REGEX REPLACE "~" " " ITR ${ITR})
      string(REGEX REPLACE "\"" "" ITR ${ITR})
      list(APPEND libs ${ITR})
    endforeach()
    set(_mkoctfile_ldirs ${libs})
    message(STATUS ${libs})
    
  else()
   separate_arguments(_mkoctfile_includedir)
  endif()

  set( OCTAVE_CXXFLAGS " ${_mkoctfile_cppflags}" )
  set( OCTAVE_LINK_FLAGS " ${_mkoctfile_ldflags} " )
  set( OCTAVE_INCLUDE_DIRS " ${_mkoctfile_includedir}")
  set( OCTAVE_LINK_DIRS ${_mkoctfile_ldirs})
  set( OCTAVE_LIBRARY ${_mkoctfile_libs})
  set( OCTAVE_LIBRARY_RELEASE " ${OCTAVE_LIBRARY} ")
  set( OCTAVE_LIBRARY_DEBUG " ${OCTAVE_LIBRARY} ")
else()
    if (OSX)
    set(libs)
    foreach(ITR ${_mkoctfile_ldirs})
      string(REGEX REPLACE "\"" "" ITR ${ITR})
      list(APPEND libs ${ITR})
    endforeach()
    set(_mkoctfile_ldirs ${libs})
    endif()	

	message(FATAL_ERROR "Unable to find mkoctfile executable")
endif()

# use octave_config
set(OCTAVE_CONFIG_EXECUTABLE OCTAVE_CONFIG_EXECUTABLE-NOTFOUND)
find_program(OCTAVE_CONFIG_EXECUTABLE 
             NAME octave-config 
             PATHS ${OCTAVE_BIN})
mark_as_advanced(OCTAVE_CONFIG_EXECUTABLE)

if(OCTAVE_CONFIG_EXECUTABLE)
  message(STATUS "Found octave-config executable")
  execute_process(
    COMMAND ${OCTAVE_CONFIG_EXECUTABLE} -v
    OUTPUT_VARIABLE OCTAVE_VERSION
    RESULT_VARIABLE _octave_config_failed)
  string(REGEX REPLACE "[\r\n]" "" OCTAVE_VERSION "${OCTAVE_VERSION}")    
  execute_process(
    COMMAND ${OCTAVE_CONFIG_EXECUTABLE} -p CANONICAL_HOST_TYPE
    OUTPUT_VARIABLE _octave_config_host_type
    RESULT_VARIABLE _octave_config_failed)
  string(REGEX REPLACE "[\r\n]" "" _octave_config_host_type "${_octave_config_host_type}")
  execute_process(
    COMMAND ${OCTAVE_CONFIG_EXECUTABLE} -p API_VERSION
    OUTPUT_VARIABLE _octave_config_api_version
    RESULT_VARIABLE _octave_config_failed)
  string(REGEX REPLACE "[\r\n]" "" _octave_config_api_version "${_octave_config_api_version}")
  execute_process(
    COMMAND ${OCTAVE_CONFIG_EXECUTABLE} -p LOCALVEROCTFILEDIR
    OUTPUT_VARIABLE _octave_config_localveroctfiledir
    RESULT_VARIABLE _octave_config_failed)
  string(REGEX REPLACE "[\r\n]" "" _octave_config_localveroctfiledir "${_octave_config_api_version}")

  set( OCTAVE_HOST_TYPE "${_octave_config_host_type}" )
  set( OCTAVE_API_VERSION "${_octave_config_api_version}" )
  set( OCTAVE_LOCALVEROCTFILEDIR "${_octave_config_localveroctfiledir}" )

else()
  message(FATAL_ERROR "Did not find octave-config executable")
endif()

IF (WIN32)
  list(APPEND OCTAVE_LINK_DIRS "$ENV{OCTAVE_ROOT}/lib/octave/${OCTAVE_VERSION}")
  message(STATUS "OCTAVE_LINK_DIRS :: ${OCTAVE_LINK_DIRS}")
  list(APPEND OCTAVE_INCLUDE_DIRS "$ENV{OCTAVE_ROOT}/include")
  list(APPEND OCTAVE_INCLUDE_DIRS "$ENV{OCTAVE_ROOT}/include/octave-${OCTAVE_VERSION}")
  list(APPEND OCTAVE_INCLUDE_DIRS "$ENV{OCTAVE_ROOT}/include/octave-${OCTAVE_VERSION}/octave")
  message(STATUS "OCTAVE_INCLUDE_DIRS :: ${OCTAVE_INCLUDE_DIRS}")
ENDIF()

FIND_LIBRARY( OCTAVE_OCTAVE_LIBRARY
			  NAMES octave liboctave
			  PATHS ${OCTAVE_LINK_DIRS}				  
			  NO_DEFAULT_PATH)
			  
FIND_LIBRARY( OCTAVE_OCTINTERP_LIBRARY
			  NAMES octinterp liboctinterp
			  PATHS ${OCTAVE_LINK_DIRS}				  
			  NO_DEFAULT_PATH)
			  
FIND_LIBRARY( OCTAVE_CRUFT_LIBRARY
			  NAMES cruft libcruft
			  PATHS ${OCTAVE_LINK_DIRS}				  
			  NO_DEFAULT_PATH)			  

SET(OCTAVE_LIBRARIES
  ${OCTAVE_OCTAVE_LIBRARY}
  ${OCTAVE_OCTINTERP_LIBRARY})

if (OCTAVE_CRUFT_LIBRARY)
    list(APPEND OCTAVE_LIBRARIES ${OCTAVE_CRUFT_LIBRARY})
endif()
				  
message(STATUS "OCTAVE_VERSION=${OCTAVE_VERSION}" )
message(STATUS "OCTAVE_OCTAVE_LIBRARY=${OCTAVE_OCTAVE_LIBRARY}")
message(STATUS "OCTAVE_OCTINTERP_LIBRARY=${OCTAVE_OCTINTERP_LIBRARY}")
message(STATUS "OCTAVE_CXXFLAGS=${_mkoctfile_cppflags}" )
message(STATUS "OCTAVE_LINK_FLAGS=${_mkoctfile_ldflags}" )
message(STATUS "OCTAVE_INCLUDE_DIRS=${_mkoctfile_includedir}")
message(STATUS "OCTAVE_LINK_DIRS=${_mkoctfile_ldirs}")
message(STATUS "OCTAVE_LIBRARY=${_mkoctfile_libs}")
message(STATUS "OCTAVE_LIBRARY_RELEASE=${OCTAVE_LIBRARY} ")
message(STATUS "OCTAVE_LIBRARY_DEBUG=${OCTAVE_LIBRARY} ")
message(STATUS "OCTAVE_LIBRARIES=${OCTAVE_LIBRARIES} ")

MARK_AS_ADVANCED(
    OCTAVE_LIBRARY_FOUND
    OCTAVE_VERSION
    OCTAVE_CXXFLAGS
    OCTAVE_LINK_FLAGS
    OCTAVE_INCLUDE_DIRS
    OCTAVE_LINK_DIRS
    OCTAVE_LIBRARY
    OCTAVE_LIBRARY_RELEASE
    OCTAVE_LIBRARY_DEBUG
)
