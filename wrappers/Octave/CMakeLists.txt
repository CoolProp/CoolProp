
cmake_minimum_required(VERSION 2.6)
set(CMAKE_USE_ABSOLUTE_PATHS ON)
set(project-name "CoolProp")
project(${project-name})

#Check for the required packaes
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
FIND_PACKAGE(Octave REQUIRED)

set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../../include")             # source reside in src/
file(GLOB_RECURSE APP_SOURCES "${CMAKE_SOURCE_DIR}/../../src/*.cpp")             # source reside in src/

SET(OCTAVE_WRAP_INCLUDE_DIRS ${INCLUDE_DIR})
find_library(OCTAVE_LIB octave PATHS ${OCTAVE_LINK_DIRS})

foreach(ITR ${OCTAVE_INCLUDE_DIRS})
  list(APPEND OCTAVE_WRAP_INCLUDE_DIRS ${ITR})
endforeach()
    
include_directories(${OCTAVE_WRAP_INCLUDE_DIRS})

message(STATUS ${OCTAVE_INCLUDE_DIRS})

set(I_FILE "CoolProp.i")

SET_SOURCE_FILES_PROPERTIES(${I_FILE} PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(CoolProp octave ${I_FILE} ${APP_SOURCES})
SWIG_LINK_LIBRARIES(CoolProp ${OCTAVE_LIB})

set_target_properties(CoolProp PROPERTIES SUFFIX ".oct" PREFIX "")

message(STATUS ${OCTAVE_LIB})

#~ get_cmake_property(_variableNames VARIABLES)
#~ foreach (_variableName ${_variableNames})
#~     message(STATUS "${_variableName}=${${_variableName}}")
#~ endforeach()

#~ find_program(MKOCTFILE_EXECUTABLE mkoctfile)
#~ if(NOT MKOCTFILE_EXECUTABLE)
#~  message(SEND_ERROR "Failed to find mkoctfile, is it in your $PATH?")
#~ endif()

#~ foreach (_variableName ${_variableNames})
#~     message(STATUS "${_variableName}=${${_variableName}}")
#~ endforeach()

#~ ## Compile the octfile
#~ set(OCT_CXXFLAGS
#~   "-O2 -Wall -Wextra -std=c++0x")

#~ set(ENV{CXXFLAGS} ${OCT_CXXFLAGS})

#~ set(OCT_FILE
#~   "bwlabeln.oct")

#~ set(OCT_SRC
#~   "${CMAKE_CURRENT_SOURCE_DIR}/bwlabeln.cc")

#~ add_custom_command(OUTPUT ${OCT_FILE}
#~   COMMAND ${MKOCTFILE_EXECUTABLE} "${OCT_SRC}"
#~   DEPENDS "bwlabeln.cc"
#~   COMMENT "Generating ${OCT_FILE}"
#~   VERBATIM)

#~ add_custom_target(bwlabeln ALL
#~   DEPENDS ${OCT_FILE})

#~ ## Install stuff
#~ file(GLOB m-files "*.m")

#~ set(data_dir "share/${project-name}")

#~ install(FILES
#~   ${m-files}
#~   DESTINATION ${data_dir})

#~ install(PROGRAMS
#~   "main.m"
#~   DESTINATION ${data_dir})

#~ install(FILES
#~   ${CMAKE_CURRENT_BINARY_DIR}/${OCT_FILE}
#~   DESTINATION ${data_dir})

#~ set(NAZ_DATA_DIR ${CMAKE_INSTALL_PREFIX}/${data_dir})
#~ configure_file(run.py.in run.py)

#~ install(PROGRAMS
#~   "${CMAKE_BINARY_DIR}/run.py"
#~   DESTINATION bin/
#~   RENAME "${project-name}")


#~ install(DIRECTORY
#~   "niak"
#~   "SB2_Release_200"
#~   "data"
#~   "Transformations"
#~   "Save_Data_with9_goodOne_Journal"
#~   DESTINATION share/${project-name}
#~   PATTERN ".hg" EXCLUDE)
