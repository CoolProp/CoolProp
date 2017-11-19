cmake_minimum_required(VERSION 3.0.0)


#######################################
#         PROJECT INFORMATION         #
#-------------------------------------#
# This CMakeLists.txt file is for the #
# CoolProp thermodynamic library.     #
#######################################

project("CoolProp" VERSION 6.1.1)
set(COOLPROP_VERSION ${PROJECT_VERSION})
string(REPLACE "." ";" COOLPROP_VERSION_LIST ${PROJECT_VERSION})
list(GET COOLPROP_VERSION_LIST 0 COOLPROP_VERSION_MAJOR)
list(GET COOLPROP_VERSION_LIST 1 COOLPROP_VERSION_MINOR)
list(GET COOLPROP_VERSION_LIST 2 COOLPROP_VERSION_PATCH)
set(COOLPROP_SOURCES_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
string(TIMESTAMP COOLPROP_CONFIGURE_TIMESTAMP "%Y-%m-%d %H:%M:%S")


#######################################
#         FIND ALL SOURCES            #
#-------------------------------------#
# The project is organised with split #
# includes and source folders, which  #
# makes it easier for developers to   #
# quickly find relevant files. This   #
# section finds all sources, headers  #
# and corresponding directories.      #
#######################################

# Get everything in the src/ directory (always), but not recursive
file(GLOB COOLPROP_MAIN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
file(GLOB COOLPROP_MAIN_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")

# These backends will be compiled in
set(COOLPROP_ENABLED_BACKENDS Cubics IF97 Helmholtz REFPROP Incompressible Tabular)
foreach (backend ${COOLPROP_ENABLED_BACKENDS})
  file(GLOB_RECURSE COOLPROP_BACKEND_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Backends/${backend}/*.cpp")
  file(GLOB_RECURSE COOLPROP_BACKEND_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/src/Backends/${backend}/*.h")
  #list(APPEND COOLPROP_MAIN_SOURCES ${COOLPROP_BACKEND_SOURCES})
  #list(APPEND COOLPROP_MAIN_HEADERS ${COOLPROP_BACKEND_HEADERS})
endforeach()

# This file is only needed for the library, normal builds do not need it.
set(COOLPROP_LIBRARY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolPropLib.cpp")
list(REMOVE_ITEM COOLPROP_MAIN_SOURCES "${COOLPROP_LIBRARY_SOURCES}")
set(COOLPROP_LIBRARY_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/CoolPropLib.h")
list(REMOVE_ITEM COOLPROP_MAIN_HEADERS "${COOLPROP_LIBRARY_HEADERS}")

# Create source groups, mainly for Visual Studio
source_group("General Sources" FILES ${COOLPROP_MAIN_SOURCES}) 
source_group("General Headers" FILES ${COOLPROP_MAIN_HEADERS})
source_group("Backend Sources" FILES ${COOLPROP_BACKEND_SOURCES}) 
source_group("Backend Headers" FILES ${COOLPROP_BACKEND_HEADERS})
source_group("Library Sources" FILES ${COOLPROP_LIBRARY_SOURCES}) 
source_group("Library Headers" FILES ${COOLPROP_LIBRARY_HEADERS})

# Add external dependencies
set(APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/Eigen")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/msgpack-c/include")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/cppformat")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")
include_directories(${APP_INCLUDE_DIRS})

# Add CMake dependencies
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/externals/cmake-modules") #https://github.com/rpavlik/cmake-modules

# Generate the version information headers, forces CMake to rerun after git commit
include(GetGitRevisionDescription)
get_git_head_revision(COOLPROP_GIT_REF COOLPROP_GIT_HASH)
git_local_changes(COOLPROP_GIT_CHG)
#message(STATUS "COOLPROP_GIT_REF: ${COOLPROP_GIT_REF}")
#message(STATUS "COOLPROP_GIT_HASH: ${COOLPROP_GIT_HASH}")
#message(STATUS "COOLPROP_GIT_CHNG: ${COOLPROP_GIT_CHNG}")
if(COOLPROP_GIT_HASH)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake/cpversion.h.in"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/cpversion.h"
      @ONLY)
else()
  message(STATUS "Could not get the git hash, leaving '${CMAKE_CURRENT_SOURCE_DIR}/include/cpversion.h' untouched.")
endif()
#
## Generate the fluid information headers
#add_custom_target(generate_fluid_headers
#                  COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/dev/generate_headers.py")

#######################################
#       CREATE DEFAULT TARGET         #
#-------------------------------------#
# All targets should be placed in     #
# subdirectories and they should also #
# have unique names. However, we can  #
# define one single default target to #
# simplify integration.               #
#######################################

set(TARGET_NAME coolprop)
add_library(${TARGET_NAME} ${COOLPROP_MAIN_SOURCES} ${COOLPROP_BACKEND_SOURCES} ${COOLPROP_LIBRARY_SOURCES} ${COOLPROP_MAIN_HEADERS} ${COOLPROP_BACKEND_HEADERS} ${COOLPROP_LIBRARY_HEADERS} ${COOLPROP_LIBRARY_EXPORTS})
if(MSVC)
  # Add postfix for debugging
  set_property(TARGET ${TARGET_NAME} PROPERTY DEBUG_POSTFIX d)
  set_property(TARGET ${TARGET_NAME} PROPERTY RELEASE_POSTFIX )
  # No lib prefix for the shared library
  set_property(TARGET ${TARGET_NAME} PROPERTY PREFIX "")
  #add_custom_command(TARGET ${TARGET_NAME} POST_BUILD
  #  COMMAND dumpbin /EXPORTS $<TARGET_FILE:${TARGET_NAME}> > ${CMAKE_CURRENT_BINARY_DIR}/exports.txt)
  #install(FILES ${CMAKE_CURRENT_BINARY_DIR}/exports.txt
  #  DESTINATION shared_library/${CMAKE_SYSTEM_NAME}/${BITNESS}bit${CONVENTION})
endif()
# For Linux
if("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
  set_property(TARGET ${TARGET_NAME} PROPERTY VERSION ${COOLPROP_VERSION})
  set_property(TARGET ${TARGET_NAME} PROPERTY SOVERSION ${COOLPROP_VERSION_MAJOR})
endif()
# For Mac systems
if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
  # see https://support.enthought.com/hc/en-us/articles/204469410-OS-X-GCC-Clang-and-Cython-in-10-9-Mavericks
  set_property(TARGET ${TARGET_NAME} APPEND_STRING PROPERTY COMPILE_FLAGS " -stdlib=libstdc++ -mmacosx-version-min=10.6")
  set_property(TARGET ${TARGET_NAME} APPEND_STRING PROPERTY LINK_FLAGS " -stdlib=libstdc++ -mmacosx-version-min=10.6")
endif()