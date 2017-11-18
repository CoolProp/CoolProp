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
  list(APPEND COOLPROP_MAIN_SOURCES ${COOLPROP_BACKEND_SOURCES})
  list(APPEND COOLPROP_MAIN_HEADERS ${COOLPROP_BACKEND_HEADERS})
endforeach()

# This file is only needed for the library, normal builds do not need it.
set(COOLPROP_LIBRARY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolPropLib.cpp")
list(REMOVE_ITEM COOLPROP_MAIN_SOURCES "${COOLPROP_LIBRARY_SOURCES}")

# Add external dependencies
set (APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}")
list (APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/Eigen")
list (APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/msgpack-c/include")
list (APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/cppformat")
list (APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
list (APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")

include_directories(${APP_INCLUDE_DIRS})

set(LIB_NAME coolprop)
add_library(${LIB_NAME} STATIC ${COOLPROP_MAIN_SOURCES} ${COOLPROP_LIBRARY_SOURCES} ${COOLPROP_LIBRARY_EXPORTS})
if(MSVC)
  set_property(TARGET ${LIB_NAME} PROPERTY DEBUG_POSTFIX d)
  set_property(TARGET ${LIB_NAME} PROPERTY RELEASE_POSTFIX )
endif()
set_property(TARGET ${LIB_NAME} PROPERTY VERSION ${COOLPROP_VERSION})
set_property(TARGET ${LIB_NAME} PROPERTY SOVERSION ${COOLPROP_VERSION_MAJOR})