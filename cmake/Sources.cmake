#######################################
#         FIND ALL SOURCES            #
#-------------------------------------#
# The project is organised with       #
# split includes and source folders   #
# this makes it easier for developers #
# to quickly find relevant includes.  #
# This section finds all sources,     #
# headers and corresponding dirs.     #
#######################################

# These backends will be compiled in
set(COOLPROP_ENABLED_BACKENDS
    Cubics
    IF97
    Helmholtz
    REFPROP
    Incompressible
    Tabular
    PCSAFT)

# Get everything in the src/ directory (always), but not recursive
file(GLOB APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

# Add the miniz source file
list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/externals/miniz-3.0.2/miniz.c")

# For each enabled backend, grab its files
foreach(backend ${COOLPROP_ENABLED_BACKENDS})
  file(GLOB_RECURSE BACKEND_SOURCES
       "${CMAKE_CURRENT_SOURCE_DIR}/src/Backends/${backend}/*.cpp")
  list(APPEND APP_SOURCES ${BACKEND_SOURCES})
endforeach()

## You can exclude this file, in case you want to run your own tests or use Catch
list(REMOVE_ITEM APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Tests/Tests.cpp")
list(REMOVE_ITEM APP_SOURCES
     "${CMAKE_CURRENT_SOURCE_DIR}/src/Tests/CoolProp-Tests.cpp")

## This file is only needed for the library, normal builds do not need it.
list(REMOVE_ITEM APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolPropLib.cpp")

# Set up include directories
set(APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/Eigen")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/msgpack-c/include")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/miniz-3.0.2")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/nlohmann-json")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/incbin")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/dev")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/boost_CoolProp")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/fmtlib/include")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/fmtlib") # should be deprecated
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")

# MSVC-specific flags for fmtlib
if(MSVC)
  # fmtlib requires that the utf-8 support be compiled in
  # TODO: add the fmt target from fmtlib directly which does this
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /utf-8 -D_CRT_SECURE_NO_WARNINGS")
endif()

## Set endianess for msgpack on ARM64 with MSVC
#if(MSVC)
#  if("${CMAKE_GENERATOR_PLATFORM}" STREQUAL "ARM64")
#    message(STATUS "Forcing msgpack-c to use little endian configuration")
#	add_compile_definitions(MSGPACK_ENDIAN_LITTLE_BYTE)
#  endif()
#endif()

include_directories(${APP_INCLUDE_DIRS})

# SWIG dependencies for wrapper modules
set(SWIG_DEPENDENCIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/DataStructures.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/CoolProp.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/AbstractState.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Configuration.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/PhaseEnvelope.h)

# Cache the sources and include dirs for use in other modules
set(COOLPROP_APP_SOURCES
    "${APP_SOURCES}"
    CACHE STRING "List of CPP sources needed for CoolProp")
set(COOLPROP_INCLUDE_DIRECTORIES
    "${APP_INCLUDE_DIRS}"
    CACHE STRING "List of include directories needed for CoolProp")
