cmake_minimum_required(VERSION 3.0.0)

set(COOLPROP_DEBUG_CMAKE ON)

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
  file(GLOB_RECURSE COOLPROP_BACKEND_SOURCES_TMP "${CMAKE_CURRENT_SOURCE_DIR}/src/Backends/${backend}/*.cpp")
  file(GLOB_RECURSE COOLPROP_BACKEND_HEADERS_TMP "${CMAKE_CURRENT_SOURCE_DIR}/src/Backends/${backend}/*.h")
  list(APPEND COOLPROP_BACKEND_SOURCES ${COOLPROP_BACKEND_SOURCES_TMP})
  list(APPEND COOLPROP_BACKEND_HEADERS ${COOLPROP_BACKEND_HEADERS_TMP})
endforeach()

# Create source groups, mainly for Visual Studio
source_group("General Sources" FILES ${COOLPROP_MAIN_SOURCES}) 
source_group("General Headers" FILES ${COOLPROP_MAIN_HEADERS})
source_group("Backend Sources" FILES ${COOLPROP_BACKEND_SOURCES}) 
source_group("Backend Headers" FILES ${COOLPROP_BACKEND_HEADERS})


#######################################
#          SET DEPENDENCIES           #
#-------------------------------------#
# We make use of a number of other    #
# libraries for various purposes.     #
#######################################

# Add external dependencies, provide hunter-based alternative later
set(APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/Eigen")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/msgpack-c/include")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/cppformat")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/Catch/single_include")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")
include_directories(${APP_INCLUDE_DIRS})

# Add CMake dependencies
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/externals/cmake-modules") #https://github.com/rpavlik/cmake-modules

# Generate the header files
add_custom_target(generate_version_headers
    COMMAND ${CMAKE_COMMAND} 
      -DCOOLPROP_DEBUG_CMAKE="${COOLPROP_DEBUG_CMAKE}" 
      -DCOOLPROP_SOURCES_ROOT="${COOLPROP_SOURCES_ROOT}" 
      -DCOOLPROP_VERSION="${COOLPROP_VERSION}" 
      -DCOOLPROP_INPUT_FILE=\"${COOLPROP_SOURCES_ROOT}/dev/cmake/cpversion.h.in\" 
      -DCOOLPROP_OUTPUT_FILE="${COOLPROP_SOURCES_ROOT}/include/cpversion.h" 
      -P "${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake/generate_version_headers.cmake"
    COMMENT "Generating ${COOLPROP_SOURCES_ROOT}/include/cpversion.h"
)
add_custom_target(generate_fluid_headers
    COMMAND ${CMAKE_COMMAND} 
      -DCOOLPROP_DEBUG_CMAKE="${COOLPROP_DEBUG_CMAKE}" 
      -DCOOLPROP_SOURCES_ROOT="${COOLPROP_SOURCES_ROOT}" 
      -P "${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake/generate_fluid_headers.cmake"
    COMMENT "Generating fluid headers from ${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake/generate_fluid_headers.cmake"
)
add_custom_target(generate_headers
    DEPENDS generate_version_headers generate_fluid_headers
  )


#######################################
#       CREATE DEFAULT TARGET         #
#-------------------------------------#
# All targets should be placed in     #
# subdirectories and they should also #
# have unique names. However, we can  #
# define one single default target    #
# that contains the commons sources   #
# for all flavours of CoolProp to     #
# simplify integration.               #
#######################################

add_library(coolprop_object_library OBJECT ${COOLPROP_MAIN_SOURCES} ${COOLPROP_BACKEND_SOURCES} ${COOLPROP_MAIN_HEADERS} ${COOLPROP_BACKEND_HEADERS})
add_dependencies(coolprop_object_library generate_headers)
set(coolprop_object_includes $<TARGET_PROPERTY:coolprop_object_library,INCLUDE_DIRECTORIES>)


#######################################
#        ADD WRAPPER TARGETS          #
#-------------------------------------#
# Add the CMakeLists.txt files for    #
# wrappers. These subprojects may     #
# include additional source files and #
# dependencies not listed here.       #
#######################################

add_subdirectory(wrappers/library)





