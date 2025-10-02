cmake_minimum_required(VERSION 3.12)

#######################################
#         CMAKE CONFIGURATION         #
#######################################

# Install prefix configuration
if(DEFINED COOLPROP_INSTALL_PREFIX)
  message(STATUS "COOLPROP_INSTALL_PREFIX=${COOLPROP_INSTALL_PREFIX}")
else()
  set(COOLPROP_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/install_root)
  message(STATUS "COOLPROP_INSTALL_PREFIX=${COOLPROP_INSTALL_PREFIX}")
endif()

set(CMAKE_INSTALL_PREFIX
    ${COOLPROP_INSTALL_PREFIX}
    CACHE PATH "default install path" FORCE)

#######################################
#         LOAD CONFIGURATION          #
#######################################

# Build options
include(cmake/Options.cmake)

# Compiler configuration
include(cmake/Compiler.cmake)

# Project version and metadata
include(cmake/Version.cmake)

# Project declaration
project(CoolProp VERSION ${COOLPROP_VERSION})
set(app_name "CoolProp")

#######################################
#         BUILD CONFIGURATION         #
#######################################

# Source files and include directories
include(cmake/Sources.cmake)

# Dependencies and external libraries
include(cmake/Dependencies.cmake)

# Core library build configuration
include(cmake/CoreLibrary.cmake)

#######################################
#         LANGUAGE WRAPPERS           #
#######################################

# All language wrapper modules (Python, Java, C#, R, Octave, etc.)
include(cmake/Wrappers.cmake)

#######################################
#         TESTING & EXAMPLES          #
#######################################

# Testing modules and examples
if(BUILD_TESTING OR COOLPROP_CATCH_MODULE OR COOLPROP_MAIN_MODULE OR
   COOLPROP_CPP_EXAMPLE_TEST OR COOLPROP_SNIPPETS OR
   COOLPROP_CLANG_ADDRESS_SANITIZER OR COOLPROP_PROFILE OR COOLPROP_COVERAGE)
  include(cmake/Testing.cmake)
endif()

#######################################
#         PROJECT COMPLETE            #
#######################################

message(STATUS "CoolProp CMake configuration complete")
