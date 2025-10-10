cmake_minimum_required(VERSION 3.15...3.27)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add cmake directory to module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Find Python
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

# Find Cython
find_package(Cython REQUIRED)

# Set root directory
set(ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")

# ============================================================================
# Generate headers and constants (must happen before building)
# ============================================================================

# Generate C++ headers from JSON files
message(STATUS "Generating CoolProp headers...")
execute_process(
    COMMAND ${Python_EXECUTABLE} "${ROOT_DIR}/dev/generate_headers.py"
    WORKING_DIRECTORY "${ROOT_DIR}/dev"
    RESULT_VARIABLE GENERATE_HEADERS_RESULT
)
if(NOT GENERATE_HEADERS_RESULT EQUAL 0)
    message(FATAL_ERROR "Failed to generate headers")
endif()

# Generate Cython constants module
message(STATUS "Generating Cython constants module...")
execute_process(
    COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/generate_constants_module.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    RESULT_VARIABLE GENERATE_CONSTANTS_RESULT
)
if(NOT GENERATE_CONSTANTS_RESULT EQUAL 0)
    message(FATAL_ERROR "Failed to generate constants module")
endif()

# ============================================================================
# Copy headers and other files needed for packaging
# ============================================================================

message(STATUS "Copying headers for packaging...")

# Remove old include directory if it exists
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/include")
    file(REMOVE_RECURSE "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/include")
endif()

# Copy main include directory
file(COPY "${ROOT_DIR}/include/"
     DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/include"
     PATTERN "*_JSON.h" EXCLUDE
     PATTERN "*_JSON_z.h" EXCLUDE
)

# Copy fmtlib headers
file(COPY "${ROOT_DIR}/externals/fmtlib/include/fmt"
     DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/include"
)

# Copy BibTeX library
file(COPY "${ROOT_DIR}/CoolPropBibTeXLibrary.bib"
     DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp"
)

# Add the main CoolProp sources
file(GLOB_RECURSE COOLPROP_SOURCES
    "${ROOT_DIR}/src/*.cpp"
)

# Remove any test files or unwanted sources
list(FILTER COOLPROP_SOURCES EXCLUDE REGEX ".*tests?/.*")

# Include directories
set(COOLPROP_INCLUDE_DIRS
    ${ROOT_DIR}
    ${ROOT_DIR}/include
    ${ROOT_DIR}/src
    ${ROOT_DIR}/dev
    ${ROOT_DIR}/externals/Eigen
    ${ROOT_DIR}/externals/fmtlib/include
    ${ROOT_DIR}/boost_CoolProp
    ${ROOT_DIR}/externals/incbin
    ${ROOT_DIR}/externals/nlohmann-json
    ${ROOT_DIR}/externals/miniz-3.0.2
    ${ROOT_DIR}/externals/msgpack-c/include
)

# Expand boost if needed
if(NOT EXISTS "${ROOT_DIR}/boost_CoolProp/boost/version.hpp")
    message(STATUS "Expanding boost sources...")
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E tar xzf "${ROOT_DIR}/dev/docker/boost_bcp_docker/boost_CoolProp.tar.xz"
        WORKING_DIRECTORY ${ROOT_DIR}
    )
endif()

# Build miniz as a static library
add_library(miniz STATIC "${ROOT_DIR}/externals/miniz-3.0.2/miniz.c")
target_include_directories(miniz PUBLIC "${ROOT_DIR}/externals/miniz-3.0.2")

# Cythonize the .pyx files
set(CYTHON_FLAGS
    --cplus
    --directive embedsignature=True
    --directive language_level=3
    --directive c_string_type=unicode
    --directive c_string_encoding=ascii
)

# Generate CoolProp module
add_custom_command(
    OUTPUT CoolProp/CoolProp.cpp
    COMMAND ${Python_EXECUTABLE} -m cython ${CYTHON_FLAGS}
        -o "${CMAKE_CURRENT_BINARY_DIR}/CoolProp/CoolProp.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/CoolProp.pyx"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/CoolProp.pyx"
    COMMENT "Cythonizing CoolProp.pyx"
)

# Generate _constants module
add_custom_command(
    OUTPUT CoolProp/_constants.cpp
    COMMAND ${Python_EXECUTABLE} -m cython ${CYTHON_FLAGS}
        -o "${CMAKE_CURRENT_BINARY_DIR}/CoolProp/_constants.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/_constants.pyx"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/_constants.pyx"
    COMMENT "Cythonizing _constants.pyx"
)

# Create CoolProp extension module
Python_add_library(CoolProp_module MODULE WITH_SOABI
    "${CMAKE_CURRENT_BINARY_DIR}/CoolProp/CoolProp.cpp"
    ${COOLPROP_SOURCES}
)

target_include_directories(CoolProp_module PRIVATE
    ${COOLPROP_INCLUDE_DIRS}
    ${Python_INCLUDE_DIRS}
)

target_link_libraries(CoolProp_module PRIVATE miniz)

# Set output name to CoolProp (not CoolProp_module)
set_target_properties(CoolProp_module PROPERTIES
    OUTPUT_NAME "CoolProp"
    PREFIX ""
)

if(WIN32)
    target_compile_options(CoolProp_module PRIVATE /utf-8 /std:c++17)
else()
    target_compile_options(CoolProp_module PRIVATE -std=c++17)
endif()

# Create _constants extension module
Python_add_library(_constants_module MODULE WITH_SOABI
    "${CMAKE_CURRENT_BINARY_DIR}/CoolProp/_constants.cpp"
)

target_include_directories(_constants_module PRIVATE
    ${COOLPROP_INCLUDE_DIRS}
    ${Python_INCLUDE_DIRS}
)

set_target_properties(_constants_module PROPERTIES
    OUTPUT_NAME "_constants"
    PREFIX ""
)

if(WIN32)
    target_compile_options(_constants_module PRIVATE /utf-8 /std:c++17)
else()
    target_compile_options(_constants_module PRIVATE -std=c++17)
endif()

# Install the modules to the CoolProp package
install(TARGETS CoolProp_module _constants_module
    LIBRARY DESTINATION CoolProp
    RUNTIME DESTINATION CoolProp
)

# Install Python package files
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/"
    DESTINATION CoolProp
    FILES_MATCHING
    PATTERN "*.py"
    PATTERN "*.pxd"
    PATTERN "*.bib"
    PATTERN "psyrc"
    PATTERN "__pycache__" EXCLUDE
    PATTERN "*.pyc" EXCLUDE
)

# Install include files if they exist
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/include")
    install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/include/"
        DESTINATION CoolProp/include
        FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
        PATTERN "*_JSON.h" EXCLUDE
        PATTERN "*_JSON_z.h" EXCLUDE
    )
endif()
