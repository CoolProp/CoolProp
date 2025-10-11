cmake_minimum_required(VERSION 3.15...3.27)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set position independent code (required for shared libraries)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Add cmake directory to module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Find Python
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

# Find Cython
find_package(Cython REQUIRED)

# Set root directory
set(ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")

# ============================================================================
# Generate headers and constants (only when dependencies change)
# ============================================================================

# Collect all JSON files that trigger header regeneration
file(GLOB FLUID_JSON_FILES "${ROOT_DIR}/dev/fluids/*.json")
file(GLOB INCOMP_JSON_FILES "${ROOT_DIR}/dev/incompressible_liquids/json/*.json")
file(GLOB CUBIC_JSON_FILES "${ROOT_DIR}/dev/cubics/*.json")
file(GLOB PCSAFT_JSON_FILES "${ROOT_DIR}/dev/pcsaft/*.json")
file(GLOB MIXTURE_JSON_FILES "${ROOT_DIR}/dev/mixtures/*.json")

# Generate C++ headers from JSON files (only when inputs change)
# The generate_headers.py script has its own timestamp checking, so this
# won't unnecessarily regenerate files
add_custom_command(
    OUTPUT
        "${ROOT_DIR}/include/all_fluids_JSON.h"
        "${ROOT_DIR}/include/all_fluids_JSON_z.h"
        "${ROOT_DIR}/include/cpversion.h"
        "${ROOT_DIR}/include/gitrevision.h"
    COMMAND ${Python_EXECUTABLE} "${ROOT_DIR}/dev/generate_headers.py"
    WORKING_DIRECTORY "${ROOT_DIR}/dev"
    DEPENDS
        "${ROOT_DIR}/dev/generate_headers.py"
        "${ROOT_DIR}/CMakeLists.txt"
        ${FLUID_JSON_FILES}
        ${INCOMP_JSON_FILES}
        ${CUBIC_JSON_FILES}
        ${PCSAFT_JSON_FILES}
        ${MIXTURE_JSON_FILES}
    COMMENT "Generating CoolProp headers from JSON files..."
)

# Create a custom target that depends on generated headers
add_custom_target(generate_headers_target
    DEPENDS
        "${ROOT_DIR}/include/all_fluids_JSON.h"
        "${ROOT_DIR}/include/all_fluids_JSON_z.h"
        "${ROOT_DIR}/include/cpversion.h"
        "${ROOT_DIR}/include/gitrevision.h"
)

# Generate Cython constants module (only when dependencies change)
add_custom_command(
    OUTPUT
        "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/_constants.pyx"
        "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/constants_header.pxd"
        "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/constants.py"
    COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/generate_constants_module.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/generate_constants_module.py"
        "${ROOT_DIR}/include/DataStructures.h"
        "${ROOT_DIR}/include/Configuration.h"
    COMMENT "Generating Cython constants module..."
)

# Create a custom target for constants generation
add_custom_target(generate_constants_target
    DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/_constants.pyx"
        "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/constants_header.pxd"
        "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/constants.py"
)

# ============================================================================
# Headers are used directly from ROOT_DIR for building (better dependency tracking)
# They will be copied to the package during the install phase
# ============================================================================

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
    DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/_constants.pyx"
        "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/constants_header.pxd"
        generate_constants_target
    COMMENT "Cythonizing _constants.pyx"
)

# Create CoolProp extension module
Python_add_library(CoolProp_module MODULE WITH_SOABI
    "${CMAKE_CURRENT_BINARY_DIR}/CoolProp/CoolProp.cpp"
    ${COOLPROP_SOURCES}
)

# Make sure generated headers exist before compiling
add_dependencies(CoolProp_module generate_headers_target)

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

# Install header files directly from ROOT_DIR (for packaging)
# These are copied at install time, not during build, to avoid
# triggering unnecessary rebuilds
install(DIRECTORY "${ROOT_DIR}/include/"
    DESTINATION CoolProp/include
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp"
    PATTERN "*_JSON.h" EXCLUDE
    PATTERN "*_JSON_z.h" EXCLUDE
)

# Install fmtlib headers (for packaging)
install(DIRECTORY "${ROOT_DIR}/externals/fmtlib/include/fmt"
    DESTINATION CoolProp/include
    FILES_MATCHING
    PATTERN "*.h"
)

# Install BibTeX library
install(FILES "${ROOT_DIR}/CoolPropBibTeXLibrary.bib"
    DESTINATION CoolProp
)
