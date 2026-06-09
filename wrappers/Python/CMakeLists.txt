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
# FindPython's Development.Module component requires CMake >= 3.18; fall back to
# the full Development component on 3.15-3.17.
if(CMAKE_VERSION VERSION_LESS 3.18)
    set(_py_dev_component Development)
else()
    set(_py_dev_component Development.Module)
endif()
find_package(Python REQUIRED COMPONENTS Interpreter ${_py_dev_component}
             OPTIONAL_COMPONENTS Development.SABIModule)

# Find Cython
find_package(Cython REQUIRED)

# Opt-in: build the nanobind core + the capsule `State` shim instead of the
# legacy Cython interface, reusing this same harness (constants generation,
# package assembly, install rules, scikit-build wheel).  Default OFF -> the
# shipped (legacy Cython) wheel is completely unchanged.  When ON, the nanobind
# core (CoolProp.<abi3>.so) replaces the legacy Cython CoolProp module and the
# link-free State capsule shim is built as CoolProp.State.
#
# A `COOLPROP_NANOBIND=ON` *environment* build flips this on via a
# scikit-build-core override (pyproject.toml) that injects
# cmake.define.COOLPROP_NANOBIND=ON; the same env var (on Py>=3.12) also turns
# on the abi3 wheel tag.  Both signals flow through ONE scikit-build-core
# override block, so they cannot diverge.  `-DCOOLPROP_NANOBIND=ON` directly
# (or via cmake.define) still works for nanobind selection; abi3-ness of the
# Cython sidecars follows scikit-build-core's own SKBUILD_SABI_VERSION below, so
# a define-only build simply yields a valid per-version nanobind wheel.
option(COOLPROP_NANOBIND "Build the nanobind-based CoolProp package (core + State capsule shim)" OFF)

# Set root directory
set(ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")

# Use CPM.cmake for dependency management (consistent with the main build)
include("${ROOT_DIR}/cmake/CPM.cmake")
include("${ROOT_DIR}/cmake/dependencies.cmake")

list(APPEND CMAKE_MODULE_PATH "${ROOT_DIR}/cmake")
include(CoolPropJSONVisibility)

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
        "${ROOT_DIR}/include/all_fluids_CBOR.h"
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
    VERBATIM
)

# Create a custom target that depends on generated headers
add_custom_target(generate_headers_target
    DEPENDS
        "${ROOT_DIR}/include/all_fluids_CBOR.h"
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
        "${ROOT_DIR}/include/CoolProp/DataStructures.h"
        "${ROOT_DIR}/include/CoolProp/Configuration.h"
    COMMENT "Generating Cython constants module..."
    VERBATIM
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
    ${Eigen_SOURCE_DIR}
    ${fmt_SOURCE_DIR}/include
    ${boost_headers_SOURCE_DIR}
    ${ROOT_DIR}/externals/incbin
    ${nlohmann_json_SOURCE_DIR}/include
    ${valijson_SOURCE_DIR}/include
    ${ROOT_DIR}/externals/miniz-3.1.1
    ${msgpack-c_SOURCE_DIR}/include
    ${IF97_SOURCE_DIR}
    ${REFPROP_headers_SOURCE_DIR}
)

# Build miniz as a static library
add_library(miniz STATIC "${ROOT_DIR}/externals/miniz-3.1.1/miniz.c")
target_include_directories(miniz PUBLIC "${ROOT_DIR}/externals/miniz-3.1.1")

# Cythonize the .pyx files
set(CYTHON_FLAGS
    --cplus
    --directive embedsignature=True
    --directive language_level=3
    --directive c_string_type=unicode
    --directive c_string_encoding=ascii
)

# Generate _constants module (common to both the legacy and nanobind interfaces)
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
    VERBATIM
)

# Extra args injected into the Cython sidecars' Python_add_library() calls.
# Empty for the legacy build (sidecars stay per-version, WITH_SOABI).  The
# nanobind branch fills this with `USE_SABI <ver>` whenever scikit-build-core
# tagged the wheel abi3 (SKBUILD_SABI_VERSION), so the State + _constants
# extensions are limited-API and the whole wheel is a single cp312-abi3 artifact.
# See the COOLPROP_NANOBIND branch below.
set(_sidecar_sabi "")

if(NOT COOLPROP_NANOBIND)

# Generate CoolProp module
add_custom_command(
    OUTPUT CoolProp/CoolProp.cpp
    COMMAND ${Python_EXECUTABLE} -m cython ${CYTHON_FLAGS}
        -o "${CMAKE_CURRENT_BINARY_DIR}/CoolProp/CoolProp.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/CoolProp.pyx"
    DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/CoolProp.pyx"
        "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/constants_header.pxd"
        generate_constants_target
    COMMENT "Cythonizing CoolProp.pyx"
    VERBATIM
)

# Create CoolProp extension module
Python_add_library(CoolProp_module MODULE WITH_SOABI
    "${CMAKE_CURRENT_BINARY_DIR}/CoolProp/CoolProp.cpp"
    ${COOLPROP_SOURCES}
)

# Make sure generated headers and constants exist before compiling
add_dependencies(CoolProp_module generate_headers_target generate_constants_target)

target_include_directories(CoolProp_module PRIVATE
    ${COOLPROP_INCLUDE_DIRS}
    ${Python_INCLUDE_DIRS}
)

target_link_libraries(CoolProp_module PRIVATE miniz)
coolprop_hide_json_symbols(CoolProp_module)

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

else()  # COOLPROP_NANOBIND: nanobind core + capsule State shim

# The nanobind core uses STABLE_ABI (abi3) on Python >= 3.12 automatically; for
# the wheel to be a single cp312-abi3 artifact the two Cython sidecars (State,
# _constants) must ALSO be limited-API.  Drive that off scikit-build-core's OWN
# signal -- SKBUILD_SABI_VERSION, which it sets to e.g. "3.12" exactly when
# wheel.py-api tagged the wheel abi3, and leaves empty otherwise.  Following the
# wheel-tag decision (rather than re-deriving abi3 from the Python version or a
# second env read) guarantees the .so ABI and the wheel tag can never disagree:
# an abi3-tagged wheel always gets abi3 sidecars, a version-specific wheel always
# gets per-version sidecars.  So Py<3.12, a define-only nanobind build (no abi3
# tag), and a direct cmake build all correctly fall back to per-version sidecars.
# USE_SABI (FindPython) defines Py_LIMITED_API, links Python::SABIModule, and
# gives the module the .abi3 SOABI suffix.
if(DEFINED SKBUILD_SABI_VERSION AND SKBUILD_SABI_VERSION)
    if(NOT TARGET Python::SABIModule)
        message(FATAL_ERROR
            "scikit-build-core requested an abi3 wheel (SABI ${SKBUILD_SABI_VERSION}) "
            "but the Development.SABIModule component (CMake >= 3.26) is missing, so "
            "the Cython sidecars cannot be built limited-API; have CMake ${CMAKE_VERSION}.")
    endif()
    set(_sidecar_sabi USE_SABI ${SKBUILD_SABI_VERSION})
endif()

# ---- nanobind core (CoolProp.<abi3>.so) replaces the legacy Cython module ----
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_DIR)
find_package(nanobind CONFIG REQUIRED)
nanobind_add_module(CoolProp_module STABLE_ABI
    ${COOLPROP_SOURCES} "${ROOT_DIR}/src/nanobind_interface.cxx")
add_dependencies(CoolProp_module generate_headers_target)
target_compile_definitions(CoolProp_module PRIVATE NANOBIND COOLPROP_NANOBIND_MODULE)
target_include_directories(CoolProp_module PRIVATE
    ${COOLPROP_INCLUDE_DIRS}
    ${Python_INCLUDE_DIRS}
)
target_link_libraries(CoolProp_module PRIVATE miniz)
set_target_properties(CoolProp_module PROPERTIES OUTPUT_NAME "CoolProp" PREFIX "")
coolprop_hide_json_symbols(CoolProp_module)

# ---- PEP 561 type stub for the nanobind core (parity with the legacy wheel) ----
# The core .so is self-contained (no Python-level imports at init), so nanobind's
# stubgen can import it directly from the build dir to emit CoolProp.pyi.
nanobind_add_stub(
    CoolProp_stub
    MODULE CoolProp
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/CoolProp.pyi"
    PYTHON_PATH "$<TARGET_FILE_DIR:CoolProp_module>"
    DEPENDS CoolProp_module
    MARKER_FILE "${CMAKE_CURRENT_BINARY_DIR}/py.typed"
    # INCLUDE_PRIVATE emits underscore-named members (e.g. the _capi capsule);
    # the pattern file replaces the auto PropsSI/HAPropsSI overloads
    # (object-typed -> overlap) with typed ones.  (AbstractState is now a public
    # type with a factory constructor -- bd CoolProp-r9sq.28 -- not a private class.)
    INCLUDE_PRIVATE
    PATTERN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/_nanobind/stub_patterns.txt"
)

# ---- capsule State shim -> CoolProp.State (link-free; forwards through _capi) ----
# No CoolProp C++ is compiled in here; the shim only needs the C-ABI struct
# header and grabs the function table from the already-imported core at runtime.
add_custom_command(
    OUTPUT CoolProp/State.cpp
    COMMAND ${Python_EXECUTABLE} -m cython ${CYTHON_FLAGS}
        -o "${CMAKE_CURRENT_BINARY_DIR}/CoolProp/State.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/_nanobind/State.pyx"
    DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/_nanobind/State.pyx"
        "${CMAKE_CURRENT_SOURCE_DIR}/_nanobind/State.pxd"
    COMMENT "Cythonizing State.pyx (capsule shim)"
    VERBATIM
)
Python_add_library(State_module MODULE WITH_SOABI ${_sidecar_sabi}
    "${CMAKE_CURRENT_BINARY_DIR}/CoolProp/State.cpp"
)
target_include_directories(State_module PRIVATE
    "${ROOT_DIR}/include"
    ${Python_INCLUDE_DIRS}
)
set_target_properties(State_module PROPERTIES OUTPUT_NAME "State" PREFIX "")
if(WIN32)
    target_compile_options(State_module PRIVATE /utf-8 /std:c++17)
else()
    target_compile_options(State_module PRIVATE -std=c++17)
endif()

endif()  # COOLPROP_NANOBIND

# Create _constants extension module.  `_sidecar_sabi` is empty for the legacy
# build (per-version) and `USE_SABI <ver>` for the nanobind abi3 build.
Python_add_library(_constants_module MODULE WITH_SOABI ${_sidecar_sabi}
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

if(NOT COOLPROP_NANOBIND)

# Install the modules to the CoolProp package
install(TARGETS CoolProp_module _constants_module
    LIBRARY DESTINATION CoolProp
    RUNTIME DESTINATION CoolProp
)

# Install Python package files.
# *.pyi and the py.typed marker ship the PEP 561 type stubs into the wheel so
# editors/type-checkers see the typed API of the compiled extension; they must
# land next to CoolProp.<soabi>.so in the installed package.
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/"
    DESTINATION CoolProp
    FILES_MATCHING
    PATTERN "*.py"
    PATTERN "*.pyi"
    PATTERN "py.typed"
    PATTERN "*.pxd"
    PATTERN "*.bib"
    PATTERN "psyrc"
    PATTERN "__pycache__" EXCLUDE
    PATTERN "*.pyc" EXCLUDE
)

else()  # COOLPROP_NANOBIND

# nanobind core + capsule State shim + generated constants module
install(TARGETS CoolProp_module State_module _constants_module
    LIBRARY DESTINATION CoolProp
    RUNTIME DESTINATION CoolProp
)
# The legacy CoolProp/*.py and *.pxd are NOT nanobind-compatible, so the package
# files are installed explicitly (not by globbing the source dir):
#  - the nanobind __init__ (from __init__.nanobind.py)
#  - the capsule shim's cimport contract, shipped as State.pxd
#  - a package __init__.pxd so `from CoolProp cimport ...` resolves downstream
#  - the generated runtime constants + cimport header (PDSim consumes these)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/_nanobind/__init__.py"
    DESTINATION CoolProp)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/_nanobind/__init__.pxd"
    DESTINATION CoolProp)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/_nanobind/State.pxd"
    DESTINATION CoolProp)
install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/constants.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/constants_header.pxd"
    DESTINATION CoolProp)
# HumidAirProp re-export shim (SI surface from the core; non-SI HAProps removed)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/_nanobind/HumidAirProp.py"
    DESTINATION CoolProp)
# PEP 561 type stub for the nanobind core + py.typed marker (generated above)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/CoolProp.pyi" "${CMAKE_CURRENT_BINARY_DIR}/py.typed"
    DESTINATION CoolProp)
# Pure-Python submodules, shipped unchanged from the source package so the v8
# import tree matches the legacy one (CoolProp.Plots/BibtexParser/GUI/tests).
# tests/ is required: it is the q2sh "existing suite runs green" done-gate.
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/Plots"
    DESTINATION CoolProp FILES_MATCHING PATTERN "*.py" PATTERN "psyrc")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/GUI"
    DESTINATION CoolProp FILES_MATCHING PATTERN "*.py")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/tests"
    DESTINATION CoolProp FILES_MATCHING PATTERN "*.py")
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/BibtexParser.py"
    DESTINATION CoolProp)

endif()  # COOLPROP_NANOBIND

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
    REGEX "detail/json\\.h$" EXCLUDE
)

# Install fmtlib headers (for packaging).  fmt is CPM-fetched (fmt_SOURCE_DIR);
# the old externals/fmtlib path was stale, so no fmt shipped -- which broke
# downstream Cython that cimports constants_header (DataStructures.h ->
# numerics.h -> detail/strings.h -> fmt/format.h, header-only via FMT_HEADER_ONLY
# unless NO_FMTLIB).  Affected the legacy wheel too; this fixes both.
install(DIRECTORY "${fmt_SOURCE_DIR}/include/fmt"
    DESTINATION CoolProp/include
    FILES_MATCHING
    PATTERN "*.h"
)

# Install BibTeX library
install(FILES "${ROOT_DIR}/CoolPropBibTeXLibrary.bib"
    DESTINATION CoolProp
)
