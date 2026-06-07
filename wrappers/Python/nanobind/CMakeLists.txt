cmake_minimum_required(VERSION 3.15)

project(CoolProp)

set(CMAKE_CXX_STANDARD 17)

set(ROOT_DIR "../../..")

# Add core CoolProp dependency to get the list of files (also brings in CPM)
add_subdirectory("${ROOT_DIR}" "_CoolProp")

# ----------------------------
# Build nanobind python module
# ----------------------------
# The binding source is src/nanobind_interface.cxx.  (It began as a 1:1
# translation of the former pybind11 interface, which has since been removed.)

# FindPython's Development.Module component requires CMake >= 3.18; fall back to
# the full Development component on 3.15-3.17 (the minimum above).
if(CMAKE_VERSION VERSION_LESS 3.18)
    set(_py_dev_component Development)
else()
    set(_py_dev_component Development.Module)
endif()
# Development.SABIModule (CPython >= 3.12, CMake >= 3.26) enables the stable ABI;
# requested optionally so the build still works on older Pythons.
find_package(Python 3.8 REQUIRED COMPONENTS Interpreter ${_py_dev_component}
             OPTIONAL_COMPONENTS Development.SABIModule)

# Locate the active interpreter's nanobind (provides nanobind-config.cmake)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_DIR)
find_package(nanobind CONFIG REQUIRED)

add_definitions(-DNANOBIND)
add_definitions(-DCOOLPROP_NANOBIND_MODULE)
include_directories(${COOLPROP_INCLUDE_DIRECTORIES})

# STABLE_ABI -> one CPython-version-independent wheel (CoolProp.abi3.so) on
# Python >= 3.12; nanobind silently falls back to a versioned extension below
# that.  This is the single-wheel-per-platform win over pybind11.
nanobind_add_module(CoolProp STABLE_ABI ${COOLPROP_APP_SOURCES} "${ROOT_DIR}/src/nanobind_interface.cxx")

# The module compiles the core sources directly, several of which include
# generated headers (gitrevision.h, cpversion.h, all_fluids_JSON*.h, ...).
# Depend on the core's generator target so a fresh checkout builds in one shot.
add_dependencies(CoolProp generate_headers)
