cmake_minimum_required(VERSION 3.0)

project(CoolProp)

# Set the standard for C++ to c++11
set(CMAKE_CXX_STANDARD 11)

set (ROOT_DIR "../../..")

# Add core CoolProp dependency as a static library
set(COOLPROP_STATIC_LIBRARY true CACHE BOOL "Use CoolProp")
set(COOLPROP_LIBRARY_NAME "_CoolProp" CACHE PATH "The name of the CoolProp static library")
add_subdirectory("${ROOT_DIR}" "_CoolProp")

# ----------------------------
# Build pybind11 python module
# ----------------------------

add_subdirectory("${ROOT_DIR}/externals/pybind11" "pybind11")
add_definitions(-DPYBIND11)
add_definitions(-DCOOLPROP_PYBIND11_MODULE)
pybind11_add_module(CoolProp ${APP_SOURCES} "${ROOT_DIR}/src/pybind11_interface.cxx")
target_link_libraries (CoolProp PUBLIC _CoolProp)