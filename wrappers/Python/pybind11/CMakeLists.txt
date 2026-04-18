cmake_minimum_required(VERSION 3.14)

project(CoolProp)

set(CMAKE_CXX_STANDARD 17)

set(ROOT_DIR "../../..")

# Add core CoolProp dependency to get the list of files (also brings in CPM)
add_subdirectory("${ROOT_DIR}" "_CoolProp")

# ----------------------------
# Build pybind11 python module
# ----------------------------

CPMAddPackage(
  NAME pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.11.0
)
add_definitions(-DPYBIND11)
add_definitions(-DCOOLPROP_PYBIND11_MODULE)
include_directories(${COOLPROP_INCLUDE_DIRECTORIES})
pybind11_add_module(CoolProp ${COOLPROP_APP_SOURCES} "${ROOT_DIR}/src/pybind11_interface.cxx")