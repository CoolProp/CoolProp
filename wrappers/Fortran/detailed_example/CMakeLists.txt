cmake_minimum_required(VERSION 3.0)
enable_language(Fortran)

project(FORTRAN_INTERFACING)

# Add CoolProp dependency
set(COOLPROP_STATIC_LIBRARY true CACHE BOOL "Use CoolProp")
set(COOLPROP_EXTERNC_LIBRARY true CACHE BOOL "Use EXTERNC functions")
add_subdirectory("../../.." "CoolProp")

# Add main code
add_executable(Main "${CMAKE_CURRENT_SOURCE_DIR}/cpinterface.f90" "${CMAKE_CURRENT_SOURCE_DIR}/cool_fortran_bind.f90" )

# Link it to its dependent libraries
target_link_libraries (Main CoolProp)