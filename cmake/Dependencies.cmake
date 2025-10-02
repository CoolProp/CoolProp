#######################################
#         REQUIRED MODULES            #
#-------------------------------------#
# CoolProp requires some standard OS  #
# features, these include:            #
# DL (CMAKE_DL_LIBS) for REFPROP      #
#######################################

# Add custom CMake module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake/Modules/")

# Find Python interpreter (needed for code generation)
message(STATUS "Looking for Python")
find_package(Python COMPONENTS Interpreter)
if(Python_Interpreter_FOUND)
  set(PYTHON_EXECUTABLE ${Python_EXECUTABLE})
endif()
if(NOT PYTHON_EXECUTABLE)
  message(WARNING "Could not find Python, be prepared for errors.")
endif()

# Find dynamic loading libraries if needed
if(CMAKE_DL_LIBS)
  find_package(${CMAKE_DL_LIBS} REQUIRED)
endif()

# Include flag manipulation functions
include(FlagFunctions) # Is found since it is in the module path.

# MSVC flag manipulation macros
macro(modify_msvc_flag_release flag_new) # Use a macro to avoid a new scope
  foreach(flag_old IN LISTS COOLPROP_MSVC_ALL)
    remove_compiler_flag_release("${flag_old} ") # add a space
    remove_compiler_flag_release(" ${flag_old}") # add a space
  endforeach()
  add_compiler_flag_release("${flag_new}")
endmacro()

macro(modify_msvc_flag_debug flag_new) # Use a macro to avoid a new scope
  foreach(flag_old IN LISTS COOLPROP_MSVC_ALL)
    remove_compiler_flag_debug("${flag_old} ") # add a space
    remove_compiler_flag_debug(" ${flag_old}") # add a space
  endforeach()
  add_compiler_flag_debug("${flag_new}")
endmacro()

macro(modify_msvc_flags flag_default) # Use a macro to avoid a new scope
  if(NOT "${COOLPROP_MSVC_REL}" STREQUAL "IGNORE")
    modify_msvc_flag_release("${COOLPROP_MSVC_REL}")
  else()
    modify_msvc_flag_release("${flag_default}")
  endif()
  if(NOT "${COOLPROP_MSVC_DBG}" STREQUAL "IGNORE")
    modify_msvc_flag_debug("${COOLPROP_MSVC_DBG}")
  else()
    modify_msvc_flag_debug("${flag_default}d")
  endif()
endmacro()

# Define COOLPROP_MSVC_ALL for the macros
set(COOLPROP_MSVC_ALL "/MTd" "/MT" "/MDd" "/MD"
    CACHE STRING "List of all MSVC runtime flags")

# Extract and decompress boost if needed
set(ZIPFN "${CMAKE_CURRENT_SOURCE_DIR}/dev/docker/boost_bcp_docker/boost_CoolProp.tar.xz")
set(OUTFN "${CMAKE_CURRENT_SOURCE_DIR}/boost_CoolProp/boost/version.hpp")
if(EXISTS ${ZIPFN})
  if(NOT EXISTS ${OUTFN})
    message(STATUS "Extracting boost...")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf ${ZIPFN}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  else()
    message(STATUS "Boost already extracted")
  endif()
endif()
