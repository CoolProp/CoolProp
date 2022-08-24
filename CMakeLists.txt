cmake_minimum_required(VERSION 2.8.11)

include(CheckIncludeFileCXX)

if(DEFINED COOLPROP_INSTALL_PREFIX)
  #set(COOLPROP_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/install_root)
  message(STATUS "COOLPROP_INSTALL_PREFIX=${COOLPROP_INSTALL_PREFIX}")
else()
  set(COOLPROP_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/install_root)
  message(STATUS "COOLPROP_INSTALL_PREFIX=${COOLPROP_INSTALL_PREFIX}")
endif()

set(CMAKE_INSTALL_PREFIX
    ${COOLPROP_INSTALL_PREFIX}
    CACHE PATH "default install path" FORCE)

#######################################
#           BUILD OPTIONS             #
#-------------------------------------#
# These options are available to be   #
# modified in the build process.      #
# packages may want to modify these   #
# to suit, or just leave as defaults. #
#######################################

option(COOLPROP_STATIC_LIBRARY "Build CoolProp as a static library (.lib, .a)"
       OFF)

option(COOLPROP_SHARED_LIBRARY "Build CoolProp as a shared library (.dll, .so)"
       OFF)

option(COOLPROP_OBJECT_LIBRARY
       "Build CoolProp objects, but do not link them (.obj, .o)" OFF)

option(COOLPROP_EES_MODULE "Build the EES module" OFF)

option(COOLPROP_WINDOWS_PACKAGE "Build the Windows installer" OFF)

option(BUILD_TESTING "Enable testing for this given builder" OFF)

option(FORCE_BITNESS_32 "Force a 32bit build regardless of the host" OFF)

option(FORCE_BITNESS_64 "Force a 64bit build regardless of the host" OFF)

option(FORCE_BITNESS_NATIVE
       "Force a native bitness build regardless of the host" OFF)

option(COOLPROP_RELEASE "Optimize the builds with the release specs" OFF)

option(COOLPROP_DEBUG "Make a debug build" OFF)

option(COOLPROP_SMATH_WORK_INPLACE "Build SMath wrapper in source directory"
       OFF)

option(
  COOLPROP_MSVC_STATIC
  "Statically link Microsoft Standard library removes dependency on MSVCRXXX.dll."
  OFF)

option(
  COOLPROP_MSVC_DYNAMIC
  "Dynamically link Microsoft Standard library to integrate with other builds."
  OFF)

option(
  COOLPROP_MSVC_DEBUG
  "Link the debug version of Microsoft Standard library to the debug builds."
  ON)

option(COOLPROP_NO_EXAMPLES
       "Do not generate example code, does only apply to some wrappers." OFF)

#option (DARWIN_USE_LIBCPP
#        "On Darwin systems, compile and link with -std=libc++ instead of the default -std=libstdc++"
#        ON)

# Force C++11 since lambdas are used in CPStrings.h
set(CMAKE_CXX_STANDARD 11)

# see
# https://stackoverflow.com/questions/52509602/cant-compile-c-program-on-a-mac-after-upgrade-to-mojave
# https://support.enthought.com/hc/en-us/articles/204469410-OS-X-GCC-Clang-and-Cython-in-10-9-Mavericks
# https://github.com/pandas-dev/pandas/pull/24274/files
# https://github.com/explosion/thinc/pull/84/files
# https://github.com/jlfaucher/builder/commit/d144d3a695949f90c5e2acff4dfd94fdcf8dcdfa
# https://github.com/CoolProp/CoolProp/issues/1778
# https://gitlab.kitware.com/cmake/cmake/issues/18396
if(DEFINED DARWIN_USE_LIBCPP)
  if(DARWIN_USE_LIBCPP)
    set(CMAKE_OSX_DEPLOYMENT_TARGET
        "10.9"
        CACHE STRING "Minimum OS X deployment version")
    set(OSX_COMPILE_FLAGS "${OSX_COMPILE_FLAGS} -stdlib=libc++")
    set(OSX_COMPILE_FLAGS "${OSX_COMPILE_FLAGS} -mmacosx-version-min=10.9")
    set(OSX_COMPILE_FLAGS "${OSX_COMPILE_FLAGS} -std=c++11")
    set(OSX_LINK_FLAGS "${OSX_LINK_FLAGS} -lc++")
    set(OSX_LINK_FLAGS "${OSX_LINK_FLAGS} -nodefaultlibs")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
  else(DARWIN_USE_LIBCPP)
    set(CMAKE_OSX_DEPLOYMENT_TARGET
        "10.5"
        CACHE STRING "Minimum OS X deployment version") # Default is 10.7
    set(OSX_COMPILE_FLAGS "${OSX_COMPILE_FLAGS} -stdlib=libstdc++")
    set(OSX_COMPILE_FLAGS "${OSX_COMPILE_FLAGS} -mmacosx-version-min=10.5")
    set(OSX_LINK_FLAGS "${OSX_LINK_FLAGS} -lstdc++")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libstdc++")
  endif(DARWIN_USE_LIBCPP)
  message(STATUS "DARWIN_USE_LIBCPP was set added some flags:")
  message(STATUS "  OSX_COMPILE_FLAGS: ${OSX_COMPILE_FLAGS}")
  message(STATUS "     OSX_LINK_FLAGS: ${OSX_LINK_FLAGS}")
else(DEFINED DARWIN_USE_LIBCPP)
  if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
    message(STATUS "OSX build detected:")
    message(
      STATUS "  You might want to pass the -DDARWIN_USE_LIBCPP=ON/OFF parameter"
    )
    message(STATUS "  to enable or disable different C++ standard libraries.")
    message(
      STATUS
        "  You can also specify the environment variable MACOSX_DEPLOYMENT_TARGET=10.9 to force clang builds."
    )
  endif("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
endif(DEFINED DARWIN_USE_LIBCPP)

#if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}${OSX_COMPILE_FLAGS}")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}${OSX_COMPILE_FLAGS}")
#endif()
# Alternatively you could use
# set_target_properties(TARGET_NAME PROPERTIES APPEND_STRING PROPERTY COMPILE_FLAGS "-mmacosx-version-min=10.5")
# set_target_properties(TARGET_NAME PROPERTIES APPEND_STRING PROPERTY LINK_FLAGS "-mmacosx-version-min=10.5")

#######################################
#         PROJECT INFORMATION         #
#-------------------------------------#
# This CMakeLists.txt file is for the #
# CoolProp thermodynamic library      #
# written by Ian Bell. The following  #
# section contains project and        #
# version information.                #
#######################################

# Project name
set(project_name "CoolProp")
set(app_name ${project_name})
project(${project_name})

# Project version
set(COOLPROP_VERSION_MAJOR 6)
set(COOLPROP_VERSION_MINOR 4)
set(COOLPROP_VERSION_PATCH 2)
set(COOLPROP_VERSION_REVISION dev)
set(COOLPROP_VERSION
    "${COOLPROP_VERSION_MAJOR}.${COOLPROP_VERSION_MINOR}.${COOLPROP_VERSION_PATCH}${COOLPROP_VERSION_REVISION}"
)
message(STATUS "CoolProp version: ${COOLPROP_VERSION}")

string(TIMESTAMP COOLPROP_YEAR 2010-%Y)
#set ( COOLPROP_YEAR "2010-2016" )
set(COOLPROP_PUBLISHER "The CoolProp developers")

# Add definitions to silence warnings in MSVC2017 related to shared ptr code.
#if (MSVC_VERSION GREATER_EQUAL 1910 AND MSVC_VERSION LESS_EQUAL 1919) # OR MSVC_TOOLSET_VERSION EQUAL 141) # This requuires CMake >= 3.7
#  add_definitions(-D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
#endif (MSVC_VERSION GREATER_EQUAL 1910 AND MSVC_VERSION LESS_EQUAL 1919)

if(MSVC
   AND NOT (MSVC_VERSION LESS 1910)
   AND NOT (MSVC_VERSION GREATER 1919))
  add_definitions(-D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
endif()

if(COOLPROP_RELEASE AND COOLPROP_DEBUG)
  message(FATAL_ERROR "You can only make a release OR and debug build.")
endif()
if(COOLPROP_RELEASE)
  set(CMAKE_BUILD_TYPE Release)
elseif(COOLPROP_DEBUG)
  set(CMAKE_BUILD_TYPE Debug)
  #ELSEIF ("${CMAKE_BUILD_TYPE}" STREQUAL "")
  #  IF("${COOLPROP_VERSION_REVISION}" STREQUAL "dev")
  #    SET(CMAKE_BUILD_TYPE Debug)
  #  ELSE ()
  #    SET(CMAKE_BUILD_TYPE Release)
  #  ENDIF()
endif()

if(COOLPROP_MSVC_STATIC)
  set(COOLPROP_MSVC_REL "/MT")
  if(COOLPROP_MSVC_DEBUG)
    set(COOLPROP_MSVC_DBG "/MTd")
  else()
    set(COOLPROP_MSVC_DBG "/MT")
  endif()
elseif(COOLPROP_MSVC_DYNAMIC)
  set(COOLPROP_MSVC_REL "/MD")
  if(COOLPROP_MSVC_DEBUG)
    set(COOLPROP_MSVC_DBG "/MDd")
  else()
    set(COOLPROP_MSVC_DBG "/MD")
  endif()
else(COOLPROP_MSVC_DYNAMIC)
  set(COOLPROP_MSVC_REL "IGNORE")
  set(COOLPROP_MSVC_DBG "IGNORE")
endif()

set(COOLPROP_MSVC_ALL "/MTd" "/MT" "/MDd" "/MD"
)# order matters "/MXd" before "/MX"

#######################################
#         FIND ALL SOURCES            #
#-------------------------------------#
# The project is organised with       #
# split includes and source folders   #
# this makes it easier for developers #
# to quickly find relevant includes.  #
# This section finds all sources,     #
# headers and corresponding dirs.     #
#######################################

# These backends will be compiled in
set(COOLPROP_ENABLED_BACKENDS
    Cubics
    IF97
    Helmholtz
    REFPROP
    Incompressible
    Tabular
    PCSAFT)

# Get everything in the src/ directory (always), but not recursive
file(GLOB APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

# For each enabled backend, grab its files
foreach(backend ${COOLPROP_ENABLED_BACKENDS})
  file(GLOB_RECURSE BACKEND_SOURCES
       "${CMAKE_CURRENT_SOURCE_DIR}/src/Backends/${backend}/*.cpp")
  list(APPEND APP_SOURCES ${BACKEND_SOURCES})
endforeach()

## You can exclude this file, in case you want to run your own tests or use Catch
list(REMOVE_ITEM APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Tests/Tests.cpp")
list(REMOVE_ITEM APP_SOURCES
     "${CMAKE_CURRENT_SOURCE_DIR}/src/Tests/CoolProp-Tests.cpp")

## This file is only needed for the library, normal builds do not need it.
list(REMOVE_ITEM APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolPropLib.cpp")

set(APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/Eigen")
list(APPEND APP_INCLUDE_DIRS
     "${CMAKE_CURRENT_SOURCE_DIR}/externals/msgpack-c/include")
list(APPEND APP_INCLUDE_DIRS
     "${CMAKE_CURRENT_SOURCE_DIR}/externals/fmtlib/include")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/fmtlib"
)# should be deprecated
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")

include_directories(${APP_INCLUDE_DIRS})

set(SWIG_DEPENDENCIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/DataStructures.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/CoolProp.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/AbstractState.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Configuration.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/PhaseEnvelope.h)

set(COOLPROP_APP_SOURCES
    ${APP_SOURCES}
    CACHE LIST "List of CPP sources needed for CoolProp")
set(COOLPROP_INCLUDE_DIRECTORIES
    ${APP_INCLUDE_DIRS}
    CACHE LIST "List of include directories needed for CoolProp")

#######################################
#         REQUIRED MODULES            #
#-------------------------------------#
# CoolProp requires some standard OS  #
# features, these include:            #
# DL (CMAKE_DL_LIBS) for REFPROP      #
#######################################
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake/Modules/")

set(Python_ADDITIONAL_VERSIONS 2.7 2.6 2.5 2.4)
find_package(PythonInterp 2.7)
if(NOT PYTHON_EXECUTABLE)
  message(STATUS "Looking for Python")
  find_package(Python COMPONENTS Interpreter)
endif()
if(NOT PYTHON_EXECUTABLE)
  message(STATUS "Looking for Python2")
  find_package(Python2 COMPONENTS Interpreter)
  if(Python2_Interpreter_FOUND)
    set(PYTHON_EXECUTABLE ${Python2_EXECUTABLE})
  endif()
endif()
if(NOT PYTHON_EXECUTABLE)
  message(STATUS "Looking for Python3")
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
  endif()
endif()
if(NOT PYTHON_EXECUTABLE)
  message(WARNING "Could not find Python, be prepared for errors.")
endif()

if(CMAKE_DL_LIBS)
  find_package(${CMAKE_DL_LIBS} REQUIRED)
endif()

include(FlagFunctions) # Is found since it is in the module path.
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

#if(MSVC)
#  add_compiler_flag_release("/EHsc")
#  add_compiler_flag_debug("/EHsc")
#endif()

#######################################
#                 BITNESS             #
#-------------------------------------#
#        Calculate if 32 or 64        #
#######################################

if(WIN32)
  if(CMAKE_CL_64)
    set(BITNESS "64")
  else()
    set(BITNESS "32")
  endif()
else()
  if(CMAKE_SIZEOF_VOID_P MATCHES "8")
    set(BITNESS "64")
  else()
    set(BITNESS "32")
  endif()
endif()

if(MSVC AND (FORCE_BITNESS_32 OR FORCE_BITNESS_64))
  message(
    STATUS
      "You cannot force a certain bitness for Visual Studio, use the generator settings for this purpose."
  )
  message(
    STATUS
      "Pass '-G \"Visual Studio 10 2010 Win64\"' to CMake to make a 64bit binary using VS2010."
  )
  message(
    STATUS
      "Pass '-G \"Visual Studio 10 2010\"' to CMake to make a 32bit binary using VS2010."
  )
  message(
    STATUS
      "Pass '-G \"Visual Studio 9 2008 Win64\"' to CMake to make a 64bit binary using VS2008."
  )
  message(
    STATUS
      "Pass '-G \"Visual Studio 9 2008\"' to CMake to make a 32bit binary using VS2008."
  )
  message(FATAL_ERROR "Fix that and try again...")
endif()

if(FORCE_BITNESS_32)
  set(BITNESS "32")
elseif(FORCE_BITNESS_64)
  set(BITNESS "64")
elseif(FORCE_BITNESS_NATIVE)
  set(BITNESS "NATIVE")
endif()

#######################################
#         SHARED POINTER              #
#-------------------------------------#
# In this section we define the       #
# flags needed to use shared_ptr      #
# reliably                            #
#######################################

include("${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake/Modules/FindSharedPtr.cmake")
find_shared_ptr()
if(NOT SHARED_PTR_FOUND)
  message(FATAL_ERROR "Must be able to find shared_ptr")
else()
  if(SHARED_PTR_TR1_MEMORY_HEADER)
    add_definitions("-DSHARED_PTR_TR1_MEMORY_HEADER")
  endif()
  if(SHARED_PTR_TR1_NAMESPACE)
    add_definitions("-DSHARED_PTR_TR1_NAMESPACE")
  endif()
endif()

#######################################
#         MAKE ARTEFACTS              #
#-------------------------------------#
# In this section we define the       #
# artefacts (exes, libs) that will be #
# made for CoolProp, these include    #
# customisation from earlier options. #
#######################################

###     FLUIDS, MIXTURES JSON       ###
add_custom_target(
  generate_headers
  COMMAND "${PYTHON_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/dev/generate_headers.py")

if(NOT COOLPROP_NO_EXAMPLES)
  add_custom_target(
    generate_examples
    COMMAND "${PYTHON_EXECUTABLE}" example_generator.py Python
            "${CMAKE_CURRENT_BINARY_DIR}/Example.py"
    COMMAND "${PYTHON_EXECUTABLE}" example_generator.py Octave
            "${CMAKE_CURRENT_BINARY_DIR}/Example.m"
    COMMAND "${PYTHON_EXECUTABLE}" example_generator.py R
            "${CMAKE_CURRENT_BINARY_DIR}/Example.R"
    #COMMAND "${PYTHON_EXECUTABLE}" example_generator.py MATLAB "${CMAKE_CURRENT_BINARY_DIR}/Example.m"
    COMMAND "${PYTHON_EXECUTABLE}" example_generator.py Java
            "${CMAKE_CURRENT_BINARY_DIR}/Example.java"
    COMMAND "${PYTHON_EXECUTABLE}" example_generator.py Csharp
            "${CMAKE_CURRENT_BINARY_DIR}/Example.cs"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dev/scripts/examples")
else()
  add_custom_target(
    generate_examples
    COMMAND
      echo
      "Example generation has been disabled with the COOLPROP_NO_EXAMPLES option."
  )
endif()

###      Library options            ###
# We already know the bitness from the earlier
# settings. Let us rely on that and only handle
# calling conventions and shared/static issues.

option(COOLPROP_STDCALL_LIBRARY
       "Build CoolProp as a 32bit shared library with stdcall" OFF)

option(COOLPROP_CDECL_LIBRARY
       "Build CoolProp as a 32bit shared library with cdecl" OFF)

option(COOLPROP_EXTERNC_LIBRARY
       "Overwrite the export settings to force extern C" OFF)

set(COOLPROP_LIBRARY_SOURCE
    "src/CoolPropLib.cpp"
    CACHE STRING "The file that contains the exported functions")

set(COOLPROP_LIBRARY_HEADER
    "include/CoolPropLib.h"
    CACHE STRING "The file that contains the export header")

set(COOLPROP_LIBRARY_NAME
    "${app_name}"
    CACHE STRING "The name of the generated library")

set(COOLPROP_LIBRARY_EXPORTS
    ""
    CACHE STRING "The file that contains the export alias list")

# Rule out cases that do not make sense
if("${BITNESS}" STREQUAL "32")
  if(COOLPROP_CDECL_LIBRARY)
    set(CONVENTION "__cdecl")
  elseif(COOLPROP_STDCALL_LIBRARY)
    set(CONVENTION "__stdcall")
  else()
    set(CONVENTION "")
  endif()
elseif("${BITNESS}" STREQUAL "64")
  if(COOLPROP_CDECL_LIBRARY)
    message(WARNING "You cannot use cdecl conventions in a 64-bit library.")
  elseif(COOLPROP_STDCALL_LIBRARY)
    message(WARNING "You cannot use stdcall conventions in a 64-bit library.")
  endif()
  set(CONVENTION "")
elseif("${BITNESS}" STREQUAL "NATIVE")
  set(CONVENTION "")
else()
  message(FATAL_ERROR "Bitness is not defined. Set it and run cmake again.")
endif()

if((COOLPROP_OBJECT_LIBRARY AND COOLPROP_STATIC_LIBRARY)
   OR (COOLPROP_OBJECT_LIBRARY AND COOLPROP_SHARED_LIBRARY)
   OR (COOLPROP_STATIC_LIBRARY AND COOLPROP_SHARED_LIBRARY))
  message(FATAL_ERROR "You can only use one of the library switches!")
endif()

if(COOLPROP_OBJECT_LIBRARY
   OR COOLPROP_STATIC_LIBRARY
   OR COOLPROP_SHARED_LIBRARY)
  # Project name
  set(LIB_NAME ${COOLPROP_LIBRARY_NAME})
  # Object, static or shared?
  if(COOLPROP_OBJECT_LIBRARY)
    add_library(${LIB_NAME} OBJECT ${APP_SOURCES})
    set(COOLPROP_LIBRARY_SOURCE "")
    set(COOLPROP_LIBRARY_HEADER "")
  elseif(COOLPROP_STATIC_LIBRARY)
    list(APPEND APP_SOURCES
         "${CMAKE_CURRENT_SOURCE_DIR}/${COOLPROP_LIBRARY_SOURCE}")
    add_library(${LIB_NAME} STATIC ${APP_SOURCES} ${COOLPROP_LIBRARY_EXPORTS})
    if(MSVC)
      # Add postfix for debugging
      set_property(TARGET ${LIB_NAME} PROPERTY DEBUG_POSTFIX d)
      set_property(TARGET ${LIB_NAME} PROPERTY RELEASE_POSTFIX)
      modify_msvc_flags("/MD")

      # Note that the default is not used if ${COOLPROP_MSVC_REL} or ${COOLPROP_MSVC_DBG} is set
    endif(MSVC)
    install(
      TARGETS ${LIB_NAME}
      DESTINATION
        static_library/${CMAKE_SYSTEM_NAME}/${BITNESS}bit_${CMAKE_CXX_COMPILER_ID}_${CMAKE_CXX_COMPILER_VERSION}
    )
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${COOLPROP_LIBRARY_HEADER}
            DESTINATION static_library)
  elseif(COOLPROP_SHARED_LIBRARY)
    list(APPEND APP_SOURCES
         "${CMAKE_CURRENT_SOURCE_DIR}/${COOLPROP_LIBRARY_SOURCE}")
    add_library(${LIB_NAME} SHARED ${APP_SOURCES} ${COOLPROP_LIBRARY_EXPORTS})
    install(
      TARGETS ${LIB_NAME}
      DESTINATION shared_library/${CMAKE_SYSTEM_NAME}/${BITNESS}bit${CONVENTION}
    )
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${COOLPROP_LIBRARY_HEADER}
            DESTINATION shared_library)
    set_property(
      TARGET ${LIB_NAME}
      APPEND_STRING
      PROPERTY COMPILE_FLAGS " -DCOOLPROP_LIB")
    # Now all the compiler specific settings for Visual Studio
    if(MSVC)
      # Add postfix for debugging
      set_property(TARGET ${LIB_NAME} PROPERTY DEBUG_POSTFIX d)
      set_property(TARGET ${LIB_NAME} PROPERTY RELEASE_POSTFIX)
      # No lib prefix for the shared library
      set_property(TARGET ${LIB_NAME} PROPERTY PREFIX "")
      modify_msvc_flags("/MT")

      # Note that the default is not used if ${COOLPROP_MSVC_REL} or ${COOLPROP_MSVC_DBG} is set
      add_custom_command(
        TARGET ${LIB_NAME}
        POST_BUILD
        COMMAND dumpbin /EXPORTS $<TARGET_FILE:${LIB_NAME}> >
                ${CMAKE_CURRENT_BINARY_DIR}/exports.txt)
      install(
        FILES ${CMAKE_CURRENT_BINARY_DIR}/exports.txt
        DESTINATION
          shared_library/${CMAKE_SYSTEM_NAME}/${BITNESS}bit${CONVENTION})
    endif()
    # For Linux
    if("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
      set_property(TARGET ${LIB_NAME} PROPERTY VERSION ${COOLPROP_VERSION})
      set_property(TARGET ${LIB_NAME} PROPERTY SOVERSION
                                               ${COOLPROP_VERSION_MAJOR})
    endif()
  else()
    message(FATAL_ERROR "You have to build a static or shared library.")
  endif()

  if(NOT COOLPROP_OBJECT_LIBRARY)
    target_link_libraries(${LIB_NAME} ${CMAKE_DL_LIBS})
  endif()

  # For windows systems, bug workaround for Eigen
  if(MSVC90)
    message(STATUS "EIGEN WORKAROUND ACTIVE!!")
    set_property(
      TARGET ${LIB_NAME}
      APPEND_STRING
      PROPERTY COMPILE_FLAGS " -DEIGEN_DONT_VECTORIZE")
  endif()

  # For mac systems, explicitly set the c++ libraries
  if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
    if(DEFINED OSX_COMPILE_FLAGS)
      set_target_properties(
        ${LIB_NAME} PROPERTIES APPEND_STRING PROPERTY COMPILE_FLAGS
                                                      "${OSX_COMPILE_FLAGS}")
    endif(DEFINED OSX_COMPILE_FLAGS)
    if(DEFINED OSX_COMPILE_FLAGS)
      set_target_properties(
        ${LIB_NAME} PROPERTIES APPEND_STRING PROPERTY LINK_FLAGS
                                                      "${OSX_LINK_FLAGS}")
    endif(DEFINED OSX_COMPILE_FLAGS)
  endif()

  # Name mangling settings
  if(COOLPROP_EXTERNC_LIBRARY)
    set_property(
      TARGET ${LIB_NAME}
      APPEND_STRING
      PROPERTY COMPILE_FLAGS " -DEXTERNC")
  endif()
  ### All options are set, we are building a library ###
  add_dependencies(${LIB_NAME} generate_headers)
  #
  if(CMAKE_VERSION VERSION_GREATER 3.0)
    # Add target include directories for easy linking with other applications
    target_include_directories(${LIB_NAME} PUBLIC ${APP_INCLUDE_DIRS})
  endif()

  # Set the bitness
  if(NOT MSVC)
    if(NOT "${BITNESS}" STREQUAL "NATIVE")
      message(STATUS "Setting bitness flag -m${BITNESS}")
      set_property(
        TARGET ${LIB_NAME}
        APPEND_STRING
        PROPERTY COMPILE_FLAGS " -m${BITNESS}")
      set_property(
        TARGET ${LIB_NAME}
        APPEND_STRING
        PROPERTY LINK_FLAGS " -m${BITNESS}")
    endif()
  endif()

  # ADD -fPIC flag if needed
  if(COOLPROP_FPIC)
    message(STATUS "Setting bitness flag -m${BITNESS}")
    set_property(
      TARGET ${LIB_NAME}
      APPEND_STRING
      PROPERTY COMPILE_FLAGS " -fPIC")
  endif()
  #
  # calling conventions
  if(("${CONVENTION}" STREQUAL "NULL") OR ("${CONVENTION}" STREQUAL ""))
    #MESSAGE(STATUS "Skipping unknown calling convention.")
  else()
    set_property(
      TARGET ${LIB_NAME}
      APPEND_STRING
      PROPERTY COMPILE_FLAGS " -DCONVENTION=${CONVENTION}")
  endif()
  #
  #set_property(SOURCE MyFile.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " -msse4.1 ")
  message(STATUS "Library compilation detected:")
  message(STATUS "Creating ${LIB_NAME}, a ${BITNESS}-bit library")
  message(STATUS "CMAKE_SYSTEM_NAME: ${CMAKE_SYSTEM_NAME}")
  message(STATUS "COOLPROP_STATIC_LIBRARY: ${COOLPROP_STATIC_LIBRARY}")
  message(STATUS "COOLPROP_SHARED_LIBRARY: ${COOLPROP_SHARED_LIBRARY}")
  message(STATUS "COOLPROP_OBJECT_LIBRARY: ${COOLPROP_OBJECT_LIBRARY}")
  message(STATUS "CONVENTION: ${CONVENTION}")
  message(STATUS "COOLPROP_LIBRARY_HEADER: ${COOLPROP_LIBRARY_HEADER}")
  message(STATUS "COOLPROP_LIBRARY_SOURCE: ${COOLPROP_LIBRARY_SOURCE}")
  #
  get_property(
    tmpVar
    TARGET ${LIB_NAME}
    PROPERTY COMPILE_FLAGS)
  message(STATUS "COMPILE_FLAGS: ${tmpVar}")
  get_property(
    tmpVar
    TARGET ${LIB_NAME}
    PROPERTY LINK_FLAGS)
  message(STATUS "LINK_FLAGS: ${tmpVar}")
  #
endif()

if(COOLPROP_IOS_TARGET)
  # Set the Base SDK (only change the SDKVER value, if for instance, you are building for iOS 5.0):
  set(SDKVER "9.2")
  set(DEVROOT
      "/Applications/Xcode.app/Contents/Developer/Platforms/iPhoneOS.platform/Developer"
  )
  set(SDKROOT "${DEVROOT}/SDKs/iPhoneOS${SDKVER}.sdk")
  if(EXISTS ${SDKROOT})
    set(CMAKE_OSX_SYSROOT "${SDKROOT}")
  else()
    message("Warning, iOS Base SDK path not found: " ${SDKROOT})
  endif()

  # Will resolve to "Standard (armv6 armv7)" on Xcode 4.0.2 and to "Standard (armv7)" on Xcode 4.2:
  set(CMAKE_OSX_ARCHITECTURES "$(ARCHS_STANDARD_32_BIT)")

  # seamless toggle between device and simulator
  set(CMAKE_XCODE_EFFECTIVE_PLATFORMS "-iphoneos;-iphonesimulator")

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
endif()

if(COOLPROP_DEBIAN_PACKAGE)
  if(NOT UNIX)
    message(
      FATAL_ERROR "COOLPROP_DEBIAN_PACKAGE can only be used on linux host")
  endif()
  list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolPropLib.cpp")
  add_library(${app_name} SHARED ${APP_SOURCES})
  set_target_properties(
    ${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB")
  set_target_properties(
    ${app_name} PROPERTIES VERSION ${COOLPROP_VERSION}
                           SOVERSION ${COOLPROP_VERSION_MAJOR})
  add_dependencies(${app_name} generate_headers)
  install(TARGETS ${app_name} DESTINATION "${CMAKE_INSTALL_PREFIX}/usr/lib")
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/CoolPropLib.h
          DESTINATION "${CMAKE_INSTALL_PREFIX}/usr/include")
endif()

if(COOLPROP_VXWORKS_MAKEFILE)

  set(INCLUDE_DIRECTORIES)
  foreach(_srcFile ${APP_INCLUDE_DIRS})
    string(CONCAT _el "-I\"" ${_srcFile} "\"")
    string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}" "$(COOLPROP_ROOT)" _el
                   "${_el}")
    list(APPEND INCLUDE_DIRECTORIES ${_el})
  endforeach()
  string(REPLACE ";" " " INCLUDE_DIRECTORIES "${INCLUDE_DIRECTORIES}")
  set(OLD_ROOT /home/ian/.wine/drive_c/)
  set(NEW_ROOT c:/)
  string(REPLACE ${OLD_ROOT} ${NEW_ROOT} INCLUDE_DIRECTORIES
                 "${INCLUDE_DIRECTORIES}")
  set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/src")
  string(REPLACE ${OLD_ROOT} ${NEW_ROOT} SRC "${SRC}")
  file(RELATIVE_PATH COOLPROP_ROOT "${CMAKE_CURRENT_BINARY_DIR}"
       "${CMAKE_CURRENT_SOURCE_DIR}")

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/Labview/vxWorks/Makefile.in"
    "vxWorksMakefile")
endif()

if(COOLPROP_VXWORKS_LIBRARY_MODULE OR COOLPROP_VXWORKS_LIBRARY)
  list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolPropLib.cpp")
  add_executable(${app_name} ${APP_SOURCES})
  set_target_properties(
    ${app_name} PROPERTIES SUFFIX ".out" COMPILE_FLAGS
                                         "${COMPILE_FLAGS} -DEXTERNC")
  add_dependencies(${app_name} generate_headers)
  install(TARGETS ${app_name}
          DESTINATION "${COOLPROP_INSTALL_PREFIX}/shared_library/VxWorks")
endif()

if(COOLPROP_PRIME_MODULE)
  if(NOT WIN32)
    message(
      FATAL_ERROR "COOLPROP_PRIME_MODULE can only be used on windows host")
  endif()
  if("${COOLPROP_PRIME_ROOT}" STREQUAL "")
    message(
      FATAL_ERROR
        "You must provide the path to MathCAD Prime Root directory using something like -DCOOLPROP_PRIME_ROOT=\"C:/Program Files/PTC/Mathcad Prime 3.1\""
    )
  else()
    message(STATUS "COOLPROP_PRIME_ROOT: ${COOLPROP_PRIME_ROOT}")
  endif()
  list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolPropLib.cpp")
  list(APPEND APP_SOURCES
       "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/MathCAD/CoolPropMathcad.cpp")
  add_library(CoolPropMathcadWrapper SHARED ${APP_SOURCES})
  include_directories("${COOLPROP_PRIME_ROOT}/Custom Functions")
  target_link_libraries(CoolPropMathcadWrapper
                        "${COOLPROP_PRIME_ROOT}/Custom Functions/mcaduser.lib")
  set_target_properties(CoolPropMathcadWrapper
                        PROPERTIES LINK_FLAGS "/ENTRY:\"DllEntryPoint\"")
  add_dependencies(CoolPropMathcadWrapper generate_headers)
  set_target_properties(CoolPropMathcadWrapper PROPERTIES SUFFIX ".dll" PREFIX
                                                                        "")
  install(TARGETS CoolPropMathcadWrapper
          DESTINATION "${COOLPROP_INSTALL_PREFIX}/MathCADPrime")
  install(
    FILES
      "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/MathCAD/Prime/CoolPropFluidProperties.mcdx"
    DESTINATION MathCADPrime)
endif()

if(COOLPROP_MATHCAD15_MODULE)
  if(NOT WIN32)
    message(
      FATAL_ERROR "COOLPROP_MATHCAD15_MODULE can only be used on windows host")
  endif()
  if("${COOLPROP_MATHCAD15_ROOT}" STREQUAL "")
    message(
      FATAL_ERROR
        "You must provide the path to MathCAD 15 Root directory using something like -DCOOLPROP_MATHCAD15_ROOT=\"C:/Program Files (x86)/Mathcad/Mathcad 15\""
    )
  else()
    message(STATUS "COOLPROP_MATHCAD15_ROOT: ${COOLPROP_MATHCAD15_ROOT}")
  endif()
  list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolPropLib.cpp")
  list(APPEND APP_SOURCES
       "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/MathCAD/CoolPropMathcad.cpp")
  add_library(CoolPropMathcadWrapper SHARED ${APP_SOURCES})
  include_directories("${COOLPROP_MATHCAD15_ROOT}/userefi/microsft/include")
  target_link_libraries(
    CoolPropMathcadWrapper
    "${COOLPROP_MATHCAD15_ROOT}/userefi/microsft/lib/mcaduser.lib")
  set_target_properties(CoolPropMathcadWrapper
                        PROPERTIES LINK_FLAGS "/ENTRY:\"DllEntryPoint\"")
  add_dependencies(CoolPropMathcadWrapper generate_headers)
  set_target_properties(CoolPropMathcadWrapper PROPERTIES SUFFIX ".dll" PREFIX
                                                                        "")
  install(TARGETS CoolPropMathcadWrapper
          DESTINATION "${COOLPROP_INSTALL_PREFIX}/MathCAD15")
  install(
    FILES
      "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/MathCAD/CoolPropFluidProperties.xmcdz"
    DESTINATION MathCAD15)
  install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/MathCAD/CoolProp_EN.xml"
          DESTINATION MathCAD15)
endif()

# EES is only compiled for 32bit Windows
if(COOLPROP_EES_MODULE)
  if(NOT "${BITNESS}" STREQUAL "32")
    message(FATAL_ERROR "You cannot build the EES wrapper as a 64-bit library.")
  endif()
  # Prepare the sources
  include_directories(${APP_INCLUDE_DIRS})
  list(APPEND APP_SOURCES "wrappers/EES/main.cpp")
  list(APPEND APP_SOURCES
       "${CMAKE_CURRENT_SOURCE_DIR}/${COOLPROP_LIBRARY_SOURCE}")
  add_library(COOLPROP_EES SHARED ${APP_SOURCES})
  # Modify the target and add dependencies
  add_dependencies(COOLPROP_EES generate_headers)
  set_target_properties(
    COOLPROP_EES
    PROPERTIES COMPILE_FLAGS
               "${COMPILE_FLAGS} -DCOOLPROP_LIB -DCONVENTION=__cdecl")
  set_target_properties(COOLPROP_EES PROPERTIES SUFFIX ".dlf" PREFIX "")
  # Creates "COOLPROP_EES.dlf"
  if(NOT MSVC)
    set_target_properties(COOLPROP_EES PROPERTIES COMPILE_FLAGS "-m32"
                                                  LINK_FLAGS "-m32")
  elseif(MSVC)
    set_target_properties(COOLPROP_EES PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                  ${CMAKE_CURRENT_BINARY_DIR})
    set_target_properties(COOLPROP_EES PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
                                                  ${CMAKE_CURRENT_BINARY_DIR})
    set_target_properties(
      COOLPROP_EES PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
                              ${CMAKE_CURRENT_BINARY_DIR})
    # etc for the other available configuration types (MinSizeRel, RelWithDebInfo)
  endif()
  # copy required files
  add_custom_command(
    TARGET COOLPROP_EES
    PRE_BUILD
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy"
      "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/EES/CoolProp.htm"
      "${CMAKE_CURRENT_BINARY_DIR}/."
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy"
      "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/EES/CoolProp.LIB"
      "${CMAKE_CURRENT_BINARY_DIR}/."
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy"
      "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/EES/CoolProp_EES_Sample.EES"
      "${CMAKE_CURRENT_BINARY_DIR}/."
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Copying the EES files to the build directory"
    VERBATIM)
  # install the generated library and the other files
  install(TARGETS COOLPROP_EES
          DESTINATION "${CMAKE_INSTALL_PREFIX}/EES/${CMAKE_SYSTEM_NAME}")
  install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/EES/CoolProp.htm"
          DESTINATION "${CMAKE_INSTALL_PREFIX}/EES/${CMAKE_SYSTEM_NAME}")
  install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/EES/CoolProp.LIB"
          DESTINATION "${CMAKE_INSTALL_PREFIX}/EES/${CMAKE_SYSTEM_NAME}")
  install(
    FILES "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/EES/CoolProp_EES_Sample.EES"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/EES/${CMAKE_SYSTEM_NAME}")
endif()

# Windows package
if(COOLPROP_WINDOWS_PACKAGE)

  message(
    STATUS "Creating Windows installer for COOLPROP_VERSION=${COOLPROP_VERSION}"
  )
  # Setting some basic build paths
  set(COOLPROP_WINDOWS_PACKAGE_32B_DIR "${CMAKE_CURRENT_BINARY_DIR}/32bitDLL")
  set(COOLPROP_WINDOWS_PACKAGE_32B_DIR_STDCALL
      "${CMAKE_CURRENT_BINARY_DIR}/32bitDLL_stdcall")
  set(COOLPROP_WINDOWS_PACKAGE_32B_DIR_CDECL
      "${CMAKE_CURRENT_BINARY_DIR}/32bitDLL_cdecl")
  set(COOLPROP_WINDOWS_PACKAGE_64B_DIR "${CMAKE_CURRENT_BINARY_DIR}/64bitDLL")
  set(COOLPROP_WINDOWS_PACKAGE_EES_DIR "${CMAKE_CURRENT_BINARY_DIR}/EES")
  set(COOLPROP_WINDOWS_PACKAGE_TMP_DIR "${CMAKE_CURRENT_BINARY_DIR}/InnoScript")
  # Pointers to the sources
  set(COOLPROP_WINDOWS_PACKAGE_EXCEL_DIR
      "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/Excel")
  set(COOLPROP_WINDOWS_PACKAGE_ISS_DIR
      "${CMAKE_CURRENT_SOURCE_DIR}/externals/ExcelAddinInstaller")
  # Generator for DLLs
  set(COOLPROP_WINDOWS_PACKAGE_DLL_GEN "${CMAKE_GENERATOR}"
  )# Use the currently selected generator, architecture is hard-coded below
  # Configure variables like version number and build year
  configure_file(
    "${COOLPROP_WINDOWS_PACKAGE_ISS_DIR}/cmake-templates/config.iss"
    "${COOLPROP_WINDOWS_PACKAGE_ISS_DIR}/config.iss")
  # Find the installer generator executable
  set(BINDIR32_ENV_NAME "ProgramFiles(x86)")
  set(BINDIR32 $ENV{${BINDIR32_ENV_NAME}})
  set(BINDIR64_ENV_NAME "ProgramFiles")
  set(BINDIR64 $ENV{${BINDIR64_ENV_NAME}})
  find_program(
    COOLPROP_WINDOWS_PACKAGE_ISS_EXE
    NAMES iscc.exe
    HINTS "${BINDIR32}/Inno Setup 6" "${BINDIR64}/Inno Setup 6")

  # ******************************************************************
  # Add the targets that prepare the build directory for the subbuilds
  # ******************************************************************
  add_custom_target(COOLPROP_WINDOWS_PACKAGE_PREPARE)
  # Prepare directories
  add_custom_command(
    TARGET COOLPROP_WINDOWS_PACKAGE_PREPARE
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "make_directory"
            "${COOLPROP_WINDOWS_PACKAGE_32B_DIR}"
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "make_directory"
            "${COOLPROP_WINDOWS_PACKAGE_32B_DIR_STDCALL}"
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "make_directory"
            "${COOLPROP_WINDOWS_PACKAGE_32B_DIR_CDECL}"
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "make_directory"
            "${COOLPROP_WINDOWS_PACKAGE_64B_DIR}"
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "make_directory"
            "${COOLPROP_WINDOWS_PACKAGE_EES_DIR}"
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "make_directory"
            "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}"
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "make_directory"
            "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/source"
    #COMMAND ${CMAKE_COMMAND} ARGS "-E" "remove_directory" "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/deploy"
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "make_directory"
            "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/deploy"
    #COMMAND ${CMAKE_COMMAND} ARGS "-E" "remove_directory" "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/bin"
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "make_directory"
            "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/bin"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Preparing the directories for the Windows installer"
    VERBATIM)

  add_custom_target(COOLPROP_WINDOWS_PACKAGE_DELETE)
  # Delete directories
  add_custom_command(
    TARGET COOLPROP_WINDOWS_PACKAGE_DELETE
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "make_directory"
            "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/source"
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "remove_directory"
            "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/deploy"
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "make_directory"
            "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/deploy"
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "remove_directory"
            "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/bin"
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "make_directory"
            "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/bin"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Removing the old build directories for the Windows installer"
    VERBATIM)

  # **************************************************************
  # Add the target for the shared libraries, 2x 32bit and 1x 64bit
  # **************************************************************
  add_custom_target(COOLPROP_WINDOWS_PACKAGE_SHARED_LIBRARIES)
  add_dependencies(COOLPROP_WINDOWS_PACKAGE_SHARED_LIBRARIES
                   COOLPROP_WINDOWS_PACKAGE_PREPARE)
  # Copy the header file
  add_custom_command(
    TARGET COOLPROP_WINDOWS_PACKAGE_SHARED_LIBRARIES
    PRE_BUILD
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/CoolPropLib.h"
      "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/source/CoolPropLib.h"
    WORKING_DIRECTORY ${COOLPROP_WINDOWS_PACKAGE_32B_DIR}
    COMMENT "Copy the header file for the CoolProp library"
    VERBATIM)
  # Build the 32bit DLLs
  add_custom_command(
    TARGET COOLPROP_WINDOWS_PACKAGE_SHARED_LIBRARIES
    PRE_BUILD
    COMMAND
      ${CMAKE_COMMAND} ARGS "-G${COOLPROP_WINDOWS_PACKAGE_DLL_GEN}" "-AWin32"
      "${CMAKE_CURRENT_SOURCE_DIR}" "-DCOOLPROP_SHARED_LIBRARY=ON"
      "-DCOOLPROP_STDCALL_LIBRARY=ON"
    COMMAND ${CMAKE_COMMAND} ARGS "--build" "." "--target" "CoolProp" "--config"
            "Release"
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy"
      "${COOLPROP_WINDOWS_PACKAGE_32B_DIR_STDCALL}/Release/CoolProp.dll"
      "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/source/CoolProp_stdcall.dll"
    WORKING_DIRECTORY ${COOLPROP_WINDOWS_PACKAGE_32B_DIR_STDCALL}
    COMMENT "Building the 32bit shared library with stdcall"
    VERBATIM)
  add_custom_command(
    TARGET COOLPROP_WINDOWS_PACKAGE_SHARED_LIBRARIES
    PRE_BUILD
    COMMAND
      ${CMAKE_COMMAND} ARGS "-G${COOLPROP_WINDOWS_PACKAGE_DLL_GEN}" "-AWin32"
      "${CMAKE_CURRENT_SOURCE_DIR}" "-DCOOLPROP_SHARED_LIBRARY=ON"
      "-DCOOLPROP_CDECL_LIBRARY=ON"
    COMMAND ${CMAKE_COMMAND} ARGS "--build" "." "--target" "CoolProp" "--config"
            "Release"
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy"
      "${COOLPROP_WINDOWS_PACKAGE_32B_DIR_CDECL}/Release/CoolProp.dll"
      "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/source/CoolProp_cdecl.dll"
    WORKING_DIRECTORY ${COOLPROP_WINDOWS_PACKAGE_32B_DIR_CDECL}
    COMMENT "Building the 32bit shared library with cdecl"
    VERBATIM)
  # Build the 64bit DLL
  add_custom_command(
    TARGET COOLPROP_WINDOWS_PACKAGE_SHARED_LIBRARIES
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} ARGS "-G${COOLPROP_WINDOWS_PACKAGE_DLL_GEN}"
            "-Ax64" "${CMAKE_CURRENT_SOURCE_DIR}" "-DCOOLPROP_SHARED_LIBRARY=ON"
    COMMAND ${CMAKE_COMMAND} ARGS "--build" "." "--target" "CoolProp" "--config"
            "Release"
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy"
      "${COOLPROP_WINDOWS_PACKAGE_64B_DIR}/Release/CoolProp.dll"
      "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/source/CoolProp_x64.dll"
    WORKING_DIRECTORY ${COOLPROP_WINDOWS_PACKAGE_64B_DIR}
    COMMENT "Building the 64bit shared library"
    VERBATIM)

  # *************************************************************
  # Add the target for EES and populate it with custom commands
  # *************************************************************
  add_custom_target(COOLPROP_WINDOWS_PACKAGE_EES)
  add_dependencies(COOLPROP_WINDOWS_PACKAGE_EES
                   COOLPROP_WINDOWS_PACKAGE_PREPARE)
  add_custom_command(
    TARGET COOLPROP_WINDOWS_PACKAGE_EES
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} ARGS "-G${COOLPROP_WINDOWS_PACKAGE_DLL_GEN}"
            "-AWin32" "${CMAKE_CURRENT_SOURCE_DIR}" "-DCOOLPROP_EES_MODULE=ON"
    COMMAND ${CMAKE_COMMAND} ARGS "--build" "." "--target" "COOLPROP_EES"
            "--config" "Release"
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy_directory"
      "${COOLPROP_WINDOWS_PACKAGE_EES_DIR}"
      "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/source/EES"
    WORKING_DIRECTORY ${COOLPROP_WINDOWS_PACKAGE_EES_DIR}
    COMMENT "Building the 32bit library for EES"
    VERBATIM)

  # *************************************************************
  # Add the target for Excel and populate it with custom commands
  # *************************************************************
  add_custom_target(COOLPROP_WINDOWS_PACKAGE_EXCEL)
  add_dependencies(
    COOLPROP_WINDOWS_PACKAGE_EXCEL COOLPROP_WINDOWS_PACKAGE_SHARED_LIBRARIES
    COOLPROP_WINDOWS_PACKAGE_PREPARE)
  # Copy the Excel files
  add_custom_command(
    TARGET COOLPROP_WINDOWS_PACKAGE_EXCEL
    PRE_BUILD
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy"
      "${COOLPROP_WINDOWS_PACKAGE_EXCEL_DIR}/CoolProp.xla"
      "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/source/"
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy"
      "${COOLPROP_WINDOWS_PACKAGE_EXCEL_DIR}/CoolProp.xlam"
      "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/source/"
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy"
      "${COOLPROP_WINDOWS_PACKAGE_EXCEL_DIR}/TestExcel.xlsx"
      "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/source/"
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "remove_directory"
            "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/bin/MicrosoftExcel/"
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "make_directory"
            "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/bin/MicrosoftExcel/"
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy"
      "${COOLPROP_WINDOWS_PACKAGE_EXCEL_DIR}/CoolProp.xla"
      "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/bin/MicrosoftExcel/"
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy"
      "${COOLPROP_WINDOWS_PACKAGE_EXCEL_DIR}/CoolProp.xlam"
      "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/bin/MicrosoftExcel/"
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy"
      "${COOLPROP_WINDOWS_PACKAGE_EXCEL_DIR}/TestExcel.xlsx"
      "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/bin/MicrosoftExcel/"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Copying the Excel files for the installer"
    VERBATIM)

  # *******************************************************************
  # Add the target for Inno Script and populate it with custom commands
  # *******************************************************************
  add_custom_target(COOLPROP_WINDOWS_PACKAGE_ISS)
  add_dependencies(
    COOLPROP_WINDOWS_PACKAGE_ISS COOLPROP_WINDOWS_PACKAGE_EXCEL
    COOLPROP_WINDOWS_PACKAGE_SHARED_LIBRARIES COOLPROP_WINDOWS_PACKAGE_PREPARE)
  # Copy the ISS files
  add_custom_command(
    TARGET COOLPROP_WINDOWS_PACKAGE_ISS
    PRE_BUILD
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy_directory"
      "${COOLPROP_WINDOWS_PACKAGE_ISS_DIR}"
      "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Copying the Inno Script files for the installer"
    VERBATIM)

  # *****************************************************************************
  # Add the target for the installer package and populate it with custom commands
  # *****************************************************************************
  add_custom_target(COOLPROP_WINDOWS_PACKAGE_INSTALLER)
  add_dependencies(COOLPROP_WINDOWS_PACKAGE_INSTALLER
                   COOLPROP_WINDOWS_PACKAGE_DELETE)
  add_dependencies(COOLPROP_WINDOWS_PACKAGE_INSTALLER
                   COOLPROP_WINDOWS_PACKAGE_PREPARE)
  add_dependencies(COOLPROP_WINDOWS_PACKAGE_INSTALLER
                   COOLPROP_WINDOWS_PACKAGE_SHARED_LIBRARIES)
  add_dependencies(COOLPROP_WINDOWS_PACKAGE_INSTALLER
                   COOLPROP_WINDOWS_PACKAGE_EES)
  add_dependencies(COOLPROP_WINDOWS_PACKAGE_INSTALLER
                   COOLPROP_WINDOWS_PACKAGE_EXCEL)
  add_dependencies(COOLPROP_WINDOWS_PACKAGE_INSTALLER
                   COOLPROP_WINDOWS_PACKAGE_ISS)
  # Build the installer and copy it to the bin directory
  add_custom_command(
    TARGET COOLPROP_WINDOWS_PACKAGE_INSTALLER
    POST_BUILD
    COMMAND ${COOLPROP_WINDOWS_PACKAGE_ISS_EXE} ARGS "addin-installer.iss"
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy_directory"
      "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/deploy"
      "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/bin/Installers/Windows"
    WORKING_DIRECTORY "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}"
    COMMENT
      "The new installer is located in '${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/bin/Installers/Windows'"
    VERBATIM)
endif()

if(COOLPROP_OCTAVE_MODULE)

  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_CXX_STANDARD 11)
  endif()

  # Must have SWIG and Octave
  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})
  find_package(Octave REQUIRED)

  # Make a src directory to deal with file permissions problem with MinGW makefile
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src)

  # Set the include folders
  set(OCTAVE_WRAP_INCLUDE_DIRS ${INCLUDE_DIR})
  foreach(ITR ${OCTAVE_INCLUDE_DIRS})
    list(APPEND OCTAVE_WRAP_INCLUDE_DIRS ${ITR})
  endforeach()
  include_directories(${OCTAVE_INCLUDE_DIRS})

  add_definitions(-DNO_ERROR_CATCHING)

  #disable internal error catching and allow swig to do the error catching itself

  set(I_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolProp.i")

  set(SWIG_OPTIONS "${COOLPROP_SWIG_OPTIONS}")
  set_source_files_properties(${I_FILE} PROPERTIES SWIG_FLAGS "${SWIG_OPTIONS}")
  set_source_files_properties(${I_FILE} PROPERTIES CPLUSPLUS ON)

  set(SWIG_MODULE_CoolProp_EXTRA_DEPS ${SWIG_DEPENDENCIES})
  swig_add_module(CoolProp octave ${I_FILE} ${APP_SOURCES})

  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # We need to see which library to link with on OSX - clang++ or stdc++
    message(STATUS "OCTAVE_OCTAVE_LIBRARY = ${OCTAVE_OCTAVE_LIBRARY}")
    if(${CMAKE_VERSION} VERSION_LESS "3.10.0")
      execute_process(COMMAND "otool -L ${OCTAVE_OCTAVE_LIBRARY} | grep libc++"
                      OUTPUT_VARIABLE COOLPROP_OCTAVE_USING_CLANG)
      message(
        STATUS "COOLPROP_OCTAVE_USING_CLANG = ${COOLPROP_OCTAVE_USING_CLANG}")
      string(STRIP "${COOLPROP_OCTAVE_USING_CLANG}" COOLPROP_OCTAVE_USING_CLANG)
    else()
      execute_process(
        COMMAND "otool -L ${OCTAVE_OCTAVE_LIBRARY}"
        COMMAND "grep libc++"
        OUTPUT_VARIABLE COOLPROP_OCTAVE_USING_CLANG
        ERROR_VARIABLE COOLPROP_OCTAVE_USING_CLANG)
      message(
        STATUS "COOLPROP_OCTAVE_USING_CLANG = ${COOLPROP_OCTAVE_USING_CLANG}")
      string(STRIP "${COOLPROP_OCTAVE_USING_CLANG}" COOLPROP_OCTAVE_USING_CLANG)
    endif()

    string(LENGTH "${COOLPROP_OCTAVE_USING_CLANG}" LEN)
    if(${LEN} GREATER 0)
      message(
        STATUS
          "Using -stdlib=libc++, this might override the settings based on DARWIN_USE_LIBCPP"
      )
      set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
    else()
      message(
        STATUS
          "Using -stdlib=libstdc++, this might override the settings based on DARWIN_USE_LIBCPP"
      )
      set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libstdc++")
    endif()
  endif()

  if(WIN32)
    include_directories($ENV{OCTAVE_ROOT}/include)
    include_directories(
      $ENV{OCTAVE_ROOT}/include/octave-${OCTAVE_VERSION}/octave)
    set_target_properties(CoolProp PROPERTIES COMPILE_FLAGS "-fpermissive")
    swig_link_libraries(CoolProp octave octinterp)
    set_target_properties(
      CoolProp
      PROPERTIES
        LINK_FLAGS
        "-L$ENV{OCTAVE_ROOT}/mingw64/lib/octave/${OCTAVE_VERSION} -L$ENV{OCTAVE_ROOT}"
    )
  else()
    swig_link_libraries(CoolProp ${OCTAVE_LIBRARIES})
  endif()

  set_target_properties(CoolProp PROPERTIES SUFFIX ".oct" PREFIX "")
  add_dependencies(${app_name} generate_headers generate_examples)

  #add_custom_command(TARGET CoolProp
  #                   POST_BUILD
  #                   COMMAND "${PYTHON_EXECUTABLE}" example_generator.py Octave "${CMAKE_CURRENT_BINARY_DIR}/Example.m"
  #                   WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dev/scripts/examples")
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/Example.m" DESTINATION Octave)
  install(
    TARGETS ${app_name}
    DESTINATION
      Octave/Octave${OCTAVE_VERSION}_${CMAKE_SYSTEM_NAME}_${BITNESS}bit)
endif()

if(COOLPROP_CSHARP_MODULE)

  # Must have SWIG and C#
  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})
  find_package(Csharp REQUIRED)

  # Make a src directory to deal with file permissions problem with MinGW makefile
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src)

  if(WIN32)
    set(MORE_SWIG_FLAGS -dllimport \"CoolProp\")
  endif()

  # Define which headers the CoolProp wrapper is dependent on
  set(SWIG_MODULE_CoolProp_EXTRA_DEPS ${SWIG_DEPENDENCIES})

  set(SWIG_OPTIONS "${COOLPROP_SWIG_OPTIONS}" "${MORE_SWIG_FLAGS}")
  string(REPLACE " " ";" SWIG_OPTIONS "${SWIG_OPTIONS}")
  message(STATUS "options passed to swig: ${SWIG_OPTIONS}")

  # Set properties before adding module
  set(I_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolProp.i")
  set_source_files_properties(${I_FILE} PROPERTIES SWIG_FLAGS "${SWIG_OPTIONS}"
                                                   CPLUSPLUS ON)

  swig_add_module(CoolProp csharp ${I_FILE} ${APP_SOURCES})

  add_definitions(-DNO_ERROR_CATCHING)

  #disable internal error catching and allow swig to do the error catching itself

  if(WIN32)
    set_target_properties(CoolProp PROPERTIES PREFIX "")
    if(MSVC)
      modify_msvc_flags("/MT")

      # Note that the default is not used if ${COOLPROP_MSVC_REL} or ${COOLPROP_MSVC_DBG} is set
    endif()
  endif()
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(CoolProp PROPERTIES PREFIX "lib")
  endif()
  if(UNIX)
    set_target_properties(CoolProp PROPERTIES PREFIX "lib")
  endif()

  add_dependencies(${app_name} generate_headers generate_examples)

  add_custom_command(
    TARGET CoolProp
    POST_BUILD
    COMMAND 7z a "${CMAKE_CURRENT_BINARY_DIR}/platform-independent.7z"
            "${CMAKE_CURRENT_BINARY_DIR}/*.cs" -x!Example.cs
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  #add_custom_command(TARGET CoolProp
  #                   POST_BUILD
  #                   COMMAND "${PYTHON_EXECUTABLE}" example_generator.py Csharp "${CMAKE_CURRENT_BINARY_DIR}/Example.cs"
  #                   WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dev/scripts/examples")
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/Example.cs" DESTINATION Csharp)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/platform-independent.7z"
          DESTINATION ${CMAKE_INSTALL_PREFIX}/Csharp)
  install(TARGETS ${app_name}
          DESTINATION Csharp/${CMAKE_SYSTEM_NAME}_${BITNESS}bit)
  enable_testing()
  if(DEFINED BUILD_TESTING)
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E make_directory
              ${CMAKE_CURRENT_SOURCE_DIR}/testing_root/Csharp${BITNESS})
    # Copy the shared object to the folder with the executable - no idea like java.library.path in C#
    install(
      TARGETS ${app_name}
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/testing_root/Csharp${BITNESS})
  endif()
  file(TO_NATIVE_PATH ${CMAKE_CURRENT_BINARY_DIR}/*.cs cp_cs_path)
  if(${BITNESS} EQUAL "32")
    set(CSHARP_PLAT "-platform:x86")
  elseif((${BITNESS} EQUAL "64"))
    set(CSHARP_PLAT "-platform:x64")
  endif()
  add_test(
    NAME Csharptestbuild
    COMMAND ${CSHARP_COMPILER} -out:Example.exe ${CSHARP_PLAT} ${cp_cs_path}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testing_root/Csharp${BITNESS})
  add_test(
    NAME Csharptestrun
    COMMAND ${CSHARP_INTERPRETER} Example.exe
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testing_root/Csharp${BITNESS})
endif()

if(COOLPROP_VBDOTNET_MODULE)

  # Must have SWIG and C#
  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})
  find_package(Csharp REQUIRED)

  # Make a src directory to deal with file permissions problem with MinGW makefile
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CoolPropVB)

  set(MORE_SWIG_FLAGS -dllimport \"CoolProp\" -namespace CoolProp)
  set(CMAKE_SWIG_OUTDIR CoolPropVB/CsharpClassLibrary)

  # Define which headers the CoolProp wrapper is dependent on
  set(SWIG_MODULE_CoolProp_EXTRA_DEPS ${SWIG_DEPENDENCIES})

  set(SWIG_OPTIONS "${MORE_SWIG_FLAGS}")
  string(REPLACE " " ";" SWIG_OPTIONS "${SWIG_OPTIONS}")
  message(STATUS "options passed to swig: ${SWIG_OPTIONS}")

  # Set properties before adding module
  set(I_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolProp.i")
  set_property(SOURCE ${I_FILE} PROPERTY CPLUSPLUS ON)
  set_property(SOURCE ${I_FILE} PROPERTY SWIG_FLAGS ${SWIG_OPTIONS})
  swig_add_module(CoolProp csharp ${I_FILE} ${APP_SOURCES})

  add_definitions(-DNO_ERROR_CATCHING)

  #disable internal error catching and allow swig to do the error catching itself

  if(WIN32)
    set_target_properties(CoolProp PROPERTIES PREFIX "")
  endif()

  add_dependencies(${app_name} generate_headers)

  add_custom_command(
    TARGET CoolProp
    PRE_BUILD
    COMMAND
      ${CMAKE_COMMAND} -E copy_directory
      ${CMAKE_CURRENT_SOURCE_DIR}/wrappers/VB.NET/CoolPropVB
      ${CMAKE_CURRENT_BINARY_DIR}/CoolPropVB
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  add_custom_command(
    TARGET CoolProp
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:CoolProp>
            ${CMAKE_CURRENT_BINARY_DIR}/CoolPropVB/CoolPropVB
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  add_custom_command(
    TARGET CoolProp
    POST_BUILD
    COMMAND
      ${CMAKE_COMMAND} -E remove
      ${CMAKE_CURRENT_BINARY_DIR}/CoolPropVB/CsharpClassLibrary/CoolPropCSHARP_wrap.cxx
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  add_custom_command(
    TARGET CoolProp
    POST_BUILD
    COMMAND 7z a "${CMAKE_CURRENT_BINARY_DIR}/VB.net_VS2012_example.7z"
            "${CMAKE_CURRENT_BINARY_DIR}/CoolPropVB"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/VB.net_VS2012_example.7z"
          DESTINATION ${CMAKE_INSTALL_PREFIX}/VB.NET)

endif()

if(COOLPROP_R_MODULE)
  if(WIN32 AND MSVC)
    message(FATAL_ERROR "Must use MinGW Makefiles generator on windows")
  endif()

  # Must have SWIG
  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})

  # Define which headers the swig wrapper is dependent on
  set(SWIG_MODULE_CoolProp_EXTRA_DEPS ${SWIG_DEPENDENCIES})

  find_package(R REQUIRED)
  include_directories(${R_INCLUDE_DIRS})

  link_directories(${R_BIN_OUT})
  if(NOT MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -m${BITNESS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -m${BITNESS}")
  endif()

  add_definitions(-DNO_ERROR_CATCHING)

  #disable internal error catching and allow swig to do the error catching itself

  # Set properties before adding module
  set(I_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolProp.i")
  set_source_files_properties(
    ${I_FILE} PROPERTIES SWIG_FLAGS "${COOLPROP_SWIG_OPTIONS}" CPLUSPLUS ON)

  swig_add_module(CoolProp r ${I_FILE} ${APP_SOURCES})
  swig_link_libraries(CoolProp "${R_LIBRARY}")

  # No lib prefix for the shared library
  set_target_properties(CoolProp PROPERTIES PREFIX "")

  add_dependencies(${app_name} generate_headers generate_examples)
  #add_custom_command(TARGET CoolProp
  #                   POST_BUILD
  #                   COMMAND "${PYTHON_EXECUTABLE}" example_generator.py R "${CMAKE_CURRENT_BINARY_DIR}/Example.R"
  #                   WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dev/scripts/examples")
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/Example.R" DESTINATION R)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/CoolProp.R" DESTINATION R)
  install(TARGETS ${app_name} DESTINATION R/${CMAKE_SYSTEM_NAME}_${BITNESS}bit)

  enable_testing()

  add_test(R_test "${R_BIN_DIR}/Rscript" Example.R)

endif()

if(COOLPROP_JAVA_MODULE)

  # Must have SWIG and Java
  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})
  find_package(Java REQUIRED)
  find_package(JNI)

  # Make a src directory to deal with file permissions problem with MinGW makefile
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src)

  message(STATUS "JAVA_INCLUDE_PATH = ${JAVA_INCLUDE_PATH}")

  include_directories(${JAVA_INCLUDE_PATH})
  include_directories(${JAVA_INCLUDE_PATH}/win32)
  include_directories(${JAVA_INCLUDE_PATH}/linux)

  set(I_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolProp.i")

  set(SWIG_OPTIONS "${COOLPROP_SWIG_OPTIONS}")
  string(REPLACE " " ";" SWIG_OPTIONS "${SWIG_OPTIONS}")
  set_source_files_properties(${I_FILE} PROPERTIES SWIG_FLAGS "${SWIG_OPTIONS}")
  set_source_files_properties(${I_FILE} PROPERTIES CPLUSPLUS ON)

  add_definitions(-DNO_ERROR_CATCHING)

  #disable internal error catching and allow swig to do the error catching itself

  set(SWIG_MODULE_CoolProp_EXTRA_DEPS ${SWIG_DEPENDENCIES})
  swig_add_module(CoolProp java ${I_FILE} ${APP_SOURCES})

  if(WIN32)
    set_target_properties(CoolProp PROPERTIES PREFIX "")
    if(MSVC)
      modify_msvc_flags("/MT")

      # Note that the default is not used if ${COOLPROP_MSVC_REL} or ${COOLPROP_MSVC_DBG} is set
    endif()
  endif()

  if(NOT MSVC)
    set_target_properties(CoolProp PROPERTIES COMPILE_FLAGS "-m${BITNESS}"
                                              LINK_FLAGS "-m${BITNESS}")
  endif()

  add_dependencies(${app_name} generate_headers generate_examples)

  add_custom_command(
    TARGET CoolProp
    POST_BUILD
    COMMAND 7z a "${CMAKE_CURRENT_BINARY_DIR}/platform-independent.7z"
            "${CMAKE_CURRENT_BINARY_DIR}/*.java" -x!Example.java
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  #add_custom_command(TARGET CoolProp
  #                   POST_BUILD
  #                   COMMAND "${PYTHON_EXECUTABLE}" example_generator.py Java "${CMAKE_CURRENT_BINARY_DIR}/Example.java"
  #                   WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dev/scripts/examples")
  #
  # Install all the generated java files
  install(
    CODE "file( GLOB _GeneratedJavaSources \"${CMAKE_CURRENT_BINARY_DIR}/*.java\" )"
    CODE "file( INSTALL \${_GeneratedJavaSources} DESTINATION ${CMAKE_INSTALL_PREFIX}/Java/platform-independent )"
  )
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/Example.java"
          DESTINATION ${CMAKE_INSTALL_PREFIX}/Java)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/platform-independent.7z"
          DESTINATION ${CMAKE_INSTALL_PREFIX}/Java)
  install(
    TARGETS ${app_name}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/Java/${CMAKE_SYSTEM_NAME}_${BITNESS}bit)
  enable_testing()
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E make_directory
            ${CMAKE_CURRENT_SOURCE_DIR}/testing_root/Java${BITNESS})
  add_test(
    NAME Javatestbuild
    COMMAND javac -d . ${CMAKE_INSTALL_PREFIX}/Java/Example.java -cp
            ${CMAKE_INSTALL_PREFIX}/Java/platform-independent
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testing_root/Java${BITNESS})
  add_test(
    NAME Javatestrun
    COMMAND
      ${Java_JAVA_EXECUTABLE}
      -Djava.library.path=${CMAKE_INSTALL_PREFIX}/Java/${CMAKE_SYSTEM_NAME}_${BITNESS}bit
      Example
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testing_root/Java${BITNESS})
endif()

# A module for Android
if(COOLPROP_ANDROID_MODULE)

  if(WIN32 AND (NOT MINGW))
    message(
      FATAL_ERROR "On windows, you must use the MinGW Makefiles generator ")
  endif()

  # For now, these must be changed manually
  set(ANDROID_MODULE_NAME "CoolProp")
  set(ANDROID_PACKAGE_NAME "CoolProp") # or blah.di.blah.CoolProp

  # Must have SWIG
  find_package(SWIG REQUIRED)

  set(I_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolProp.i")

  list(APPEND APP_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/jni/CoolProp_wrap.cxx)
  string(REPLACE ";" " " APP_INCLUDE_DIRS "${APP_INCLUDE_DIRS}")
  string(REPLACE ";" " " APP_SOURCES "${APP_SOURCES}")

  file(MAKE_DIRECTORY jni)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/Android/Android.mk.template"
    "${CMAKE_CURRENT_BINARY_DIR}/jni/Android.mk")
  file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/Android/Application.mk"
       DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/jni")
  string(REPLACE "." "/" ANDROID_PACKAGE_PATH "${ANDROID_PACKAGE_NAME}")
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${ANDROID_PACKAGE_PATH}")

  message(STATUS "WORKING_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}")
  get_filename_component(NDK_BUILD_PATH "${NDK_PATH}/ndk-build" ABSOLUTE)
  get_filename_component(SRC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src" ABSOLUTE)
  get_filename_component(INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include"
                         ABSOLUTE)

  add_custom_target(
    CoolProp ALL
    COMMAND ${NDK_BUILD_PATH}
    DEPENDS jni/CoolProp_wrap.cxx
    VERBATIM)

  add_custom_command(
    OUTPUT jni/CoolProp_wrap.cxx
    COMMAND
      ${SWIG_EXECUTABLE} -v -c++ -java -I${SRC_PATH} -I${INCLUDE_PATH} -o
      ${CMAKE_CURRENT_BINARY_DIR}/jni/CoolProp_wrap.cxx -package
      ${ANDROID_PACKAGE_NAME} -outdir
      ${CMAKE_CURRENT_BINARY_DIR}/${ANDROID_PACKAGE_PATH} ${I_FILE}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    VERBATIM)

  add_dependencies(CoolProp generate_headers)

endif()

if(COOLPROP_PHP_MODULE)

  # Must have SWIG
  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})

  execute_process(
    COMMAND php-config --includes
    OUTPUT_VARIABLE php_config_includes
    RESULT_VARIABLE php_config_failed)
  if(php_config_failed)
    message(FATAL_ERROR "calling \"php-config --includes\" failed; message:"
                        ${php_config_includes})
  endif()
  string(STRIP "${php_config_includes}" php_config_includes)
  string(REPLACE "-I" "" PHP_INCLUDES "${php_config_includes}")
  separate_arguments(PHP_INCLUDES)

  message(STATUS "php includes=${PHP_INCLUDES}")
  include_directories(${PHP_INCLUDES})

  add_definitions(-DNO_ERROR_CATCHING)

  #disable internal error catching and allow swig to do the error catching itself

  set(I_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolProp.i")

  set(SWIG_OPTIONS "${COOLPROP_SWIG_OPTIONS}")
  string(REPLACE " " ";" SWIG_OPTIONS "${SWIG_OPTIONS}")
  set_source_files_properties(${I_FILE} PROPERTIES SWIG_FLAGS "${SWIG_OPTIONS}")
  set_source_files_properties(${I_FILE} PROPERTIES CPLUSPLUS ON)

  set(SWIG_MODULE_CoolProp_EXTRA_DEPS ${SWIG_DEPENDENCIES})
  swig_add_module(CoolProp php ${I_FILE} ${APP_SOURCES})

  if(WIN32)
    set_target_properties(CoolProp PROPERTIES PREFIX "")
  endif()

  if(NOT MSVC)
    set_target_properties(CoolProp PROPERTIES COMPILE_FLAGS "-m${BITNESS}"
                                              LINK_FLAGS "-m${BITNESS}")
  endif()
  add_dependencies(CoolProp generate_headers)

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CoolProp.php
          DESTINATION ${CMAKE_INSTALL_PREFIX}/PHP/cross-platform)
  install(TARGETS ${app_name}
          DESTINATION ${CMAKE_INSTALL_PREFIX}/PHP/${CMAKE_SYSTEM_NAME})

endif()

function(JOIN VALUES GLUE OUTPUT)
  string(REGEX REPLACE "([^\\]|^);" "\\1${GLUE}" _TMP_STR "${VALUES}")
  string(REGEX REPLACE "[\\](.)" "\\1" _TMP_STR "${_TMP_STR}") #fixes escaping
  set(${OUTPUT}
      "${_TMP_STR}"
      PARENT_SCOPE)
endfunction()

if(COOLPROP_PYTHON_BINARIES)
  if(WIN32)
    set(COOLPROP_PYTHON_BINARY_VERSIONS
        bdist_wheel --dist-dir ${CMAKE_INSTALL_PREFIX}/Python bdist_wininst
        --dist-dir ${CMAKE_INSTALL_PREFIX}/Python)
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(COOLPROP_PYTHON_BINARY_VERSIONS bdist_wheel --dist-dir
                                        ${CMAKE_INSTALL_PREFIX}/Python)
  endif()

  add_custom_target(
    CoolProp
    COMMAND python setup.py ${COOLPROP_PYTHON_BINARY_VERSIONS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/wrappers/Python)

endif()

if(COOLPROP_PYTHON_PYPI)

  add_custom_target(
    CoolProp
    COMMAND python prepare_pypi.py --dist-dir=${CMAKE_INSTALL_PREFIX}/Python
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/wrappers/Python/pypi)
endif()

if(COOLPROP_LIBREOFFICE_MODULE)

  if("${LO_PROGRAM_PATH}" STREQUAL "")
    message(
      FATAL_ERROR
        "You must provide the path to LibreOffice programs, something like -DLO_PROGRAM_PATH=/usr/lib/libreoffice/program"
    )
  else()
    message(STATUS "LO_PROGRAM_PATH: ${LO_PROGRAM_PATH}")
  endif()

  if("${LO_SDK_PATH}" STREQUAL "")
    message(
      FATAL_ERROR
        "You must provide the path to LibreOffice SDK, something like -DLO_SDK_PATH=/usr/lib/libreoffice/sdk"
    )
  else()
    message(STATUS "LO_SDK_PATH: ${LO_SDK_PATH}")
  endif()

  add_custom_target(CoolPropLibreOfficeAddin)

  # set paths for LibreOffice tools
  set(LO_IDLC "${LO_SDK_PATH}/bin/idlc")
  set(LO_REGMERGE "${LO_PROGRAM_PATH}/regmerge")
  set(COOLPROP_LIBREOFFICE_TMP_DIR "${CMAKE_CURRENT_BINARY_DIR}/LibreOffice")

  # set version strings for LibreOffice extension
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/LibreOffice/src/description.xml.in"
    "${COOLPROP_LIBREOFFICE_TMP_DIR}/src/description.xml")
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/LibreOffice/src/scripts/scripts.py.in"
    "${COOLPROP_LIBREOFFICE_TMP_DIR}/src/scripts/scripts.py")

  add_custom_command(
    TARGET CoolPropLibreOfficeAddin
    # copy source files to build directory
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy_directory"
      "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/LibreOffice/src"
      "${COOLPROP_LIBREOFFICE_TMP_DIR}/src"
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "remove"
      "${COOLPROP_LIBREOFFICE_TMP_DIR}/src/description.xml.in"
      "${COOLPROP_LIBREOFFICE_TMP_DIR}/src/scripts/scripts.py.in"
    # rebuild the registry database file (rdb)
    COMMAND ${LO_IDLC} XCoolProp.idl -I. -I${LO_SDK_PATH}/idl -O. XCoolProp.idl
    COMMAND ${LO_REGMERGE} XCoolProp.rdb /UCR XCoolProp.urd
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "remove" XCoolProp.urd
    # download and bundle latest Python pip package (py2.py3, platform independent)
    COMMAND pip download pip -d pythonpath
    COMMAND 7z x "./pythonpath/pip-*.whl" -y -opythonpath
    # download and bundle latest Python certifi package (py2.py3, platform independent)
    COMMAND pip download certifi -d pythonpath
    COMMAND 7z x "./pythonpath/certifi-*.whl" -y -opythonpath
    # add license file
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "make_directory"
            "${COOLPROP_LIBREOFFICE_TMP_DIR}/src/license"
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy" "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE"
      "${COOLPROP_LIBREOFFICE_TMP_DIR}/src/license/."
    # package complete folder to extension
    COMMAND 7z a -tzip "../CoolProp.oxt"
    # copy example spreadsheet file
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy"
      "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/LibreOffice/TestLibreOffice.ods"
      "${COOLPROP_LIBREOFFICE_TMP_DIR}/."
    WORKING_DIRECTORY ${COOLPROP_LIBREOFFICE_TMP_DIR}/src
    COMMENT "Building LibreOffice wrapper"
    VERBATIM)
endif()

if(COOLPROP_JAVASCRIPT_MODULE)
  # cmake -DCOOLPROP_JAVASCRIPT_MODULE=ON
  #       -DCMAKE_TOOLCHAIN_FILE=${EMSCRIPTEN}/cmake/Platform/Emscripten.cmake
  #       ../..

  # Toolchain MUST be defined in the call to CMake

  if(MSVC)
    message(
      FATAL_ERROR
        "Cannot use visual studio, use MinGW Makefiles generator on windows")
  endif()

  add_definitions(-sDISABLE_EXCEPTION_CATCHING=0)
  # If you want a monolithic file with no async memory loading, define EMSCRIPTEN_NO_MEMORY_INIT_FILE
  if(EMSCRIPTEN_NO_MEMORY_INIT_FILE)
    set(EMSCRIPTEN_INIT_FLAG "--memory-init-file 0")
  else()
    set(EMSCRIPTEN_INIT_FLAG "--memory-init-file 1")
  endif()

  set(CMAKE_EXE_LINKER_FLAGS
      "--bind ${EMSCRIPTEN_INIT_FLAG} -s ASSERTIONS=1 -s DISABLE_EXCEPTION_CATCHING=0"
  )
  set(CMAKE_BUILD_TYPE Release)

  list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolPropLib.cpp"
       "${CMAKE_CURRENT_SOURCE_DIR}/src/emscripten_interface.cxx")
  include_directories(${APP_INCLUDE_DIRS})
  add_executable(coolprop ${APP_SOURCES})
  add_dependencies(coolprop generate_headers)
  set_target_properties(coolprop PROPERTIES PREFIX "" SUFFIX .js)
  #install (TARGETS coolprop DESTINATION ${CMAKE_INSTALL_PREFIX}/Javascript)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/coolprop.js"
          DESTINATION ${CMAKE_INSTALL_PREFIX}/Javascript)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/coolprop.wasm"
          DESTINATION ${CMAKE_INSTALL_PREFIX}/Javascript)
  #install (FILES "${CMAKE_CURRENT_BINARY_DIR}/install_manifest.txt" DESTINATION ${CMAKE_INSTALL_PREFIX}/Javascript)
  install(
    FILES
      "${CMAKE_CURRENT_SOURCE_DIR}/Web/coolprop/wrappers/Javascript/index.html"
    DESTINATION ${CMAKE_INSTALL_PREFIX}/Javascript)
endif()

if(COOLPROP_MATHEMATICA_MODULE)

  if(MSVC)
    modify_msvc_flags("/MT")

    # Note that the default is not used if ${COOLPROP_MSVC_REL} or ${COOLPROP_MSVC_DBG} is set
  endif()

  set(CMAKE_MODULE_PATH
      ${CMAKE_MODULE_PATH}
      "${CMAKE_CURRENT_SOURCE_DIR}/externals/FindMathematica/CMake/Mathematica/"
  )
  find_package(Mathematica COMPONENTS WolframLibrary)
  message(
    STATUS
      "Mathematica_WolframLibrary_FOUND=${Mathematica_WolframLibrary_FOUND}")
  message(
    STATUS
      "Mathematica_WolframLibrary_INCLUDE_DIR=${Mathematica_WolframLibrary_INCLUDE_DIR}"
  )
  message(STATUS "Mathematica_USERBASE_DIR=${Mathematica_USERBASE_DIR}")

  list(
    APPEND APP_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/Mathematica/CoolPropMathematica.cpp")
  list(APPEND APP_INCLUDE_DIRS "${Mathematica_WolframLibrary_INCLUDE_DIR}")
  include_directories(${APP_INCLUDE_DIRS})
  add_library(CoolProp SHARED ${APP_SOURCES})
  add_dependencies(CoolProp generate_headers)

  if(MSVC)
    add_custom_command(
      TARGET ${app_name}
      POST_BUILD
      COMMAND dumpbin /EXPORTS $<TARGET_FILE:CoolProp> >
              ${CMAKE_CURRENT_BINARY_DIR}/exports.txt)
  endif()

  install(FILES $<TARGET_FILE:CoolProp>
          DESTINATION Mathematica/${CMAKE_SYSTEM_NAME})
  install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/Mathematica/example.nb"
          DESTINATION Mathematica)
endif()

if(COOLPROP_SMATH_MODULE)
  if(COOLPROP_SMATH_WORK_INPLACE)
    set(COOLPROP_WORK_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
  else()
    set(COOLPROP_WORK_BASE_DIR ${CMAKE_CURRENT_BINARY_DIR})
  endif()
  set(COOLPROP_VERSION
      ${COOLPROP_VERSION_MAJOR}.${COOLPROP_VERSION_MINOR}.${COOLPROP_VERSION_PATCH}.0
  )
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/SMath/coolprop_wrapper/Properties/AssemblyInfo.cs.template"
    "${COOLPROP_WORK_BASE_DIR}/wrappers/SMath/coolprop_wrapper/Properties/AssemblyInfo.cs"
  )
  message(
    STATUS
      "Generated ${COOLPROP_WORK_BASE_DIR}/wrappers/SMath/coolprop_wrapper/Properties/AssemblyInfo.cs"
  )
  file(WRITE "${COOLPROP_WORK_BASE_DIR}/wrappers/SMath/config.ini"
       "${COOLPROP_VERSION}")
  message(
    STATUS "Generated ${COOLPROP_WORK_BASE_DIR}/wrappers/SMath/config.ini")
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/SMath/coolprop_wrapper/install.bat.template"
    "${COOLPROP_WORK_BASE_DIR}/wrappers/SMath/coolprop_wrapper/install.bat")
  message(
    STATUS
      "Generated ${COOLPROP_WORK_BASE_DIR}/wrappers/SMath/coolprop_wrapper/install.bat"
  )
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/SMath/coolprop_wrapper/build_zip.bat.template"
    "${COOLPROP_WORK_BASE_DIR}/wrappers/SMath/coolprop_wrapper/build_zip.bat")
  message(
    STATUS
      "Generated ${COOLPROP_WORK_BASE_DIR}/wrappers/SMath/coolprop_wrapper/build_zip.bat"
  )
  file(TO_NATIVE_PATH
       "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/SMath/coolprop_wrapper"
       DOS_STYLE_SOURCE_DIR)
  file(TO_NATIVE_PATH
       "${COOLPROP_WORK_BASE_DIR}/wrappers/SMath/coolprop_wrapper"
       DOS_STYLE_TARGET_DIR)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/SMath/coolprop_wrapper/coolprop_wrapper.csproj.template"
    "${COOLPROP_WORK_BASE_DIR}/wrappers/SMath/coolprop_wrapper/coolprop_wrapper.csproj"
  )
  message(
    STATUS
      "Generated ${COOLPROP_WORK_BASE_DIR}/wrappers/SMath/coolprop_wrapper/coolprop_wrapper.csproj"
  )
  include_external_msproject(
    CoolPropWrapper
    ${COOLPROP_WORK_BASE_DIR}/wrappers/SMath/coolprop_wrapper/coolprop_wrapper.csproj
    TYPE FAE04EC0-301F-11D3-BF4B-00C04F79EFBC
    PLATFORM AnyCPU)
  message(
    STATUS
      "C# project ${COOLPROP_WORK_BASE_DIR}/wrappers/SMath/coolprop_wrapper/coolprop_wrapper.csproj included"
  )

endif()

# Use like cmake ..\CoolProp.git -DCOOLPROP_MY_MAIN=dev/coverity/main.cxx
if(COOLPROP_MY_MAIN)
  list(APPEND APP_SOURCES "${COOLPROP_MY_MAIN}")
  add_executable(Main ${APP_SOURCES})
  add_dependencies(Main generate_headers)
  if(UNIX)
    target_link_libraries(Main ${CMAKE_DL_LIBS})
  endif()
endif()

if(COOLPROP_MAIN_MODULE)
  # Allow you to independently add back the testing CPP files
  if(COOLPROP_TEST)
    list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Tests/Tests.cpp")
    list(APPEND APP_SOURCES
         "${CMAKE_CURRENT_SOURCE_DIR}/src/Tests/CoolProp-Tests.cpp")
  endif()
  list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cxx")
  add_executable(Main ${APP_SOURCES})
  add_dependencies(Main generate_headers)
  if(COOLPROP_TEST)
    set_target_properties(Main PROPERTIES COMPILE_FLAGS
                                          "${COMPILE_FLAGS} -DENABLE_CATCH")
  endif()

  if(COOLPROP_IWYU)
    find_program(iwyu_path NAMES include-what-you-use iwyu)
    if(NOT iwyu_path)
      message(FATAL_ERROR "Could not find the program include-what-you-use")
    endif()
    set_property(TARGET Main PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${iwyu_path})
  endif()

  if(UNIX)
    target_link_libraries(Main ${CMAKE_DL_LIBS})
  endif()
endif()

###      COOLPROP TESTING APP       ###
if(COOLPROP_CATCH_MODULE)

  add_subdirectory("externals/Catch2")
  enable_testing()

  list(APPEND APP_SOURCES
       "${CMAKE_CURRENT_SOURCE_DIR}/src/Tests/CoolProp-Tests.cpp")

  # CATCH TEST, compile everything with catch and set test entry point
  add_executable(CatchTestRunner ${APP_SOURCES})

  target_link_libraries(CatchTestRunner PRIVATE Catch2::Catch2WithMain)
  target_compile_definitions(CatchTestRunner PRIVATE ENABLE_CATCH)
  if(UNIX)
    target_link_libraries(CatchTestRunner ${CMAKE_DL_LIBS})
  endif()

  include(CTest)
  include(${CMAKE_CURRENT_SOURCE_DIR}/externals/Catch2/extras/Catch.cmake)
  catch_discover_tests(CatchTestRunner)

  if(COOLPROP_IWYU)
    find_program(iwyu_path NAMES include-what-you-use iwyu)
    if(NOT iwyu_path)
      message(FATAL_ERROR "Could not find the program include-what-you-use")
    endif()
    set_property(TARGET CatchTestRunner PROPERTY CXX_INCLUDE_WHAT_YOU_USE
                                                 ${iwyu_path})
  endif()
endif()

if(COOLPROP_CPP_EXAMPLE_TEST)
  # C++ Documentation Test
  add_executable(docuTest.exe "Web/examples/C++/Example.cpp")
  add_dependencies(docuTest.exe ${app_name})
  target_link_libraries(docuTest.exe ${app_name})
  if(UNIX)
    target_link_libraries(docuTest.exe ${CMAKE_DL_LIBS})
  endif()
  add_test(DocumentationTest docuTest.exe)
endif()

if(COOLPROP_SNIPPETS)
  list(APPEND APP_SOURCES
       "${CMAKE_CURRENT_SOURCE_DIR}/${COOLPROP_LIBRARY_SOURCE}")
  # Make the static library with which the snippets will be linked
  add_library(${app_name} STATIC ${APP_SOURCES})
  add_dependencies(${app_name} generate_headers)
  set_property(
    TARGET ${app_name}
    APPEND_STRING
    PROPERTY COMPILE_FLAGS " -DEXTERNC")

  # Collect all the snippets
  file(GLOB_RECURSE snippets
       "${CMAKE_CURRENT_SOURCE_DIR}/Web/coolprop/snippets/*.cxx")

  message(STATUS "snippets found = ${snippets}")
  foreach(snippet ${snippets})

    get_filename_component(snippet_name ${snippet} NAME)
    get_filename_component(snippet_exe ${snippet} NAME_WE)
    message(STATUS "snippet_name = ${snippet_name}")

    add_executable(${snippet_exe} ${snippet})
    add_dependencies(${snippet_exe} CoolProp)
    target_link_libraries(${snippet_exe} CoolProp)
    if(UNIX)
      target_link_libraries(${snippet_exe} ${CMAKE_DL_LIBS})
    endif()

    if(MSVC)
      set_target_properties(
        ${snippet_exe} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                  ${CMAKE_CURRENT_BINARY_DIR}/bin)
      set_target_properties(
        ${snippet_exe} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
                                  ${CMAKE_CURRENT_BINARY_DIR}/bin)
      set_target_properties(
        ${snippet_exe} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
                                  ${CMAKE_CURRENT_BINARY_DIR}/bin)
      # etc for the other available configuration types (MinSizeRel, RelWithDebInfo)
      set(BIN_PATH "${CMAKE_CURRENT_BINARY_DIR}/bin")
    else()
      set(BIN_PATH "${CMAKE_CURRENT_BINARY_DIR}")
    endif()

    set_property(
      TARGET ${snippet_exe}
      APPEND_STRING
      PROPERTY COMPILE_FLAGS " -DEXTERNC")

    # Run it and save the output to a file with .output appended
    add_custom_command(
      TARGET ${snippet_exe}
      POST_BUILD
      COMMAND
        ${BIN_PATH}/${snippet_exe} >
        ${CMAKE_CURRENT_SOURCE_DIR}/Web/coolprop/snippets/${snippet_name}.output
    )

  endforeach()

endif()

if(COOLPROP_CLANG_ADDRESS_SANITIZER)

  set(CMAKE_CXX_FLAGS "-fsanitize=address -g")
  list(APPEND APP_SOURCES
       "${CMAKE_CURRENT_SOURCE_DIR}/src/Tests/catch_always_return_success.cxx")
  # CATCH TEST, compile everything with catch and set test entry point
  add_executable(CatchTestRunner ${APP_SOURCES})
  add_dependencies(CatchTestRunner generate_headers)
  set_target_properties(
    CatchTestRunner PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DENABLE_CATCH")
  set(CMAKE_CXX_FLAGS "-O1")
  set(CMAKE_EXE_LINKER_FLAGS
      "-fsanitize=address -fno-omit-frame-pointer -lstdc++")
  if(UNIX)
    target_link_libraries(CatchTestRunner ${CMAKE_DL_LIBS})
  endif()
  add_custom_command(
    TARGET CatchTestRunner
    POST_BUILD
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/CatchTestRunner)
endif()

if(COOLPROP_PROFILE)
  if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-g -O2")
    set(CMAKE_C_FLAGS "-g -O2")
  endif()
endif()

if(COOLPROP_COVERAGE)
  if(CMAKE_COMPILER_IS_GNUCXX)
    # See also http://stackoverflow.com/a/16536401 (detailed guide on using gcov with cmake)
    include(CodeCoverage)
    set(CMAKE_CXX_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
    set(CMAKE_C_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
    setup_target_for_coverage(CoolProp_coverage Main coverage)
  endif()
endif()

# TODO: check relevance of http://www.cmake.org/Wiki/BuildingWinDLL

#include_directories("${CMAKE_CURRENT_SOURCE_DIR}/CoolProp")
#FILE(GLOB coolprop_files "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/*.cpp")
#add_library(coolprop STATIC ${coolprop_files})
