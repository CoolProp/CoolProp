cmake_minimum_required(VERSION 2.8)


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
set (CoolProp_VERSION_MAJOR 5)
set (CoolProp_VERSION_MINOR 0)
set (CoolProp_VERSION_PATCH 0)
set (CoolProp_VERSION ${CoolProp_VERSION_MAJOR}.${CoolProp_VERSION_MINOR}.${CoolProp_VERSION_PATCH})

#######################################
#           BUILD OPTIONS             #
#-------------------------------------#
# These options are available to be   #
# modified in the build process.      #
# packages may want to modify these   #
# to suit, or just leave as defaults. #
#######################################

option (COOLPROP_STATIC_LIBRARY
        "Build CoolProp as a static library (.lib, .a)"
        OFF)

option (COOLPROP_SHARED_LIBRARY
        "Build CoolProp as a shared library (.dll, .so)"
        OFF)
        
option (COOLPROP_EES_MODULE
        "Build the EES module"
        OFF)
        
option (COOLPROP_TESTING 
        "Build with test flags and include main() from catch" 
        OFF)

#######################################
#         FIND ALL SOURCES            #
#-------------------------------------#
# The project is organised with       #
# split includes and source folders   #
# this makes it easier for developers #
# to quickly find relevant includes.  #
# This section finds all sources,     #
# headers and corrosponding dirs.     #
#######################################

file(GLOB_RECURSE APP_SOURCES "src/*.cpp")
file(GLOB_RECURSE APP_HEADERS "include/*.h" "src/*.h") # "externals/*.hpp")
list(REMOVE_ITEM  APP_SOURCES "${CMAKE_SOURCE_DIR}/src/Tests/Tests.cpp")

set (APP_INCLUDE_DIRS "")
foreach (_headerFile ${APP_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND APP_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES APP_INCLUDE_DIRS)

include_directories(${APP_INCLUDE_DIRS})


#######################################
#         REQUIRED MODULES            #
#-------------------------------------#
# CoolProp requires some standard OS  #
# features, these include:            #
# DL (CMAKE_DL_LIBS) for REFPROP      #
#######################################
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/dev/cmake/Modules/")

set(Python_ADDITIONAL_VERSIONS 2.7 2.6 2.5 2.4)
find_package (PythonInterp 2.7 REQUIRED)
if(UNIX)
    find_package (${CMAKE_DL_LIBS} REQUIRED)
endif()


#######################################
#         MAKE ARTIFACTS              #
#-------------------------------------#
# In this section we define the       #
# artifacts (exes, libs) that will be #
# made for coolprop, these include    #
# customisation from earier options.  #
#######################################

###     FLUIDS, MIXTURES JSON       ###
add_custom_target(generate_headers
                  COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_SOURCE_DIR}/dev/generate_headers.py")

###      COOLPROP LIB or DLL        ###
if (COOLPROP_STATIC_LIBRARY)
    add_library(${app_name} STATIC ${APP_SOURCES})
    add_dependencies (${app_name} generate_headers)
endif()

##### COOLPROP SHARED LIBRARY ######
if (COOLPROP_SHARED_LIBRARY)
    add_library(${app_name} SHARED ${APP_SOURCES})
    set_target_properties (${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB")
    add_dependencies (${app_name} generate_headers)
endif()

if (COOLPROP_64BIT_SHARED_LIBRARY)
    add_library(${app_name} SHARED ${APP_SOURCES})
    set_target_properties (${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB")
    if (!MSVC)
      set_target_properties(${app_name} PROPERTIES COMPILE_FLAGS "-m64" LINK_FLAGS "-m64")
    endif()
    add_dependencies (${app_name} generate_headers)
    set_target_properties(${app_name} PROPERTIES PREFIX "")
endif()

if (COOLPROP_32BIT_CDECL_SHARED_LIBRARY)
    add_library(${app_name} SHARED ${APP_SOURCES})
    set_target_properties (${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB -DCONVENTION=__cdecl")
    if (!MSVC)
      set_target_properties(${app_name} PROPERTIES COMPILE_FLAGS "-m32" LINK_FLAGS "-m32")
    endif()
    add_dependencies (${app_name} generate_headers)
    set_target_properties(${app_name} PROPERTIES PREFIX "")
endif()

if (COOLPROP_32BIT_STDCALL_SHARED_LIBRARY)
    add_library(${app_name} SHARED ${APP_SOURCES})
    set_target_properties (${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB -DCONVENTION=__stdcall")
    if (!MSVC)
      set_target_properties(${app_name} PROPERTIES COMPILE_FLAGS "-m32" LINK_FLAGS "-m32")
    endif()
    add_dependencies (${app_name} generate_headers)
    set_target_properties(${app_name} PROPERTIES PREFIX "")
endif()

if (COOLPROP_EES_MODULE)
  list (APPEND APP_SOURCES "wrappers/EES/main.cpp")
  include_directories(${APP_INCLUDE_DIRS})
  add_library(COOLPROP_EES SHARED ${APP_SOURCES})
  set_target_properties (COOLPROP_EES PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB -DCONVENTION=__cdecl")
  add_dependencies (COOLPROP_EES generate_headers)
  set_target_properties(COOLPROP_EES PROPERTIES SUFFIX ".dlf" PREFIX "")
  if (!MSVC)
    set_target_properties(COOLPROP_EES PROPERTIES COMPILE_FLAGS "-m32" LINK_FLAGS "-m32")
  endif()
endif()

if (COOLPROP_OCTAVE_MODULE)

  # Must have SWIG and Octave
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE(${SWIG_USE_FILE})
  FIND_PACKAGE(Octave REQUIRED)
  
  # Find the required libraries
  find_library(OCTAVE_LIB octave PATHS ${OCTAVE_LINK_DIRS})
  find_library(OCTINTERP_LIB octinterp PATHS ${OCTAVE_LINK_DIRS})
  find_library(CRUFT_LIB cruft PATHS ${OCTAVE_LINK_DIRS})
  
  # Set the include folders
  SET(OCTAVE_WRAP_INCLUDE_DIRS ${INCLUDE_DIR})
  foreach(ITR ${OCTAVE_INCLUDE_DIRS})
    list(APPEND OCTAVE_WRAP_INCLUDE_DIRS ${ITR})
  endforeach()
  include_directories(${OCTAVE_WRAP_INCLUDE_DIRS})

  set(I_FILE "${CMAKE_SOURCE_DIR}/src/CoolProp.i")

  SET_SOURCE_FILES_PROPERTIES(${I_FILE} PROPERTIES CPLUSPLUS ON)
  SWIG_ADD_MODULE(CoolProp octave ${I_FILE} ${APP_SOURCES})
  SWIG_LINK_LIBRARIES(CoolProp ${OCTAVE_LIB} ${OCTINTERP_LIB} ${CRUFT_LIB})

  set_target_properties(CoolProp PROPERTIES SUFFIX ".oct" PREFIX "")
  
endif()

# NOT WORKING!
if (COOLPROP_MATHEMATICA_MODULE)

  # copy "C:\Program Files\Wolfram Research\Mathematica\9.0\SystemFiles\IncludeFiles\C\WolframLibrary.h" .
  # REM folder obtained from FileNameJoin[{$BaseDirectory, "SystemFiles", "LibraryResources", $SystemID}] in Mathematica
  # copy CoolProp.dll "C:\ProgramData\Mathematica\SystemFiles\LibraryResources\Windows-x86-64\"

  list (APPEND APP_SOURCES "wrappers/Mathematica/CoolPropMathematica.cpp")
  include_directories(${APP_INCLUDE_DIRS})
  add_library(COOLProp SHARED ${APP_SOURCES})
  set_target_properties (COOLPROP_EES PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB -DCONVENTION=__cdecl")
  add_dependencies (COOLPROP_EES generate_headers)
  set_target_properties(COOLPROP_EES PROPERTIES SUFFIX ".dlf" PREFIX "")
endif()

###      COOLPROP TESTING APP       ###
if (COOLPROP_TESTING)
  list (APPEND APP_INCLUDE_DIRS "externals/Catch/include")
  include_directories(${APP_INCLUDE_DIRS})
  enable_testing()

  # CATCH TEST, compile everything with catch and set test entry point
  add_executable        (testRunner.exe ${APP_SOURCES})
  add_dependencies      (testRunner.exe generate_headers)
  set_target_properties (testRunner.exe PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DENABLE_CATCH")
  set_source_files_properties(src/PolyMath.cpp PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCATCH_CONFIG_MAIN")
  if(UNIX)
    target_link_libraries (testRunner.exe ${CMAKE_DL_LIBS})
  endif()
  add_test(ProceedureTests testRunner.exe)

  # C++ Documentation Test
  add_executable        (docuTest.exe "Web/examples/C++/Example.cpp")
  add_dependencies      (docuTest.exe ${app_name})
  target_link_libraries (docuTest.exe ${app_name})
  if(UNIX)
    target_link_libraries (docuTest.exe ${CMAKE_DL_LIBS})
  endif()
  add_test(DocumentationTest docuTest.exe)
endif()


# TODO: check relevance of http://www.cmake.org/Wiki/BuildingWinDLL

#include_directories("${CMAKE_CURRENT_SOURCE_DIR}/CoolProp")
#FILE(GLOB coolprop_files "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/*.cpp")
#add_library(coolprop STATIC ${coolprop_files})
