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
        "Build and install CoolProp as a STATIC library (.lib, .a) as opposed to SHARED (.dll, .so)"
        ON)

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

set (APP_INCLUDE_DIRS "")
foreach (_headerFile ${APP_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND APP_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES APP_INCLUDE_DIRS)
include_directories(${APP_INCLUDE_DIRS})


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
                  COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/dev/generate_headers.py")


###      COOLPROP LIB or DLL        ###
if (COOLPROP_STATIC_LIBRARY)
    add_library(${app_name} STATIC ${APP_SOURCES})
else()
    add_library(${app_name} SHARED ${APP_SOURCES})
endif()
add_dependencies (${app_name} generate_headers)


###      COOLPROP TESTING APP       ###
if (COOLPROP_TESTING)
  list (APPEND APP_INCLUDE_DIRS "externals/Catch/include")
  include_directories(${APP_INCLUDE_DIRS})
  enable_testing()

  # CATCH TEST, compile everything with catch and set test entry point
  add_executable        (testRunner.exe ${APP_SOURCES})
  add_dependencies      (testRunner.exe generate_headers)
  set_target_properties (testRunner.exe PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DENABLE_CATCH")
  set_source_files_properties(src/Helmholtz.cpp PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCATCH_CONFIG_MAIN")
  target_link_libraries (testRunner.exe ${CMAKE_DL_LIBS})
  add_test(ProceedureTests testRunner.exe)


  # C++ Documentation Test
  add_executable        (docuTest.exe "Web/examples/C++/Example.cpp")
  add_dependencies      (docuTest.exe ${app_name})
  target_link_libraries (docuTest.exe ${app_name})
  target_link_libraries (docuTest.exe ${CMAKE_DL_LIBS})
  add_test(DocumentationTest docuTest.exe)
endif()


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
find_package (${CMAKE_DL_LIBS} REQUIRED)

# TODO: check relevance of http://www.cmake.org/Wiki/BuildingWinDLL

#include_directories("${CMAKE_CURRENT_SOURCE_DIR}/CoolProp")
#FILE(GLOB coolprop_files "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/*.cpp")
#add_library(coolprop STATIC ${coolprop_files})
