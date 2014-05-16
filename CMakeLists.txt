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
        

#######################################
#         FIND ALL SOURCES            #
#-------------------------------------#
# The project is organised with       #
# split includes and source folders   #
# this makes it easier for developers #
# to quickly find relevant includes.  #
#######################################

file(GLOB_RECURSE APP_SOURCES "src/*.cpp")             # source reside in src/
file(GLOB_RECURSE APP_HEADERS "include/*.h" "src/*.h") # includes reside in include/ and src/

set (APP_INCLUDE_DIRS "")
foreach (_headerFile ${APP_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND APP_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES APP_INCLUDE_DIRS)
include_directories(${APP_INCLUDE_DIRS})

# Build CoolProp as a shared or static library depending on build option
if (COOLPROP_STATIC_LIBRARY)
    add_library(${app_name} STATIC ${APP_SOURCES})
else()
    add_library(${app_name} SHARED ${APP_SOURCES})
endif()

target_link_libraries(${app_name} ${CMAKE_DL_LIBS})

#######################################
#         REQUIRED MODULES            #
#-------------------------------------#
# CoolProp requires some standard OS  #
# features, these include:            #
# DL (CMAKE_DL_LIBS) for REFPROP      #
#######################################
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/dev/cmake/modules/")

set(Python_ADDITIONAL_VERSIONS 2.7 2.6 2.5 2.4)
find_package (PythonInterp REQUIRED)

ADD_CUSTOM_COMMAND (TARGET ${app_name} PRE_BUILD 
                    COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/dev/JSON_to_CPP.py" "--path" ${CMAKE_CURRENT_SOURCE_DIR}
                    COMMENT "Information string for prebuild execution"                     
) 



# TODO: check relevance of http://www.cmake.org/Wiki/BuildingWinDLL



#include_directories("${CMAKE_CURRENT_SOURCE_DIR}/CoolProp")
#FILE(GLOB coolprop_files "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/*.cpp")
#add_library(coolprop STATIC ${coolprop_files})
