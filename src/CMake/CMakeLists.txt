cmake_minimum_required (VERSION 2.8)
project (CoolProp)

# Offer the user the choice of overriding the directories
set (PROJECT_SRC_DIR src CACHE PATH "Project directory for source files")
set (PROJECT_BIN_DIR bin CACHE PATH "Build directory for executables")
set (PROJECT_LIB_DIR lib CACHE PATH "Build directory for libraries")
set (PROJECT_INC_DIR include CACHE PATH "Build directory for header files")

# The version number.
set (CoolProp_VERSION_MAJOR 5)
set (CoolProp_VERSION_MINOR 0)
set (CoolProp_VERSION_PATCH 0)
set (CoolProp_VERSION ${CoolProp_VERSION_MAJOR}.${CoolProp_VERSION_MINOR}.${CoolProp_VERSION_PATCH})

FIND_PACKAGE(Git)
EXECUTE_PROCESS(
     COMMAND ${GIT_EXECUTABLE} describe --match=CatchMeIfYouCan --always --abbrev=40
     OUTPUT_VARIABLE CoolProp_REVISION
     OUTPUT_STRIP_TRAILING_WHITESPACE
 )

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SRC_DIR}/CoolPropConfig.h.in"
  "${PROJECT_SRC_DIR}/CoolPropConfig.h"
  )
 
# add the binary tree to the search path for include files
# so that we will find CoolPropConfig.h
include_directories("${PROJECT_SRC_DIR}" "${PROJECT_BIN_DIR}" "${PROJECT_LIB_DIR}" "${PROJECT_INC_DIR}")

# add the executable
add_executable(CoolProp ${PROJECT_SRC_DIR}/CoolPropMain.cpp)
