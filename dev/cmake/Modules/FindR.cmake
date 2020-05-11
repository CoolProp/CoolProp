# - Find R
# This module finds the libraries corresponding to the R program
# It is based on the code from http://github.com/usnistgov/REFPROP-wrappers which is in the public domain

# This code sets the following variables:
#
#  R_LIBRARY             - path to the R shared library
#  R_INCLUDE_DIRS        - include directory(ies)
#  R_BIN_DIR             - directory that contains the R executable being used

SET(R_FOUND 0)

IF( "${R_BIN}" STREQUAL "")
    MESSAGE(FATAL_ERROR "You must pass the CMake variable R_BIN pointing to the directory that contains your executable")
ENDIF()

find_program(R_EXEC
             NAMES R
             PATHS ${R_BIN}
             NO_SYSTEM_ENVIRONMENT_PATH
             NO_CMAKE_SYSTEM_PATH
             NO_DEFAULT_PATH
             NO_CMAKE_PATH
             )
MESSAGE(STATUS "R_EXEC= ${R_EXEC}")
# Parse the output of the R path command, removing whitespace
FUNCTION(chomp arg1 arg2)
    string(REPLACE "\n" ";" arg1list ${arg1})
    list(GET arg1list 1 arg1line)
    string(LENGTH ${arg1line} len)
    MATH(EXPR length ${len}-6)
    string(SUBSTRING ${arg1line} 5 ${length} arg1path)
    # see http://www.cmake.org/pipermail/cmake/2008-November/025423.html
    set(${arg2} ${arg1path} PARENT_SCOPE)
ENDFUNCTION(chomp)

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/get_home.R R.home())
execute_process(
    COMMAND ${R_EXEC} --quiet -f ${CMAKE_CURRENT_BINARY_DIR}/get_home.R OUTPUT_VARIABLE R_HOME_TEXT RESULT_VARIABLE R_HOME_RESULT
)
file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/get_home.R)
MESSAGE(STATUS "R_HOME_TEXT = ${R_HOME_TEXT} w/ RESULT=${R_HOME_RESULT}")
chomp(${R_HOME_TEXT} R_HOME_TEXT)

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/get_include.R R.home(component=\"include\"))
execute_process(
    COMMAND ${R_EXEC} --quiet -f ${CMAKE_CURRENT_BINARY_DIR}/get_include.R OUTPUT_VARIABLE R_INCLUDE_TEXT RESULT_VARIABLE R_INCLUDE_RESULT
)
file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/get_include.R)
chomp(${R_INCLUDE_TEXT} R_INCLUDE_DIRS)
MESSAGE(STATUS "R_INCLUDE = ${R_INCLUDE_DIRS} w/ RESULT=${R_INCLUDE_RESULT}")

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/get_bin.R R.home(component=\"bin\"))
execute_process(
    COMMAND ${R_EXEC} --quiet -f ${CMAKE_CURRENT_BINARY_DIR}/get_bin.R OUTPUT_VARIABLE R_BIN_TEXT RESULT_VARIABLE R_BIN_RESULT
)
file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/get_bin.R)
chomp(${R_BIN_TEXT} R_BIN_DIR)
MESSAGE(STATUS "R_BIN_TEXT = ${R_BIN_DIR} w/ RESULT=${R_BIN_RESULT}")

if (WIN32)
# Bug in cmake 3.17, need to search for file rather than library
find_file(
          R_LIBRARY
          NAMES R.dll
          PATHS ${R_BIN_DIR}
          NO_DEFAULT_PATH
)
else()
set(R_PATH_LOC "${R_HOME_TEXT}/lib")
MESSAGE(STATUS "R_PATH_LOC = ${R_PATH_LOC}")
find_library(
          R_LIBRARY
          NAMES R
          PATHS ${R_PATH_LOC}
          NO_DEFAULT_PATH
          #NO_CMAKE_PATH
          )
endif()
MESSAGE(STATUS "R_LIBRARY = ${R_LIBRARY}")

find_package_handle_standard_args(R DEFAULT_MSG
                                  R_LIBRARY R_INCLUDE_DIRS)

MARK_AS_ADVANCED(
  R_BIN_DIR
  R_LIBRARY
  R_INCLUDE_DIRS
)