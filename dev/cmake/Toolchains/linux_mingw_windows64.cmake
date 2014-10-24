# run with cmake -DCMAKE_TOOLCHAIN_FILE=...
# the name of the target operating system
SET(CMAKE_SYSTEM_NAME Windows)

#######################################
# Choose an appropriate compiler prefix
#######################################
#
# for classical mingw32
# see http://www.mingw.org/
#set(COMPILER_PREFIX "i586-mingw32msvc")
#
# for 32 or 64 bits mingw-w64
# see http://mingw-w64.sourceforge.net/
#set(COMPILER_PREFIX   "i686-w64-mingw32")
#set(COMPILER_PREFIX "x86_64-w64-mingw32")
#
set(COMPILER_PREFIX "x86_64-w64-mingw32")
#
#######################################
# which compilers to use for C and C++
#######################################
#
# Either search for the compiler ...
#find_program(CMAKE_RC_COMPILER NAMES ${COMPILER_PREFIX}-windres)
# 
# ... or rely on the system path
#SET(CMAKE_RC_COMPILER ${COMPILER_PREFIX}-windres)
#
#find_program(CMAKE_RC_COMPILER NAMES ${COMPILER_PREFIX}-windres)
#find_program(CMAKE_C_COMPILER NAMES ${COMPILER_PREFIX}-gcc)
#find_program(CMAKE_CXX_COMPILER NAMES ${COMPILER_PREFIX}-g++)
SET(CMAKE_RC_COMPILER ${COMPILER_PREFIX}-windres)
SET(CMAKE_C_COMPILER ${COMPILER_PREFIX}-gcc)
SET(CMAKE_CXX_COMPILER ${COMPILER_PREFIX}-g++)
#
#
#######################################
# Where to look for target binaries
#######################################
#
# here is the target environment located
set(CMAKE_FIND_ROOT_PATH /usr/${COMPILER_PREFIX})
#
# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search 
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
#