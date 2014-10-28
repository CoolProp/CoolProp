cmake_minimum_required(VERSION 2.8)


if (DEFINED COOLPROP_INSTALL_PREFIX)
    #set(COOLPROP_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install_root)
    message(STATUS "COOLPROP_INSTALL_PREFIX=${COOLPROP_INSTALL_PREFIX}")
else()
    set(COOLPROP_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install_root)
    message(STATUS "COOLPROP_INSTALL_PREFIX=${COOLPROP_INSTALL_PREFIX}")
endif()

set(CMAKE_INSTALL_PREFIX ${COOLPROP_INSTALL_PREFIX} CACHE PATH "default install path" FORCE)  
  
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
set (COOLPROP_VERSION_MAJOR 5)
set (COOLPROP_VERSION_MINOR 0)
set (COOLPROP_VERSION_PATCH 0)
set (COOLPROP_VERSION ${COOLPROP_VERSION_MAJOR}.${COOLPROP_VERSION_MINOR}.${COOLPROP_VERSION_PATCH})

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
        
option (BUILD_TESTING 
        "Enable testing for this given builder" 
        OFF)

option (FORCE_BITNESS_32
       "Force a 32bit build regardless of the host"
       OFF)

option (FORCE_BITNESS_64
       "Force a 64bit build regardless of the host"
       OFF)
        
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

file(GLOB_RECURSE APP_SOURCES "src/*.cpp")
file(GLOB_RECURSE APP_HEADERS "include/*.h" "src/*.h")

## You can exclude this file, in case you want to run your own tests or use Catch
list(REMOVE_ITEM  APP_SOURCES "${CMAKE_SOURCE_DIR}/src/Tests/Tests.cpp")
list(REMOVE_ITEM  APP_SOURCES "${CMAKE_SOURCE_DIR}/src/Tests/CoolProp-Tests.cpp")

## This file is only needed for the library, normal builds do not need it.
list(REMOVE_ITEM  APP_SOURCES "${CMAKE_SOURCE_DIR}/src/CoolPropLib.cpp")


set (APP_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/externals/Eigen")
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
#                 BITNESS             #
#-------------------------------------#
#        Calculate if 32 or 64        #
#######################################

if(WIN32)
  if (CMAKE_CL_64)
    SET(BITNESS "64")
  else()
    SET(BITNESS "32")
  endif()
else()
  if (CMAKE_SIZEOF_VOID_P MATCHES "8")
    SET(BITNESS "64")
  else()
    SET(BITNESS "32")
  endif()
endif()

IF(MSVC AND (FORCE_BITNESS_32 OR FORCE_BITNESS_64))
  MESSAGE(STATUS "You cannot force a certain bitness for Visual Studio, use the generator settings for this purpose." )
  MESSAGE(STATUS "Pass '-G \"Visual Studio 10 2010 Win64\"' to CMake to make a 64bit binary using VS2010." )
  MESSAGE(STATUS "Pass '-G \"Visual Studio 10 2010\"' to CMake to make a 32bit binary using VS2010." )
  MESSAGE(STATUS "Pass '-G \"Visual Studio 9 2008 Win64\"' to CMake to make a 64bit binary using VS2008." )
  MESSAGE(STATUS "Pass '-G \"Visual Studio 9 2008\"' to CMake to make a 32bit binary using VS2008." )
  MESSAGE(FATAL_ERROR "Fix that and try again...")  
ENDIF()

IF(FORCE_BITNESS_32)
  SET(BITNESS "32")
ELSEIF(FORCE_BITNESS_64)
  SET(BITNESS "64")
ENDIF()

if (!MSVC)
  #set_target_properties(${app_name} PROPERTIES COMPILE_FLAGS "-m${BITNESS}" LINK_FLAGS "-m${BITNESS}")
  #set_target_properties(${app_name} PROPERTIES                        COMPILE_FLAGS "-m${BITNESS}" LINK_FLAGS "-m${BITNESS}")
  set_target_properties(${app_name} PROPERTIES APPEND_STRING PROPERTY COMPILE_FLAGS "-m${BITNESS}")
  set_target_properties(${app_name} PROPERTIES APPEND_STRING PROPERTY    LINK_FLAGS "-m${BITNESS}") 
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
if (COOLPROP_SHARED_LIBRARY)
    list(APPEND APP_SOURCES "${CMAKE_SOURCE_DIR}/src/CoolPropLib.cpp")
endif()

if (COOLPROP_STATIC_LIBRARY_MODULE OR COOLPROP_STATIC_LIBRARY OR COOLPROP_EXTERNC_STATIC_LIBRARY)
    list(APPEND APP_SOURCES "${CMAKE_SOURCE_DIR}/src/CoolPropLib.cpp")
    add_library(${app_name} STATIC ${APP_SOURCES})
    add_dependencies (${app_name} generate_headers)
    install (TARGETS ${app_name} DESTINATION static_library/${CMAKE_SYSTEM_NAME})
endif()

###      Process some more options  ###

if (COOLPROP_EXTERNC_STATIC_LIBRARY)
  set_target_properties (${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DEXTERNC")
endif()


if (COOLPROP_64BIT_SHARED_LIBRARY_MODULE OR COOLPROP_64BIT_SHARED_LIBRARY)
    list(APPEND APP_SOURCES "${CMAKE_SOURCE_DIR}/src/CoolPropLib.cpp")
    add_library(${app_name} SHARED ${APP_SOURCES})
    set_target_properties (${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB")
    if (!MSVC)
      set_target_properties(${app_name} PROPERTIES COMPILE_FLAGS "-m64" LINK_FLAGS "-m64")
    endif()
    add_dependencies (${app_name} generate_headers)
    
    if (MSVC)
      
      # No lib prefix for the shared library
      set_target_properties(${app_name} PROPERTIES PREFIX "")
      add_custom_command(TARGET ${app_name}
                         POST_BUILD
                         COMMAND dumpbin /EXPORTS $<TARGET_FILE:CoolProp> > ${CMAKE_CURRENT_BINARY_DIR}/exports.txt)
      install (FILES ${CMAKE_CURRENT_BINARY_DIR}/exports.txt DESTINATION shared_library/${CMAKE_SYSTEM_NAME}/64bit)
    endif()
    install (TARGETS ${app_name} DESTINATION shared_library/${CMAKE_SYSTEM_NAME}/64bit)
    install (FILES ${CMAKE_SOURCE_DIR}/include/CoolPropLib.h DESTINATION shared_library)
endif()

if (COOLPROP_DEBIAN_PACKAGE)
    if(!UNIX)
      message(FATAL_ERROR "COOLPROP_DEBIAN_PACKAGE can only be used on linux host")
    endif()
    list(APPEND APP_SOURCES "${CMAKE_SOURCE_DIR}/src/CoolPropLib.cpp")
    add_library(${app_name} SHARED ${APP_SOURCES})
    set_target_properties (${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB")
    set_target_properties (${app_name} PROPERTIES VERSION ${COOLPROP_VERSION} SOVERSION ${COOLPROP_VERSION_MAJOR})
    add_dependencies (${app_name} generate_headers)
    install (TARGETS ${app_name} DESTINATION "${CMAKE_INSTALL_PREFIX}/usr/lib")
    install (FILES ${CMAKE_SOURCE_DIR}/include/CoolPropLib.h DESTINATION "${CMAKE_INSTALL_PREFIX}/usr/include")
endif()

if (COOLPROP_32BIT_SHARED_LIBRARY_LINUX_MODULE OR COOLPROP_32BIT_SHARED_LIBRARY_LINUX)
    if(!UNIX)
      message(FATAL_ERROR "COOLPROP_32BIT_SHARED_LIBRARY_LINUX_MODULE can only be used on linux host")
    endif()
    list(APPEND APP_SOURCES "${CMAKE_SOURCE_DIR}/src/CoolPropLib.cpp")
    add_library(${app_name} SHARED ${APP_SOURCES})
    set_target_properties (${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DEXTERNC -DCONVENTION= -m32" LINK_FLAGS "-m32")
    add_dependencies (${app_name} generate_headers)
    install (TARGETS ${app_name} DESTINATION shared_library/${CMAKE_SYSTEM_NAME}/32bit__cdecl_calling_convention)
    install (FILES ${CMAKE_SOURCE_DIR}/include/CoolPropLib.h DESTINATION shared_library)
endif()

if (COOLPROP_32BIT_CDECL_SHARED_LIBRARY_MODULE OR COOLPROP_32BIT_CDECL_SHARED_LIBRARY)
    list(APPEND APP_SOURCES "${CMAKE_SOURCE_DIR}/src/CoolPropLib.cpp")
    add_library(${app_name} SHARED ${APP_SOURCES})
    set_target_properties (${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB -DCONVENTION=__cdecl")
    if (!MSVC)
      set_target_properties(${app_name} PROPERTIES COMPILE_FLAGS "-m32" LINK_FLAGS "-m32")
    endif()
    add_dependencies (${app_name} generate_headers)
    set_target_properties(${app_name} PROPERTIES PREFIX "")
    if (MSVC)
      add_custom_command(TARGET ${app_name}
                         POST_BUILD
                         COMMAND dumpbin /EXPORTS $<TARGET_FILE:CoolProp> > ${CMAKE_CURRENT_BINARY_DIR}/exports.txt)
      install (FILES ${CMAKE_CURRENT_BINARY_DIR}/exports.txt DESTINATION shared_library/${CMAKE_SYSTEM_NAME}/32bit__cdecl_calling_convention)
    endif()
    install (TARGETS ${app_name} DESTINATION shared_library/${CMAKE_SYSTEM_NAME}/32bit__cdecl_calling_convention)
    install (FILES ${CMAKE_SOURCE_DIR}/include/CoolPropLib.h DESTINATION shared_library)
endif()

if (COOLPROP_32BIT_STDCALL_SHARED_LIBRARY_MODULE OR COOLPROP_32BIT_STDCALL_SHARED_LIBRARY)
    list(APPEND APP_SOURCES "${CMAKE_SOURCE_DIR}/src/CoolPropLib.cpp")
    add_library(${app_name} SHARED ${APP_SOURCES})
    set_target_properties (${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB -DCONVENTION=__stdcall")
    if (!MSVC)
      set_target_properties(${app_name} PROPERTIES COMPILE_FLAGS "-m32" LINK_FLAGS "-m32")
    endif()
    add_dependencies (${app_name} generate_headers)
    set_target_properties(${app_name} PROPERTIES PREFIX "")
    if (MSVC)
      add_custom_command(TARGET ${app_name}
                         POST_BUILD
                         COMMAND dumpbin /EXPORTS $<TARGET_FILE:CoolProp> > ${CMAKE_CURRENT_BINARY_DIR}/exports.txt)
      install (FILES ${CMAKE_CURRENT_BINARY_DIR}/exports.txt DESTINATION shared_library/${CMAKE_SYSTEM_NAME}/32bit__stdcall_calling_convention)
    endif()
    install (TARGETS ${app_name} DESTINATION shared_library/${CMAKE_SYSTEM_NAME}/32bit__stdcall_calling_convention)
    install (FILES ${CMAKE_SOURCE_DIR}/include/CoolPropLib.h DESTINATION shared_library)
endif()

# EES is only compiled for windows
if (COOLPROP_EES_MODULE)
  list (APPEND APP_SOURCES "wrappers/EES/main.cpp")
  list(APPEND APP_SOURCES "${CMAKE_SOURCE_DIR}/src/CoolPropLib.cpp")
  include_directories(${APP_INCLUDE_DIRS})
  add_library(COOLPROP_EES SHARED ${APP_SOURCES})
  set_target_properties (COOLPROP_EES PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB -DCONVENTION=__cdecl")
  add_dependencies (COOLPROP_EES generate_headers)
  set_target_properties(COOLPROP_EES PROPERTIES SUFFIX ".dlf" PREFIX "")
  message(STATUS "Injecting the version COOLPROP_VERSION=${COOLPROP_VERSION}")
  # Put the version into the InnoInstaller setup file
  CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/wrappers/EES/BuildInnoInstaller.iss.in" "${CMAKE_CURRENT_BINARY_DIR}/BuildInnoInstaller.iss")
  if (!MSVC)
    set_target_properties(COOLPROP_EES PROPERTIES COMPILE_FLAGS "-m32" LINK_FLAGS "-m32")
  endif()
  add_dependencies (COOLPROP_EES generate_headers)
  
  if ( MSVC )
    set_target_properties( COOLPROP_EES PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
    set_target_properties( COOLPROP_EES PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_BINARY_DIR})
    set_target_properties( COOLPROP_EES PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR})
    # etc for the other available configuration types (MinSizeRel, RelWithDebInfo)
  endif ()
  # copy required files
  set(REQUIRED_FILES
    CoolProp.htm
    CoolProp.LIB
    CoolProp_EES_Sample.EES
    )
    foreach (_FILE ${REQUIRED_FILES})
        add_custom_command(TARGET COOLPROP_EES
                           PRE_BUILD
                           COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/wrappers/EES/${_FILE} ${CMAKE_CURRENT_BINARY_DIR}/.)
    endforeach()
  
  # Run InnoSetup to make the installer (must be on your path)
  add_custom_command(TARGET COOLPROP_EES
                     POST_BUILD
                     COMMAND iscc /cc BuildInnoInstaller.iss
                     WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  # install the generated executable installer from InnoSetup
  install(
    CODE "file( INSTALL ${CMAKE_CURRENT_BINARY_DIR}/Output/SetupCOOLPROP_EES.exe DESTINATION ${CMAKE_INSTALL_PREFIX}/EES/${CMAKE_SYSTEM_NAME} )"
  )
endif()

if (COOLPROP_OCTAVE_MODULE)
  
  # Must have SWIG and Octave
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE(${SWIG_USE_FILE})
  FIND_PACKAGE(Octave REQUIRED)
  
  # Make a src directory to deal with file permissions problem with MinGW makefile
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src)
  
  # Set the include folders
  SET(OCTAVE_WRAP_INCLUDE_DIRS ${INCLUDE_DIR})
  foreach(ITR ${OCTAVE_INCLUDE_DIRS})
    list(APPEND OCTAVE_WRAP_INCLUDE_DIRS ${ITR})
  endforeach()
  include_directories(${OCTAVE_INCLUDE_DIRS})

  set(I_FILE "${CMAKE_SOURCE_DIR}/src/CoolProp.i")

  SET_SOURCE_FILES_PROPERTIES(${I_FILE} PROPERTIES CPLUSPLUS ON)
  
  SET(SWIG_MODULE_CoolProp_EXTRA_DEPS
    ${CMAKE_SOURCE_DIR}/include/DataStructures.h
    ${CMAKE_SOURCE_DIR}/include/CoolProp.h
    ${CMAKE_SOURCE_DIR}/include/AbstractState.h
  ) 
  SWIG_ADD_MODULE(CoolProp octave ${I_FILE} ${APP_SOURCES})
  
  if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # We need to see which library to link with on OSX - clang++ or stdc++
    execute_process(COMMAND "otool -L ${OCTAVE_OCTAVE_LIB} | grep libc++" OUTPUT_VARIABLE COOLPROP_OCTAVE_USING_CLANG)
    if ("${COOLPROP_OCTAVE_USING_CLANG} " STREQUAL " ")
      message(STATUS "Using -stdlib=libc++")
      set_target_properties(CoolProp PROPERTIES LINK_FLAGS "-stdlib=libc++")
    else()
      message(STATUS "Using -stdlib=libstdc++")
    endif()
  endif()
  
  if (WIN32)
    include_directories($ENV{OCTAVE_ROOT}/include)
    include_directories($ENV{OCTAVE_ROOT}/include/octave-${OCTAVE_VERSION}/octave)
    set_target_properties(CoolProp PROPERTIES COMPILE_FLAGS "-fpermissive")
    SWIG_LINK_LIBRARIES(CoolProp octave octinterp)
    set_target_properties(CoolProp PROPERTIES LINK_FLAGS "-L$ENV{OCTAVE_ROOT}/lib/octave/${OCTAVE_VERSION}")
  else()
    SWIG_LINK_LIBRARIES(CoolProp ${OCTAVE_LIBRARIES})  
  endif()

  set_target_properties(CoolProp PROPERTIES SUFFIX ".oct" PREFIX "")
  add_dependencies (${app_name} generate_headers)  
  
  install (FILES ${CMAKE_SOURCE_DIR}/wrappers/Octave/Example.m DESTINATION Octave)
  install (TARGETS ${app_name} DESTINATION Octave/Octave${OCTAVE_VERSION}_${CMAKE_SYSTEM_NAME}_${BITNESS}bit)
endif()

if (COOLPROP_CSHARP_MODULE)

  # Must have SWIG and C#
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE(${SWIG_USE_FILE})
  FIND_PACKAGE(Csharp REQUIRED)
  
  # Make a src directory to deal with file permissions problem with MinGW makefile
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src)

  set(I_FILE "${CMAKE_SOURCE_DIR}/src/CoolProp.i")
	
  SET(SWIG_MODULE_CoolProp_EXTRA_DEPS
    ${CMAKE_SOURCE_DIR}/include/DataStructures.h
    ${CMAKE_SOURCE_DIR}/include/CoolProp.h
    ${CMAKE_SOURCE_DIR}/include/AbstractState.h
  )     
  # Set properties before adding module
  SET_SOURCE_FILES_PROPERTIES(${I_FILE} PROPERTIES CPLUSPLUS ON)
  if (WIN32)
    SET(CMAKE_SWIG_FLAGS -dllimport \"CoolProp\")
    #SET_SOURCE_FILES_PROPERTIES(${I_FILE} PROPERTIES SWIG_FLAGS "-dllimport CoolProp")
  endif()
  SWIG_ADD_MODULE(CoolProp csharp ${I_FILE} ${APP_SOURCES})

  if (WIN32)
    set_target_properties(CoolProp PROPERTIES PREFIX "")
  endif()
  
  add_dependencies (${app_name} generate_headers)  
  
  add_custom_command(TARGET CoolProp
                     POST_BUILD
                     COMMAND 7z a "${CMAKE_CURRENT_BINARY_DIR}/platform-independent.7z" "${CMAKE_CURRENT_BINARY_DIR}/*.cs"
                     WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  
  # Install all the generated cs files
  install(
    CODE "file( GLOB _GeneratedCsharpSources \"${CMAKE_CURRENT_BINARY_DIR}/*.cs\" )"
    CODE "file( INSTALL \${_GeneratedCsharpSources} DESTINATION ${CMAKE_INSTALL_PREFIX}/Csharp/platform-independent )"
  )
  install (FILES ${CMAKE_SOURCE_DIR}/Web/coolprop/wrappers/Csharp/Example.cs DESTINATION Csharp)
  install (FILES "${CMAKE_CURRENT_BINARY_DIR}/platform-independent.7z" DESTINATION ${CMAKE_INSTALL_PREFIX}/Csharp)
  install (TARGETS ${app_name} DESTINATION Csharp/${CMAKE_SYSTEM_NAME}_${BITNESS}bit)
  enable_testing()
  if (DEFINED BUILD_TESTING)
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/testing_root/Csharp${BITNESS})
    # Copy the shared object to the folder with the executable - no idea like java.library.path in C#
    install (TARGETS ${app_name} DESTINATION ${CMAKE_SOURCE_DIR}/testing_root/Csharp${BITNESS})
  endif()
  file(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX}/Csharp/platform-independent/*.cs cp_cs_path)
  file(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX}/Csharp/Example.cs cp_example_path)
  if (${BITNESS} EQUAL "32")
    set(CSHARP_PLAT "-platform:x86")
  elseif((${BITNESS} EQUAL "64"))
    set(CSHARP_PLAT "-platform:x64")
  endif()
  add_test(NAME Csharptestbuild
           COMMAND ${CSHARP_COMPILER} -out:Example.exe ${CSHARP_PLAT} ${cp_cs_path} ${cp_example_path}
           WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing_root/Csharp${BITNESS})
  add_test(NAME Csharptestrun
           COMMAND ${CSHARP_INTERPRETER} Example.exe
           WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing_root/Csharp${BITNESS})
endif()

if (COOLPROP_JAVA_MODULE)

  # Must have SWIG and Java
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE(${SWIG_USE_FILE})
  find_package(Java REQUIRED)
  find_package(JNI)

  # Make a src directory to deal with file permissions problem with MinGW makefile
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src)
  
  message(STATUS "JAVA_INCLUDE_PATH = ${JAVA_INCLUDE_PATH}")

  include_directories(${JAVA_INCLUDE_PATH})
  include_directories(${JAVA_INCLUDE_PATH}/win32)
  include_directories(${JAVA_INCLUDE_PATH}/linux)

  set(I_FILE "${CMAKE_SOURCE_DIR}/src/CoolProp.i")

  SET_SOURCE_FILES_PROPERTIES(${I_FILE} PROPERTIES CPLUSPLUS ON)

  SET(SWIG_MODULE_CoolProp_EXTRA_DEPS
    ${CMAKE_SOURCE_DIR}/include/DataStructures.h
    ${CMAKE_SOURCE_DIR}/include/CoolProp.h
    ${CMAKE_SOURCE_DIR}/include/AbstractState.h
  ) 
  SWIG_ADD_MODULE(CoolProp java ${I_FILE} ${APP_SOURCES})

  if (WIN32)
    set_target_properties(CoolProp PROPERTIES PREFIX "")
  endif()
  
  if (!MSVC)
    set_target_properties(COOLPROP_EES PROPERTIES COMPILE_FLAGS "-m64" LINK_FLAGS "-m64")
  endif()
  add_dependencies (${app_name} generate_headers)
  
  add_custom_command(TARGET CoolProp
                     POST_BUILD
                     COMMAND 7z a "${CMAKE_CURRENT_BINARY_DIR}/platform-independent.7z" "${CMAKE_CURRENT_BINARY_DIR}/*.java"
                     WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  
  # Install all the generated java files
  install(
    CODE "file( GLOB _GeneratedJavaSources \"${CMAKE_CURRENT_BINARY_DIR}/*.java\" )"
    CODE "file( INSTALL \${_GeneratedJavaSources} DESTINATION ${CMAKE_INSTALL_PREFIX}/Java/platform-independent )"
  )
  install (FILES ${CMAKE_SOURCE_DIR}/Web/coolprop/wrappers/Java/Example.java DESTINATION ${CMAKE_INSTALL_PREFIX}/Java)
  install (FILES "${CMAKE_CURRENT_BINARY_DIR}/platform-independent.7z" DESTINATION ${CMAKE_INSTALL_PREFIX}/Java)
  install (TARGETS ${app_name} DESTINATION ${CMAKE_INSTALL_PREFIX}/Java/${CMAKE_SYSTEM_NAME}_${BITNESS}bit)
  enable_testing()
  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/testing_root/Java${BITNESS})
  add_test(NAME Javatestbuild
           COMMAND javac -d . ${CMAKE_INSTALL_PREFIX}/Java/Example.java -cp ${CMAKE_INSTALL_PREFIX}/Java/platform-independent
           WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing_root/Java${BITNESS})
  add_test(NAME Javatestrun
           COMMAND ${Java_JAVA_EXECUTABLE} -Djava.library.path=${CMAKE_INSTALL_PREFIX}/Java/${CMAKE_SYSTEM_NAME}_${BITNESS}bit Example
           WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing_root/Java${BITNESS})
endif()

if (COOLPROP_PHP_MODULE)

  # Must have SWIG
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE(${SWIG_USE_FILE})

  execute_process(
    COMMAND php-config --includes
    OUTPUT_VARIABLE php_config_includes
    RESULT_VARIABLE php_config_failed
    )
  string(STRIP ${php_config_includes} php_config_includes)
  string(REPLACE "-I" "" PHP_INCLUDES ${php_config_includes})
  SEPARATE_ARGUMENTS(PHP_INCLUDES)

  message(STATUS "php includes=${PHP_INCLUDES}")
  include_directories(${PHP_INCLUDES})

  set(I_FILE "${CMAKE_SOURCE_DIR}/src/CoolProp.i")

  SET_SOURCE_FILES_PROPERTIES(${I_FILE} PROPERTIES CPLUSPLUS ON)

  SET(SWIG_MODULE_CoolProp_EXTRA_DEPS
    ${CMAKE_SOURCE_DIR}/include/DataStructures.h
    ${CMAKE_SOURCE_DIR}/include/CoolProp.h
    ${CMAKE_SOURCE_DIR}/include/AbstractState.h
  ) 
  SWIG_ADD_MODULE(CoolProp php ${I_FILE} ${APP_SOURCES})

  if (WIN32)
    set_target_properties(CoolProp PROPERTIES PREFIX "")
  endif()
  
  if (!MSVC)
    set_target_properties(CoolProp PROPERTIES COMPILE_FLAGS "-m64" LINK_FLAGS "-m64")
  endif()
  add_dependencies (CoolProp generate_headers)
  
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/CoolProp.php DESTINATION ${CMAKE_INSTALL_PREFIX}/PHP/cross-platform)
  install (TARGETS ${app_name} DESTINATION ${CMAKE_INSTALL_PREFIX}/PHP/${CMAKE_SYSTEM_NAME})

endif()

function(JOIN VALUES GLUE OUTPUT)
  string (REGEX REPLACE "([^\\]|^);" "\\1${GLUE}" _TMP_STR "${VALUES}")
  string (REGEX REPLACE "[\\](.)" "\\1" _TMP_STR "${_TMP_STR}") #fixes escaping
  set (${OUTPUT} "${_TMP_STR}" PARENT_SCOPE)
endfunction()

if (COOLPROP_SCILAB_SWIG_MODULE)
  
  # Must have SWIG
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE(${SWIG_USE_FILE})

  IF( "$ENV{SCILAB_ROOT}" STREQUAL "" )
      MESSAGE(STATUS "SCILAB_ROOT environment variable not set." )
      MESSAGE(STATUS "In Linux this can be done in your user .bashrc file by appending the corresponding line, e.g:" )
      MESSAGE(STATUS "export SCILAB_ROOT=/usr/local/SCILAB" )
      MESSAGE(STATUS "In Windows this can be done by adding system variable, e.g:" )
      MESSAGE(STATUS "SCILAB_ROOT=D:\\Program Files\\SCILAB" )
  	MESSAGE(FATAL_ERROR "Quitting.")
  endif()

  SET(SWIG_MODULE_CoolProp_EXTRA_DEPS
    ${CMAKE_SOURCE_DIR}/include/DataStructures.h
    ${CMAKE_SOURCE_DIR}/include/CoolProp.h
    ${CMAKE_SOURCE_DIR}/include/AbstractState.h
  ) 
  set(SCILAB_INCLUDE_PATH "$ENV{SCILAB_ROOT}/include")
    	
  set(I_FILE "${CMAKE_SOURCE_DIR}/src/CoolProp.i")

  include_directories(${SCILAB_INCLUDE_PATH})  
  include_directories(${SCILAB_INCLUDE_PATH}/scilab)

  set(SWIG_OPTIONS "-nobuilder")
  SET_SOURCE_FILES_PROPERTIES(${I_FILE} PROPERTIES SWIG_FLAGS "${SWIG_OPTIONS}")
  SET_SOURCE_FILES_PROPERTIES(${I_FILE} PROPERTIES CPLUSPLUS ON)

  list(APPEND APP_SOURCES gw_CoolProp.c)
  set_source_files_properties(gw_CoolProp.c PROPERTIES GENERATED TRUE)

  SWIG_ADD_MODULE(CoolProp scilab ${I_FILE} ${APP_SOURCES})

  add_dependencies (CoolProp generate_headers)

endif()

if (COOLPROP_MATLAB_SWIG_MODULE)
  
  # Must have SWIG
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE(${SWIG_USE_FILE})

  find_package(Matlab REQUIRED)

  IF(MATLAB_FOUND)
      message(STATUS "MATLAB Found, MATLAB MEX will be compiled.")
  ELSE(MATLAB_FOUND)
      MESSAGE("MATLAB not found...nothing will be built.")
  ENDIF(MATLAB_FOUND)
  
  set(I_FILE "${CMAKE_SOURCE_DIR}/src/CoolProp.i")

  list (APPEND APP_SOURCES ${CMAKE_SOURCE_DIR}/wrappers/MATLAB/Matlabdef.def) # To export mexFunction

  SET_SOURCE_FILES_PROPERTIES(${I_FILE} PROPERTIES CPLUSPLUS ON)

  SWIG_ADD_MODULE(CoolPropMATLAB_wrap matlab ${I_FILE} ${APP_SOURCES})
  SWIG_LINK_LIBRARIES(CoolPropMATLAB_wrap ${MATLAB_LIBRARIES})
  
  add_definitions(/DMATLAB_MEX_FILE) #define matlab macros
  add_definitions(/DMX_COMPAT_32)

  if(WIN32) # 32-bit or 64-bit mex
    if (CMAKE_CL_64)
        SET_TARGET_PROPERTIES(CoolPropMATLAB_wrap PROPERTIES PREFIX "" SUFFIX .mexw64)
    else()
        SET_TARGET_PROPERTIES(CoolPropMATLAB_wrap PROPERTIES SUFFIX .mexw32)
    endif()
  else()
    if (CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET_TARGET_PROPERTIES(CoolPropMATLAB_wrap PROPERTIES PREFIX "" SUFFIX .mexa64 PREFIX "")
    else()
        SET_TARGET_PROPERTIES(CoolPropMATLAB_wrap PROPERTIES PREFIX "" SUFFIX .mexglx PREFIX "")
    endif()
  endif()
  add_dependencies (CoolPropMATLAB_wrap generate_headers)
endif()

if (COOLPROP_MATLAB_MODULE)

  add_definitions(/DMATLAB_MEX_FILE) #define matlab macros
  add_definitions(/DMX_COMPAT_32)

  find_package(Matlab REQUIRED)

  IF(MATLAB_FOUND)
      message(STATUS "MATLAB Found, MATLAB MEX will be compiled.")
  ELSE(MATLAB_FOUND)
      MESSAGE("MATLAB not found...nothing will be built.")
  ENDIF(MATLAB_FOUND)

  # set up matlab libraries
  INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIR})
  add_library(PropsSI SHARED ${APP_SOURCES} ${CMAKE_SOURCE_DIR}/wrappers/MATLAB/Matlabdef.def ${CMAKE_SOURCE_DIR}/wrappers/MATLAB/PropsSI.cpp)
  add_library(HAPropsSI SHARED ${APP_SOURCES} ${CMAKE_SOURCE_DIR}/wrappers/MATLAB/Matlabdef.def ${CMAKE_SOURCE_DIR}/wrappers/MATLAB/HAProps.cpp)
  target_link_libraries(PropsSI ${MATLAB_LIBRARIES})
  target_link_libraries(HAPropsSI ${MATLAB_LIBRARIES})
 
  set_target_properties(PropsSI PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}) #Put .lib in this directory so it won't get installed
  set_target_properties(HAPropsSI PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}) #Put .lib in this directory so it won't get installed

  if(WIN32) # 32-bit or 64-bit mex
    if (CMAKE_CL_64)
        SET_TARGET_PROPERTIES(PropsSI PROPERTIES PREFIX "" SUFFIX .mexw64)
        SET_TARGET_PROPERTIES(HAPropsSI PROPERTIES PREFIX "" SUFFIX .mexw64)
    else()
        SET_TARGET_PROPERTIES(PropsSI PROPERTIES SUFFIX .mexw32)
        SET_TARGET_PROPERTIES(HAPropsSI PROPERTIES SUFFIX .mexw32)
    endif()
  else()
  	if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        if (${BITNESS} EQUAL "32")
		    SET_TARGET_PROPERTIES(PropsSI PROPERTIES PREFIX "" SUFFIX .mexmaci32 PREFIX "")
    	    SET_TARGET_PROPERTIES(HAPropsSI PROPERTIES PREFIX "" SUFFIX .mexmaci32 PREFIX "")        
  		elseif((${BITNESS} EQUAL "64"))
	        SET_TARGET_PROPERTIES(PropsSI PROPERTIES PREFIX "" SUFFIX .mexmaci64 PREFIX "")
    	    SET_TARGET_PROPERTIES(HAPropsSI PROPERTIES PREFIX "" SUFFIX .mexmaci64 PREFIX "")        
  		endif()
    else()
    	if (CMAKE_SIZEOF_VOID_P MATCHES "8")
        	SET_TARGET_PROPERTIES(PropsSI PROPERTIES PREFIX "" SUFFIX .mexa64 PREFIX "")
        	SET_TARGET_PROPERTIES(HAPropsSI PROPERTIES PREFIX "" SUFFIX .mexa64 PREFIX "")
    	else()
        	SET_TARGET_PROPERTIES(PropsSI PROPERTIES PREFIX "" SUFFIX .mexglx PREFIX "")
        	SET_TARGET_PROPERTIES(HAPropsSI PROPERTIES PREFIX "" SUFFIX .mexglx PREFIX "")
    	endif()
    endif()
  endif()
  add_dependencies (PropsSI generate_headers)
  add_dependencies (HAPropsSI generate_headers)
  
  install (TARGETS PropsSI HAPropsSI DESTINATION ${CMAKE_INSTALL_PREFIX}/MATLAB/${CMAKE_SYSTEM_NAME}_${BITNESS}bit)
endif()

if (COOLPROP_PYTHON_BINARIES)
  if (WIN32)
    set(COOLPROP_PYTHON_BINARY_VERSIONS  bdist_wheel --dist-dir ${CMAKE_INSTALL_PREFIX}/Python bdist_wininst --dist-dir ${CMAKE_INSTALL_PREFIX}/Python)
  elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(COOLPROP_PYTHON_BINARY_VERSIONS  bdist_wheel --dist-dir ${CMAKE_INSTALL_PREFIX}/Python)
  endif()
  
  add_custom_target(CoolProp
                    COMMAND python setup.py ${COOLPROP_PYTHON_BINARY_VERSIONS} 
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/wrappers/Python
                    )
  
endif()

if (COOLPROP_PYTHON_PYPI)
    
  add_custom_target(CoolProp
                    COMMAND python prepare_pypi.py --dist-dir=${CMAKE_INSTALL_PREFIX}/Python
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/wrappers/Python/pypi
                    )
endif()

if (COOLPROP_JAVASCRIPT_MODULE)
  # cmake -DCOOLPROP_JAVASCRIPT_MODULE=ON 
  #       -DCMAKE_TOOLCHAIN_FILE=${EMSCRIPTEN}/cmake/Platform/Emscripten.cmake 
  #       ../..

  # Toolchain MUST be defined in the call to CMake

  add_definitions(-s DISABLE_EXCEPTION_CATCHING=0)
  add_definitions(-s ASSERTIONS=1)
  add_definitions(-DEXTERNC)
  set(CMAKE_EXE_LINKER_FLAGS "-DEXTERNC -s EXPORTED_FUNCTIONS=\"['_main','_F2K','_Props1SI','_PropsSI','_get_global_param_string','_HAProps']\"")
  set(CMAKE_BUILD_TYPE Release)
  
  list(APPEND APP_SOURCES "${CMAKE_SOURCE_DIR}/src/CoolPropLib.cpp")
  include_directories(${APP_INCLUDE_DIRS})
  add_executable(coolprop ${APP_SOURCES})
  add_dependencies (coolprop generate_headers)
  SET_TARGET_PROPERTIES(coolprop PROPERTIES PREFIX "" SUFFIX .js)
  install (TARGETS coolprop DESTINATION ${CMAKE_INSTALL_PREFIX}/Javascript)
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

if (COOLPROP_MAIN_MODULE)
  list(APPEND APP_SOURCES "${CMAKE_SOURCE_DIR}/src/main.cxx")
  add_executable        (Main ${APP_SOURCES})
  add_dependencies      (Main generate_headers)
if(UNIX)
    target_link_libraries (Main ${CMAKE_DL_LIBS})
  endif()
endif()

###      COOLPROP TESTING APP       ###
if (COOLPROP_CATCH_MODULE)
  enable_testing()

  list(APPEND APP_SOURCES "${CMAKE_SOURCE_DIR}/src/Tests/test_main.cxx")
  list(APPEND APP_SOURCES "${CMAKE_SOURCE_DIR}/src/Tests/CoolProp-Tests.cpp")
  
  # CATCH TEST, compile everything with catch and set test entry point
  add_executable        (CatchTestRunner ${APP_SOURCES})
  add_dependencies      (CatchTestRunner generate_headers)
  set_target_properties (CatchTestRunner PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DENABLE_CATCH")
  if(UNIX)
    target_link_libraries (CatchTestRunner ${CMAKE_DL_LIBS})
  endif()
  add_test(ProcedureTests CatchTestRunner)
endif()

if (COOLPROP_CPP_EXAMPLE_TEST)
  # C++ Documentation Test
  add_executable        (docuTest.exe "Web/examples/C++/Example.cpp")
  add_dependencies      (docuTest.exe ${app_name})
  target_link_libraries (docuTest.exe ${app_name})
  if(UNIX)
    target_link_libraries (docuTest.exe ${CMAKE_DL_LIBS})
  endif()
  add_test(DocumentationTest docuTest.exe)
endif()

if (COOLPROP_SNIPPETS)

  # Make the static library with which the snippets will be linked
  add_library(${app_name} STATIC ${APP_SOURCES})
  add_dependencies (${app_name} generate_headers)
  set_target_properties (${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB -DCONVENTION=__cdecl")
  
  # Collect all the snippets
  file(GLOB_RECURSE snippets "${CMAKE_SOURCE_DIR}/Web/coolprop/snippets/*.cxx")
  
  message(STATUS "snippets found = ${snippets}")
  foreach (snippet ${snippets})
    
    get_filename_component(snippet_name ${snippet} NAME)
    message(STATUS "snippet_name = ${snippet_name}")
    
    add_executable        (${snippet_name} ${snippet})
    add_dependencies      (${snippet_name} CoolProp)
    target_link_libraries (${snippet_name} CoolProp)
    if(UNIX)
      target_link_libraries (${snippet_name} ${CMAKE_DL_LIBS})
    endif()
    
    if ( MSVC )
        set_target_properties( ${snippet_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin )
        set_target_properties( ${snippet_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_BINARY_DIR}/bin )
        set_target_properties( ${snippet_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/bin )
        # etc for the other available configuration types (MinSizeRel, RelWithDebInfo)
    endif ()
    # Run it and save the output to a file with .output appended
    message(STATUS "${CMAKE_CURRENT_BINARY_DIR}/bin/${snippet_name} > ${CMAKE_CURRENT_BINARY_DIR}/bin/${snippet_name}.output")
    
    add_custom_command(TARGET ${snippet_name}
                       POST_BUILD
                       COMMAND ${CMAKE_CURRENT_BINARY_DIR}/bin/${snippet_name} > ${CMAKE_SOURCE_DIR}/Web/coolprop/snippets/${snippet_name}.output)
  
  endforeach()
  
endif()

if (COOLPROP_CLANG_ADDRESS_SANITIZER)

  SET(CMAKE_CXX_FLAGS "-fsanitize=address -g")
  list(APPEND APP_SOURCES "${CMAKE_SOURCE_DIR}/src/Tests/test_main.cxx")
  # CATCH TEST, compile everything with catch and set test entry point
  add_executable        (CatchTestRunner ${APP_SOURCES})
  add_dependencies      (CatchTestRunner generate_headers)
  set_target_properties (CatchTestRunner PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DENABLE_CATCH")
  set(CMAKE_EXE_LINKER_FLAGS "-fsanitize=address -lstdc++")
  if(UNIX)
    target_link_libraries (CatchTestRunner ${CMAKE_DL_LIBS})
  endif()
  add_custom_command(TARGET CatchTestRunner
                     POST_BUILD
                     COMMAND ${CMAKE_CURRENT_BINARY_DIR}/CatchTestRunner)
endif()



if (COOLPROP_PROFILE)
  if(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "-g -O2")
    SET(CMAKE_C_FLAGS "-g -O2")
  endif()
endif()


if (COOLPROP_COVERAGE)
  if(CMAKE_COMPILER_IS_GNUCXX)
    # See also http://stackoverflow.com/a/16536401 (detailed guide on using gcov with cmake)
    include(CodeCoverage)
    SET(CMAKE_CXX_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
    SET(CMAKE_C_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
    setup_target_for_coverage(CoolProp_coverage Main coverage)
  endif()
endif()


# TODO: check relevance of http://www.cmake.org/Wiki/BuildingWinDLL

#include_directories("${CMAKE_CURRENT_SOURCE_DIR}/CoolProp")
#FILE(GLOB coolprop_files "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/*.cpp")
#add_library(coolprop STATIC ${coolprop_files})
