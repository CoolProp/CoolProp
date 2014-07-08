cmake_minimum_required(VERSION 2.8)


if (DEFINED COOLPROP_INSTALL_PREFIX)
    set(COOLPROP_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install_root)
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
file(GLOB_RECURSE APP_HEADERS "include/*.h" "src/*.h")
list(REMOVE_ITEM  APP_SOURCES "${CMAKE_SOURCE_DIR}/src/Tests/Tests.cpp")
## You can exclude this file as well, in case you want to run your own tests

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
    set_target_properties (${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB -DCONVENTION=__cdecl")
    install (TARGETS ${app_name} DESTINATION static_library/${CMAKE_SYSTEM_NAME})
endif()

if (COOLPROP_64BIT_SHARED_LIBRARY)
    add_library(${app_name} SHARED ${APP_SOURCES})
    set_target_properties (${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB")
    if (!MSVC)
      set_target_properties(${app_name} PROPERTIES COMPILE_FLAGS "-m64" LINK_FLAGS "-m64")
    endif()
    add_dependencies (${app_name} generate_headers)
    set_target_properties(${app_name} PROPERTIES PREFIX "")
    install (TARGETS ${app_name} DESTINATION shared_library/${CMAKE_SYSTEM_NAME}/64bit)
endif()

if (COOLPROP_32BIT_CDECL_SHARED_LIBRARY)
    add_library(${app_name} SHARED ${APP_SOURCES})
    set_target_properties (${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB -DCONVENTION=__cdecl")
    if (!MSVC)
      set_target_properties(${app_name} PROPERTIES COMPILE_FLAGS "-m32" LINK_FLAGS "-m32")
    endif()
    add_dependencies (${app_name} generate_headers)
    set_target_properties(${app_name} PROPERTIES PREFIX "")
    install (TARGETS ${app_name} DESTINATION shared_library/${CMAKE_SYSTEM_NAME}/32bit__cdecl_calling_convention)
endif()

if (COOLPROP_32BIT_STDCALL_SHARED_LIBRARY)
    add_library(${app_name} SHARED ${APP_SOURCES})
    set_target_properties (${app_name} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB -DCONVENTION=__stdcall")
    if (!MSVC)
      set_target_properties(${app_name} PROPERTIES COMPILE_FLAGS "-m32" LINK_FLAGS "-m32")
    endif()
    add_dependencies (${app_name} generate_headers)
    set_target_properties(${app_name} PROPERTIES PREFIX "")
    install (TARGETS ${app_name} DESTINATION shared_library/${CMAKE_SYSTEM_NAME}/32bit__stdcall_calling_convention)
endif()

# EES is only compiled for windows
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
  add_dependencies (${app_name} generate_headers)  
  install (TARGETS ${app_name} DESTINATION EES/${CMAKE_SYSTEM_NAME})
endif()

if (COOLPROP_OCTAVE_MODULE)

  # Must have SWIG and Octave
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE(${SWIG_USE_FILE})
  FIND_PACKAGE(Octave REQUIRED)
  
  # Make a src directory to deal with file permissions problem with MinGW makefile
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src)
  
  # Find the required libraries
  find_library(OCTAVE_LIB liboctave.a PATHS ${OCTAVE_LINK_DIRS})
  find_library(OCTINTERP_LIB liboctinterp.a PATHS ${OCTAVE_LINK_DIRS})
  find_library(CRUFT_LIB libcruft.a PATHS ${OCTAVE_LINK_DIRS})
  
  message(STATUS "OCTAVE_LIB=${OCTAVE_LIB}")
  message(STATUS "OCTINTERP_LIB=${OCTINTERP_LIB}")
  message(STATUS "CRUFT_LIB=${CRUFT_LIB}")
  
  # Set the include folders
  SET(OCTAVE_WRAP_INCLUDE_DIRS ${INCLUDE_DIR})
  foreach(ITR ${OCTAVE_INCLUDE_DIRS})
    list(APPEND OCTAVE_WRAP_INCLUDE_DIRS ${ITR})
  endforeach()
  include_directories(${OCTAVE_WRAP_INCLUDE_DIRS})

  set(I_FILE "${CMAKE_SOURCE_DIR}/src/CoolProp.i")

  SET_SOURCE_FILES_PROPERTIES(${I_FILE} PROPERTIES CPLUSPLUS ON)
  SWIG_ADD_MODULE(CoolProp octave ${I_FILE} ${APP_SOURCES})
  if (!OSX)
    SWIG_LINK_LIBRARIES(CoolProp ${OCTAVE_LIB} ${OCTINTERP_LIB} ${CRUFT_LIB})  
  else()
    SWIG_LINK_LIBRARIES(CoolProp ${OCTAVE_LIB} ${OCTINTERP_LIB})
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
  
  # Install all the generated cs files
  install(
    CODE "file( GLOB _GeneratedCsharpSources \"${CMAKE_CURRENT_BINARY_DIR}/*.cs\" )"
    CODE "file( INSTALL \${_GeneratedCsharpSources} DESTINATION ${CMAKE_INSTALL_PREFIX}/Csharp/platform-independent )"
  )
  install (FILES ${CMAKE_SOURCE_DIR}/wrappers/C\#/Example.cs DESTINATION Csharp)
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
  elseif((${BITNESS} EQUAL "32"))
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

  set(I_FILE "${CMAKE_SOURCE_DIR}/src/CoolProp.i")

  SET_SOURCE_FILES_PROPERTIES(${I_FILE} PROPERTIES CPLUSPLUS ON)

  SWIG_ADD_MODULE(CoolProp java ${I_FILE} ${APP_SOURCES})

  if (WIN32)
    set_target_properties(CoolProp PROPERTIES PREFIX "")
  endif()
  
  if (!MSVC)
    set_target_properties(COOLPROP_EES PROPERTIES COMPILE_FLAGS "-m64" LINK_FLAGS "-m64")
  endif()
  add_dependencies (${app_name} generate_headers)
  
  # Install all the generated java files
  install(
    CODE "file( GLOB _GeneratedJavaSources \"${CMAKE_CURRENT_BINARY_DIR}/*.java\" )"
    CODE "file( INSTALL \${_GeneratedJavaSources} DESTINATION ${CMAKE_INSTALL_PREFIX}/Java/platform-independent )"
  )
  install (FILES ${CMAKE_SOURCE_DIR}/wrappers/Java/Example.java DESTINATION ${CMAKE_INSTALL_PREFIX}/Java)
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

  if(WIN32) # 32-bit or 64-bit mex
    if (CMAKE_CL_64)
        SET_TARGET_PROPERTIES(PropsSI PROPERTIES PREFIX "" SUFFIX .mexw64)
        SET_TARGET_PROPERTIES(HAPropsSI PROPERTIES PREFIX "" SUFFIX .mexw64)
    else()
        SET_TARGET_PROPERTIES(PropsSI PROPERTIES SUFFIX .mexw32)
        SET_TARGET_PROPERTIES(HAPropsSI PROPERTIES SUFFIX .mexw32)
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
  add_dependencies (PropsSI generate_headers)
  add_dependencies (HAPropsSI generate_headers)
endif()

if (COOLPROP_PYTHON_LOCAL_INSTALL)
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/Python/setup.py.in")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/module/__init__.py")
    set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

    configure_file(${SETUP_PY_IN} ${SETUP_PY})

    add_custom_command(OUTPUT ${OUTPUT}
                       COMMAND ${PYTHON} ${SETUP_PY} build
                       COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                       DEPENDS ${DEPS})

    add_custom_target(target ALL DEPENDS ${OUTPUT})

    install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install)")
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
if (COOLPROP_CATCH_MODULE)
  enable_testing()

  list(APPEND APP_SOURCES "${CMAKE_SOURCE_DIR}/src/Tests/test_main.cxx")
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


# TODO: check relevance of http://www.cmake.org/Wiki/BuildingWinDLL

#include_directories("${CMAKE_CURRENT_SOURCE_DIR}/CoolProp")
#FILE(GLOB coolprop_files "${CMAKE_CURRENT_SOURCE_DIR}/CoolProp/*.cpp")
#add_library(coolprop STATIC ${coolprop_files})
