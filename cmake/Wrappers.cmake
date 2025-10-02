if(COOLPROP_PRIME_MODULE)
  if(NOT WIN32)
    message(
      FATAL_ERROR "COOLPROP_PRIME_MODULE can only be used on windows host")
  endif()
  if("${COOLPROP_PRIME_ROOT}" STREQUAL "")
    message(
      FATAL_ERROR
        "You must provide the path to Mathcad Prime Root directory using something like -DCOOLPROP_PRIME_ROOT=\"C:/Program Files/PTC/Mathcad Prime 3.1\""
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
          DESTINATION "${COOLPROP_INSTALL_PREFIX}/MathcadPrime")
  install(
    FILES
      "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/MathCAD/Prime/CoolPropFluidProperties.mcdx"
    DESTINATION MathcadPrime)
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
  set(COOLPROP_WINDOWS_PACKAGE_ARM64_DIR "${CMAKE_CURRENT_BINARY_DIR}/arm64bitDLL")
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
            "${COOLPROP_WINDOWS_PACKAGE_ARM64_DIR}"
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
    COMMENT "Building the 64bit shared library with x86_64 architecture"
    VERBATIM)

  # Build the 64bit DLL
  add_custom_command(
    TARGET COOLPROP_WINDOWS_PACKAGE_SHARED_LIBRARIES
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} ARGS "-G${COOLPROP_WINDOWS_PACKAGE_DLL_GEN}"
            "-Aarm64" "${CMAKE_CURRENT_SOURCE_DIR}" "-DCOOLPROP_SHARED_LIBRARY=ON"
    COMMAND ${CMAKE_COMMAND} ARGS "--build" "." "--target" "CoolProp" "--config"
            "Release"
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy"
      "${COOLPROP_WINDOWS_PACKAGE_ARM64_DIR}/Release/CoolProp.dll"
      "${COOLPROP_WINDOWS_PACKAGE_TMP_DIR}/source/CoolProp_arm64.dll"
    WORKING_DIRECTORY ${COOLPROP_WINDOWS_PACKAGE_ARM64_DIR}
    COMMENT "Building the 64bit shared library with arm64 architecture"
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
    # https://stackoverflow.com/questions/7035734/how-do-i-get-the-parent-directory-path-of-a-path-in-cmake
    get_filename_component(PARENT_DIR ${ITR} DIRECTORY)
    list(APPEND OCTAVE_WRAP_INCLUDE_DIRS ${PARENT_DIR})
    message(STATUS "PARENT_DIR: ${PARENT_DIR}")
  endforeach()
  include_directories(${OCTAVE_WRAP_INCLUDE_DIRS})
  message(STATUS "OCTAVE_WRAP_INCLUDE_DIRS: ${OCTAVE_WRAP_INCLUDE_DIRS}")

  # Disable internal error catching and allow swig to do the error catching itself
  add_definitions(-DNO_ERROR_CATCHING)

  set(I_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolProp.i")

  set(SWIG_OPTIONS "${COOLPROP_SWIG_OPTIONS}")
  set_source_files_properties(${I_FILE} PROPERTIES SWIG_FLAGS "${SWIG_OPTIONS}")
  set_source_files_properties(${I_FILE} PROPERTIES CPLUSPLUS ON)

  set(SWIG_MODULE_CoolProp_EXTRA_DEPS ${SWIG_DEPENDENCIES})
  swig_add_library(
    CoolProp
    LANGUAGE octave 
    SOURCES ${I_FILE} ${APP_SOURCES}
  )

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
  include_directories(${JAVA_INCLUDE_PATH}/darwin)

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

  # set paths for LibreOffice tools
  set(LO_UNOIDL_WRITE "${LO_SDK_PATH}/bin/unoidl-write")
  set(COOLPROP_LIBREOFFICE_TMP_DIR "${CMAKE_CURRENT_BINARY_DIR}/LibreOffice")

  # set version strings for LibreOffice extension
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/LibreOffice/src/description.xml.in"
    "${COOLPROP_LIBREOFFICE_TMP_DIR}/src/description.xml")
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/LibreOffice/src/scripts/scripts.py.in"
    "${COOLPROP_LIBREOFFICE_TMP_DIR}/src/scripts/scripts.py")

  add_custom_target(CoolPropLibreOfficeAddin ALL DEPENDS CoolProp.oxt)
    
  add_custom_command(
    OUTPUT CoolProp.oxt
    # copy source files to build directory
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "copy_directory"
      "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/LibreOffice/src"
      "${COOLPROP_LIBREOFFICE_TMP_DIR}/src"
    COMMAND
      ${CMAKE_COMMAND} ARGS "-E" "remove"
      "${COOLPROP_LIBREOFFICE_TMP_DIR}/src/description.xml.in"
      "${COOLPROP_LIBREOFFICE_TMP_DIR}/src/scripts/scripts.py.in"
    # build the registry database file (rdb)
    COMMAND
      ${LO_UNOIDL_WRITE}
      ${LO_PROGRAM_PATH}/types.rdb ${LO_PROGRAM_PATH}/types/offapi.rdb
      XCoolProp.idl XCoolProp.rdb
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
    WORKING_DIRECTORY ${COOLPROP_LIBREOFFICE_TMP_DIR}/src
    COMMENT "Building LibreOffice wrapper"
    VERBATIM)

   # install LibreOffice extension and example spreadsheet file
   install(FILES "${COOLPROP_LIBREOFFICE_TMP_DIR}/CoolProp.oxt"
           DESTINATION "${COOLPROP_INSTALL_PREFIX}/LibreOffice")
   install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/LibreOffice/TestLibreOffice.ods"
           DESTINATION "${COOLPROP_INSTALL_PREFIX}/LibreOffice")
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
  add_definitions(-DCOOLPROP_NO_INCBIN)
  # # If you want a monolithic file with no async memory loading, define EMSCRIPTEN_NO_MEMORY_INIT_FILE
  # if(EMSCRIPTEN_NO_MEMORY_INIT_FILE)
  #   set(EMSCRIPTEN_INIT_FLAG "--memory-init-file 0")
  # else()
  #   set(EMSCRIPTEN_INIT_FLAG "--memory-init-file 1")
  # endif()

  set(CMAKE_EXE_LINKER_FLAGS
      "-lembind ${EMSCRIPTEN_INIT_FLAG} -s ASSERTIONS=1 -s DISABLE_EXCEPTION_CATCHING=0 -sALLOW_MEMORY_GROWTH=1 -s EXPORT_ES6=1 -s MODULARIZE=1"
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

