#######################################
#         CORE LIBRARY BUILD          #
#-------------------------------------#
# Build the main CoolProp library in  #
# static, shared, or object form      #
#######################################

# Include needed CMake modules
include(CheckIncludeFileCXX)

#######################################
#         BITNESS DETECTION           #
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
  message(STATUS "You cannot force a certain bitness for Visual Studio, use the generator settings for this purpose.")
  message(STATUS "Pass '-G \"Visual Studio 10 2010 Win64\"' to CMake to make a 64bit binary using VS2010.")
  message(STATUS "Pass '-G \"Visual Studio 10 2010\"' to CMake to make a 32bit binary using VS2010.")
  message(STATUS "Pass '-G \"Visual Studio 9 2008 Win64\"' to CMake to make a 64bit binary using VS2008.")
  message(STATUS "Pass '-G \"Visual Studio 9 2008\"' to CMake to make a 32bit binary using VS2008.")
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
#      CODE GENERATION TARGETS        #
#######################################

# Generate headers from fluid/mixture definitions
add_custom_target(
  generate_headers
  COMMAND "${PYTHON_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/dev/generate_headers.py")

# Generate examples for various languages
if(NOT COOLPROP_NO_EXAMPLES)
  add_custom_target(
    generate_examples
    COMMAND "${PYTHON_EXECUTABLE}" example_generator.py Python
            "${CMAKE_CURRENT_BINARY_DIR}/Example.py"
    COMMAND "${PYTHON_EXECUTABLE}" example_generator.py Octave
            "${CMAKE_CURRENT_BINARY_DIR}/Example.m"
    COMMAND "${PYTHON_EXECUTABLE}" example_generator.py R
            "${CMAKE_CURRENT_BINARY_DIR}/Example.R"
    COMMAND "${PYTHON_EXECUTABLE}" example_generator.py Java
            "${CMAKE_CURRENT_BINARY_DIR}/Example.java"
    COMMAND "${PYTHON_EXECUTABLE}" example_generator.py Csharp
            "${CMAKE_CURRENT_BINARY_DIR}/Example.cs"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dev/scripts/examples")
else()
  add_custom_target(
    generate_examples
    COMMAND echo "Example generation has been disabled with the COOLPROP_NO_EXAMPLES option.")
endif()

#######################################
#      LIBRARY CONFIGURATION          #
#######################################

set(COOLPROP_LIBRARY_SOURCE
    "src/CoolPropLib.cpp"
    CACHE STRING "The file that contains the exported functions")

set(COOLPROP_LIBRARY_HEADER
    "include/CoolPropLib.h"
    CACHE STRING "The file that contains the export header")

set(COOLPROP_LIBRARY_NAME
    "CoolProp"
    CACHE STRING "The name of the generated library")

set(COOLPROP_LIBRARY_EXPORTS
    ""
    CACHE STRING "The file that contains the export alias list")

# Determine calling convention based on bitness
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

# Validate library type options
if((COOLPROP_OBJECT_LIBRARY AND COOLPROP_STATIC_LIBRARY)
   OR (COOLPROP_OBJECT_LIBRARY AND COOLPROP_SHARED_LIBRARY)
   OR (COOLPROP_STATIC_LIBRARY AND COOLPROP_SHARED_LIBRARY))
  message(FATAL_ERROR "You can only use one of the library switches!")
endif()

#######################################
#      BUILD LIBRARY TARGET           #
#######################################

if(COOLPROP_OBJECT_LIBRARY OR COOLPROP_STATIC_LIBRARY OR COOLPROP_SHARED_LIBRARY)
  set(LIB_NAME ${COOLPROP_LIBRARY_NAME})

  # Create library target based on type
  if(COOLPROP_OBJECT_LIBRARY)
    add_library(${LIB_NAME} OBJECT ${APP_SOURCES})
    set(COOLPROP_LIBRARY_SOURCE "")
    set(COOLPROP_LIBRARY_HEADER "")

  elseif(COOLPROP_STATIC_LIBRARY)
    list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${COOLPROP_LIBRARY_SOURCE}")
    add_library(${LIB_NAME} STATIC ${APP_SOURCES} ${COOLPROP_LIBRARY_EXPORTS})

    if(MSVC)
      set_property(TARGET ${LIB_NAME} PROPERTY DEBUG_POSTFIX d)
      set_property(TARGET ${LIB_NAME} PROPERTY RELEASE_POSTFIX)
      modify_msvc_flags("/MD")
    endif()

    install(TARGETS ${LIB_NAME}
            DESTINATION static_library/${CMAKE_SYSTEM_NAME}/${BITNESS}bit_${CMAKE_CXX_COMPILER_ID}_${CMAKE_CXX_COMPILER_VERSION})
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${COOLPROP_LIBRARY_HEADER}
            DESTINATION static_library)

  elseif(COOLPROP_SHARED_LIBRARY)
    list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${COOLPROP_LIBRARY_SOURCE}")
    add_library(${LIB_NAME} SHARED ${APP_SOURCES} ${COOLPROP_LIBRARY_EXPORTS})

    # Determine output folder (special handling for ARM64)
    if(MSVC AND ("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "ARM64" OR "${CMAKE_VS_PLATFORM_NAME}" STREQUAL "arm64"))
      set(OUTPUT_FOLDER "shared_library/${CMAKE_SYSTEM_NAME}/${BITNESS}bit__arm64")
    else()
      set(OUTPUT_FOLDER "shared_library/${CMAKE_SYSTEM_NAME}/${BITNESS}bit${CONVENTION}")
    endif()

    install(TARGETS ${LIB_NAME} DESTINATION ${OUTPUT_FOLDER})
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${COOLPROP_LIBRARY_HEADER}
            DESTINATION shared_library)

    set_property(TARGET ${LIB_NAME} APPEND_STRING PROPERTY COMPILE_FLAGS " -DCOOLPROP_LIB")

    # MSVC-specific settings
    if(MSVC)
      set_property(TARGET ${LIB_NAME} PROPERTY DEBUG_POSTFIX d)
      set_property(TARGET ${LIB_NAME} PROPERTY RELEASE_POSTFIX)
      set_property(TARGET ${LIB_NAME} PROPERTY PREFIX "")
      modify_msvc_flags("/MT")

      # Generate exports and headers dumps
      add_custom_command(TARGET ${LIB_NAME} POST_BUILD
                         COMMAND dumpbin /EXPORTS $<TARGET_FILE:${LIB_NAME}> >
                                 ${CMAKE_CURRENT_BINARY_DIR}/exports.txt)
      install(FILES ${CMAKE_CURRENT_BINARY_DIR}/exports.txt DESTINATION ${OUTPUT_FOLDER})

      add_custom_command(TARGET ${LIB_NAME} POST_BUILD
                         COMMAND dumpbin /HEADERS $<TARGET_FILE:${LIB_NAME}> >
                                 ${CMAKE_CURRENT_BINARY_DIR}/headers.txt)
      install(FILES ${CMAKE_CURRENT_BINARY_DIR}/headers.txt DESTINATION ${OUTPUT_FOLDER})
    endif()

    # Linux-specific settings
    if("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
      set_property(TARGET ${LIB_NAME} PROPERTY VERSION ${COOLPROP_VERSION})
      set_property(TARGET ${LIB_NAME} PROPERTY SOVERSION ${COOLPROP_VERSION_MAJOR})
    endif()
  else()
    message(FATAL_ERROR "You have to build a static or shared library.")
  endif()

  # Link with dynamic loading libraries
  if(NOT COOLPROP_OBJECT_LIBRARY)
    target_link_libraries(${LIB_NAME} ${CMAKE_DL_LIBS})
  endif()

  # Eigen workaround for MSVC 2008
  if(MSVC90)
    message(STATUS "EIGEN WORKAROUND ACTIVE!!")
    set_property(TARGET ${LIB_NAME} APPEND_STRING PROPERTY COMPILE_FLAGS " -DEIGEN_DONT_VECTORIZE")
  endif()

  # macOS-specific settings
  if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
    if(DEFINED OSX_COMPILE_FLAGS)
      set_target_properties(${LIB_NAME} PROPERTIES APPEND_STRING PROPERTY COMPILE_FLAGS "${OSX_COMPILE_FLAGS}")
    endif()
    if(DEFINED OSX_LINK_FLAGS)
      set_target_properties(${LIB_NAME} PROPERTIES APPEND_STRING PROPERTY LINK_FLAGS "${OSX_LINK_FLAGS}")
    endif()
  endif()

  # Name mangling settings
  if(COOLPROP_EXTERNC_LIBRARY)
    set_property(TARGET ${LIB_NAME} APPEND_STRING PROPERTY COMPILE_FLAGS " -DEXTERNC")
  endif()

  # Dependencies
  add_dependencies(${LIB_NAME} generate_headers)

  # Include directories
  if(CMAKE_VERSION VERSION_GREATER 3.0)
    target_include_directories(${LIB_NAME} PUBLIC ${APP_INCLUDE_DIRS})
  endif()

  # Set bitness flags
  if(NOT MSVC)
    if(NOT "${BITNESS}" STREQUAL "NATIVE")
      message(STATUS "Setting bitness flag -m${BITNESS}")
      set_property(TARGET ${LIB_NAME} APPEND_STRING PROPERTY COMPILE_FLAGS " -m${BITNESS}")
      set_property(TARGET ${LIB_NAME} APPEND_STRING PROPERTY LINK_FLAGS " -m${BITNESS}")
    endif()
  endif()

  # Position-independent code flag
  if(COOLPROP_FPIC)
    message(STATUS "Setting fPIC flag")
    set_property(TARGET ${LIB_NAME} APPEND_STRING PROPERTY COMPILE_FLAGS " -fPIC")
  endif()

  # Calling conventions
  if(NOT ("${CONVENTION}" STREQUAL "NULL" OR "${CONVENTION}" STREQUAL ""))
    set_property(TARGET ${LIB_NAME} APPEND_STRING PROPERTY COMPILE_FLAGS " -DCONVENTION=${CONVENTION}")
  endif()

  # Status messages
  message(STATUS "Library compilation detected:")
  message(STATUS "Creating ${LIB_NAME}, a ${BITNESS}-bit library")
  message(STATUS "CMAKE_SYSTEM_NAME: ${CMAKE_SYSTEM_NAME}")
  message(STATUS "COOLPROP_STATIC_LIBRARY: ${COOLPROP_STATIC_LIBRARY}")
  message(STATUS "COOLPROP_SHARED_LIBRARY: ${COOLPROP_SHARED_LIBRARY}")
  message(STATUS "COOLPROP_OBJECT_LIBRARY: ${COOLPROP_OBJECT_LIBRARY}")
  message(STATUS "CONVENTION: ${CONVENTION}")
  message(STATUS "COOLPROP_LIBRARY_HEADER: ${COOLPROP_LIBRARY_HEADER}")
  message(STATUS "COOLPROP_LIBRARY_SOURCE: ${COOLPROP_LIBRARY_SOURCE}")

  get_property(tmpVar TARGET ${LIB_NAME} PROPERTY COMPILE_FLAGS)
  message(STATUS "COMPILE_FLAGS: ${tmpVar}")
  get_property(tmpVar TARGET ${LIB_NAME} PROPERTY LINK_FLAGS)
  message(STATUS "LINK_FLAGS: ${tmpVar}")
endif()

#######################################
#      PLATFORM-SPECIFIC BUILDS       #
#######################################

# iOS target
if(COOLPROP_IOS_TARGET)
  set(SDKVER "9.2")
  set(DEVROOT "/Applications/Xcode.app/Contents/Developer/Platforms/iPhoneOS.platform/Developer")
  set(SDKROOT "${DEVROOT}/SDKs/iPhoneOS${SDKVER}.sdk")
  if(EXISTS ${SDKROOT})
    set(CMAKE_OSX_SYSROOT "${SDKROOT}")
  else()
    message("Warning, iOS Base SDK path not found: " ${SDKROOT})
  endif()
  set(CMAKE_OSX_ARCHITECTURES "$(ARCHS_STANDARD_32_BIT)")
  set(CMAKE_XCODE_EFFECTIVE_PLATFORMS "-iphoneos;-iphonesimulator")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
endif()

# Debian package
if(COOLPROP_DEBIAN_PACKAGE)
  if(NOT UNIX)
    message(FATAL_ERROR "COOLPROP_DEBIAN_PACKAGE can only be used on linux host")
  endif()
  list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolPropLib.cpp")
  add_library(CoolProp SHARED ${APP_SOURCES})
  set_target_properties(CoolProp PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB")
  set_target_properties(CoolProp PROPERTIES VERSION ${COOLPROP_VERSION}
                                            SOVERSION ${COOLPROP_VERSION_MAJOR})
  add_dependencies(CoolProp generate_headers)
  install(TARGETS CoolProp DESTINATION "${CMAKE_INSTALL_PREFIX}/usr/lib")
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/CoolPropLib.h
          DESTINATION "${CMAKE_INSTALL_PREFIX}/usr/include")
endif()

# VxWorks makefile generation
if(COOLPROP_VXWORKS_MAKEFILE)
  set(INCLUDE_DIRECTORIES)
  foreach(_srcFile ${APP_INCLUDE_DIRS})
    string(CONCAT _el "-I\"" ${_srcFile} "\"")
    string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}" "$(COOLPROP_ROOT)" _el "${_el}")
    list(APPEND INCLUDE_DIRECTORIES ${_el})
  endforeach()
  string(REPLACE ";" " " INCLUDE_DIRECTORIES "${INCLUDE_DIRECTORIES}")
  set(OLD_ROOT /home/ian/.wine/drive_c/)
  set(NEW_ROOT c:/)
  string(REPLACE ${OLD_ROOT} ${NEW_ROOT} INCLUDE_DIRECTORIES "${INCLUDE_DIRECTORIES}")
  set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/src")
  string(REPLACE ${OLD_ROOT} ${NEW_ROOT} SRC "${SRC}")
  file(RELATIVE_PATH COOLPROP_ROOT "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/wrappers/Labview/vxWorks/Makefile.in" "vxWorksMakefile")
endif()

# VxWorks library
if(COOLPROP_VXWORKS_LIBRARY_MODULE OR COOLPROP_VXWORKS_LIBRARY)
  list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/CoolPropLib.cpp")
  add_executable(CoolProp ${APP_SOURCES})
  set_target_properties(CoolProp PROPERTIES SUFFIX ".out"
                                            COMPILE_FLAGS "${COMPILE_FLAGS} -DEXTERNC")
  add_dependencies(CoolProp generate_headers)
  install(TARGETS CoolProp DESTINATION "${COOLPROP_INSTALL_PREFIX}/shared_library/VxWorks")
endif()
