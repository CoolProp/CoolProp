include_guard(GLOBAL)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(COOLPROP_STDCALL_LIBRARY
       "Build CoolProp as a 32-bit shared library with stdcall" OFF)
option(COOLPROP_CDECL_LIBRARY
       "Build CoolProp as a 32-bit shared library with cdecl" OFF)
option(COOLPROP_EXTERNC_LIBRARY
       "Force C linkage for the CoolProp C API" OFF)

set(COOLPROP_LIBRARY_SOURCE
    "src/CoolPropLib.cpp"
    CACHE STRING "The file that contains the exported functions")
set(COOLPROP_LIBRARY_HEADER
    "include/CoolProp/CoolPropLib.h"
    CACHE STRING "The file that contains the export header")
set(COOLPROP_LIBRARY_NAME
    "${app_name}"
    CACHE STRING "The legacy name of the generated library target")
set(COOLPROP_LIBRARY_EXPORTS
    ""
    CACHE STRING "The file that contains the export alias list")

set(COOLPROP_DEFAULT_LIBRARY
    "SHARED"
    CACHE STRING
          "Variant selected by CoolProp::CoolProp when static and shared are built")
set_property(CACHE COOLPROP_DEFAULT_LIBRARY PROPERTY STRINGS STATIC SHARED)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(_coolprop_top_level_install_default ON)
else()
  set(_coolprop_top_level_install_default OFF)
endif()
option(COOLPROP_INSTALL_CMAKE_PACKAGE
       "Install the standard CoolProp CMake package and development headers"
       ${_coolprop_top_level_install_default})
option(COOLPROP_INSTALL_LEGACY_LAYOUT
       "Install the historical static_library/shared_library release layout"
       ${_coolprop_top_level_install_default})
set(COOLPROP_INSTALL_CMAKEDIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/CoolProp"
    CACHE STRING "CoolProp CMake package installation directory")

function(_coolprop_set_msvc_runtime target linkage)
  if(NOT MSVC)
    return()
  endif()

  if(COOLPROP_MSVC_STATIC)
    set(_runtime "MultiThreaded")
    set(_runtime_flag "/MT")
  elseif(COOLPROP_MSVC_DYNAMIC)
    set(_runtime "MultiThreadedDLL")
    set(_runtime_flag "/MD")
  elseif(linkage STREQUAL "STATIC")
    # Preserve the existing defaults: static CoolProp uses the DLL CRT.
    set(_runtime "MultiThreadedDLL")
    set(_runtime_flag "/MD")
  else()
    # Preserve the existing defaults: shared CoolProp embeds the static CRT.
    set(_runtime "MultiThreaded")
    set(_runtime_flag "/MT")
  endif()

  # Debug configurations must always use the matching debug CRT.  Keeping a
  # release CRT in a d-postfixed library causes _ITERATOR_DEBUG_LEVEL and heap
  # ownership mismatches in otherwise ordinary Debug consumers.
  if(_runtime STREQUAL "MultiThreadedDLL")
    set(_runtime "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
  else()
    set(_runtime "MultiThreaded$<$<CONFIG:Debug>:Debug>")
  endif()

  set_property(TARGET "${target}" PROPERTY MSVC_RUNTIME_LIBRARY "${_runtime}")

  # A parent with cmake_minimum_required(VERSION 3.14) may have enabled MSVC
  # languages while CMP0091 was OLD. In that case CMake leaves this internal
  # default empty and ignores MSVC_RUNTIME_LIBRARY. Append an explicit
  # target-local option so add_subdirectory consumers still get the requested
  # runtime without rewriting the parent's global CMAKE_CXX_FLAGS_* values.
  if(_coolprop_needs_msvc_runtime_fallback)
    target_compile_options(
      "${target}"
      PRIVATE "$<$<CONFIG:Debug>:${_runtime_flag}d>"
              "$<$<NOT:$<CONFIG:Debug>>:${_runtime_flag}>")
  endif()
endfunction()

# Keep executables linked through the canonical target on the same CRT as the
# concrete library variant selected by the producer.
function(_coolprop_set_canonical_msvc_runtime target)
  if(COOLPROP_STATIC_LIBRARY AND COOLPROP_SHARED_LIBRARY)
    string(TOUPPER "${COOLPROP_DEFAULT_LIBRARY}" _linkage)
  elseif(COOLPROP_STATIC_LIBRARY)
    set(_linkage STATIC)
  elseif(COOLPROP_SHARED_LIBRARY)
    set(_linkage SHARED)
  else()
    return()
  endif()

  _coolprop_set_msvc_runtime("${target}" "${_linkage}")
endfunction()

function(_coolprop_configure_core_target target linkage)
  # CMake applies cxx_* features only to C++ translation units. This promotes
  # C++ consumers to the standard required by the public headers without
  # imposing a C++ compiler on projects which link only the shared C API.
  target_compile_features("${target}" PUBLIC cxx_std_17)
  target_include_directories(
    "${target}"
    PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    PRIVATE ${APP_INCLUDE_DIRS})

  target_compile_definitions("${target}" PRIVATE VALIJSON_USE_EXCEPTIONS=1)
  target_link_libraries("${target}" PUBLIC coolprop_eigen_headers
                                           coolprop_fmt_headers)

  if(NOT linkage STREQUAL "OBJECT")
    target_link_libraries("${target}" PRIVATE Threads::Threads)
    if(CMAKE_DL_LIBS)
      target_link_libraries("${target}" PRIVATE ${CMAKE_DL_LIBS})
    endif()
  endif()

  add_dependencies("${target}" generate_headers)

  if(MSVC90)
    target_compile_definitions("${target}" PRIVATE EIGEN_DONT_VECTORIZE)
  endif()

  if(APPLE)
    if(DEFINED OSX_COMPILE_FLAGS)
      set_property(TARGET "${target}" APPEND_STRING PROPERTY COMPILE_FLAGS
                                                               "${OSX_COMPILE_FLAGS}")
    endif()
    if(DEFINED OSX_LINK_FLAGS AND NOT linkage STREQUAL "OBJECT")
      set_property(TARGET "${target}" APPEND_STRING PROPERTY LINK_FLAGS
                                                            "${OSX_LINK_FLAGS}")
    endif()
  endif()

  if(COOLPROP_EXTERNC_LIBRARY)
    target_compile_definitions("${target}" PUBLIC EXTERNC)
  endif()

  if(NOT MSVC AND NOT BITNESS STREQUAL "NATIVE")
    target_compile_options("${target}" PRIVATE "-m${BITNESS}")
    if(NOT linkage STREQUAL "OBJECT")
      target_link_options("${target}" PRIVATE "-m${BITNESS}")
    endif()
  endif()

  if(COOLPROP_FPIC)
    set_property(TARGET "${target}" PROPERTY POSITION_INDEPENDENT_CODE ON)
  endif()

  if(NOT CONVENTION STREQUAL "")
    target_compile_definitions("${target}" PUBLIC "CONVENTION=${CONVENTION}")
  endif()

  if(linkage STREQUAL "STATIC")
    _coolprop_set_msvc_runtime("${target}" STATIC)
  elseif(linkage STREQUAL "SHARED")
    target_compile_definitions("${target}"
                               PRIVATE COOLPROP_SHARED_LIBRARY_BUILD
                               INTERFACE COOLPROP_SHARED_LIBRARY_USE)
    _coolprop_set_msvc_runtime("${target}" SHARED)
    coolprop_hide_json_symbols("${target}")

    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
      set_target_properties(
        "${target}"
        PROPERTIES VERSION "${COOLPROP_VERSION}"
                   SOVERSION "${COOLPROP_VERSION_MAJOR}")
    endif()
  endif()
endfunction()

function(_coolprop_install_legacy_headers directory)
  install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/${COOLPROP_LIBRARY_HEADER}"
          DESTINATION "${directory}")
  install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/CoolProp"
    DESTINATION "${directory}/include"
    FILES_MATCHING
    PATTERN "*.h"
    REGEX "detail/(json|msgpack)\\.h$" EXCLUDE)
  install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
    DESTINATION "${directory}/include"
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "CoolProp" EXCLUDE
    PATTERN "*_JSON*.h" EXCLUDE
    PATTERN "*_CBOR*.h" EXCLUDE
    PATTERN "CPmsgpack.h" EXCLUDE
    PATTERN "gitrevision.h" EXCLUDE
    PATTERN "cpversion.h" EXCLUDE)
endfunction()

function(_coolprop_install_standard_headers)
  install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/CoolProp"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING
    PATTERN "*.h"
    REGEX "detail/(json|msgpack)\\.h$" EXCLUDE)
  install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "CoolProp" EXCLUDE
    PATTERN "*_JSON*.h" EXCLUDE
    PATTERN "*_CBOR*.h" EXCLUDE
    PATTERN "CPmsgpack.h" EXCLUDE
    PATTERN "gitrevision.h" EXCLUDE
    PATTERN "cpversion.h" EXCLUDE)

  # Eigen and fmt are header-only implementation dependencies which also occur
  # in installed public headers.  Install private copies under the CoolProp
  # include tree so the exported package is self-contained and relocatable.
  set(_third_party_dir "${CMAKE_INSTALL_INCLUDEDIR}/CoolProp/third_party")
  install(DIRECTORY "${Eigen_SOURCE_DIR}/Eigen"
          DESTINATION "${_third_party_dir}/eigen")
  install(FILES "${Eigen_SOURCE_DIR}/unsupported/Eigen/Polynomials"
          DESTINATION "${_third_party_dir}/eigen/unsupported/Eigen")
  install(DIRECTORY "${Eigen_SOURCE_DIR}/unsupported/Eigen/src/Polynomials"
          DESTINATION "${_third_party_dir}/eigen/unsupported/Eigen/src")
  install(DIRECTORY "${fmt_SOURCE_DIR}/include/fmt"
          DESTINATION "${_third_party_dir}/fmt")

  install(FILES "${Eigen_SOURCE_DIR}/COPYING.MPL2"
          DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/licenses/CoolProp/Eigen")
  install(FILES "${fmt_SOURCE_DIR}/LICENSE"
          DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/licenses/CoolProp/fmt")
endfunction()

function(coolprop_add_library_targets)
  if(COOLPROP_CDECL_LIBRARY AND COOLPROP_STDCALL_LIBRARY)
    message(
      FATAL_ERROR
        "COOLPROP_CDECL_LIBRARY and COOLPROP_STDCALL_LIBRARY are mutually exclusive"
    )
  endif()

  if(COOLPROP_MSVC_STATIC AND COOLPROP_MSVC_DYNAMIC)
    message(
      FATAL_ERROR
        "COOLPROP_MSVC_STATIC and COOLPROP_MSVC_DYNAMIC are mutually exclusive"
    )
  endif()

  if(COOLPROP_OBJECT_LIBRARY
     AND (COOLPROP_STATIC_LIBRARY OR COOLPROP_SHARED_LIBRARY))
    message(
      FATAL_ERROR
        "COOLPROP_OBJECT_LIBRARY cannot be combined with COOLPROP_STATIC_LIBRARY or COOLPROP_SHARED_LIBRARY"
    )
  endif()

  string(TOUPPER "${COOLPROP_DEFAULT_LIBRARY}" _default_library)
  if(NOT _default_library STREQUAL "STATIC"
     AND NOT _default_library STREQUAL "SHARED")
    message(FATAL_ERROR
            "COOLPROP_DEFAULT_LIBRARY must be STATIC or SHARED")
  endif()

  if("${BITNESS}" STREQUAL "32")
    if(COOLPROP_CDECL_LIBRARY)
      set(CONVENTION "__cdecl" PARENT_SCOPE)
      set(CONVENTION "__cdecl")
    elseif(COOLPROP_STDCALL_LIBRARY)
      set(CONVENTION "__stdcall" PARENT_SCOPE)
      set(CONVENTION "__stdcall")
    else()
      set(CONVENTION "" PARENT_SCOPE)
      set(CONVENTION "")
    endif()
  elseif("${BITNESS}" STREQUAL "64" OR "${BITNESS}" STREQUAL "NATIVE")
    if(COOLPROP_CDECL_LIBRARY OR COOLPROP_STDCALL_LIBRARY)
      message(WARNING
              "Explicit x86 calling conventions are ignored for ${BITNESS}-bit builds")
    endif()
    set(CONVENTION "" PARENT_SCOPE)
    set(CONVENTION "")
  else()
    message(FATAL_ERROR "Bitness is not defined. Set it and run CMake again.")
  endif()

  if(NOT (COOLPROP_OBJECT_LIBRARY OR COOLPROP_STATIC_LIBRARY
          OR COOLPROP_SHARED_LIBRARY))
    return()
  endif()

  find_package(Threads REQUIRED)

  add_library(coolprop_eigen_headers INTERFACE)
  set_property(TARGET coolprop_eigen_headers PROPERTY EXPORT_NAME
                                                        _EigenHeaders)
  target_include_directories(
    coolprop_eigen_headers SYSTEM
    INTERFACE "$<BUILD_INTERFACE:${Eigen_SOURCE_DIR}>"
              "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/CoolProp/third_party/eigen>"
  )

  add_library(coolprop_fmt_headers INTERFACE)
  set_property(TARGET coolprop_fmt_headers PROPERTY EXPORT_NAME _FmtHeaders)
  target_include_directories(
    coolprop_fmt_headers SYSTEM
    INTERFACE "$<BUILD_INTERFACE:${fmt_SOURCE_DIR}/include>"
              "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/CoolProp/third_party/fmt>"
  )
  if(MSVC)
    target_compile_options(coolprop_fmt_headers
                           INTERFACE "$<$<COMPILE_LANGUAGE:CXX>:/utf-8>")
  endif()

  set(_core_sources ${APP_SOURCES})
  if(NOT COOLPROP_LIBRARY_SOURCE STREQUAL "")
    list(APPEND _core_sources
         "${CMAKE_CURRENT_SOURCE_DIR}/${COOLPROP_LIBRARY_SOURCE}")
  endif()

  set(_shared_sources ${_core_sources})
  if((MSVC OR MINGW) AND COOLPROP_SHARED_LIBRARY)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/dev/CoolProp.rc.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/CoolProp.rc" @ONLY)
    list(APPEND _shared_sources "${CMAKE_CURRENT_BINARY_DIR}/CoolProp.rc")
  endif()

  set(_static_target "")
  set(_shared_target "")
  set(_canonical_target "")
  set(_export_targets "")

  if(COOLPROP_OBJECT_LIBRARY)
    add_library("${COOLPROP_LIBRARY_NAME}" OBJECT ${APP_SOURCES})
    set(_canonical_target "${COOLPROP_LIBRARY_NAME}")
    _coolprop_configure_core_target("${_canonical_target}" OBJECT)
    if(NOT TARGET CoolProp::Objects)
      add_library(CoolProp::Objects ALIAS "${_canonical_target}")
    endif()
  elseif(COOLPROP_STATIC_LIBRARY AND COOLPROP_SHARED_LIBRARY)
    set(_static_target "${COOLPROP_LIBRARY_NAME}_static")
    set(_shared_target "${COOLPROP_LIBRARY_NAME}_shared")
    set(_canonical_target "${COOLPROP_LIBRARY_NAME}")

    add_library("${_static_target}" STATIC ${_core_sources})
    add_library("${_shared_target}" SHARED ${_shared_sources}
                                           ${COOLPROP_LIBRARY_EXPORTS})
    add_library("${_canonical_target}" INTERFACE)

    if(_default_library STREQUAL "STATIC")
      target_link_libraries("${_canonical_target}" INTERFACE
                            "${_static_target}")
    else()
      target_link_libraries("${_canonical_target}" INTERFACE
                            "${_shared_target}")
    endif()

    set_target_properties("${_canonical_target}" PROPERTIES EXPORT_NAME
                                                             CoolProp)
    set_target_properties("${_static_target}" PROPERTIES EXPORT_NAME Static)
    set_target_properties("${_shared_target}" PROPERTIES EXPORT_NAME Shared
                                                         OUTPUT_NAME
                                                         "${COOLPROP_LIBRARY_NAME}")
    if(WIN32)
      set_target_properties("${_static_target}" PROPERTIES OUTPUT_NAME
                                                           "${COOLPROP_LIBRARY_NAME}_static")
    else()
      set_target_properties("${_static_target}" PROPERTIES OUTPUT_NAME
                                                           "${COOLPROP_LIBRARY_NAME}")
    endif()

    list(APPEND _export_targets "${_canonical_target}" "${_static_target}"
         "${_shared_target}")
  elseif(COOLPROP_STATIC_LIBRARY)
    set(_static_target "${COOLPROP_LIBRARY_NAME}")
    set(_canonical_target "${_static_target}")
    add_library("${_static_target}" STATIC ${_core_sources})
    set_target_properties("${_static_target}" PROPERTIES EXPORT_NAME CoolProp)
    list(APPEND _export_targets "${_static_target}")
  elseif(COOLPROP_SHARED_LIBRARY)
    set(_shared_target "${COOLPROP_LIBRARY_NAME}")
    set(_canonical_target "${_shared_target}")
    add_library("${_shared_target}" SHARED ${_shared_sources}
                                           ${COOLPROP_LIBRARY_EXPORTS})
    set_target_properties("${_shared_target}" PROPERTIES EXPORT_NAME CoolProp)
    list(APPEND _export_targets "${_shared_target}")
  endif()

  if(_static_target)
    _coolprop_configure_core_target("${_static_target}" STATIC)
    if(MSVC)
      set_target_properties("${_static_target}" PROPERTIES DEBUG_POSTFIX d)
    endif()
  endif()

  if(_shared_target)
    _coolprop_configure_core_target("${_shared_target}" SHARED)

    if(MSVC OR (MINGW AND DEFINED ENV{MSYSTEM}))
      set_target_properties("${_shared_target}" PROPERTIES DEBUG_POSTFIX d
                                                         PREFIX "")
    endif()

    if(MINGW AND DEFINED ENV{MSYSTEM})
      set_target_properties("${_shared_target}" PROPERTIES IMPORT_PREFIX ""
                                                         IMPORT_SUFFIX ".a")
      target_link_options("${_shared_target}" PRIVATE -static-libgcc
                                                      -static-libstdc++)
      target_link_libraries("${_shared_target}" PRIVATE -Wl,-Bstatic
                                                        -lwinpthread
                                                        -Wl,-Bdynamic)
    endif()
  endif()

  if(_canonical_target AND NOT TARGET CoolProp::CoolProp)
    add_library(CoolProp::CoolProp ALIAS "${_canonical_target}")
  endif()
  if(_static_target AND NOT TARGET CoolProp::Static)
    add_library(CoolProp::Static ALIAS "${_static_target}")
  endif()
  if(_shared_target AND NOT TARGET CoolProp::Shared)
    add_library(CoolProp::Shared ALIAS "${_shared_target}")
  endif()

  if(COOLPROP_INSTALL_LEGACY_LAYOUT AND _static_target)
    set(_legacy_static_dir
        "static_library/${CMAKE_SYSTEM_NAME}/${BITNESS}bit_${CMAKE_CXX_COMPILER_ID}_${CMAKE_CXX_COMPILER_VERSION}"
    )
    install(TARGETS "${_static_target}" DESTINATION "${_legacy_static_dir}")
    _coolprop_install_legacy_headers(static_library)
  endif()

  if(COOLPROP_INSTALL_LEGACY_LAYOUT AND _shared_target)
    if(MSVC AND (CMAKE_VS_PLATFORM_NAME STREQUAL "ARM64"
                 OR CMAKE_VS_PLATFORM_NAME STREQUAL "arm64"))
      set(_legacy_shared_dir
          "shared_library/${CMAKE_SYSTEM_NAME}/${BITNESS}bit__arm64")
    else()
      set(_legacy_shared_dir
          "shared_library/${CMAKE_SYSTEM_NAME}/${BITNESS}bit${CONVENTION}")
    endif()
    install(TARGETS "${_shared_target}" DESTINATION "${_legacy_shared_dir}")
    _coolprop_install_legacy_headers(shared_library)

    if(MSVC)
      add_custom_command(
        TARGET "${_shared_target}"
        POST_BUILD
        COMMAND dumpbin /EXPORTS "$<TARGET_FILE:${_shared_target}>" >
                "${CMAKE_CURRENT_BINARY_DIR}/exports.txt"
        VERBATIM)
      add_custom_command(
        TARGET "${_shared_target}"
        POST_BUILD
        COMMAND dumpbin /HEADERS "$<TARGET_FILE:${_shared_target}>" >
                "${CMAKE_CURRENT_BINARY_DIR}/headers.txt"
        VERBATIM)
      install(FILES "${CMAKE_CURRENT_BINARY_DIR}/exports.txt"
                    "${CMAKE_CURRENT_BINARY_DIR}/headers.txt"
              DESTINATION "${_legacy_shared_dir}")
    endif()
  endif()

  if(COOLPROP_INSTALL_CMAKE_PACKAGE AND _export_targets)
    install(
      TARGETS ${_export_targets} coolprop_eigen_headers coolprop_fmt_headers
      EXPORT CoolPropTargets
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")
    _coolprop_install_standard_headers()

    install(
      EXPORT CoolPropTargets
      FILE CoolPropTargets.cmake
      NAMESPACE CoolProp::
      DESTINATION "${COOLPROP_INSTALL_CMAKEDIR}")

    set(COOLPROP_PACKAGE_HAS_STATIC OFF)
    set(COOLPROP_PACKAGE_HAS_SHARED OFF)
    if(_static_target)
      set(COOLPROP_PACKAGE_HAS_STATIC ON)
    endif()
    if(_shared_target)
      set(COOLPROP_PACKAGE_HAS_SHARED ON)
    endif()
    configure_package_config_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/CoolPropConfig.cmake.in"
      "${CMAKE_CURRENT_BINARY_DIR}/CoolPropConfig.cmake"
      INSTALL_DESTINATION "${COOLPROP_INSTALL_CMAKEDIR}"
      PATH_VARS CMAKE_INSTALL_INCLUDEDIR)
    write_basic_package_version_file(
      "${CMAKE_CURRENT_BINARY_DIR}/CoolPropConfigVersion.cmake"
      VERSION "${COOLPROP_VERSION_NUMERIC}"
      COMPATIBILITY SameMajorVersion)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/CoolPropConfig.cmake"
                  "${CMAKE_CURRENT_BINARY_DIR}/CoolPropConfigVersion.cmake"
            DESTINATION "${COOLPROP_INSTALL_CMAKEDIR}")
  endif()

  message(STATUS "CoolProp library targets:")
  message(STATUS "  static: ${_static_target}")
  message(STATUS "  shared: ${_shared_target}")
  message(STATUS "  canonical: ${_canonical_target}")
  message(STATUS "  bitness: ${BITNESS}")
  if(CONVENTION STREQUAL "")
    message(STATUS "  calling convention: default")
  else()
    message(STATUS "  calling convention: ${CONVENTION}")
  endif()
endfunction()
