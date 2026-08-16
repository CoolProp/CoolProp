cmake_minimum_required(VERSION 3.14)

option(COOLPROP_C_ONLY "Configure a consumer with no C++ compiler enabled" OFF)
if(COOLPROP_C_ONLY)
  project(CoolPropCMakeConsumerSmoke LANGUAGES C)
else()
  project(CoolPropCMakeConsumerSmoke LANGUAGES C CXX)
endif()

# Keep the consumer deliberately below CoolProp's public C++ requirement.  A
# correctly exported CoolProp target must promote the C++ smoke executable to
# C++17 through its INTERFACE_COMPILE_FEATURES usage requirement.
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
if(NOT COOLPROP_C_ONLY)
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

set(
  COOLPROP_EXPECTED_LINKAGE
  ""
  CACHE STRING "Expected installed CoolProp linkage: STATIC, SHARED, or empty")
set_property(CACHE COOLPROP_EXPECTED_LINKAGE PROPERTY STRINGS "" STATIC SHARED)
string(TOUPPER "${COOLPROP_EXPECTED_LINKAGE}" _coolprop_expected_linkage)

set(COOLPROP_SOURCE_DIR "" CACHE PATH
    "Optional CoolProp checkout to consume with add_subdirectory")

if(NOT _coolprop_expected_linkage STREQUAL ""
   AND NOT _coolprop_expected_linkage STREQUAL "STATIC"
   AND NOT _coolprop_expected_linkage STREQUAL "SHARED")
  message(
    FATAL_ERROR
      "COOLPROP_EXPECTED_LINKAGE must be STATIC, SHARED, or empty; got '${COOLPROP_EXPECTED_LINKAGE}'"
  )
endif()

if(COOLPROP_C_ONLY AND NOT _coolprop_expected_linkage STREQUAL "SHARED")
  message(
    FATAL_ERROR
      "COOLPROP_C_ONLY exercises the portable shared-library C ABI and requires COOLPROP_EXPECTED_LINKAGE=SHARED"
  )
endif()

if(COOLPROP_SOURCE_DIR)
  get_filename_component(_coolprop_source_dir "${COOLPROP_SOURCE_DIR}" ABSOLUTE)
  set(_coolprop_prefix_sentinel
      "${CMAKE_CURRENT_BINARY_DIR}/parent-install-prefix")
  set(CMAKE_INSTALL_PREFIX
      "${_coolprop_prefix_sentinel}"
      CACHE PATH "Consumer install prefix" FORCE)
  set(COOLPROP_STATIC_LIBRARY ON CACHE BOOL "" FORCE)
  set(COOLPROP_SHARED_LIBRARY OFF CACHE BOOL "" FORCE)
  set(COOLPROP_OBJECT_LIBRARY OFF CACHE BOOL "" FORCE)
  set(COOLPROP_INSTALL_CMAKE_PACKAGE OFF CACHE BOOL "" FORCE)
  add_subdirectory("${_coolprop_source_dir}"
                   "${CMAKE_CURRENT_BINARY_DIR}/coolprop-build")

  if(NOT "${CMAKE_INSTALL_PREFIX}" STREQUAL "${_coolprop_prefix_sentinel}")
    message(
      FATAL_ERROR
        "CoolProp changed its parent project's CMAKE_INSTALL_PREFIX from '${_coolprop_prefix_sentinel}' to '${CMAKE_INSTALL_PREFIX}'"
    )
  endif()
  if(COOLPROP_INSTALL_LEGACY_LAYOUT)
    message(
      FATAL_ERROR
        "CoolProp enabled its legacy release install layout inside a parent project"
    )
  endif()

  if(MSVC AND NOT CMAKE_MSVC_RUNTIME_LIBRARY_DEFAULT)
    get_target_property(_coolprop_source_compile_options CoolProp
                        COMPILE_OPTIONS)
    string(FIND "${_coolprop_source_compile_options}" "/MD"
                _coolprop_msvc_runtime_index)
    if(_coolprop_msvc_runtime_index EQUAL -1)
      message(
        FATAL_ERROR
          "CoolProp did not provide a target-local MSVC runtime fallback for a CMP0091-OLD parent"
      )
    endif()
  endif()

  set(_coolprop_required_variant_targets CoolProp::Static)
else()
  if(_coolprop_expected_linkage STREQUAL "STATIC")
    find_package(CoolProp 8 CONFIG REQUIRED COMPONENTS Static)
    set(_coolprop_required_variant_targets CoolProp::Static)
  elseif(_coolprop_expected_linkage STREQUAL "SHARED")
    find_package(CoolProp 8 CONFIG REQUIRED COMPONENTS Shared)
    set(_coolprop_required_variant_targets CoolProp::Shared)
  else()
    find_package(CoolProp 8 CONFIG REQUIRED COMPONENTS Static Shared)
    set(_coolprop_required_variant_targets CoolProp::Static CoolProp::Shared)
  endif()
endif()

if(NOT TARGET CoolProp::CoolProp)
  message(
    FATAL_ERROR
      "The installed CoolProp package did not define the canonical CoolProp::CoolProp target"
  )
endif()

foreach(_coolprop_required_target IN LISTS _coolprop_required_variant_targets)
  if(NOT TARGET "${_coolprop_required_target}")
    message(
      FATAL_ERROR
        "The CoolProp package reported the requested component as available but did not define ${_coolprop_required_target}"
    )
  endif()
endforeach()

# CoolProp::CoolProp is the required public API. Link through all requested
# linkage-specific targets too, so the documented component targets cannot
# silently disappear.
set(_coolprop_provider_targets CoolProp::CoolProp
                               ${_coolprop_required_variant_targets})
list(REMOVE_DUPLICATES _coolprop_provider_targets)

enable_testing()

foreach(_coolprop_provider_target IN LISTS _coolprop_provider_targets)
  get_target_property(_coolprop_target_type "${_coolprop_provider_target}" TYPE)
  set(_coolprop_effective_linkage "")

  if(_coolprop_provider_target STREQUAL "CoolProp::Static")
    set(_coolprop_effective_linkage STATIC)
    if(NOT _coolprop_target_type STREQUAL "STATIC_LIBRARY"
       AND NOT _coolprop_target_type STREQUAL "INTERFACE_LIBRARY")
      message(
        FATAL_ERROR
          "CoolProp::Static exists but has incompatible type ${_coolprop_target_type}"
      )
    endif()
  elseif(_coolprop_provider_target STREQUAL "CoolProp::Shared")
    set(_coolprop_effective_linkage SHARED)
    if(NOT _coolprop_target_type STREQUAL "SHARED_LIBRARY"
       AND NOT _coolprop_target_type STREQUAL "INTERFACE_LIBRARY")
      message(
        FATAL_ERROR
          "CoolProp::Shared exists but has incompatible type ${_coolprop_target_type}"
      )
    endif()
  elseif(_coolprop_target_type STREQUAL "STATIC_LIBRARY")
    set(_coolprop_effective_linkage STATIC)
  elseif(_coolprop_target_type STREQUAL "SHARED_LIBRARY")
    set(_coolprop_effective_linkage SHARED)
  elseif(_coolprop_target_type STREQUAL "INTERFACE_LIBRARY")
    set(_coolprop_effective_linkage "${_coolprop_expected_linkage}")
  else()
    message(
      FATAL_ERROR
        "${_coolprop_provider_target} has unsupported type ${_coolprop_target_type}"
    )
  endif()

  # A canonical selector or single-variant compatibility target may be an
  # INTERFACE wrapper. Resolve one level to find its actual provider type and,
  # on Windows, the DLL that must be copied beside the smoke executable.
  set(_coolprop_runtime_target "${_coolprop_provider_target}")
  if(_coolprop_target_type STREQUAL "INTERFACE_LIBRARY")
    get_target_property(
      _coolprop_interface_links "${_coolprop_provider_target}"
      INTERFACE_LINK_LIBRARIES)
    foreach(_coolprop_link_candidate IN LISTS _coolprop_interface_links)
      if(TARGET "${_coolprop_link_candidate}")
        get_target_property(_coolprop_link_candidate_type
                            "${_coolprop_link_candidate}" TYPE)
        if(NOT _coolprop_link_candidate_type STREQUAL "INTERFACE_LIBRARY")
          set(_coolprop_runtime_target "${_coolprop_link_candidate}")
          break()
        endif()
      endif()
    endforeach()
  endif()

  get_target_property(_coolprop_runtime_target_type
                      "${_coolprop_runtime_target}" TYPE)
  if(_coolprop_effective_linkage STREQUAL "")
    if(_coolprop_runtime_target_type STREQUAL "STATIC_LIBRARY")
      set(_coolprop_effective_linkage STATIC)
    elseif(_coolprop_runtime_target_type STREQUAL "SHARED_LIBRARY")
      set(_coolprop_effective_linkage SHARED)
    endif()
  endif()

  if(_coolprop_provider_target STREQUAL "CoolProp::CoolProp"
     AND _coolprop_expected_linkage STREQUAL "STATIC"
     AND NOT _coolprop_effective_linkage STREQUAL "STATIC")
    message(
      FATAL_ERROR
        "Expected a static CoolProp package, but CoolProp::CoolProp has type ${_coolprop_target_type}"
    )
  endif()
  if(_coolprop_provider_target STREQUAL "CoolProp::CoolProp"
     AND _coolprop_expected_linkage STREQUAL "SHARED"
     AND NOT _coolprop_effective_linkage STREQUAL "SHARED")
    message(
      FATAL_ERROR
        "Expected a shared CoolProp package, but CoolProp::CoolProp has type ${_coolprop_target_type}"
    )
  endif()

  string(REPLACE "::" "_" _coolprop_target_id "${_coolprop_provider_target}")
  string(TOLOWER "${_coolprop_target_id}" _coolprop_target_id)

  set(_coolprop_c_executable "coolprop_c_smoke_${_coolprop_target_id}")
  if(_coolprop_effective_linkage STREQUAL "STATIC")
    set(_coolprop_c_api_source c_api_cpp.cpp)
  else()
    set(_coolprop_c_api_source c_api.c)
  endif()
  add_executable("${_coolprop_c_executable}" "${_coolprop_c_api_source}")
  target_link_libraries("${_coolprop_c_executable}" PRIVATE
                        "${_coolprop_provider_target}")

  set(_coolprop_smoke_executables "${_coolprop_c_executable}")
  set(_coolprop_smoke_tests "${_coolprop_target_id}.c-api")
  add_test(NAME "${_coolprop_target_id}.c-api"
           COMMAND "${_coolprop_c_executable}")

  # The shared library's portable contract is the C ABI in CoolPropLib.h.
  # Its C++ symbols are not exported on every platform (notably Windows), so
  # exercise the documented C++ API only for an installed static archive.
  if(_coolprop_effective_linkage STREQUAL "STATIC" AND NOT COOLPROP_C_ONLY)
    set(_coolprop_cxx_executable "coolprop_cxx_smoke_${_coolprop_target_id}")
    add_executable("${_coolprop_cxx_executable}" cpp_api.cpp)
    target_link_libraries("${_coolprop_cxx_executable}" PRIVATE
                          "${_coolprop_provider_target}")
    list(APPEND _coolprop_smoke_executables "${_coolprop_cxx_executable}")
    list(APPEND _coolprop_smoke_tests "${_coolprop_target_id}.cxx-api")
    add_test(NAME "${_coolprop_target_id}.cxx-api"
             COMMAND "${_coolprop_cxx_executable}")
  endif()

  # Windows does not search the linked DLL's install directory at runtime.
  # Copying TARGET_FILE is CMake-3.14-compatible.
  if(WIN32 AND _coolprop_runtime_target_type STREQUAL "SHARED_LIBRARY")
    foreach(_coolprop_executable IN LISTS _coolprop_smoke_executables)
      add_custom_command(
        TARGET "${_coolprop_executable}"
        POST_BUILD
        COMMAND
          "${CMAKE_COMMAND}" -E copy_if_different
          "$<TARGET_FILE:${_coolprop_runtime_target}>"
          "$<TARGET_FILE_DIR:${_coolprop_executable}>"
        VERBATIM)
    endforeach()
  endif()

  set_tests_properties(${_coolprop_smoke_tests}
                       PROPERTIES LABELS "cmake-consumer")

  message(
    STATUS
      "Testing ${_coolprop_provider_target} (${_coolprop_target_type}, ${_coolprop_effective_linkage})"
  )
endforeach()
