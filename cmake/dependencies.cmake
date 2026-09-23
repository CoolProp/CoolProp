# Dependencies for CoolProp managed via CPM.cmake
# CPM_SOURCE_CACHE can be overridden via environment variable (e.g. ~/.cache/CPM)
# to share the download cache across git worktrees and build directories.
# Without a stable cache location, FetchContent re-runs on every cmake configure,
# touching header timestamps and forcing a complete C++ rebuild each time.
if(NOT DEFINED CPM_SOURCE_CACHE AND "$ENV{CPM_SOURCE_CACHE}" STREQUAL "")
  set(CPM_SOURCE_CACHE "${CMAKE_CURRENT_LIST_DIR}/../.cpm_cache"
      CACHE PATH "CPM source cache")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/CPM.cmake")

option(COOLPROP_VENDOR_THIRD_PARTY
       "Bundle Eigen and fmt headers with the installed CoolProp package" ON)

# ── Core header-only deps ──────────────────────────────────────────────────

# Use tarball instead of git clone: gitlab.com cloning over CMake's
# FetchContent custom-build step is flaky on Windows runners — when the
# clone retries internally MSBuild still propagates the first failure
# (see PR #2890/#2905 for similar wheel-build flake mitigations).
if(COOLPROP_VENDOR_THIRD_PARTY)
  # Vendoring requires source headers and licenses, even when CPM normally
  # prefers installed packages. OFF below is the system-package alternative.
  if(CPM_Eigen_SOURCE)
    get_filename_component(Eigen_SOURCE_DIR "${CPM_Eigen_SOURCE}" ABSOLUTE)
  else()
    CPMAddPackage(
      NAME Eigen
      VERSION 5.0.1
      URL https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz
      URL_HASH SHA256=e9c326dc8c05cd1e044c71f30f1b2e34a6161a3b6ecf445d56b53ff1669e3dec
      DOWNLOAD_ONLY YES
      FORCE YES
    )
  endif()
  set(COOLPROP_EIGEN_INCLUDE_DIRS "${Eigen_SOURCE_DIR}")
else()
  find_package(Eigen3 CONFIG REQUIRED)
  if(Eigen3_VERSION VERSION_LESS 3.4)
    message(FATAL_ERROR "CoolProp requires Eigen 3.4 or newer; found ${Eigen3_VERSION}.")
  endif()
  get_target_property(COOLPROP_EIGEN_INCLUDE_DIRS Eigen3::Eigen
                      INTERFACE_INCLUDE_DIRECTORIES)
endif()

CPMAddPackage(
  NAME msgpack-c
  GIT_REPOSITORY https://github.com/msgpack/msgpack-c
  GIT_TAG        919908742b4fdbc575e77fe1a8657e70c9573c44
  DOWNLOAD_ONLY  YES
)

# nlohmann/json — replacement for rapidjson (GH: RapidJSON→nlohmann migration).
# Header-only; included via the hidden-visibility wrapper include/CoolProp/detail/json.h.
CPMAddPackage(
  NAME nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.12.0
  DOWNLOAD_ONLY  YES   # header-only; we only need the include dir
)

# Valijson — header-only JSON-Schema (draft-7) validator that validates an
# nlohmann::json instance directly via its bundled adapter. Used for runtime
# validation of user-supplied PC-SAFT / cubic fluids.
#
# Use the release tarball, NOT a git clone: valijson's repo carries test-only
# submodules (googletest, yaml-cpp, nlohmann-json, rapidjson, …) whose deeply-
# nested test-fixture filenames exceed Windows' MAX_PATH and break a recursive
# git checkout (notably the Tauri GUI build, which configures from an already-
# deep path).  GIT_SUBMODULES "" does NOT help here because CMP0097 defaults to
# OLD under this project's cmake_minimum_required, where empty means "all".  The
# GitHub source archive contains no submodule contents, so the tarball sidesteps
# both that and the Windows git-clone flakiness (same rationale as Eigen above).
CPMAddPackage(
  NAME valijson
  VERSION 1.0.6
  URL https://github.com/tristanpenman/valijson/archive/refs/tags/v1.0.6.tar.gz
  URL_HASH SHA256=bf0839de19510ff7792d8a8aca94ea11a288775726b36c4c9a2662651870f8da
  DOWNLOAD_ONLY YES   # header-only; we only need include/
)

CPMAddPackage(
  NAME IF97
  GIT_REPOSITORY https://github.com/CoolProp/IF97
  GIT_TAG        7aaced024a702f0985474bf293cdaae9c8d06521
  DOWNLOAD_ONLY  YES
)

CPMAddPackage(
  NAME REFPROP_headers
  GIT_REPOSITORY https://github.com/CoolProp/REFPROP-headers.git
  GIT_TAG        b4faab1b73911c32c4b69c526c7e92f74edb67de
  DOWNLOAD_ONLY  YES
)

CPMAddPackage(
  NAME boost_headers
  GIT_REPOSITORY https://github.com/CoolProp/boost-headers.git
  GIT_TAG        c68104660ca4bd80d0d5cb34c4eba0cf5bab3f73
  DOWNLOAD_ONLY  YES
)

CPMAddPackage(
  NAME multicomplex
  GIT_REPOSITORY https://github.com/usnistgov/multicomplex
  GIT_TAG        39bf9ca52c7882ff0788bb9087c7548ebd8fba4c
  DOWNLOAD_ONLY  YES
)

# ── fmt (header-only use; disable fmt's own tests/docs) ───────────────────

if(COOLPROP_VENDOR_THIRD_PARTY)
  if(CPM_fmt_SOURCE)
    get_filename_component(fmt_SOURCE_DIR "${CPM_fmt_SOURCE}" ABSOLUTE)
  else()
    CPMAddPackage(
      NAME fmt
      GIT_REPOSITORY https://github.com/fmtlib/fmt.git
      GIT_TAG        12.0.0
      DOWNLOAD_ONLY YES
      FORCE YES
    )
  endif()
  set(COOLPROP_FMT_INCLUDE_DIRS "${fmt_SOURCE_DIR}/include")

  # Diagnose unavailable CPM sources here, not halfway through installation.
  foreach(_header_or_license IN ITEMS
      "${Eigen_SOURCE_DIR}/Eigen/Core"
      "${Eigen_SOURCE_DIR}/unsupported/Eigen/Polynomials"
      "${Eigen_SOURCE_DIR}/unsupported/Eigen/src/Polynomials"
      "${Eigen_SOURCE_DIR}/COPYING.README"
      "${Eigen_SOURCE_DIR}/COPYING.MPL2"
      "${Eigen_SOURCE_DIR}/COPYING.APACHE"
      "${Eigen_SOURCE_DIR}/COPYING.BSD"
      "${Eigen_SOURCE_DIR}/COPYING.MINPACK"
      "${fmt_SOURCE_DIR}/include/fmt/format.h"
      "${fmt_SOURCE_DIR}/LICENSE")
    if(NOT EXISTS "${_header_or_license}")
      message(FATAL_ERROR
              "Cannot vendor Eigen/fmt: '${_header_or_license}' is missing. Provide their CPM sources or set COOLPROP_VENDOR_THIRD_PARTY=OFF to use installed packages.")
    endif()
  endforeach()
else()
  find_package(fmt CONFIG REQUIRED)
  get_target_property(COOLPROP_FMT_INCLUDE_DIRS fmt::fmt-header-only
                      INTERFACE_INCLUDE_DIRECTORIES)
endif()

# ── Catch2 (testing only) ──────────────────────────────────────────────────
# Fetched and add_subdirectory'd only when COOLPROP_CATCH_MODULE is ON.

if(COOLPROP_CATCH_MODULE)
  CPMAddPackage(
    NAME Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2
    GIT_TAG        v3.8.0
  )
endif()

# ── Windows packaging helpers (optional) ──────────────────────────────────

if(COOLPROP_WINDOWS_PACKAGE)
  # This revision carries the EesUserLib64 task, which installs the 64-bit EES
  # library into the Userlib64 folder.
  CPMAddPackage(
    NAME ExcelAddinInstaller
    GIT_REPOSITORY https://github.com/CoolProp/ExcelAddinInstaller.git
    GIT_TAG        7fba5c452ed2490d830b516fe453c6411d2066ad
    DOWNLOAD_ONLY  YES
  )
endif()

# ── Mathematica (optional) ─────────────────────────────────────────────────

if(COOLPROP_MATHEMATICA_MODULE)
  CPMAddPackage(
    NAME FindMathematica
    GIT_REPOSITORY https://github.com/sakra/FindMathematica
    GIT_TAG        4.2.0
    DOWNLOAD_ONLY  YES
  )
endif()
