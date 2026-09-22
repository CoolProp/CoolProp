# Dependencies for CoolProp managed via CPM.cmake
# CPM_SOURCE_CACHE can be overridden via environment variable (e.g. ~/.cache/CPM)
# to share the download cache across git worktrees and build directories.
# Without a stable cache location, FetchContent re-runs on every cmake configure,
# touching header timestamps and forcing a complete C++ rebuild each time.
if(NOT DEFINED CPM_SOURCE_CACHE AND "$ENV{CPM_SOURCE_CACHE}" STREQUAL "")
  set(CPM_SOURCE_CACHE "${CMAKE_CURRENT_LIST_DIR}/../.cpm_cache" CACHE PATH "CPM source cache")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/CPM.cmake")

# --- Offline (vendored) dependency sources, GH #3388 ----------------------
#
# Every distribution build system (sbuild, mock, OBS, makepkg in a clean
# chroot) builds inside a sandbox with no network at all, so CPM cannot
# download anything there.  The release tarball produced by
# dev/packaging/make-release-tarball.sh therefore ships a pre-unpacked copy of
# each dependency under externals/cpm/<name>/.
#
# Where such a directory exists we hand it to CPM as CPM_<name>_SOURCE, which
# makes CPMAddPackage use that directory verbatim and skip the download
# entirely.  A normal developer checkout has no externals/cpm/, so nothing
# changes for day-to-day work: CPM keeps fetching and caching as before.
set(COOLPROP_VENDORED_DEPS_DIR
    "${CMAKE_CURRENT_LIST_DIR}/../externals/cpm"
    CACHE PATH
    "Directory of pre-unpacked dependency sources used for offline builds")

# Normalised copy, used for the path comparison in the offline gate at the
# bottom of this file.  The default above contains a "/../" segment, and CPM
# echoes whatever it was handed straight back into CPM_PACKAGE_*_SOURCE_DIR.
get_filename_component(COOLPROP_VENDORED_DEPS_REAL
                       "${COOLPROP_VENDORED_DEPS_DIR}" REALPATH)

# Fail-closed switch for packaging builds.  With this ON, a dependency that is
# not vendored aborts the configure step instead of quietly reaching for the
# network, which is what would otherwise happen the moment somebody adds a
# dependency and forgets to re-run dev/packaging/vendor-deps.sh.
option(COOLPROP_REQUIRE_VENDORED_DEPS
       "Abort the configure step if any dependency is not vendored (offline packaging builds)"
       OFF)

# Point CPM at the vendored copy of one dependency, if the tarball shipped one.
# A macro rather than a function, so that the CPM_<name>_SOURCE variable it sets
# lands in the scope where CPMAddPackage is called a few lines below.
macro(coolprop_vendor_dependency _dep_name)
  if(NOT "${CPM_${_dep_name}_SOURCE}" STREQUAL "")
    # Somebody pointed CPM at a local checkout by hand; leave that alone.
    # Note this does NOT exempt it from the offline gate at the bottom of this
    # file: with COOLPROP_REQUIRE_VENDORED_DEPS=ON the gate still demands that
    # every resolved package sit under COOLPROP_VENDORED_DEPS_DIR, so a hand
    # override outside that directory is refused there.  That is deliberate -
    # an offline packaging build must not quietly pick up a developer's working
    # checkout - but it does mean the two options are not combinable.
    message(STATUS "CPM: ${_dep_name} overridden by CPM_${_dep_name}_SOURCE")
  elseif(IS_DIRECTORY "${COOLPROP_VENDORED_DEPS_DIR}/${_dep_name}")
    set(CPM_${_dep_name}_SOURCE "${COOLPROP_VENDORED_DEPS_DIR}/${_dep_name}")
    message(STATUS "CPM: ${_dep_name} from vendored source (offline)")
  elseif(COOLPROP_REQUIRE_VENDORED_DEPS)
    message(
      FATAL_ERROR
        "COOLPROP_REQUIRE_VENDORED_DEPS is ON but '${COOLPROP_VENDORED_DEPS_DIR}/${_dep_name}' "
        "does not exist, so ${_dep_name} would be downloaded.  Run "
        "dev/packaging/vendor-deps.sh to populate it, or build from a release tarball.")
  endif()
endmacro()

# Every dependency fetched below.  The optional Windows- and Mathematica-only
# packages are deliberately absent: they are not part of a Linux package build.
set(COOLPROP_CPM_DEPENDENCIES
    Eigen
    msgpack-c
    nlohmann_json
    valijson
    IF97
    REFPROP_headers
    boost_headers
    multicomplex
    fmt)
if(COOLPROP_CATCH_MODULE)
  list(APPEND COOLPROP_CPM_DEPENDENCIES Catch2)
endif()

foreach(_dep IN LISTS COOLPROP_CPM_DEPENDENCIES)
  coolprop_vendor_dependency(${_dep})
endforeach()

# ── Core header-only deps ──────────────────────────────────────────────────

# Use tarball instead of git clone: gitlab.com cloning over CMake's
# FetchContent custom-build step is flaky on Windows runners — when the
# clone retries internally MSBuild still propagates the first failure
# (see PR #2890/#2905 for similar wheel-build flake mitigations).
CPMAddPackage(
  NAME Eigen
  VERSION 5.0.1
  URL https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz
  URL_HASH SHA256=e9c326dc8c05cd1e044c71f30f1b2e34a6161a3b6ecf445d56b53ff1669e3dec
  DOWNLOAD_ONLY YES   # header-only; skip Eigen's own CMake targets
)

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

CPMAddPackage(
  NAME fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG        12.0.0
  OPTIONS
    "FMT_INSTALL OFF"
    "FMT_TEST OFF"
    "FMT_DOC OFF"
  DOWNLOAD_ONLY  YES   # CoolProp uses fmt in header-only mode via FMT_HEADER_ONLY
)

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

# --- Offline gate, GH #3388 -----------------------------------------------
#
# coolprop_vendor_dependency() above only knows the names listed in
# COOLPROP_CPM_DEPENDENCIES, so on its own it would let a dependency that
# somebody adds later download from the network while still reporting an
# offline build.  This check does not depend on that list: CPM_PACKAGES holds
# every package CPMAddPackage actually resolved, whatever its name, and each one
# has to have come out of the vendored directory.
if(COOLPROP_REQUIRE_VENDORED_DEPS)
  # Refuse a vendor root that cannot discriminate.
  #
  # Test the RAW variable for emptiness, not the resolved one: measured on
  # CMake 3.28, get_filename_component(REALPATH) turns an empty path into
  # CMAKE_CURRENT_SOURCE_DIR rather than leaving it empty, so a check on the
  # resolved value never fires.  With the root silently equal to the source
  # tree, the prefix test below accepts an ordinary in-tree build/_deps/...
  # download as "vendored" (verified: it matches at position 0) and the gate
  # passes on exactly the thing it exists to catch.
  #
  # The source and binary directories are rejected for the same reason even
  # when spelled out in full, and "/" because every path starts with it.  The
  # macro above happens to abort first in today's code, but a gate must not
  # depend on another check firing.
  # Reject ANCESTRY, not equality.  The check below is a path-prefix test, so a
  # root that merely sits ABOVE the build tree makes every in-tree download
  # match it -- /tmp, or the parent of the source tree, and so on.  Testing for
  # equality with the source and binary directories closes two spellings and
  # leaves the whole class open (measured: with the root at any ancestor of the
  # build directory, a package downloaded into <build>/_deps/ is accepted and
  # the gate reports "Offline build verified").
  #
  # Both sides are REALPATH-resolved before comparing, because
  # COOLPROP_VENDORED_DEPS_REAL is resolved while CMAKE_CURRENT_*_DIR keeps
  # whatever spelling the caller used; a source tree reached through a symlink
  # would otherwise never compare equal.  "/" is kept as its own case purely for
  # the message: it resolves so that the prefix test below rejects every
  # package anyway (fail-closed), but saying "at or above the source tree" up
  # front beats ten confusing per-package errors.
  get_filename_component(_cp_src_real "${CMAKE_CURRENT_SOURCE_DIR}" REALPATH)
  get_filename_component(_cp_bin_real "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)
  string(FIND "${_cp_src_real}/" "${COOLPROP_VENDORED_DEPS_REAL}/" _cp_src_under)
  string(FIND "${_cp_bin_real}/" "${COOLPROP_VENDORED_DEPS_REAL}/" _cp_bin_under)
  if("${COOLPROP_VENDORED_DEPS_DIR}" STREQUAL ""
     OR "${COOLPROP_VENDORED_DEPS_REAL}" STREQUAL ""
     OR "${COOLPROP_VENDORED_DEPS_REAL}" STREQUAL "/"
     OR _cp_src_under EQUAL 0
     OR _cp_bin_under EQUAL 0)
    message(
      FATAL_ERROR
        "Offline build requested but COOLPROP_VENDORED_DEPS_DIR ('${COOLPROP_VENDORED_DEPS_DIR}') "
        "resolves to '${COOLPROP_VENDORED_DEPS_REAL}', which is at or above the source or build "
        "directory and therefore cannot be told apart from an ordinary in-tree download.  Set it "
        "to the directory holding the vendored sources (externals/cpm in a release tarball).")
  endif()
  if(NOT CPM_PACKAGES)
    message(
      FATAL_ERROR
        "Offline build requested but CPM reported no packages at all, so this check "
        "cannot confirm anything.  CPM.cmake may have changed how it records packages; "
        "fix this gate before shipping a package built from this tree.")
  endif()
  foreach(_pkg IN LISTS CPM_PACKAGES)
    get_filename_component(_pkg_dir "${CPM_PACKAGE_${_pkg}_SOURCE_DIR}" REALPATH)
    string(FIND "${_pkg_dir}" "${COOLPROP_VENDORED_DEPS_REAL}/" _vendor_pos)
    if(NOT _vendor_pos EQUAL 0)
      message(
        FATAL_ERROR
          "Offline build requested (COOLPROP_REQUIRE_VENDORED_DEPS=ON) but package "
          "'${_pkg}' was resolved from '${_pkg_dir}', which is not inside "
          "'${COOLPROP_VENDORED_DEPS_REAL}'.  Run dev/packaging/vendor-deps.sh, or "
          "build from a release tarball made by dev/packaging/make-release-tarball.sh.")
    endif()
  endforeach()
  message(STATUS "Offline build verified: all ${CPM_PACKAGES} came from vendored sources")
endif()
