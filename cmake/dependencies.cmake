# Dependencies for CoolProp managed via CPM.cmake
# Set CPM_SOURCE_CACHE in your environment (e.g. ~/.cache/CPM) to share the
# download cache across git worktrees and build directories.

include("${CMAKE_CURRENT_LIST_DIR}/CPM.cmake")

# ── Core header-only deps ──────────────────────────────────────────────────

CPMAddPackage(
  NAME Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG        3.4.0
  DOWNLOAD_ONLY  YES   # header-only; skip Eigen's own CMake targets
)

CPMAddPackage(
  NAME msgpack-c
  GIT_REPOSITORY https://github.com/msgpack/msgpack-c
  GIT_TAG        919908742b4fdbc575e77fe1a8657e70c9573c44
  DOWNLOAD_ONLY  YES
)

CPMAddPackage(
  NAME rapidjson
  GIT_REPOSITORY https://github.com/Tencent/rapidjson.git
  GIT_TAG        24b5e7a8b27f42fa16b96fc70aade9106cf7102f
  DOWNLOAD_ONLY  YES
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
  CPMAddPackage(
    NAME ExcelAddinInstaller
    GIT_REPOSITORY https://github.com/CoolProp/ExcelAddinInstaller.git
    GIT_TAG        db8ce41cdb02079a2d9242ea08f3633e8a1d38b0
    DOWNLOAD_ONLY  YES
  )
endif()

# ── Mathematica (optional) ─────────────────────────────────────────────────

if(COOLPROP_MATHEMATICA_MODULE)
  CPMAddPackage(
    NAME FindMathematica
    GIT_REPOSITORY https://github.com/sakra/FindMathematica
    GIT_TAG        99d8c95eed1a517f4a06511ea40b2cee6477a8c1
    DOWNLOAD_ONLY  YES
  )
endif()
