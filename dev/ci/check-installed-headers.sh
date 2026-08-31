#!/usr/bin/env bash
# Fail if forbidden internal headers are shipped in the installed tree.
#
# Three assertions:
#
#  1. detail/json.h must NOT ship.  That header #includes nlohmann/json.hpp
#     and valijson; it is internal-only and is excluded by CMake's install
#     rules.  Its presence in the installed tree is a hard failure.
#
#  2. No shipped header may #include nlohmann or valijson.  This is the
#     install-side companion to the symbol-leak gate (check-json-symbols.sh):
#     even if detail/json.h is correctly excluded, a newly added public header
#     that transitively pulls nlohmann/valijson would be a regression.
#
# Generic positive control: at least ONE installed header matching
# */detail/*.h (other than detail/json.h) must be present.  This proves the
# CMake EXCLUDE is specific to json.h and is not over-broad (e.g. wiping the
# whole detail/ subtree).  Candidates that satisfy this today:
# detail/configuration_keys.h, detail/strings.h, detail/tools.h,
# detail/CachedElement.h, detail/filepaths.h,
# detail/PlatformDetermination.h, detail/state_capi.h, detail/atomic_write.h.
# (detail/json.h and detail/msgpack.h are NOT candidates -- both are excluded.)
# This control does NOT name detail/rapidjson.h (that header was deleted in
# Phase Final); any of the above is sufficient.
#
# Fail-closed: a failed `cmake --install`, an empty install tree, or a
# violated assertion is a FAILURE, never a vacuous pass (same discipline as
# check-json-symbols.sh — no `|| true` that masks a real failure).
set -euo pipefail

BUILD_DIR="${1:-build_shared}"
if [ ! -d "$BUILD_DIR" ]; then
    echo "FAIL: build dir '$BUILD_DIR' not found (configure+build a shared build first)" >&2
    exit 2
fi

PREFIX="$(mktemp -d)"
trap 'rm -rf "$PREFIX"' EXIT

if ! cmake --install "$BUILD_DIR" --prefix "$PREFIX" >/tmp/install-headers-check.log 2>&1; then
    echo "FAIL: cmake --install '$BUILD_DIR' failed (see /tmp/install-headers-check.log)" >&2
    exit 1
fi

# Guard against a vacuous pass: a real install ships many headers.
NHEADERS="$(find "$PREFIX" -name '*.h' | grep -c . || true)"
if [ "$NHEADERS" -eq 0 ]; then
    echo "FAIL: install tree under $PREFIX has no headers — wrong build dir or empty install? Cannot validate." >&2
    exit 1
fi

# Generic positive control: at least one detail/*.h other than detail/json.h
# must ship, proving the CMake EXCLUDE is specific to json.h, not over-broad.
DETAIL_OTHER="$(find "$PREFIX" -path '*/detail/*.h' ! -name 'json.h' | head -1)"
if [ -z "$DETAIL_OTHER" ]; then
    echo "FAIL: no detail/*.h (other than json.h) found in the install tree — CMake EXCLUDE too broad or install layout changed." >&2
    exit 1
fi

# Assertion 1: the internal-only detail/json.h and detail/msgpack.h must NOT
# ship.  detail/json.h pulls nlohmann/json.hpp + valijson; detail/msgpack.h
# pulls msgpack.hpp.  None of those third-party headers are installed, so
# shipping either makes the installed tree non-self-contained for downstream.
LEAKED="$(find "$PREFIX" \( -path '*/detail/json.h' -o -path '*/detail/msgpack.h' -o -name 'CPmsgpack.h' \) || true)"
if [ -n "$LEAKED" ]; then
    echo "FAIL: an internal-only header (detail/json.h, detail/msgpack.h, or the CPmsgpack.h shim) is shipped in the installed headers:" >&2
    printf '%s\n' "$LEAKED" >&2
    exit 1
fi

# Assertion 2: no shipped header may #include nlohmann, valijson, msgpack or
# boost (the internal-only third-party libraries that are NOT installed).
# grep -El exits 1 on no match (the PASS case) and 0 on match (the FAIL case).
# The `|| true` prevents set -e from aborting on the no-match exit code; the
# actual gate is the `-n` test below, which is fail-closed.
LEAKING_INCLUDES="$(grep -rEl '#[[:space:]]*include[[:space:]]*[<"]((nlohmann|valijson|boost)/|msgpack)' \
    "$PREFIX" --include='*.h' --include='*.hpp' || true)"
if [ -n "$LEAKING_INCLUDES" ]; then
    echo "FAIL: installed header(s) #include nlohmann/valijson/msgpack/boost (must stay internal):" >&2
    printf '%s\n' "$LEAKING_INCLUDES" >&2
    exit 1
fi

# Assertion 3: every installed CoolProp-owned header must compile STANDALONE.
# Vendored third-party headers are intentionally excluded: their internal
# headers are not part of CoolProp's self-containedness contract.  This catches
# the general non-self-contained case the grep above cannot -- e.g. a shipped
# header that #includes a non-installed header by some other path (the old
# detail/msgpack.h pulling msgpack.hpp slipped past the nlohmann/valijson grep).
# Compile each CoolProp-owned *.h as its own translation unit.  That surface
# depends only on Eigen (the fluids/numerics/superancillary tiers) and on fmt
# unless NO_FMTLIB is defined; resolve Eigen from the build's CPM cache and
# define NO_FMTLIB so fmt is not required.  Boost is deliberately NOT on the
# path: the superancillary rootfinder routes through an out-of-line helper
# (src/superancillary.cpp), so a regression that re-introduces a boost include
# into an installed header fails this compile (and assertion 2 above).
CXX_BIN="${CXX:-c++}"
CACHE="$BUILD_DIR/CMakeCache.txt"
if [ ! -f "$CACHE" ]; then
    echo "FAIL: $CACHE not found -- '$BUILD_DIR' is not a configured CMake build dir; cannot resolve Eigen for the self-containedness check." >&2
    exit 1
fi
EIGEN_DIR="$(sed -n 's/^CPM_PACKAGE_Eigen_SOURCE_DIR:INTERNAL=//p' "$CACHE" | head -1)"
if [ -z "$EIGEN_DIR" ] || [ ! -d "$EIGEN_DIR" ]; then
    echo "FAIL: could not resolve the Eigen include dir from $CACHE -- cannot run the self-containedness check (fail-closed)." >&2
    exit 1
fi
INSTALL_INCLUDEDIR="$(sed -n 's/^CMAKE_INSTALL_INCLUDEDIR:[^=]*=//p' "$CACHE" | head -1)"
if [ -z "$INSTALL_INCLUDEDIR" ]; then
    echo "FAIL: could not resolve CMAKE_INSTALL_INCLUDEDIR from $CACHE -- cannot select the canonical package include root." >&2
    exit 1
fi
INC_ROOT="$PREFIX/$INSTALL_INCLUDEDIR"
if [ ! -f "$INC_ROOT/CoolProp/CoolProp.h" ]; then
    echo "FAIL: canonical installed header $INC_ROOT/CoolProp/CoolProp.h is missing -- cannot validate self-containedness." >&2
    exit 1
fi
SC_LOG=/tmp/installed-headers-selfcontained.log
: >"$SC_LOG"
SC_FAIL=0
SC_HEADERS=0
while IFS= read -r hdr; do
    SC_HEADERS=$((SC_HEADERS + 1))
    rel="${hdr#"$INC_ROOT"/}"
    if ! printf '#include "%s"\n' "$rel" | "$CXX_BIN" -std=c++17 -fsyntax-only \
            -DNO_FMTLIB -DCOOLPROP_NO_DEPRECATED_HEADER_WARNINGS \
            -I"$INC_ROOT" -I"$EIGEN_DIR" -x c++ - >>"$SC_LOG" 2>&1; then
        echo "FAIL: installed header is not self-contained: $rel" >&2
        SC_FAIL=$((SC_FAIL + 1))
    fi
done < <(find "$INC_ROOT" -path "$INC_ROOT/CoolProp/third_party" -prune -o -type f -name '*.h' -print)
if [ "$SC_HEADERS" -eq 0 ]; then
    echo "FAIL: canonical include root $INC_ROOT contains no CoolProp-owned headers -- cannot validate self-containedness." >&2
    exit 1
fi
if [ "$SC_FAIL" -ne 0 ]; then
    echo "FAIL: $SC_FAIL installed header(s) do not compile standalone (see $SC_LOG)." >&2
    echo "      A shipped header that #includes a non-installed header (internal or third-party)" >&2
    echo "      is broken for downstream consumers." >&2
    exit 1
fi

echo "OK: internal json/msgpack headers not installed; no shipped header pulls nlohmann/valijson/msgpack/boost; all ${SC_HEADERS} CoolProp-owned installed headers compile standalone (-DNO_FMTLIB, Eigen-only on the path) from ${BUILD_DIR}"
