#!/usr/bin/env bash
# Fail if forbidden internal headers are shipped in the installed tree.
#
# Two assertions:
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
# detail/CachedElement.h, detail/filepaths.h, detail/msgpack.h,
# detail/PlatformDetermination.h, detail/state_capi.h, detail/atomic_write.h.
# This control does NOT name detail/rapidjson.h because Phase Final will
# delete that header; any of the above is sufficient.
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

# Assertion 1: detail/json.h must NOT ship.
LEAKED="$(find "$PREFIX" -path '*/detail/json.h' || true)"
if [ -n "$LEAKED" ]; then
    echo "FAIL: detail/json.h is shipped in the installed headers (it pulls nlohmann/json.hpp + valijson):" >&2
    printf '%s\n' "$LEAKED" >&2
    exit 1
fi

# Assertion 2: no shipped header may #include nlohmann or valijson.
# grep -El exits 1 on no match (the PASS case) and 0 on match (the FAIL case).
# The `|| true` prevents set -e from aborting on the no-match exit code; the
# actual gate is the `-n` test below, which is fail-closed.
LEAKING_INCLUDES="$(grep -rEl '#[[:space:]]*include[[:space:]]*[<"]((nlohmann|valijson)/)' \
    "$PREFIX" --include='*.h' --include='*.hpp' || true)"
if [ -n "$LEAKING_INCLUDES" ]; then
    echo "FAIL: installed header(s) #include nlohmann/valijson (must stay internal):" >&2
    printf '%s\n' "$LEAKING_INCLUDES" >&2
    exit 1
fi

echo "OK: detail/json.h not installed; no shipped header pulls nlohmann/valijson; detail/ control header present (${NHEADERS} headers from ${BUILD_DIR})"
