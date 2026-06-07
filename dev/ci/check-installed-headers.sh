#!/usr/bin/env bash
# Fail if forbidden internal headers are shipped in the installed tree.
#
# Enforced today: include/CoolProp/detail/json.h (which #includes
# nlohmann/json.hpp + valijson) must NOT be installed.  detail/rapidjson.h MUST
# still be installed (positive control — it is reachable via Configuration.h and
# is only removed at Phase Final; its presence proves the exclude is specific,
# not over-broad).
#
# Fail-closed: a failed `cmake --install`, an empty install tree, or a missing
# positive-control header is a FAILURE, never a vacuous pass (same discipline as
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

# Positive control: detail/rapidjson.h MUST ship (exclude must not be over-broad).
if ! find "$PREFIX" -path '*/detail/rapidjson.h' | grep -q .; then
    echo "FAIL: detail/rapidjson.h is absent from the install tree — exclude too broad or install layout changed." >&2
    exit 1
fi

# The assertion: detail/json.h must NOT ship.
LEAKED="$(find "$PREFIX" -path '*/detail/json.h' || true)"
if [ -n "$LEAKED" ]; then
    echo "FAIL: detail/json.h is shipped in the installed headers (it pulls nlohmann/json.hpp + valijson):" >&2
    printf '%s\n' "$LEAKED" >&2
    exit 1
fi

echo "OK: detail/json.h not installed; detail/rapidjson.h present (${NHEADERS} headers from ${BUILD_DIR})"

# TODO(superancillary de-leak): once superancillary.h no longer #includes
# nlohmann/json.hpp, broaden this to assert that NO shipped header pulls
# nlohmann/valijson (grep the installed tree), making it the install-side
# companion to the symbol-leak gate.
