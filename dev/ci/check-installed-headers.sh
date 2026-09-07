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

# Assertion 3: every installed header must compile STANDALONE.  This catches the
# general non-self-contained case the grep above cannot -- e.g. a shipped header
# that #includes a non-installed header by some other path (the old
# detail/msgpack.h pulling msgpack.hpp slipped past the nlohmann/valijson grep).
# Compile each installed *.h as its own translation unit.  The installed surface
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
INC_ROOT="$(find "$PREFIX" -path '*/include/CoolProp/CoolProp.h' | head -1)"
INC_ROOT="${INC_ROOT%/CoolProp/CoolProp.h}"
if [ -z "$INC_ROOT" ] || [ ! -d "$INC_ROOT" ]; then
    echo "FAIL: could not locate the installed include root (no CoolProp/CoolProp.h) -- cannot validate self-containedness." >&2
    exit 1
fi
SC_LOG=/tmp/installed-headers-selfcontained.log
: >"$SC_LOG"

# The sweep is ~93 INDEPENDENT -fsyntax-only invocations, so it parallelises
# cleanly across cores.  Serially it was 47.7 s -- the single largest item in
# ./dev/ci/preflight.sh's 122 s floor, and ~1000x the 0.05 s that `cmake
# --install` itself costs (bd CoolProp-5vun).
SC_JOBS="${COOLPROP_HEADER_CHECK_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"
# Clamp the auto-detected value.  getconf reports HOST cores, not a container's
# cgroup quota, and a syntax-only TU peaks around 180 MB here -- so on a large
# machine an unbounded -P invites the OOM killer, which takes cc1plus rather than
# the bash child.  The child would then see a non-zero compiler status and report
# "installed header is not self-contained", i.e. a false regression rather than
# an infrastructure error.  8 keeps nearly all the speedup (the sweep is only ~93
# items); an explicit COOLPROP_HEADER_CHECK_JOBS is honoured as given -- including
# 0, which xargs reads as "unbounded" and which therefore re-opens the very OOM
# case the clamp exists to prevent.  That is the point of an override, but set it
# deliberately.
if [ -z "${COOLPROP_HEADER_CHECK_JOBS:-}" ] && [ "$SC_JOBS" -gt 8 ] 2>/dev/null; then
    SC_JOBS=8
fi
SC_TMP="$(mktemp -d)"
# Supersede the earlier trap so the scratch dir is cleaned up too.  Both paths
# are quoted; PREFIX and SC_TMP are mktemp output, but the quoting is what keeps
# a future edit from turning this into an unbounded rm.
trap 'rm -rf "$PREFIX" "$SC_TMP"' EXIT

SC_TOTAL="$(find "$INC_ROOT" -name '*.h' | grep -c . || true)"
[ -n "$SC_TOTAL" ] || SC_TOTAL=0
if [ "$SC_TOTAL" -eq 0 ]; then
    echo "FAIL: no *.h found under the installed include root ($INC_ROOT) -- the self-containedness sweep would verify nothing." >&2
    exit 1
fi

export CXX_BIN INC_ROOT EIGEN_DIR SC_TMP

# One header per child.  Each child writes to its OWN files: two parallel
# writers appending to a single log interleave mid-line, which would corrupt the
# very diagnostics this gate prints.  A child records `ok.` or `fail.` either
# way and exits 0, so xargs' own status stays reserved for infrastructure failure
# (couldn't spawn, bad interpreter) rather than a compile result.  The one
# deviation is safe: if the marker path itself cannot be written (ENAMETOOLONG on
# a pathologically deep header) the child exits non-zero and surfaces through
# xargs as exactly that kind of infrastructure failure, and the count guard below
# catches it independently.
#
# -print0/-0 keeps a path containing spaces as one argument.  The command is
# passed inline rather than via `export -f`: exported bash functions do not
# survive a version mismatch between the parent shell and the `bash` on PATH.
SC_XARGS_RC=0
find "$INC_ROOT" -name '*.h' -print0 \
    | xargs -0 -P "$SC_JOBS" -n1 bash -c '
        # SHELLOPTS is not exported, so this child does NOT inherit the parent'"'"'s
        # pipefail.  Without it a failing printf would hand the compiler an empty
        # translation unit, which exits 0 -- recording a header as self-contained
        # without ever having compiled it.
        set -o pipefail
        hdr="$1"
        rel="${hdr#"$INC_ROOT"/}"
        # Percent-encoding, and the % must be escaped BEFORE / is mapped onto
        # it -- that ordering is what makes this reversible, so two distinct
        # headers cannot land on one marker file and silently lose a result.
        # Mapping every separator onto _ is not injective (a plain tr collides
        # "a/b.h" with "a_b.h"); neither is doubling _ first, which still
        # collides "a_/b.h" with "a/_b.h" -- both become "a___b.h".
        slug="$(printf "%s" "$rel" | sed "s/%/%25/g; s#/#%2F#g")"
        if printf "#include \"%s\"\n" "$rel" | "$CXX_BIN" -std=c++17 -fsyntax-only \
                -DNO_FMTLIB -DCOOLPROP_NO_DEPRECATED_HEADER_WARNINGS \
                -I"$INC_ROOT" -I"$EIGEN_DIR" -x c++ - >"$SC_TMP/out.$slug" 2>&1; then
            printf "%s\n" "$rel" >"$SC_TMP/ok.$slug"
        else
            printf "%s\n" "$rel" >"$SC_TMP/fail.$slug"
        fi
    ' _ || SC_XARGS_RC=$?
if [ "$SC_XARGS_RC" -ne 0 ]; then
    echo "FAIL: the parallel header sweep could not run (xargs exit $SC_XARGS_RC) -- treating as a failure rather than a pass." >&2
    exit 1
fi

# Collect the per-header compiler output in a stable, deterministic order.  It is
# sorted by encoded header path -- NOT the directory order find happened to yield
# when this ran serially -- so the log is reproducible run to run.  LC_ALL=C is
# what makes that true ACROSS environments: the default collation ignores
# punctuation, so a C-locale CI runner and a UTF-8 local shell would otherwise
# order the same slugs differently.
find "$SC_TMP" -name 'out.*' | LC_ALL=C sort | while IFS= read -r out; do
    cat "$out" >>"$SC_LOG"
done

# Parentheses around the -o alternation are belt-and-braces: with only -name
# predicates here both forms behave identically, but they keep the grouping
# explicit so adding another test (say -type f) later cannot silently bind to
# one branch only.
SC_CHECKED="$(find "$SC_TMP" \( -name 'ok.*' -o -name 'fail.*' \) | grep -c . || true)"
[ -n "$SC_CHECKED" ] || SC_CHECKED=0
SC_FAIL="$(find "$SC_TMP" -name 'fail.*' | grep -c . || true)"
[ -n "$SC_FAIL" ] || SC_FAIL=0

# Fail-closed guard with no serial equivalent: a child killed before it wrote
# either marker (OOM, SIGKILL) leaves no failure behind, so counting only
# `fail.` markers would report a clean pass over a sweep that never finished.
# Every header found must have produced exactly one result.
if [ "$SC_CHECKED" -ne "$SC_TOTAL" ]; then
    echo "FAIL: the header sweep recorded $SC_CHECKED result(s) for $SC_TOTAL header(s) -- some check did not complete, so this is not a pass." >&2
    exit 1
fi

# The success line below names shared_library/CoolPropLib.h as *the* header
# outside the sweep.  Assert that rather than trusting it: if an install rule
# ever ships another header outside INC_ROOT, the message would understate the
# gap -- the same overstatement this change just fixed in the other direction.
#
# Compare the IDENTITY of the unswept set, not merely its size.  A count-only
# check ("exactly one header is unswept") passes unchanged if CoolPropLib.h stops
# being installed and some different header takes its place outside the sweep:
# the total still differs by one, the gate still goes green, and the success line
# then names a file that is no longer the exemption.  Enumerating the paths is
# what makes the claim true rather than merely arithmetically consistent.
#
# A mismatch is a hard failure, which also means the gate stops asserting a stale
# exemption once the sweep widens to cover CoolPropLib.h (bd CoolProp-1n5g).
# Identity is asserted on the FILE NAME, not the full installed path.  What the
# gate needs to know is "the one header we do not compile is still CoolPropLib.h,
# not some other header that quietly stopped being swept".  Pinning the whole
# relative path would additionally hard-code one install layout, so a change to
# the install prefix or directory structure would fail the gate for a reason that
# has nothing to do with header coverage.  The full path is printed either way,
# so a layout change is still visible in the log.
SC_UNSWEPT_EXPECTED="CoolPropLib.h"
find "$PREFIX" -name '*.h' ! -path "$INC_ROOT/*" | LC_ALL=C sort >"$SC_TMP/unswept"
# -printf '%P' would be tidier but is GNU-only; strip the prefix in the shell so
# this keeps working under BSD/macOS find.
SC_UNSWEPT_REL="$(while IFS= read -r h; do printf '%s\n' "${h#"$PREFIX"/}"; done <"$SC_TMP/unswept")"
SC_UNSWEPT_N="$(printf '%s' "$SC_UNSWEPT_REL" | grep -c . || true)"
[ -n "$SC_UNSWEPT_N" ] || SC_UNSWEPT_N=0
if [ "$SC_UNSWEPT_N" -ne 1 ] || [ "${SC_UNSWEPT_REL##*/}" != "$SC_UNSWEPT_EXPECTED" ]; then
    echo "FAIL: expected exactly one installed *.h outside the include root, named '$SC_UNSWEPT_EXPECTED'; found $SC_UNSWEPT_N:" >&2
    printf '%s\n' "${SC_UNSWEPT_REL:-(none)}" | sed 's/^/        /' >&2
    echo "      Either a header is newly shipped outside the swept tree -- in which case it is going unchecked -- or the sweep has widened and this assertion plus the success message need updating (bd CoolProp-1n5g)." >&2
    exit 1
fi
SC_UNSWEPT_SHOWN="$SC_UNSWEPT_REL"

if [ "$SC_FAIL" -ne 0 ]; then
    find "$SC_TMP" -name 'fail.*' | LC_ALL=C sort | while IFS= read -r marker; do
        echo "FAIL: installed header is not self-contained: $(cat "$marker")" >&2
    done
    echo "FAIL: $SC_FAIL installed header(s) do not compile standalone (see $SC_LOG)." >&2
    echo "      A shipped header that #includes a non-installed header (internal or third-party)" >&2
    echo "      is broken for downstream consumers." >&2
    exit 1
fi

# ${SC_TOTAL}, not ${NHEADERS}: the sweep compiles the headers under the include
# root, while NHEADERS counts every *.h in the install prefix.  The difference is
# exactly one file -- the top-level shared_library/CoolPropLib.h, the single-
# include C/DLL API header, which sits outside INC_ROOT and is therefore NOT
# swept (bd CoolProp-1n5g).  Claiming NHEADERS here overstated the
# check by that one header.
echo "OK: internal json/msgpack headers not installed; no shipped header pulls nlohmann/valijson/msgpack/boost; all ${SC_TOTAL} headers under the include root compile standalone (of ${NHEADERS} installed *.h; ${SC_UNSWEPT_SHOWN} is not swept) with -DNO_FMTLIB, Eigen-only on the path, from ${BUILD_DIR}"
