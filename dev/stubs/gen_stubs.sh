#!/usr/bin/env bash
#
# Generate the committed CoolProp.pyi type stub for the Cython wrapper.
#
# This is a DEV/CI tool, not part of the package build.  It needs a modern
# Cython (>= 3.1) + stubgen-pyx, unrelated to the Cython that *compiles*
# CoolProp (0.29 era).  To keep the build untouched we use an isolated,
# version-pinned virtualenv under dev/stubs/.venv-stubgen.
#
# CoolProp.pyx ``include``s HumidAirProp.pyx and AbstractState.pyx, so a single
# stubgen run over CoolProp.pyx yields the whole CoolProp.CoolProp module stub
# (AbstractState, State, PropsSI, ... all in one file) — which matches how the
# extension is actually built (one CoolProp.<soabi>.so).
#
# Output is deterministic so CI can assert it is in sync:
#
#     dev/stubs/gen_stubs.sh
#     git diff --exit-code -- 'wrappers/Python/CoolProp/*.pyi'
#
# A nonzero diff means a .pyx changed without regenerating the stub.
#
set -euo pipefail

# Pin generator versions: the drift gate compares byte-for-byte, so the tool
# versions must be reproducible.  Bump deliberately (and regenerate) together.
CYTHON_VERSION="${CYTHON_VERSION:-3.2.5}"
STUBGEN_PYX_VERSION="${STUBGEN_PYX_VERSION:-0.2.11}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"
PKG_DIR="$REPO_ROOT/wrappers/Python/CoolProp"
VENV="$HERE/.venv-stubgen"
PYTHON="${PYTHON:-python3}"

# --- isolated, pinned toolchain ------------------------------------------------
# Reuse the venv only if it already holds the pinned versions.  Gating on mere
# existence would silently reuse a stale toolchain after a version bump and
# defeat the byte-for-byte determinism the drift gate relies on.
STAMP="$VENV/.versions"
WANT="cython==$CYTHON_VERSION stubgen-pyx==$STUBGEN_PYX_VERSION"
if [ ! -x "$VENV/bin/stubgen-pyx" ] || [ "$(cat "$STAMP" 2>/dev/null)" != "$WANT" ]; then
    echo ">> creating stubgen venv at $VENV (cython==$CYTHON_VERSION, stubgen-pyx==$STUBGEN_PYX_VERSION)"
    rm -rf "$VENV"
    "$PYTHON" -m venv "$VENV"
    "$VENV/bin/pip" install -q --upgrade pip
    "$VENV/bin/pip" install -q "cython==$CYTHON_VERSION" "stubgen-pyx==$STUBGEN_PYX_VERSION"
    printf '%s' "$WANT" > "$STAMP"
fi
STUBGEN="$VENV/bin/stubgen-pyx"
VPY="$VENV/bin/python"

# --- generate into a throwaway copy of the package -----------------------------
# stubgen needs the package laid out as ``CoolProp/`` (so ``from CoolProp.foo
# cimport ...`` resolves) and CoolProp.pyx must parse.  We copy the sources to a
# temp dir and strip trailing ``;`` (a C-ism that the modern Cython parser
# rejects) WITHOUT touching the committed sources.
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/CoolProp"
cp -f "$PKG_DIR"/*.pyx "$PKG_DIR"/*.pxd "$TMP/CoolProp/" 2>/dev/null || true
for f in "$TMP/CoolProp"/*.pyx; do
    sed -i.bak 's/;[[:space:]]*$//' "$f" && rm -f "$f.bak"
done

RAW="$TMP/CoolProp.raw.pyi"
echo ">> stubgen-pyx  CoolProp/CoolProp.pyx (expands includes)"
( cd "$TMP" && "$STUBGEN" --file CoolProp/CoolProp.pyx --output-file "$RAW" --continue-on-error . )

echo ">> postprocess (+ hand-written overloads) -> CoolProp.pyi"
"$VPY" "$HERE/postprocess.py" "$RAW" "$PKG_DIR/CoolProp.pyi" "$HERE/overloads_coolprop.pyi"

echo ">> done"
