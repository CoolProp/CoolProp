#!/usr/bin/env bash
# Wrapper for the pre-commit clang-tidy hook (CoolProp-2uw.6).
#
# Skips gracefully when:
#   - clang-tidy is not on PATH (contributor without LLVM installed)
#   - $COOLPROP_BUILD_DIR/compile_commands.json is missing (no build/ yet)
#
# Otherwise runs `clang-tidy -p $BUILD_DIR <files>` and inherits all check
# selection + WarningsAsErrors from the repo-root .clang-tidy config.
#
# Override the build directory:
#   COOLPROP_BUILD_DIR=build_debug pre-commit run --hook-stage manual clang-tidy

set -euo pipefail

BUILD_DIR="${COOLPROP_BUILD_DIR:-build}"
COMPDB="$BUILD_DIR/compile_commands.json"

if ! command -v clang-tidy >/dev/null 2>&1; then
  echo "warning: clang-tidy not on PATH; skipping (install LLVM 18+ for parity with CI)" >&2
  exit 0
fi

if [ ! -f "$COMPDB" ]; then
  echo "warning: $COMPDB not found; skipping clang-tidy" >&2
  echo "         configure cmake first: cmake -G Ninja -B $BUILD_DIR -S ." >&2
  echo "         (or set COOLPROP_BUILD_DIR=<path> to point at your existing build dir)" >&2
  exit 0
fi

exec clang-tidy -p "$BUILD_DIR" "$@"
