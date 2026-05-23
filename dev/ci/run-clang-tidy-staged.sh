#!/usr/bin/env bash
# Wrapper for the pre-commit clang-tidy hook (CoolProp-2uw.6).
#
# Skips gracefully when:
#   - clang-tidy is not on PATH AND not findable via Homebrew LLVM paths
#   - $COOLPROP_BUILD_DIR/compile_commands.json is missing (no build/ yet)
#
# Otherwise runs `clang-tidy -p $BUILD_DIR <files>` and inherits all check
# selection + WarningsAsErrors from the repo-root .clang-tidy config.
#
# macOS note (per issue #2926): Homebrew clang-tidy needs explicit
# -isysroot + libc++ include args, otherwise standard headers don't
# resolve and analysis silently degrades (bugprone-infinite-loop,
# bugprone-virtual-near-miss false positives).  We auto-add those here.
#
# Override the build directory:
#   COOLPROP_BUILD_DIR=build_debug pre-commit run --hook-stage manual clang-tidy

set -euo pipefail

BUILD_DIR="${COOLPROP_BUILD_DIR:-build}"
COMPDB="$BUILD_DIR/compile_commands.json"

# Locate clang-tidy.  Prefer PATH (Linux/CI), fall back to common
# Homebrew install locations on macOS (per issue #2926 reproduction
# notes).  Pin to llvm@18 by default for CI parity; allow override via
# COOLPROP_CLANG_TIDY.
if [ -n "${COOLPROP_CLANG_TIDY:-}" ]; then
  CLANG_TIDY="$COOLPROP_CLANG_TIDY"
elif command -v clang-tidy >/dev/null 2>&1; then
  CLANG_TIDY="$(command -v clang-tidy)"
elif [ "$(uname -s)" = "Darwin" ]; then
  # On macOS, prefer the newest available Homebrew llvm.  Apple's libc++
  # uses C++23 builtins (__builtin_clzg, __builtin_ctzg, __builtin_addcb,
  # ...) that clang-tidy 18 doesn't know about, so parsing <bitset> /
  # <charconv> emits a wave of bogus clang-diagnostic-error noise.
  # clang-tidy 21+ handles them.  CI runs on Linux with matching headers,
  # so the 18 pin there is fine.
  if [ -x "/opt/homebrew/opt/llvm@21/bin/clang-tidy" ]; then
    CLANG_TIDY="/opt/homebrew/opt/llvm@21/bin/clang-tidy"
  elif [ -x "/opt/homebrew/opt/llvm/bin/clang-tidy" ]; then
    CLANG_TIDY="/opt/homebrew/opt/llvm/bin/clang-tidy"
  elif [ -x "/opt/homebrew/opt/llvm@18/bin/clang-tidy" ]; then
    CLANG_TIDY="/opt/homebrew/opt/llvm@18/bin/clang-tidy"
  fi
fi
if [ -z "${CLANG_TIDY:-}" ]; then
  echo "warning: clang-tidy not on PATH and not installed via Homebrew; skipping" >&2
  echo "         install with: brew install llvm   (macOS)" >&2
  echo "         (or set COOLPROP_CLANG_TIDY=/path/to/clang-tidy to point at a specific binary)" >&2
  exit 0
fi

if [ ! -f "$COMPDB" ]; then
  echo "warning: $COMPDB not found; skipping clang-tidy" >&2
  echo "         configure cmake first: cmake -G Ninja -B $BUILD_DIR -S ." >&2
  echo "         (or set COOLPROP_BUILD_DIR=<path> to point at your existing build dir)" >&2
  exit 0
fi

# Per issue #2926: macOS Homebrew clang-tidy can't find <vector>, <string>,
# etc. without explicit sysroot + libc++ include hints.  CI on Linux is
# unaffected (system headers resolve naturally).  Detect macOS via `uname`
# and add the args only there to keep Linux invocations unchanged.
EXTRA_ARGS=()
if [ "$(uname -s)" = "Darwin" ] && command -v xcrun >/dev/null 2>&1; then
  SDK_PATH="$(xcrun --show-sdk-path 2>/dev/null || true)"
  if [ -n "$SDK_PATH" ]; then
    EXTRA_ARGS+=(
      "--extra-arg=-isysroot$SDK_PATH"
      "--extra-arg=-stdlib=libc++"
      "--extra-arg=-isystem$SDK_PATH/usr/include/c++/v1"
    )
  fi
fi

exec "$CLANG_TIDY" -p "$BUILD_DIR" "${EXTRA_ARGS[@]}" "$@"
