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
# notes); allow override via COOLPROP_CLANG_TIDY.  Version 19+ is required
# -- see the check after this block.
if [ -n "${COOLPROP_CLANG_TIDY:-}" ]; then
  CLANG_TIDY="$COOLPROP_CLANG_TIDY"
elif command -v clang-tidy >/dev/null 2>&1; then
  CLANG_TIDY="$(command -v clang-tidy)"
elif [ "$(uname -s)" = "Darwin" ]; then
  # On macOS, prefer the newest available Homebrew llvm.  Apple's libc++
  # uses C++23 builtins (__builtin_clzg, __builtin_ctzg, __builtin_addcb,
  # ...) that clang-tidy 18 doesn't know about, so parsing <bitset> /
  # <charconv> emits a wave of bogus clang-diagnostic-error noise.
  # clang-tidy 21+ handles them.
  if [ -x "/opt/homebrew/opt/llvm@21/bin/clang-tidy" ]; then
    CLANG_TIDY="/opt/homebrew/opt/llvm@21/bin/clang-tidy"
  elif [ -x "/opt/homebrew/opt/llvm/bin/clang-tidy" ]; then
    CLANG_TIDY="/opt/homebrew/opt/llvm/bin/clang-tidy"
  elif [ -x "/opt/homebrew/opt/llvm@18/bin/clang-tidy" ]; then
    # Found only so the version check below rejects it by name, rather
    # than reporting "not found" and skipping.
    CLANG_TIDY="/opt/homebrew/opt/llvm@18/bin/clang-tidy"
  fi
fi
if [ -z "${CLANG_TIDY:-}" ]; then
  echo "warning: clang-tidy not on PATH and not installed via Homebrew; skipping" >&2
  echo "         install with: brew install llvm   (macOS)" >&2
  echo "         (or set COOLPROP_CLANG_TIDY=/path/to/clang-tidy to point at a specific binary)" >&2
  exit 0
fi

# .clang-tidy uses ExcludeHeaderFilterRegex, which clang-tidy 19 introduced.
# clang-tidy 18 does not skip that unknown key -- it prints a parse error,
# discards the whole config and runs its DEFAULT checks, still exiting 0, so
# an old binary would "pass" against the wrong rule set.  Fail rather than skip: the message
# starts with "error: " so preflight counts it as a finding.
# `|| true` is deliberate: a binary that will not run leaves CT_MAJOR empty
# and lands in the error branch below; without it, pipefail + set -e would
# exit here with no message, which preflight's log grep would read as clean.
CT_MAJOR="$("$CLANG_TIDY" --version 2>/dev/null | sed -nE 's/.*version ([0-9]+)\..*/\1/p' | head -1 || true)"
if [ -z "$CT_MAJOR" ] || [ "$CT_MAJOR" -lt 19 ]; then
  echo "error: $CLANG_TIDY is version '${CT_MAJOR:-unknown}'; clang-tidy >= 19 is required (.clang-tidy uses ExcludeHeaderFilterRegex)" >&2
  echo "       install with: brew install llvm   (macOS), or set COOLPROP_CLANG_TIDY=/path/to/clang-tidy-19+" >&2
  exit 1
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

# Header scope.  .clang-tidy reports findings in project headers; that is
# what a changed-lines run (-line-filter, as CI's clang-tidy-diff does) wants,
# since the line filter confines them to lines the change touched.  A
# whole-file run has no such bound: one src/*.cpp would drag in ~600
# pre-existing findings from the headers it includes and fail every push.  So
# without a line filter, keep the scope to the file itself (a regex that
# matches no path), which is what the old, broken '*' filter did in effect.
HAS_LINE_FILTER=0
for arg in "$@"; do
  case "$arg" in -line-filter | --line-filter | -line-filter=* | --line-filter=*) HAS_LINE_FILTER=1 ;; esac
done
if [ "$HAS_LINE_FILTER" = 0 ]; then
  EXTRA_ARGS+=("--header-filter=^\$")
fi

# ${EXTRA_ARGS[@]+...}: macOS's /bin/bash 3.2 aborts on an empty array under
# set -u ("unbound variable"), and that abort carries no "error: " line.
exec "$CLANG_TIDY" -p "$BUILD_DIR" ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"} "$@"
