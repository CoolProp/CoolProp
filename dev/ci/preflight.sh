#!/usr/bin/env bash
#
# preflight.sh — local pre-push quality gate for CoolProp.
#
# Runs the same checks CI runs on the diff between the current HEAD and
# its upstream / origin/master, so a passing preflight predicts a green
# CI.  Closes CoolProp-6r6.
#
# Designed to be the body of a pre-push git hook (which `git commit
# --no-verify` does NOT bypass — only `git push --no-verify` does).
# Invoke directly to spot-check at any time:
#
#   ./dev/ci/preflight.sh                # check against origin/master
#   ./dev/ci/preflight.sh --base=HEAD~1  # check against an earlier ref
#   ./dev/ci/preflight.sh --skip=cppcheck,clang-tidy   # subset
#   ./dev/ci/preflight.sh --skip=json-symbols          # subset
#
# Tools resolved at runtime:
#   - clang-format     : uvx clang-format@<version-from-.pre-commit-config>
#   - cppcheck         : system binary (graceful skip if missing)
#   - clang-tidy       : delegates to dev/ci/run-clang-tidy-staged.sh
#                         (which already graceful-skips when clang-tidy or
#                          compile_commands.json isn't around)
#   - semgrep          : uvx semgrep with p/cpp + p/security-audit rulesets
#   - Catch2 tests     : ./build_catch/CatchTestRunner with tag scope
#                        auto-selected from the changed paths
#
# Exit codes: 0 = all checks passed (or were skipped intentionally),
# non-zero = at least one check failed.  When invoked as a pre-push hook
# a non-zero exit will block the push; agents using `--no-verify` to skip
# the hook are doing so deliberately and should run preflight separately.

set -euo pipefail

# ---------- arg parsing ----------------------------------------------

BASE_REF="origin/master"
SKIP_CHECKS=""
for arg in "$@"; do
    case "$arg" in
        --base=*) BASE_REF="${arg#*=}" ;;
        --skip=*) SKIP_CHECKS="${arg#*=}" ;;
        --help|-h)
            sed -n '2,30p' "$0"
            exit 0
            ;;
        *)
            echo "preflight: unknown arg '$arg' (use --base=<ref> or --skip=<csv>)" >&2
            exit 2
            ;;
    esac
done

skip_check() {
    [[ ",$SKIP_CHECKS," == *",$1,"* ]]
}

# ---------- locate repo + cd to root ---------------------------------

REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

# ---------- resolve diff set -----------------------------------------

# Files changed on this branch vs base, restricted to extensions the
# C-family checks care about.  Filters out deleted files (.cpp.cpp at
# the same path would otherwise be checked even though it no longer
# exists on disk).
CHANGED_CPP="$(git diff --name-only --diff-filter=ACMR "$BASE_REF"...HEAD -- '*.cpp' '*.h' '*.hpp' '*.cc' '*.cxx' || true)"
# Also pick up uncommitted changes in the working tree — preflight is
# meant to gate pushes, but agents often run it mid-edit too.
UNSTAGED_CPP="$(git diff --name-only --diff-filter=ACMR -- '*.cpp' '*.h' '*.hpp' '*.cc' '*.cxx' || true)"
ALL_CPP="$(printf '%s\n%s\n' "$CHANGED_CPP" "$UNSTAGED_CPP" | sort -u | grep -v '^$' || true)"

# All paths changed (any extension) — used for tag auto-selection.
ALL_PATHS="$(git diff --name-only "$BASE_REF"...HEAD; git diff --name-only)"
ALL_PATHS="$(printf '%s\n' "$ALL_PATHS" | sort -u | grep -v '^$' || true)"

# ---------- pretty helpers -------------------------------------------

PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0
FAIL_NAMES=()

step() {
    printf '\n\033[1m== %s ==\033[0m\n' "$1"
}

ok() {
    printf '\033[32m✓ %s\033[0m\n' "$1"
    PASS_COUNT=$((PASS_COUNT + 1))
}

fail() {
    printf '\033[31m✗ %s\033[0m\n' "$1"
    FAIL_COUNT=$((FAIL_COUNT + 1))
    FAIL_NAMES+=("$1")
}

skip() {
    printf '\033[33m- %s (skipped: %s)\033[0m\n' "$1" "$2"
    SKIP_COUNT=$((SKIP_COUNT + 1))
}

# ---------- check 1: clang-format ------------------------------------

step "clang-format dry-run vs $BASE_REF"
if skip_check clang-format; then
    skip "clang-format" "--skip=clang-format"
elif [ -z "$ALL_CPP" ]; then
    skip "clang-format" "no C/C++ files in diff"
else
    # Read the pinned version from .pre-commit-config.yaml so CI + hook +
    # preflight all use the same binary.
    CF_VER="$(awk '/mirrors-clang-format/{f=1} f && /^[[:space:]]*rev:/{print $2; exit}' .pre-commit-config.yaml | tr -d 'v"' || true)"
    if [ -z "$CF_VER" ]; then
        CF_VER="18.1.8"
    fi
    # uvx caches the binary; first invocation downloads, subsequent are
    # ~instant.  --dry-run --Werror returns non-zero on any reformatting.
    if uvx clang-format@"$CF_VER" --dry-run --Werror $ALL_CPP 2>&1 | grep -q .; then
        fail "clang-format (run: uvx clang-format@$CF_VER -i <files>)"
    else
        ok "clang-format ($CF_VER, $(printf '%s\n' "$ALL_CPP" | wc -l | tr -d ' ') file(s))"
    fi
fi

# ---------- check 2: build CatchTestRunner ---------------------------

step "build CatchTestRunner"
if skip_check build; then
    skip "build" "--skip=build"
elif [ ! -d build_catch ]; then
    # Auto-configure on first run.  Same flags as CI.
    cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON >/dev/null 2>&1 || {
        fail "build (cmake configure failed; run cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON manually)"
    }
fi
if ! skip_check build && [ -d build_catch ]; then
    if cmake --build build_catch --target CatchTestRunner -j8 2>&1 | tee /tmp/preflight-build.log | tail -5 | grep -qE "error:|FAILED"; then
        fail "build (see /tmp/preflight-build.log)"
    else
        ok "build CatchTestRunner"
    fi
fi

# ---------- check 2b: JSON symbol-leak gate (shared library) ---------
#
# Headline enforced invariant of the RapidJSON->nlohmann migration: no
# nlohmann/valijson symbol may be exported from CoolProp's shared
# products.  Visibility attributes only take effect once linked into a
# shared object, so this MUST inspect a .so/.dylib — never the static
# .a (the Catch runner links the archive).  preflight builds the Catch
# runner against a static/object lib, so we maintain a dedicated
# build_shared dir for this gate.  RapidJSON is intentionally NOT in the
# default pattern (it is still exported today, by design, mid-migration).
step "JSON symbol leak (shared library)"
if skip_check json-symbols; then
    skip "json-symbols" "--skip=json-symbols"
elif skip_check build; then
    skip "json-symbols" "--skip=build (shared build needed)"
else
    if [ ! -d build_shared ]; then
        cmake -B build_shared -S . -DCOOLPROP_SHARED_LIBRARY=ON -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1 || true
    fi
    if [ -d build_shared ]; then
        cmake --build build_shared -j8 >/tmp/preflight-shared-build.log 2>&1 || true
    fi
    SHARED_LIB="$(find build_shared \( -name 'libCoolProp.so' -o -name 'libCoolProp.dylib' \) 2>/dev/null | head -1 || true)"
    if [ -z "$SHARED_LIB" ] || [ ! -f "$SHARED_LIB" ]; then
        skip "json-symbols" "no shared CoolProp library found (see /tmp/preflight-shared-build.log)"
    elif ./dev/ci/check-json-symbols.sh "$SHARED_LIB"; then
        ok "json-symbols (no nlohmann/valijson exported from $SHARED_LIB)"
    else
        fail "json-symbols (nlohmann/valijson symbols exported from $SHARED_LIB)"
    fi
fi

# ---------- check 3: Catch2 tests with auto-selected tag scope -------

step "Catch2 tests"
if skip_check tests; then
    skip "tests" "--skip=tests"
elif [ ! -x ./build_catch/CatchTestRunner ]; then
    skip "tests" "CatchTestRunner not built"
else
    # Tag scope selection.  Path -> tag mapping mirrors how CI's broad
    # workflow runs the full suite, but skips the expensive `[slow]`
    # tests by default for fast local feedback.  Pass --slow to include.
    TAG_FILTER=""
    if printf '%s\n' "$ALL_PATHS" | grep -qE "^(src/SBTL/|include/CoolProp/sbtl/|src/Backends/SVDSBTL/|src/Region/|src/SVD/|include/CoolProp/region/|include/CoolProp/svd/)"; then
        # SBTL/SVDSBTL surface area touched — run the umbrella tags.
        # [SBTL] catches the adapter-layer tests (serializer round-trip,
        # multi-fluid PH preset) that [SVDSBTL] alone misses.
        TAG_FILTER="[SBTL],[SVDSBTL],[SVDComponents],[region],[!benchmark]"
    elif printf '%s\n' "$ALL_PATHS" | grep -qE "^(src/Backends/Helmholtz/|src/Backends/REFPROP/)"; then
        # HEOS / REFPROP path touched — broader sweep including transport
        # and flash routines.
        TAG_FILTER="[Helmholtz],[REFPROP],[!benchmark]"
    else
        # Default: run everything fast (skip the [slow] long tests).
        TAG_FILTER="[!slow][!benchmark]"
    fi
    echo "  tag filter: $TAG_FILTER"
    if ./build_catch/CatchTestRunner $TAG_FILTER 2>&1 | tee /tmp/preflight-tests.log | tail -3 | grep -qE "failed|Errors:"; then
        fail "tests (see /tmp/preflight-tests.log)"
    else
        ok "tests ($TAG_FILTER)"
    fi
fi

# ---------- check 4: cppcheck ----------------------------------------

step "cppcheck"
if skip_check cppcheck; then
    skip "cppcheck" "--skip=cppcheck"
elif ! command -v cppcheck >/dev/null 2>&1; then
    skip "cppcheck" "cppcheck not on PATH (brew install cppcheck)"
elif [ -z "$ALL_CPP" ]; then
    skip "cppcheck" "no C/C++ files in diff"
else
    # --error-exitcode=1 surfaces any warning/style/error as a hard fail.
    # Same rules CI's informational cppcheck job uses.
    # --language=c++ + --std=c++17 force the C++ parser on .h files
    # (cppcheck otherwise picks C and rejects `namespace`).  Matches
    # the CI cppcheck workflow's invocation.
    #
    # --enable=warning (NOT style/performance/portability): style
    # findings are opinion and the CI cppcheck workflow runs in
    # "informational" mode so they don't block PRs.  Preflight mirrors
    # that — warnings are real-bug-class (uninit vars, null deref,
    # buffer overflow) and worth blocking.
    if ! cppcheck --enable=warning --error-exitcode=1 --quiet --inline-suppr --language=c++ --std=c++17 \
                  --suppress=missingIncludeSystem --suppress=unknownMacro $ALL_CPP 2>/tmp/preflight-cppcheck.log; then
        tail -30 /tmp/preflight-cppcheck.log
        fail "cppcheck (see /tmp/preflight-cppcheck.log)"
    else
        ok "cppcheck ($(printf '%s\n' "$ALL_CPP" | wc -l | tr -d ' ') file(s))"
    fi
fi

# ---------- check 5: clang-tidy diff-only ----------------------------

step "clang-tidy (diff-only, signal-filtered)"
if skip_check clang-tidy; then
    skip "clang-tidy" "--skip=clang-tidy"
elif [ -z "$ALL_CPP" ]; then
    skip "clang-tidy" "no C/C++ files in diff"
elif [ ! -f build_catch/compile_commands.json ]; then
    skip "clang-tidy" "build_catch/compile_commands.json missing (cmake configure with -DCMAKE_EXPORT_COMPILE_COMMANDS=ON)"
else
    # Noise filter: clang-tidy checks that CI explicitly elected NOT to
    # gate on, per issue #2926's "Filtered as noise" section.  These
    # dominate the raw output (Catch2 macro expansions, REFPROP C-API
    # surface, identifier-reserved patterns from numerical-derivative
    # naming) without delivering signal worth blocking on.  Findings
    # matching ANY of these check names are subtracted from the gating
    # count; preflight passes if the remaining (signal) count is zero.
    #
    # Sourced from #2926 — keep in sync if that triage report updates.
    # New noise classes that recur across PRs without yielding action
    # belong here too.
    CLANG_TIDY_NOISE_CHECKS=(
        cppcoreguidelines-avoid-do-while
        cert-err58-cpp
        modernize-avoid-c-arrays
        cppcoreguidelines-init-variables
        cppcoreguidelines-pro-bounds-pointer-arithmetic
        cert-dcl37-c
        cert-dcl51-cpp
        bugprone-reserved-identifier
        cert-msc32-c
        cert-msc51-cpp
        clang-analyzer-optin.core.EnumCastOutOfRange
        # AbstractState::AbstractState() calls the (virtual) clear() to
        # initialize members; the base impl is independent of any
        # overrides.  Refactoring around the warning would mean splitting
        # clear() into virtual + non-virtual halves across the whole
        # backend hierarchy.  Cppcheck classifies the same finding as
        # `style` (not warning), and CI's clang-tidy job runs
        # clang-tidy-diff (changed lines only) so it never reports this.
        # Keeping the suppression scoped to preflight to match.
        clang-analyzer-optin.cplusplus.VirtualCall
    )
    NOISE_PATTERN="$(IFS='|'; echo "${CLANG_TIDY_NOISE_CHECKS[*]}")"

    CPP_ONLY="$(printf '%s\n' "$ALL_CPP" | grep -E '\.(cpp|cc|cxx)$' || true)"
    if [ -z "$CPP_ONLY" ]; then
        skip "clang-tidy" "no .cpp files in diff (headers covered transitively)"
    else
        COOLPROP_BUILD_DIR=build_catch ./dev/ci/run-clang-tidy-staged.sh $CPP_ONLY >/tmp/preflight-clang-tidy.log 2>&1 || true
        if grep -q "^warning:.*skipping" /tmp/preflight-clang-tidy.log; then
            skip "clang-tidy" "$(grep -m1 '^warning:' /tmp/preflight-clang-tidy.log | sed 's/^warning: //')"
        else
            RAW="$(grep -cE 'warning: |error: ' /tmp/preflight-clang-tidy.log 2>/dev/null | head -1 || echo 0)"
            # Each finding line ends with `[<check-name>,-warnings-as-errors]`
            # or `[<check-name>]`.  Match the bracketed check name and
            # exclude any line whose name is in NOISE_PATTERN.
            SIGNAL_LINES="$(grep -E 'warning: |error: ' /tmp/preflight-clang-tidy.log 2>/dev/null \
                | grep -vE "\\[($NOISE_PATTERN)(,|\\])" || true)"
            SIGNAL_COUNT="$(printf '%s\n' "$SIGNAL_LINES" | grep -c . || echo 0)"
            if [ "$SIGNAL_COUNT" -gt 0 ]; then
                printf '\n--- signal findings (noise-filtered, see #2926) ---\n'
                printf '%s\n' "$SIGNAL_LINES" | head -30
                printf '%s\n' "$SIGNAL_LINES" > /tmp/preflight-clang-tidy-signal.log
                fail "clang-tidy ($SIGNAL_COUNT signal / $RAW raw findings; see /tmp/preflight-clang-tidy-signal.log)"
            else
                ok "clang-tidy ($(printf '%s\n' "$CPP_ONLY" | wc -l | tr -d ' ') .cpp file(s); 0 signal / $RAW raw findings)"
            fi
        fi
    fi
fi

# ---------- check 6: semgrep (CodeQL-class catches) ------------------

step "semgrep (cpp + security-audit)"
if skip_check semgrep; then
    skip "semgrep" "--skip=semgrep"
elif [ -z "$ALL_CPP" ]; then
    skip "semgrep" "no C/C++ files in diff"
else
    # uvx-resolved semgrep with the p/security-audit ruleset.  Pin
    # Python 3.12 so semgrep's opentelemetry dep doesn't trip on the
    # missing pkg_resources in Python 3.9.  (p/cpp returns 404 on
    # semgrep registry as of 2026; security-audit catches the major
    # CodeQL-class issues that have slipped through previous PRs.
    # Custom rules for any pattern not in p/security-audit can be
    # added under .semgrep/ and configured here.)
    SEMGREP_CONFIG="--config p/security-audit"
    if [ -d ".semgrep" ]; then
        SEMGREP_CONFIG="$SEMGREP_CONFIG --config .semgrep/"
    fi
    if ! uvx --python 3.12 semgrep $SEMGREP_CONFIG --error --quiet $ALL_CPP 2>/tmp/preflight-semgrep.log; then
        tail -30 /tmp/preflight-semgrep.log
        fail "semgrep (see /tmp/preflight-semgrep.log)"
    else
        ok "semgrep ($(printf '%s\n' "$ALL_CPP" | wc -l | tr -d ' ') file(s))"
    fi
fi

# ---------- summary --------------------------------------------------

echo
echo "──────────────────────────────────────────────────────"
echo "preflight summary: $PASS_COUNT passed / $FAIL_COUNT failed / $SKIP_COUNT skipped"
echo "──────────────────────────────────────────────────────"

if [ $FAIL_COUNT -gt 0 ]; then
    printf '\033[31mFAIL:\033[0m\n'
    for n in "${FAIL_NAMES[@]}"; do printf '  - %s\n' "$n"; done
    exit 1
fi

# ---------- pre-PR code-reviewer reminder ----------------------------
#
# Pre-push shell hooks can't mechanically invoke a Claude Code subagent
# (subagents are an in-conversation construct, not a CLI).  Print a
# loud reminder so the next step before `gh pr create` is clear.  See
# CLAUDE.md "Pre-PR adversarial review" for the canonical invocation.
#
# This banner ALWAYS prints when preflight passes — agents and humans
# both see it.  Skip the banner if --skip=banner is passed (useful for
# successive iteration runs where the reviewer was already run).
if ! skip_check banner; then
    printf '\n\033[1;36m┌─────────────────────────────────────────────────────────┐\033[0m\n'
    printf '\033[1;36m│  REMINDER: before `gh pr create`, run code-reviewer.   │\033[0m\n'
    printf '\033[1;36m│  See CLAUDE.md "Pre-PR adversarial review" for the     │\033[0m\n'
    printf '\033[1;36m│  exact Agent({subagent_type: ...}) invocation.         │\033[0m\n'
    printf '\033[1;36m└─────────────────────────────────────────────────────────┘\033[0m\n'
fi

exit 0
