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
#   ./dev/ci/preflight.sh --skip=install-headers        # subset
#   ./dev/ci/preflight.sh --skip=incomp-sanity          # subset
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
        # Append rather than assign: a repeated --skip= used to overwrite the
        # earlier one, so `--skip=a --skip=b` silently skipped only b and ran a.
        # Both forms now work -- CSV in one flag, or the flag repeated.
        --skip=*) SKIP_CHECKS="${SKIP_CHECKS:+$SKIP_CHECKS,}${arg#*=}" ;;
        --help|-h)
            # Print the header comment block, stopping at the first
            # non-comment line. A hardcoded end line silently truncated this
            # mid-sentence every time a usage line was added to the header.
            sed -n '2,${/^#/!q;p;}' "$0"
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
    # ~instant.  --dry-run --Werror prints one diagnostic per violation AND
    # exits non-zero.
    #
    # Capture the output and test that, rather than piping into `grep -q`:
    # under `set -o pipefail` (line 35) the pipeline inherits clang-format's
    # non-zero exit, so `if <pipeline>` was false exactly when violations
    # existed and the gate reported PASS.  With no violations clang-format
    # exits 0 but grep finds nothing and exits 1 -- also false, also PASS.
    # Both branches led to PASS, so the gate could never fail.
    # ALL_CPP is a newline-separated scalar; read it into an array so paths
    # containing spaces or glob characters survive as single arguments.
    CF_FILES=()
    while IFS= read -r _f; do
        [ -n "$_f" ] && CF_FILES+=("$_f")
    done <<< "$ALL_CPP"

    # Branch on clang-format's EXIT STATUS, not on whether it printed
    # anything.  uvx writes progress to stderr on a cold cache
    # ("Downloading clang-format (1.3MiB)"), which 2>&1 captures, so a
    # non-empty-output test reports a formatting failure on a clean tree
    # the first time preflight runs on any machine.
    CF_OUT="$(uvx clang-format@"$CF_VER" --dry-run --Werror "${CF_FILES[@]}" 2>&1)" && CF_RC=0 || CF_RC=$?
    if [ "$CF_RC" -ne 0 ]; then
        printf '%s\n' "$CF_OUT" | head -20
        fail "clang-format (run: uvx clang-format@$CF_VER -i <files>)"
    else
        ok "clang-format ($CF_VER, ${#CF_FILES[@]} file(s))"
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
    # Test cmake's own exit status.  The previous form piped the build log
    # through `tee | tail -5 | grep -qE "error:|FAILED"` under `set -o
    # pipefail`: a failing build made the pipeline non-zero, so `if
    # <pipeline>` was false and the gate reported the build as PASSING.
    # Grepping the last 5 lines was independently unreliable -- a link
    # error, an OOM kill or a cmake usage error need not print "error:".
    if cmake --build build_catch --target CatchTestRunner -j8 >/tmp/preflight-build.log 2>&1; then
        ok "build CatchTestRunner"
    else
        tail -20 /tmp/preflight-build.log
        fail "build (see /tmp/preflight-build.log)"
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
# build_shared dir for this gate.  The gate's default pattern is
# `nlohmann|valijson|rapidjson`; RapidJSON has been removed, so its
# symbols must not be exported either.
step "JSON symbol leak (shared library)"
if skip_check json-symbols; then
    skip "json-symbols" "--skip=json-symbols"
elif skip_check build; then
    skip "json-symbols" "--skip=build (shared build needed)"
else
    # A leak gate that silently skips on a broken shared build is fail-open —
    # surface configure/build failures as a hard fail. Intentional skips go
    # through --skip=json-symbols (handled above), not through swallowed errors.
    SHARED_OK=1
    if [ ! -d build_shared ]; then
        if ! cmake -B build_shared -S . -DCOOLPROP_SHARED_LIBRARY=ON -DCMAKE_BUILD_TYPE=Release >/tmp/preflight-shared-build.log 2>&1; then
            fail "json-symbols (shared configure failed; see /tmp/preflight-shared-build.log)"
            SHARED_OK=0
        fi
    fi
    if [ "$SHARED_OK" = 1 ] && ! cmake --build build_shared -j8 >/tmp/preflight-shared-build.log 2>&1; then
        fail "json-symbols (shared build failed; see /tmp/preflight-shared-build.log)"
        SHARED_OK=0
    fi
    SHARED_LIB=""
    if [ "$SHARED_OK" = 1 ]; then
        SHARED_LIB="$(find build_shared \( -name 'libCoolProp.so' -o -name 'libCoolProp.dylib' \) 2>/dev/null | head -1 || true)"
    fi
    if [ "$SHARED_OK" != 1 ]; then
        : # configure/build failure already reported as a fail above
    elif [ -z "$SHARED_LIB" ] || [ ! -f "$SHARED_LIB" ]; then
        fail "json-symbols (shared build succeeded but no libCoolProp.so/.dylib found)"
    elif ./dev/ci/check-json-symbols.sh "$SHARED_LIB"; then
        ok "json-symbols (no nlohmann/valijson exported from $SHARED_LIB)"
    else
        fail "json-symbols (nlohmann/valijson symbols exported from $SHARED_LIB)"
    fi
fi

# ---------- check 2c: installed-header hygiene -----------------------
#
# Companion to the symbol-leak gate on the install side: assert that
# detail/json.h (which pulls nlohmann/json.hpp + valijson) is not shipped
# in the installed headers.  Reuses build_shared (built by the json-symbols
# step above).  Fail-closed lives in the script.
step "installed-header hygiene"
if skip_check install-headers; then
    skip "install-headers" "--skip=install-headers"
elif skip_check build; then
    skip "install-headers" "--skip=build (shared build needed)"
elif [ ! -d build_shared ]; then
    skip "install-headers" "build_shared not available (run without --skip=json-symbols)"
elif ./dev/ci/check-installed-headers.sh build_shared; then
    ok "install-headers (detail/json.h not shipped; no installed header pulls nlohmann/valijson)"
else
    fail "install-headers (detail/json.h shipped, or a header pulls nlohmann/valijson, or install failed; see /tmp/install-headers-check.log)"
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
    # tests by default for fast local feedback.  (There is no --slow flag;
    # run `./build_catch/CatchTestRunner "[slow]"` directly for those.)
    # Catch2 filter syntax, since three of these were wrong before:
    #   ~[tag]   EXCLUDES a tag.  `[!slow]` does NOT exclude -- it selects a
    #            literal tag named "!slow", which no test carries, so
    #            `[!slow][!benchmark]` matched 0 test cases and this gate
    #            passed while running NOTHING.
    #   ,        inside one spec is OR.
    #   [!benchmark] is a real Catch2 tag, but benchmarks are HIDDEN from the
    #            default set already (`~[!benchmark]` and no filter both list
    #            468).  So appending `,[!benchmark]` to an OR-list ADDED the
    #            benchmarks instead of excluding them.  No benchmark term is
    #            needed; --benchmark-samples stays a CI concern.
    # Separate argv specs are AND-ed (intersected), not OR-ed, so an OR-list
    # must be one comma-separated argument.
    TAG_FILTER=""
    if printf '%s\n' "$ALL_PATHS" | grep -qE "^(src/SBTL/|include/CoolProp/sbtl/|src/Backends/SVDSBTL/|src/Region/|src/SVD/|include/CoolProp/region/|include/CoolProp/svd/)"; then
        # SBTL/SVDSBTL surface area touched — run the umbrella tags.
        # [SBTL] catches the adapter-layer tests (serializer round-trip,
        # multi-fluid PH preset) that [SVDSBTL] alone misses.
        TAG_FILTER="[SBTL],[SVDSBTL],[SVDComponents],[region]"
    elif printf '%s\n' "$ALL_PATHS" | grep -qE "^(src/Backends/Helmholtz/|src/Backends/REFPROP/)"; then
        # HEOS / REFPROP path touched — broader sweep including transport
        # and flash routines.
        TAG_FILTER="[Helmholtz],[REFPROP]"
    elif printf '%s\n' "$ALL_PATHS" | grep -qE "^src/Backends/Cubics/"; then
        # Cubic backends touched.  [cubic] is the umbrella (every [cubic_*]
        # test now carries it too — Catch2 tags are exact-match, not prefix).
        # [mixture_derivs2] is the composition-derivative gate: DerivativeFixture
        # is instantiated for PengRobinsonBackend and SRKBackend and finite-
        # differences the whole fugacity chain, so it is where a wrong
        # composition derivative shows up.
        TAG_FILTER="[cubic],[volume_translation],[mixture_derivs2],[michelsen]"
    else
        # Default: run everything fast (skip the [slow] long tests).
        TAG_FILTER="~[slow]"
    fi
    # The expression-DSL surface is ORTHOGONAL to the branches above, so it has to
    # be OR-ed in rather than being another elif.  Its parse branch and dispatch arm
    # live under src/Backends/Helmholtz/, so any change to them selects
    # "[Helmholtz],[REFPROP]" -- which contains ZERO [expression] test cases.  Before
    # this, a branch touching the DSL, its tests and shipped fluid data reported a
    # green preflight having run none of them.
    if printf '%s\n' "$ALL_PATHS" | grep -qE "^(src/expression/|include/CoolProp/expression/|src/Backends/Helmholtz/|src/Tests/CoolProp-Tests-Expression\.cpp|dev/fluids/)"; then
        case "$TAG_FILTER" in
            "~[slow]") : ;;  # already runs everything except [slow]
            *) TAG_FILTER="${TAG_FILTER},[expression]" ;;
        esac
    fi
    echo "  tag filter: $TAG_FILTER"
    # Gate on the runner's EXIT CODE, not on grepping its output.  The old
    # form piped into `grep -qE "failed|Errors:"`, so the `if` saw grep's
    # status and the runner's was discarded -- a zero-match run (exit 2,
    # "No tests ran") contains neither word and was reported as a pass.
    # Also require a non-zero test count, so a filter that stops matching
    # after a rename fails loudly instead of silently testing nothing.  That
    # count alone only catches a TOTALLY stale filter, though: rename one tag
    # out of the comma-separated OR-lists below and the rest still match, so
    # the gate would pass while testing less.  `--warn UnmatchedTestSpec` on
    # the run closes that -- Catch2 then exits 3 (UnmatchedTestSpecExitCode)
    # if any single term matched nothing.  Verified: "[cbor],[NoSuchTag]"
    # exits 3 on a run.  Note it does NOT work on --list-tests (exits 0
    # there), which is why it is on the run and not the listing.
    #
    # Count the cases via `--list-tests --verbosity quiet`, which prints one
    # test name per line and nothing else.  Deliberately NOT parsing the
    # human-readable "N matching test cases" summary: that string is a
    # presentation detail that a Catch2 upgrade can reword, and if it ever
    # stopped matching, the count would silently read 0.  Line counting also
    # lets the listing's own exit status stay meaningful -- a non-zero exit
    # here means the listing itself failed (missing/broken runner), which is
    # distinct from a filter that legitimately matches nothing (exit 0, no
    # lines).  No `2>/dev/null` and no `|| echo 0`: swallowing either the
    # stderr or the status is what lets a gate fail open.
    test_log="$(mktemp "${TMPDIR:-/tmp}/preflight-tests.XXXXXX")"
    if ! listed_tests=$(./build_catch/CatchTestRunner "$TAG_FILTER" \
                            --list-tests --verbosity quiet); then
        fail "tests (could not list cases for filter '$TAG_FILTER' -- is the runner intact?)"
    else
        matched=$(printf '%s\n' "$listed_tests" | awk 'NF { c++ } END { print c + 0 }')
        if [ "$matched" -eq 0 ]; then
            fail "tests (filter '$TAG_FILTER' matched 0 test cases -- filter is stale, not a pass)"
        # tee to the log but NOT to the terminal: streaming all ~410 cases here
        # would bury the cppcheck/clang-tidy/semgrep results and the summary
        # below it.  `>/dev/null` does not cost the exit status -- under
        # pipefail the pipeline still reports the runner's non-zero status, not
        # tee's (verified).  On failure the tail is echoed so there is context
        # without having to open the log.
        elif ./build_catch/CatchTestRunner "$TAG_FILTER" \
                 --warn UnmatchedTestSpec 2>&1 | tee "$test_log" >/dev/null; then
            ok "tests ($TAG_FILTER, $matched cases listed)"
        else
            # `|| true` guards the DISPLAY only: without it a tail failure
            # would abort the script under `set -e` before `fail` records the
            # result.  It cannot mask the gate -- `fail` runs unconditionally.
            tail -15 "$test_log" || true
            fail "tests ($TAG_FILTER; full log: $test_log)"
        fi
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
            # `|| echo 0` appended a second line (grep -c prints 0 then exits
            # 1).  Harmless here, but the same construct on SIGNAL_COUNT below
            # fed a numeric test and errored on every clean run.
            RAW="$(grep -cE 'warning: |error: ' /tmp/preflight-clang-tidy.log 2>/dev/null | head -1 || true)"
            [ -n "$RAW" ] || RAW=0
            # Each finding line ends with `[<check-name>,-warnings-as-errors]`
            # or `[<check-name>]`.  Match the bracketed check name and
            # exclude any line whose name is in NOISE_PATTERN.
            SIGNAL_LINES="$(grep -E 'warning: |error: ' /tmp/preflight-clang-tidy.log 2>/dev/null \
                | grep -vE "\\[($NOISE_PATTERN)(,|\\])" || true)"
            SIGNAL_COUNT="$(printf '%s\n' "$SIGNAL_LINES" | grep -c . || true)"
            [ -n "$SIGNAL_COUNT" ] || SIGNAL_COUNT=0
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

# ---------- check 7: fluid-data schema validation --------------------
#
# Build-time correctness gate for embedded fluid data (RapidJSON->nlohmann
# migration spec, section 5): validate the source JSON data files under dev/
# against their committed JSON schemas before they're compiled into headers
# and embedded.  Scoped to runs where a dev/*.json file changed.  The
# validator is pure Python and runs independently of the C++ build; resolve
# the jsonschema dependency through uvx the same way clang-format/semgrep are
# resolved (with a graceful fallback to a system jsonschema if uvx is absent).
step "fluid-data schema validation"
if skip_check schema-validate; then
    skip "schema-validate" "--skip=schema-validate"
elif ! printf '%s\n' "$ALL_PATHS" | grep -qE '^dev/(pcsaft|cubics|mixtures)/.*\.json$'; then
    skip "schema-validate" "no dev/{pcsaft,cubics,mixtures}/*.json files in diff"
else
    SCHEMA_LOG=/tmp/preflight-schema-validate.log
    SCHEMA_RC=0
    if command -v uvx >/dev/null 2>&1; then
        uvx --from jsonschema python dev/validate_fluid_schemas.py >"$SCHEMA_LOG" 2>&1 || SCHEMA_RC=$?
    elif command -v python3 >/dev/null 2>&1; then
        python3 dev/validate_fluid_schemas.py >"$SCHEMA_LOG" 2>&1 || SCHEMA_RC=$?
    else
        SCHEMA_RC=127
        echo "no uvx or python3 on PATH" >"$SCHEMA_LOG"
    fi
    if [ "$SCHEMA_RC" -eq 0 ]; then
        SCHEMA_N="$(grep -c '^OK' "$SCHEMA_LOG" 2>/dev/null || true)"
        [ -n "$SCHEMA_N" ] || SCHEMA_N=0
        ok "schema-validate ($SCHEMA_N data file(s) validated)"
    else
        tail -30 "$SCHEMA_LOG"
        fail "schema-validate (see $SCHEMA_LOG)"
    fi
fi

# ---------- check 8: incompressible JSON sanity -----------------------
#
# Guards the committed json/*.json against unfitted placeholders, all-zero
# templates, non-numeric or non-finite values, and blocks the C++ loader would
# reject; also the grid-axis ordering contract and the golden-master refit.
# Nothing ran any of it before this check.  Scoped to runs touching the
# incompressible data or its writer.  The pytest path runs the whole directory;
# the fallback below runs only test_json_sanity.py, the one module that needs
# neither numpy nor scipy.
step "incompressible JSON sanity"
if skip_check incomp-sanity; then
    skip "incomp-sanity" "--skip=incomp-sanity"
elif ! printf '%s\n' "$ALL_PATHS" | grep -qE '^dev/incompressible_liquids/'; then
    skip "incomp-sanity" "no dev/incompressible_liquids/ files in diff"
else
    INCOMP_LOG=/tmp/preflight-incomp-sanity.log
    INCOMP_RC=0
    if ! command -v python3 >/dev/null 2>&1; then
        INCOMP_RC=127
        echo "no python3 on PATH" >"$INCOMP_LOG"
    elif python3 -c 'import pytest' >/dev/null 2>&1; then
        # --color=no is load-bearing: with PY_COLORS/FORCE_COLOR set pytest
        # emits ANSI even when redirected, so the count grep below scores 0 and
        # a passing run is reported as "verified nothing".
        python3 -m pytest dev/incompressible_liquids/ -q --color=no >"$INCOMP_LOG" 2>&1 || INCOMP_RC=$?
    else
        # pytest is not required: the checks are plain asserts, so call them
        # directly rather than skip the gate.  Exiting non-zero on an empty or
        # renamed module matters, else it would report a clean pass.
        python3 - >"$INCOMP_LOG" 2>&1 <<'PY' || INCOMP_RC=$?
import importlib.util, inspect, pathlib, sys

path = pathlib.Path("dev/incompressible_liquids/test_json_sanity.py")
spec = importlib.util.spec_from_file_location("test_json_sanity", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
names = sorted(n for n in dir(module) if n.startswith("test_"))
if not names:
    sys.exit("no test_* functions found in {0}".format(path))
for name in names:
    func = getattr(module, name)
    # Calling a generator function only builds a generator; no assert runs.
    # pytest errors on yield-tests, so match that instead of passing green.
    if inspect.isgeneratorfunction(func):
        sys.exit("{0} is a generator function; its asserts would never run".format(name))
    result = func()
    if result is not None:
        sys.exit("{0} returned {1!r}, expected None".format(name, result))
    print("OK", name)
PY
    fi
    if [ "$INCOMP_RC" -eq 0 ]; then
        # grep -c prints 0 and returns 1, so `|| true` (not `|| echo 0`) keeps
        # one line.  The count gates the pass: pytest exits 0 when every test is
        # skipped, and a green "0 check group(s)" would be a fail-open.
        INCOMP_N="$(grep -cE '^(OK|[0-9]+ passed)' "$INCOMP_LOG" 2>/dev/null || true)"
        [ -n "$INCOMP_N" ] || INCOMP_N=0
        if [ "$INCOMP_N" -gt 0 ]; then
            ok "incomp-sanity ($INCOMP_N check group(s))"
        else
            tail -30 "$INCOMP_LOG"
            fail "incomp-sanity (ran but verified nothing; all tests skipped?)"
        fi
    else
        tail -30 "$INCOMP_LOG"
        fail "incomp-sanity (see $INCOMP_LOG)"
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
