#!/usr/bin/env bash
#
# preflight.sh — local pre-push quality gate for CoolProp.
#
# Runs the same checks CI runs on the diff between the merge-base with
# origin/master (or --base) and the working tree -- committed, staged and
# unstaged changes -- so a passing preflight predicts a green CI.  Lint
# findings count on changed lines only; stage logs go to a per-run
# directory printed at the start.  Closes CoolProp-6r6.
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
#   ./dev/ci/preflight.sh --jobs=4       # cap parallelism (default: all cores)
#
# Tools resolved at runtime:
#   - clang-format     : uvx clang-format@<version-from-.pre-commit-config>
#   - cppcheck         : system binary (graceful skip if missing)
#   - clang-tidy       : delegates to dev/ci/run-clang-tidy-staged.sh
#                         (which already graceful-skips when clang-tidy or
#                          compile_commands.json isn't around)
#   - semgrep          : uvx semgrep with p/cpp + p/security-audit rulesets
#   - Catch2 tests     : ./build_catch/CatchTestRunner with tag scope
#                        auto-selected from the changed paths, sharded
#                        across --jobs cores by dev/ci/run-catch-sharded.sh;
#                        BENCHMARK bodies run once, as in CI
#
# Exit codes: 0 = all checks passed (or were skipped intentionally),
# non-zero = at least one check failed.  When invoked as a pre-push hook
# a non-zero exit will block the push; agents using `--no-verify` to skip
# the hook are doing so deliberately and should run preflight separately.

set -euo pipefail

# ---------- arg parsing ----------------------------------------------

BASE_REF="origin/master"
SKIP_CHECKS=""
JOBS=""
JOBS_GIVEN=0
for arg in "$@"; do
    case "$arg" in
        --base=*) BASE_REF="${arg#*=}" ;;
        # Append rather than assign: a repeated --skip= used to overwrite the
        # earlier one, so `--skip=a --skip=b` silently skipped only b and ran a.
        # Both forms now work -- CSV in one flag, or the flag repeated.
        --skip=*) SKIP_CHECKS="${SKIP_CHECKS:+$SKIP_CHECKS,}${arg#*=}" ;;
        --jobs=*) JOBS="${arg#*=}"; JOBS_GIVEN=1 ;;
        --help|-h)
            # Print the header comment block, stopping at the first
            # non-comment line. A hardcoded end line silently truncated this
            # mid-sentence every time a usage line was added to the header.
            sed -n '2,${/^#/!q;p;}' "$0"
            exit 0
            ;;
        *)
            echo "preflight: unknown arg '$arg' (use --base=<ref>, --skip=<csv> or --jobs=<n>)" >&2
            exit 2
            ;;
    esac
done

skip_check() {
    [[ ",$SKIP_CHECKS," == *",$1,"* ]]
}

# Parallelism for the builds, the sharded test run and clang-tidy.
# The all-cores default applies only when --jobs is absent; an explicit
# `--jobs=` falls through to the validation below and is rejected.
if [ "$JOBS_GIVEN" = 0 ]; then
    JOBS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)"
fi
# Leading zeros are rejected, not stripped: bash arithmetic reads 010 as
# octal 8 and 08 as an error, while xargs -P reads them as decimal.
case "$JOBS" in '' | *[!0-9]* | 0*)
    echo "preflight: --jobs must be a positive integer without leading zeros, got '$JOBS'" >&2
    exit 2
    ;;
esac

# Share CPM dependency downloads across worktrees.  cmake/dependencies.cmake
# already defaults the cache to <worktree>/.cpm_cache, so build dirs WITHIN a
# worktree share it -- but every new worktree still re-clones all of the
# dependencies on its first configure, which was most of the 19 s that
# configure took.  CPM stores the value as a cache variable per build dir, so
# this only affects build dirs configured from here on; existing ones keep
# what they have.  CPM takes a file lock per package, so concurrent worktrees
# can share one cache.  `-` rather than `:-`: an explicitly empty
# CPM_SOURCE_CACHE is an opt-out and is respected.
export CPM_SOURCE_CACHE="${CPM_SOURCE_CACHE-$HOME/.cache/CPM}"

# REFPROP for the [REFPROP] cases.  Every REFPROP-backed case calls
# Skip_if_No_REFPROP(), so without a REFPROP the gate still prints green --
# with ~30 cases silently skipped.  The backend looks at
# $COOLPROP_REFPROP_ROOT, then /opt/refprop; nothing ever set the variable,
# so on a machine with REFPROP installed anywhere else those cases never ran.
# Probe the usual install locations and export the first one that looks like
# a REFPROP root (shared library + fluid files).  The test summary reports
# the skipped count either way.
case "$(uname -s)" in
    Darwin) REFPROP_LIBS="librefprop.dylib" ;;
    Linux) REFPROP_LIBS="librefprop.so" ;;
    *) REFPROP_LIBS="REFPRP64.dll REFPROP.dll" ;;
esac
refprop_root_ok() {
    [ -d "$1" ] || return 1
    [ -d "$1/FLUIDS" ] || [ -d "$1/fluids" ] || return 1
    local lib
    for lib in $REFPROP_LIBS; do
        [ -e "$1/$lib" ] && return 0
    done
    return 1
}
if [ -n "${COOLPROP_REFPROP_ROOT+set}" ] && [ -z "$COOLPROP_REFPROP_ROOT" ]; then
    # Explicitly empty: opt out.  Unset it so the backend does not treat ""
    # as a root, and skip the probe.
    unset COOLPROP_REFPROP_ROOT
    REFPROP_NOTE="REFPROP: disabled (COOLPROP_REFPROP_ROOT set empty); [REFPROP] cases will SKIP"
elif [ -n "${COOLPROP_REFPROP_ROOT:-}" ]; then
    REFPROP_NOTE="REFPROP: \$COOLPROP_REFPROP_ROOT=$COOLPROP_REFPROP_ROOT"
elif refprop_root_ok /opt/refprop; then
    REFPROP_NOTE="REFPROP: /opt/refprop (backend default)"
else
    REFPROP_NOTE="REFPROP: not found; [REFPROP] cases will SKIP (set COOLPROP_REFPROP_ROOT; set it empty to opt out of the probe)"
    for rp in "$HOME/REFPROP10" "$HOME/REFPROP" "$HOME/refprop" /Applications/REFPROP; do
        if refprop_root_ok "$rp"; then
            export COOLPROP_REFPROP_ROOT="$rp"
            REFPROP_NOTE="REFPROP: $rp (auto-detected; set COOLPROP_REFPROP_ROOT to override)"
            break
        fi
    done
fi

# Ninja when available: it is what CI uses, and it schedules the build better
# than Unix Makefiles.  Only applied when preflight configures a build dir
# itself (the dir does not exist yet), so it can never clash with the
# generator an existing CMakeCache.txt records.
CMAKE_GEN_ARGS=()
if command -v ninja >/dev/null 2>&1; then
    CMAKE_GEN_ARGS=(-G Ninja)
fi

# ---------- locate repo + cd to root ---------------------------------

REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

# ---------- per-run log directory ------------------------------------
#
# Every stage log lives in a fresh directory for THIS run.  They used to be
# fixed /tmp/preflight-*.log paths shared by every worktree on the machine,
# so two concurrent runs overwrote each other's logs: one run reported 48
# clang-tidy findings that all belonged to another worktree's run.  That is
# fail-open as well as fail-closed -- a gate that greps a log can read the
# other run's clean log.  Kept after the run so the paths in the messages
# stay valid.
PF_LOGDIR="$(mktemp -d "${TMPDIR:-/tmp}/preflight.XXXXXX")"
echo "preflight logs: $PF_LOGDIR"

# ---------- build-dir validation --------------------------------------
#
# ensure_build_dir <dir> <log> <cmake -D args...>
#
# Returns 0 when <dir> holds a usable configure of THIS checkout (configuring
# it first if needed), non-zero if configuring failed (see <log>).
#
# The old guard was `[ ! -d dir ]`: an interrupted configure left the dir
# behind with a truncated CMakeCache.txt, the guard then skipped configuring
# forever, and every later run failed the build or json-symbols with a
# misleading message.  A dir is usable only if its cache names this project,
# points at THIS source tree (a dir copied from another worktree does not),
# and its generator file exists.  A dir that fails that is one preflight
# owns (build_catch / build_shared, both gitignored), so it is removed and
# configured from scratch.
build_dir_problem() {
    local dir="$1" cache="$1/CMakeCache.txt" home
    [ -f "$cache" ] || { echo "no CMakeCache.txt"; return; }
    grep -qx 'CMAKE_PROJECT_NAME:STATIC=CoolProp' "$cache" || { echo "CMakeCache.txt does not name project CoolProp (truncated?)"; return; }
    home="$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' "$cache")"
    if [ -z "$home" ] || [ "$(cd "$home" 2>/dev/null && pwd -P)" != "$(pwd -P)" ]; then
        echo "CMakeCache.txt belongs to a different source tree (${home:-unknown})"
        return
    fi
    [ -f "$dir/build.ninja" ] || [ -f "$dir/Makefile" ] || { echo "no build.ninja or Makefile (configure did not finish)"; return; }
}
ensure_build_dir() {
    local dir="$1" log="$2" problem
    shift 2
    if [ -d "$dir" ]; then
        problem="$(build_dir_problem "$dir")"
        [ -z "$problem" ] && return 0
        echo "  $dir is not a usable build dir ($problem); reconfiguring from scratch"
        rm -rf "$dir"
    fi
    cmake -B "$dir" -S . ${CMAKE_GEN_ARGS[@]+"${CMAKE_GEN_ARGS[@]}"} "$@" >"$log" 2>&1
}

# ---------- resolve diff set -----------------------------------------

# Everything is diffed against the merge-base with the WORKING TREE, which
# covers committed, staged and unstaged changes in one diff.  The old form
# unioned `BASE...HEAD` with `git diff` (unstaged only), so staged-but-
# uncommitted edits were never checked, and its line numbers would have
# referred to two different versions of a file.
#
# An unresolvable base is a hard error.  The old `git diff ... || true`
# turned a typo in --base into an EMPTY diff, and every check then reported
# "skipped: no C/C++ files in diff" -- preflight passed having checked
# nothing.
if ! MERGE_BASE="$(git merge-base "$BASE_REF" HEAD 2>/dev/null)"; then
    echo "preflight: cannot find a merge-base between '$BASE_REF' and HEAD (git fetch origin, or pass --base=<ref>)" >&2
    exit 2
fi
CPP_GLOBS=('*.cpp' '*.h' '*.hpp' '*.cc' '*.cxx')

# One git diff invocation for everything below, with every output-shaping
# setting pinned so user config cannot change what the parsers read:
# diff.mnemonicPrefix / srcPrefix / dstPrefix / noprefix change the "b/"
# prefix, core.quotePath quotes non-ASCII paths, diff.relative rewrites paths
# against the cwd, and external diff / textconv replace the output entirely.
# --no-renames: a moved file shows as delete + add, so ALL_PATHS carries the
# OLD path too and a file moved out of a rule's directory still triggers it.
git_diff() {
    git -c core.quotePath=false -c diff.noprefix=false -c diff.mnemonicPrefix=false \
        diff --no-color --no-ext-diff --no-textconv --no-renames --no-relative \
        --src-prefix=a/ --dst-prefix=b/ "$@"
}

# C-family files that still exist (deleted files are filtered out: there is
# nothing on disk to check).
ALL_CPP="$(git_diff --name-only --diff-filter=ACM "$MERGE_BASE" -- "${CPP_GLOBS[@]}")"
ALL_CPP="$(printf '%s\n' "$ALL_CPP" | sort -u | grep -v '^$' || true)"
# The same list as an array, for passing to tools: expanding $ALL_CPP
# unquoted split a path containing a space into several nonexistent files,
# and cppcheck then skipped the real one without failing.
ALL_CPP_ARR=()
while IFS= read -r _f; do
    [ -n "$_f" ] && ALL_CPP_ARR+=("$_f")
done <<< "$ALL_CPP"

# All paths changed (any extension, deletions included) — used for tag
# auto-selection.
ALL_PATHS="$(git_diff --name-only "$MERGE_BASE")"
ALL_PATHS="$(printf '%s\n' "$ALL_PATHS" | sort -u | grep -v '^$' || true)"

# Added/modified line ranges of the C-family files, one "path<TAB>first<TAB>last"
# per hunk, with line numbers in the working-tree file -- what cppcheck and
# clang-tidy see.  Pure-deletion hunks (+c,0) add no lines and are dropped.
# CI's clang-tidy job scopes to changed lines the same way (clang-tidy-diff),
# so this also stops preflight failing on findings that predate the branch.
#
# The file name is taken only from the "+++ " line of a file's HEADER (between
# "diff --git" and its first "@@"): inside a hunk, an added source line that
# starts with "++ " also reads "+++ ...".  git appends a TAB to that line when
# the path contains a space, so one trailing TAB is stripped.
CHANGED_RANGES="$(git_diff -U0 --diff-filter=ACM "$MERGE_BASE" -- "${CPP_GLOBS[@]}" | awk '
    /^diff --git / { hdr = 1; f = ""; next }
    hdr && /^\+\+\+ / {
        f = substr($0, 5); sub(/\t$/, "", f)
        f = (substr(f, 1, 2) == "b/") ? substr(f, 3) : ""
        next
    }
    /^@@ / {
        hdr = 0
        if (f != "" && match($0, /\+[0-9]+(,[0-9]+)?/)) {
            n = split(substr($0, RSTART + 1, RLENGTH - 1), p, ",")
            cnt = (n > 1) ? p[2] + 0 : 1
            if (cnt > 0) printf "%s\t%d\t%d\n", f, p[1], p[1] + cnt - 1
        }
    }')"

# Completeness: every C-family file that gained lines must have ranges.  Both
# lint gates DROP findings outside CHANGED_RANGES, so a parse that loses a
# file (a path shape or git setting the parser did not anticipate) would make
# them pass having checked nothing.  --numstat -z is a separate, NUL-safe
# listing of the same diff; if the two disagree, stop.
RANGE_FILES="$(printf '%s\n' "$CHANGED_RANGES" | awk -F'\t' 'NF == 3 { print $1 }' | sort -u)"
MISSING_RANGES=""
# Written to a file first so a failing git aborts here under set -e; inside
# the process substitution its status would be lost and the check would pass
# on an empty listing.
git_diff --numstat -z --diff-filter=ACM "$MERGE_BASE" -- "${CPP_GLOBS[@]}" >"$PF_LOGDIR/numstat.z"
while IFS=$'\t' read -r -d '' ns_added _ns_deleted ns_path; do
    case "$ns_added" in '' | - | 0) continue ;; esac   # binary, or deletions only
    if ! printf '%s\n' "$RANGE_FILES" | grep -qxF -- "$ns_path"; then
        MISSING_RANGES="$MISSING_RANGES $ns_path"
    fi
done <"$PF_LOGDIR/numstat.z"
if [ -n "$MISSING_RANGES" ]; then
    echo "preflight: could not map changed lines for:$MISSING_RANGES" >&2
    echo "           (the cppcheck/clang-tidy line filters would skip these files; refusing to run)" >&2
    exit 2
fi

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
elif ! ensure_build_dir build_catch "$PF_LOGDIR/configure-build_catch.log" \
        -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release; then
    # Release is not optional: CoolProp sets no default CMAKE_BUILD_TYPE, so
    # without it the runner is built unoptimized and every test stage after
    # this is several times slower.
    tail -20 "$PF_LOGDIR/configure-build_catch.log" || true
    fail "build (cmake configure of build_catch failed; see $PF_LOGDIR/configure-build_catch.log)"
else
    # Test cmake's own exit status.  The previous form piped the build log
    # through `tee | tail -5 | grep -qE "error:|FAILED"` under `set -o
    # pipefail`: a failing build made the pipeline non-zero, so `if
    # <pipeline>` was false and the gate reported the build as PASSING.
    # Grepping the last 5 lines was independently unreliable -- a link
    # error, an OOM kill or a cmake usage error need not print "error:".
    if cmake --build build_catch --target CatchTestRunner -j"$JOBS" >"$PF_LOGDIR/build.log" 2>&1; then
        ok "build CatchTestRunner"
    else
        tail -20 "$PF_LOGDIR/build.log"
        fail "build (see $PF_LOGDIR/build.log)"
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
    if ! ensure_build_dir build_shared "$PF_LOGDIR/configure-build_shared.log" \
            -DCOOLPROP_SHARED_LIBRARY=ON -DCMAKE_BUILD_TYPE=Release; then
        fail "json-symbols (shared configure failed; see $PF_LOGDIR/configure-build_shared.log)"
        SHARED_OK=0
    fi
    if [ "$SHARED_OK" = 1 ] && ! cmake --build build_shared -j"$JOBS" >"$PF_LOGDIR/shared-build.log" 2>&1; then
        fail "json-symbols (shared build failed; see $PF_LOGDIR/shared-build.log)"
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
elif mkdir "$PF_LOGDIR/install-headers" \
     && COOLPROP_CI_LOGDIR="$PF_LOGDIR/install-headers" ./dev/ci/check-installed-headers.sh build_shared; then
    ok "install-headers (detail/json.h not shipped; no installed header pulls nlohmann/valijson)"
else
    fail "install-headers (detail/json.h shipped, or a header pulls nlohmann/valijson, or install failed; see $PF_LOGDIR/install-headers)"
fi

# ---------- check 3: Catch2 tests with auto-selected tag scope -------

step "Catch2 tests"
if skip_check tests; then
    skip "tests" "--skip=tests"
elif [ ! -x ./build_catch/CatchTestRunner ]; then
    skip "tests" "CatchTestRunner not built"
else
    # Tag scope selection: ~[slow] ALWAYS, plus the tags of EVERY path rule
    # the diff matches.  (There is no --slow flag; run
    # `./build_catch/CatchTestRunner "[slow]"` directly for the rest.)
    #
    # This used to be a first-match if/elif chain whose branches REPLACED the
    # default.  Because ~[slow] (the else branch) is the broad set and every
    # rule's list was narrower, matching a rule shrank the sweep: a Helmholtz
    # diff ran 69 cases instead of ~560, and a diff touching both SBTL and
    # Helmholtz paths never saw [flash].  A dev/fluids change selected only
    # the SBTL tags, so the [melting] suite -- which reads that data -- was not
    # run, and a red [melting] shipped green through preflight (#3153).  Now
    # the rules can only ADD: since ~[slow] already holds every non-slow case,
    # what a rule contributes is its [slow]-tagged cases (today: SBTL +57,
    # Helmholtz/REFPROP +6; the cubic, expression and melting rules add none,
    # and stay so a future [slow] tag in those suites is picked up).  That is
    # affordable only because the run is sharded across cores.
    #
    # Catch2 filter syntax, since three of these were wrong before:
    #   ~[tag]   EXCLUDES a tag.  `[!slow]` does NOT exclude -- it selects a
    #            literal tag named "!slow", which no test carries, so
    #            `[!slow][!benchmark]` matched 0 test cases and this gate
    #            passed while running NOTHING.
    #   ,        inside one spec is OR, so "~[slow],[SBTL]" is every non-slow
    #            case plus every [SBTL] case (slow or not).
    #   [!benchmark] is a real Catch2 tag, but benchmarks are HIDDEN from the
    #            default set already (`~[!benchmark]` and no filter both list
    #            468).  So appending `,[!benchmark]` to an OR-list ADDED the
    #            benchmarks instead of excluding them.
    # Separate argv specs are AND-ed (intersected), not OR-ed, so an OR-list
    # must be one comma-separated argument.
    EXTRA_TAGS=""
    add_tags() {
        local t
        for t in "$@"; do
            case ",$EXTRA_TAGS," in
                *",$t,"*) ;;
                *) EXTRA_TAGS="${EXTRA_TAGS:+$EXTRA_TAGS,}$t" ;;
            esac
        done
    }
    diff_touches() {
        printf '%s\n' "$ALL_PATHS" | grep -qE "$1"
    }
    if diff_touches "^(src/SBTL/|include/CoolProp/sbtl/|src/Backends/SVDSBTL/|src/Region/|src/SVD/|include/CoolProp/region/|include/CoolProp/svd/|dev/fluids/|dev/mixtures/)"; then
        # SBTL/SVDSBTL surface area.  [SBTL] catches the adapter-layer tests
        # (serializer round-trip, multi-fluid PH preset) that [SVDSBTL] alone
        # misses.
        #
        # dev/fluids and dev/mixtures are in this list because the SVD tables are
        # SAMPLED from that data: change a fluid and the cached table for it is
        # stale, but nothing on the load path notices (the cache filename hashes
        # build options, not fluid data, and the serializer kRevision check only
        # sees format changes).  The tests that compare a table against HEOS are
        # tagged [slow], so without this a fluid-data-only change never runs
        # them.  PR #3352 shipped a new R-32 viscosity and CI caught the 8 %
        # table/HEOS mismatch that preflight had missed.
        add_tags "[SBTL]" "[SVDSBTL]" "[SVDComponents]" "[region]"
    fi
    if diff_touches "^(src/Backends/Helmholtz/|src/Backends/REFPROP/)"; then
        # HEOS / REFPROP.  [flash],[mixture] matter: the PT-flash two-phase
        # residual test and the all_deltaonly agreement test sat red on master
        # until #3323 because a Helmholtz diff never reached them (bd
        # CoolProp-n2qs).  Under the ~[slow] floor they run anyway; the rule's
        # remaining contribution is the 6 [slow] cases in these tags.
        add_tags "[Helmholtz]" "[REFPROP]" "[flash]" "[mixture]"
    fi
    if diff_touches "^src/Backends/Cubics/"; then
        # Cubic backends.  [cubic] is the umbrella (every [cubic_*] test carries
        # it -- Catch2 tags are exact-match, not prefix).  [mixture_derivs2]
        # finite-differences the whole fugacity chain for PengRobinsonBackend
        # and SRKBackend.  [helmholtz] is here because
        # HelmholtzConsistencyFixture builds a ResidualHelmholtzGeneralizedCubic
        # over SRK and PengRobinson and finite-differences all 14 derivatives,
        # including the third- and fourth-order delta terms that feed the
        # critical-point and stability routines -- nothing else reaches them.
        # [GERG] and [json_validation] cover the cubic JSON payload and
        # change_EOS.
        add_tags "[cubic]" "[volume_translation]" "[mixture_derivs2]" "[michelsen]" "[change_EOS]" "[GERG]" "[json_validation]" "[helmholtz]"
    fi
    if diff_touches "^(src/expression/|include/CoolProp/expression/|src/Backends/Helmholtz/|src/Tests/CoolProp-Tests-Expression\.cpp|dev/fluids/)"; then
        # The expression DSL; its parse branch and dispatch arm live under
        # src/Backends/Helmholtz/ and its data under dev/fluids/.
        add_tags "[expression]"
    fi
    if diff_touches "^(dev/fluids/|src/Backends/Helmholtz/|src/Tests/CoolProp-Tests-[A-Za-z]*Melting\.cpp)"; then
        # Melting lines: the fluid JSON carries them, and the Helmholtz backend
        # (MeltingCaloric.cpp, FluidLibrary, Ancillaries, the flash routines)
        # evaluates them.  r1w7.4.
        add_tags "[melting]"
    fi
    TAG_FILTER="~[slow]${EXTRA_TAGS:+,$EXTRA_TAGS}"
    echo "  tag filter: $TAG_FILTER"
    echo "  $REFPROP_NOTE"
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
    test_logdir="$PF_LOGDIR/tests"
    mkdir "$test_logdir"
    if ! listed_tests=$(./build_catch/CatchTestRunner "$TAG_FILTER" \
                            --list-tests --verbosity quiet); then
        fail "tests (could not list cases for filter '$TAG_FILTER' -- is the runner intact?)"
    else
        matched=$(printf '%s\n' "$listed_tests" | awk 'NF { c++ } END { print c + 0 }')
        if [ "$matched" -eq 0 ]; then
            fail "tests (filter '$TAG_FILTER' matched 0 test cases -- filter is stale, not a pass)"
        # Sharded across $JOBS cores (see dev/ci/run-catch-sharded.sh for the
        # load-balancing and fail-closed argument).  The runner's own summary
        # goes to the log dir, not the terminal: streaming ~560 cases here
        # would bury the cppcheck/clang-tidy/semgrep results and the summary
        # below it.  Beyond its exit status, the sharded runner checks that
        # the cases that ran add up to $matched, so a shard that crashed or
        # silently ran nothing fails rather than shrinking the gate.
        #
        # BENCHMARK flags: the same one-sample settings CI uses.  Catch2's
        # defaults take 100 samples plus analysis, which made
        # [superanc],[caching] alone 25.9 s.  --skip-benchmarks would be
        # faster still but is NOT safe: a BENCHMARK body that throws fails its
        # test case, and some bodies are the only place a code path runs
        # ("Performance regression for TS; on/off" [2438] is the only
        # update(SmolarT_INPUTS) call in its case; the #2773 update_with_guesses
        # paths in [caching] run only inside BENCHMARK).  Running each body
        # once costs ~1 s over skipping them once the run is sharded.
        elif mkdir "$test_logdir/shards" \
             && ./dev/ci/run-catch-sharded.sh ./build_catch/CatchTestRunner "$TAG_FILTER" "$matched" "$JOBS" "$test_logdir/shards" \
                 --warn UnmatchedTestSpec \
                 --benchmark-samples 1 --benchmark-no-analysis --benchmark-warmup-time 0 \
                 >"$test_logdir/summary.txt" 2>&1; then
            # The runner's last line carries the skipped count; surface it in
            # the verdict so a green gate cannot hide cases that never ran.
            ok "tests ($TAG_FILTER, $matched cases; $(tail -1 "$test_logdir/summary.txt" | sed -n 's/.*(\([0-9]* skipped\)).*/\1/p'))"
        else
            # `|| true` guards the DISPLAY only: without it a failed cat
            # would abort the script under `set -e` before `fail` records the
            # result.  It cannot mask the gate -- `fail` runs unconditionally.
            cat "$test_logdir/summary.txt" || true
            fail "tests ($TAG_FILTER; per-shard logs: $test_logdir/shards)"
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
    #
    # Scoped to CHANGED LINES, like CI's diff-only lint jobs.  The whole-file
    # form failed on findings that predate the branch, so any diff touching a
    # file with old findings could never pass and pushes routinely needed
    # --no-verify -- which also skips every other check.  cppcheck has no
    # line filter, so it analyses the whole files and its findings are
    # filtered afterwards against CHANGED_RANGES.
    #
    # Fail-closed rules for the filter:
    #   - exit status must be 0 or 1; anything else (crash, bad option) fails.
    #   - exit 1 with an EMPTY findings log fails: cppcheck also exits 1 for
    #     usage errors (printed on stdout, not in this log), which must not
    #     read as "findings, all filtered away".
    #   - analysis-failure ids (syntaxError, internal*, preprocessorError-
    #     Directive, cppcheckError) and findings at line 0 are kept whatever
    #     their line: they mean part of the file was not analysed at all.
    #     (unknownMacro, the other such id, stays suppressed as before.)
    #   - a finding line that does not parse is kept, not dropped.
    printf '%s\n' "$CHANGED_RANGES" >"$PF_LOGDIR/changed-ranges.tsv"
    CPPCHECK_RC=0
    cppcheck --enable=warning --error-exitcode=1 --quiet --inline-suppr --language=c++ --std=c++17 \
             --suppress=missingIncludeSystem --suppress=unknownMacro \
             --template='{file}\t{line}\t{severity}\t{id}\t{message}' \
             "${ALL_CPP_ARR[@]}" 2>"$PF_LOGDIR/cppcheck.log" || CPPCHECK_RC=$?
    if [ "$CPPCHECK_RC" -ne 0 ] && [ "$CPPCHECK_RC" -ne 1 ]; then
        tail -30 "$PF_LOGDIR/cppcheck.log" || true
        fail "cppcheck (exited $CPPCHECK_RC -- crashed or rejected its arguments; see $PF_LOGDIR/cppcheck.log)"
    else
        awk -F'\t' '
            FILENAME == ARGV[1] { if (NF == 3) { n++; rf[n] = $1; rs[n] = $2 + 0; re[n] = $3 + 0 } next }
            NF < 5 || $2 !~ /^[0-9]+$/ { print; next }        # unparseable: keep
            $4 ~ /^(syntaxError|internalAstError|internalError|cppcheckError|preprocessorErrorDirective)$/ || $2 == 0 { print; next }
            { for (i = 1; i <= n; i++) if (rf[i] == $1 && $2 >= rs[i] && $2 <= re[i]) { print; next } }
        ' "$PF_LOGDIR/changed-ranges.tsv" "$PF_LOGDIR/cppcheck.log" >"$PF_LOGDIR/cppcheck-changed-lines.log"
        CPPCHECK_ALL="$(awk 'NF { c++ } END { print c + 0 }' "$PF_LOGDIR/cppcheck.log")"
        CPPCHECK_KEPT="$(awk 'NF { c++ } END { print c + 0 }' "$PF_LOGDIR/cppcheck-changed-lines.log")"
        if [ "$CPPCHECK_RC" -eq 1 ] && [ "$CPPCHECK_ALL" -eq 0 ]; then
            fail "cppcheck (exited 1 but reported no findings -- usage error?; see $PF_LOGDIR/cppcheck.log)"
        elif [ "$CPPCHECK_KEPT" -gt 0 ]; then
            head -30 "$PF_LOGDIR/cppcheck-changed-lines.log" || true
            fail "cppcheck ($CPPCHECK_KEPT finding(s) on changed lines or not attributable to a line; see $PF_LOGDIR/cppcheck-changed-lines.log)"
        else
            ok "cppcheck ($(printf '%s\n' "$ALL_CPP" | wc -l | tr -d ' ') file(s); 0 on changed lines, $CPPCHECK_ALL elsewhere in those files)"
        fi
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
        # Scoped to CHANGED LINES via clang-tidy's -line-filter, the same
        # scoping CI's clang-tidy-diff job uses.  The whole-file form failed on
        # findings that predate the branch (VLERoutines.cpp alone has 31), so
        # touching such a file made the gate unpassable.  One JSON filter lists
        # the changed ranges of every changed C-family file, headers included,
        # so a finding reported in a changed header line is kept too.  Verified:
        # the filter is a path-suffix match, and clang-tidy's exit status counts
        # only the findings that survive it.
        CT_FILTER="$(printf '%s\n' "$CHANGED_RANGES" | awk -F'\t' '
            NF == 3 {
                f = $1; gsub(/\\/, "\\\\", f); gsub(/"/, "\\\"", f)
                if (!(f in seen)) { seen[f] = 1; order[++n] = f }
                r[f] = r[f] (r[f] == "" ? "" : ",") "[" $2 "," $3 "]"
            }
            END {
                printf "["
                for (i = 1; i <= n; i++) printf "%s{\"name\":\"%s\",\"lines\":[%s]}", (i > 1 ? "," : ""), order[i], r[order[i]]
                printf "]"
            }')"
        # Only .cpp files with at least one added line: a file whose hunks are
        # all deletions has nothing for the filter to keep.  Named, not dropped.
        CT_FILES="$(printf '%s\n' "$CHANGED_RANGES" | awk -F'\t' 'NF == 3 && $1 ~ /\.(cpp|cc|cxx)$/ && !s[$1]++ { print $1 }')"
        CT_NO_ADDED="$(printf '%s\n' "$CPP_ONLY" | grep -vxF -f <(printf '%s\n' "$CT_FILES") || true)"
        if [ -n "$CT_NO_ADDED" ]; then
            echo "  no added lines (deletions only), not analysed: $(printf '%s ' $CT_NO_ADDED)"
        fi
    fi
    if [ -n "$CPP_ONLY" ] && [ -z "$CT_FILES" ]; then
        skip "clang-tidy" "changed .cpp files have no added lines"
    elif [ -n "$CPP_ONLY" ]; then
        # One clang-tidy per file, $JOBS at a time (it was one process over
        # every file, ~47 s each, serially).  Each file logs to its own
        # zero-padded file and records its OWN exit status; the logs are joined
        # in order afterwards so parallel writers never interleave lines.
        CT_LOGDIR="$PF_LOGDIR/clang-tidy"
        mkdir "$CT_LOGDIR"
        ct_i=0
        while IFS= read -r ct_f; do
            [ -n "$ct_f" ] || continue
            printf '%04d\0%s\0' "$ct_i" "$ct_f"
            ct_i=$((ct_i + 1))
        done <<< "$CT_FILES" \
            | xargs -0 -n 2 -P "$JOBS" bash -c \
                'COOLPROP_BUILD_DIR=build_catch ./dev/ci/run-clang-tidy-staged.sh "-line-filter=$1" "$3" >"$0/$2.log" 2>&1 && rc=0 || rc=$?; echo "$rc" >"$0/$2.rc"' \
                "$CT_LOGDIR" "$CT_FILTER"
        cat "$CT_LOGDIR"/*.log >"$PF_LOGDIR/clang-tidy.log"
        # Per-file verdicts first.  The verdict used to come from grepping the
        # combined log alone, behind a `|| true` on the run: a clang-tidy that
        # crashed or died before printing anything left no "warning:" line
        # and passed.  Each file's exit status now has to agree with its log:
        # 0 is fine, 1 must come with at least one finding (it is what
        # WarningsAsErrors produces), anything else -- a signal, a usage
        # error -- fails.  The wrapper's graceful skip (no clang-tidy binary,
        # no compile_commands.json) prints a "skipping" warning and exits 0;
        # the stage is reported as skipped only when EVERY file says so, so
        # one file's marker cannot hide another file's crash.
        CT_BAD=()
        CT_SKIPPED=0
        CT_TOTAL=0
        ct_i=0
        while IFS= read -r ct_f; do
            [ -n "$ct_f" ] || continue
            ct_id="$(printf '%04d' "$ct_i")"
            ct_i=$((ct_i + 1))
            CT_TOTAL=$((CT_TOTAL + 1))
            if [ ! -f "$CT_LOGDIR/$ct_id.rc" ]; then
                CT_BAD+=("$ct_f (no exit status recorded)")
                continue
            fi
            ct_rc="$(cat "$CT_LOGDIR/$ct_id.rc")"
            if [ "$ct_rc" = 0 ] && grep -aq '^warning:.*skipping' "$CT_LOGDIR/$ct_id.log"; then
                CT_SKIPPED=$((CT_SKIPPED + 1))
                continue
            fi
            ct_n="$(grep -acE 'warning: |error: ' "$CT_LOGDIR/$ct_id.log" || true)"
            if [ "$ct_rc" = 1 ] && [ "${ct_n:-0}" -eq 0 ]; then
                CT_BAD+=("$ct_f (exit 1 with no findings)")
            elif [ "$ct_rc" != 0 ] && [ "$ct_rc" != 1 ]; then
                CT_BAD+=("$ct_f (exit $ct_rc)")
            fi
        done <<< "$CT_FILES"
        if [ "${#CT_BAD[@]}" -eq 0 ] && [ "$CT_TOTAL" -gt 0 ] && [ "$CT_SKIPPED" -eq "$CT_TOTAL" ]; then
            skip "clang-tidy" "$(grep -a -m1 '^warning:' "$PF_LOGDIR/clang-tidy.log" | sed 's/^warning: //')"
        elif [ "$CT_SKIPPED" -gt 0 ]; then
            [ "${#CT_BAD[@]}" -eq 0 ] || printf '  %s\n' "${CT_BAD[@]}"
            fail "clang-tidy (skipped $CT_SKIPPED of $CT_TOTAL file(s) -- a partial skip is not a pass; logs in $CT_LOGDIR)"
        else
            # -a on every grep: with the binary-content heuristic, grep prints
            # "Binary file ... matches" instead of the lines.  One observed run
            # reported "1 signal / 1434 raw" where the true signal count was
            # 132; had that one line matched the noise list, the gate would
            # have passed with every finding hidden.  The per-file check above
            # backs this up: a file that exited 1 must show a finding line.
            RAW="$(grep -acE 'warning: |error: ' "$PF_LOGDIR/clang-tidy.log" || true)"
            [ -n "$RAW" ] || RAW=0
            # Each finding line ends with `[<check-name>,-warnings-as-errors]`
            # or `[<check-name>]`.  Match the bracketed check name and
            # exclude any line whose name is in NOISE_PATTERN.
            SIGNAL_LINES="$(grep -aE 'warning: |error: ' "$PF_LOGDIR/clang-tidy.log" \
                | grep -avE "\\[($NOISE_PATTERN)(,|\\])" || true)"
            SIGNAL_COUNT="$(printf '%s\n' "$SIGNAL_LINES" | grep -ac . || true)"
            [ -n "$SIGNAL_COUNT" ] || SIGNAL_COUNT=0
            if [ "${#CT_BAD[@]}" -gt 0 ]; then
                printf '  %s\n' "${CT_BAD[@]}"
                fail "clang-tidy (${#CT_BAD[@]} file(s) did not complete normally; logs in $CT_LOGDIR)"
            elif [ "$SIGNAL_COUNT" -gt 0 ]; then
                printf '\n--- signal findings on changed lines (noise-filtered, see #2926) ---\n'
                printf '%s\n' "$SIGNAL_LINES" | head -30
                printf '%s\n' "$SIGNAL_LINES" > "$PF_LOGDIR/clang-tidy-signal.log"
                fail "clang-tidy ($SIGNAL_COUNT signal / $RAW raw findings on changed lines; see $PF_LOGDIR/clang-tidy-signal.log)"
            else
                ok "clang-tidy ($(printf '%s\n' "$CT_FILES" | wc -l | tr -d ' ') .cpp file(s); 0 signal / $RAW raw findings on changed lines)"
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
    if ! uvx --python 3.12 semgrep $SEMGREP_CONFIG --error --quiet "${ALL_CPP_ARR[@]}" 2>"$PF_LOGDIR/semgrep.log"; then
        tail -30 "$PF_LOGDIR/semgrep.log"
        fail "semgrep (see $PF_LOGDIR/semgrep.log)"
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
    SCHEMA_LOG="$PF_LOGDIR/schema-validate.log"
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
    INCOMP_LOG="$PF_LOGDIR/incomp-sanity.log"
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
