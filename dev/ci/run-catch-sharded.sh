#!/usr/bin/env bash
# Run a Catch2 test filter split across CPU cores.
#
#   run-catch-sharded.sh <runner> <filter> <expected-cases> <jobs> <logdir> [extra runner args...]
#
# Exit status: 0 = every case ran and passed; 1 = a test failed, a shard
# crashed, or the cases that ran do not add up to <expected-cases>;
# 2 = usage / infrastructure error.
#
# Why sharding: the suite is ~560 cases in one single-threaded process, and
# its cost is concentrated (the ten slowest cases are ~60 % of the time).
# Catch2's --shard-count splits the case list into CONTIGUOUS chunks, and the
# heavy cases sit next to each other in the source files, so N shards on N
# cores only buys ~2x: two chunks inherit all the heavy cases.  Cutting the
# list into ~4x more shards than cores and feeding them through a work queue
# (xargs -P) lets the cheap shards drain around the expensive ones; measured
# on ~[slow], 72 s of case time -> ~21 s simulated makespan on 8 workers.
#
# Fail-closed properties (each one answers "what makes this pass when it
# should fail?"):
#   - Every shard records its OWN exit status to <logdir>/<i>.rc.  A missing
#     .rc (the child never ran, or was killed before writing it) is a
#     failure, not a pass.
#   - A shard with zero cases exits 2 in Catch2 ("No tests ran"), so the
#     shard count is capped at <expected-cases>; every shard is non-empty.
#   - Catch2 exits 4 when every case in a run was SKIPped.  A serial run over
#     the whole filter only does that if the entire selection skipped, but a
#     shard can easily hold nothing except, say, one REFPROP case that skips
#     here.  So exit 4 is accepted per shard, and the serial semantics are
#     restored over the union: if no case in ANY shard ran to a pass or fail
#     (everything skipped), the run fails.
#   - Completeness: each shard also writes Catch2's XML reporter output, and
#     the per-shard case totals (successes+failures+expectedFailures+skips
#     from <OverallResultsCases>) must sum to exactly <expected-cases>.  That
#     catches a shard that exited 0 but ran fewer cases than it was given.
#     The same element's failures= must be 0 in every shard, independently
#     of the exit status (defence in depth; Catch2 already exits non-zero).
#   - <logdir> must be empty, so a stale .rc/.xml from an earlier run can
#     never stand in for a shard that did not run this time.
#     It reads a machine-readable attribute, not the console summary, which
#     a Catch2 upgrade could reword.
#   - The caller keeps --warn UnmatchedTestSpec in <extra args>; Catch2
#     evaluates it against the whole filter before sharding, so a stale tag
#     still exits 3 in every shard (verified on Catch2 3.8.0).

set -euo pipefail

if [ "$#" -lt 5 ]; then
    echo "usage: $0 <runner> <filter> <expected-cases> <jobs> <logdir> [extra runner args...]" >&2
    exit 2
fi

RUNNER="$1"
FILTER="$2"
EXPECTED="$3"
JOBS="$4"
LOGDIR="$5"
shift 5

case "$EXPECTED" in '' | *[!0-9]*)
    echo "run-catch-sharded: expected-cases must be a non-negative integer, got '$EXPECTED'" >&2
    exit 2
    ;;
esac
case "$JOBS" in '' | *[!0-9]* | 0*)
    echo "run-catch-sharded: jobs must be a positive integer without leading zeros, got '$JOBS'" >&2
    exit 2
    ;;
esac
if [ "$EXPECTED" -eq 0 ]; then
    # The caller already treats a zero-case filter as a stale filter; never
    # turn it into a silent pass here.
    echo "run-catch-sharded: expected-cases is 0 -- refusing to run an empty selection" >&2
    exit 2
fi
if [ ! -x "$RUNNER" ]; then
    echo "run-catch-sharded: runner '$RUNNER' is not executable" >&2
    exit 2
fi
if [ ! -d "$LOGDIR" ]; then
    echo "run-catch-sharded: log directory '$LOGDIR' does not exist" >&2
    exit 2
fi
if [ -n "$(ls -A "$LOGDIR")" ]; then
    echo "run-catch-sharded: log directory '$LOGDIR' is not empty" >&2
    exit 2
fi

SHARDS=$((JOBS * 4))
if [ "$SHARDS" -gt "$EXPECTED" ]; then
    SHARDS="$EXPECTED"
fi

# One child per shard.  Arguments reach the child positionally ("$@"), never
# through string interpolation, so filter characters such as ~ [ ] , survive
# unmangled.  The child's own status is always 0 once it has written the .rc
# file, which keeps xargs' exit status reserved for infrastructure failure.
# shellcheck disable=SC2016  # $1..$@ are expanded by the child shell
CHILD='
    i="$1"; runner="$2"; filter="$3"; shards="$4"; logdir="$5"; shift 5
    "$runner" "$filter" "$@" \
        --shard-count "$shards" --shard-index "$i" \
        --reporter "console::out=$logdir/$i.log" \
        --reporter "xml::out=$logdir/$i.xml" \
        >"$logdir/$i.stdout" 2>&1 && rc=0 || rc=$?
    echo "$rc" >"$logdir/$i.rc"
'

XARGS_RC=0
i=0
while [ "$i" -lt "$SHARDS" ]; do
    printf '%s\n' "$i"
    i=$((i + 1))
done | xargs -P "$JOBS" -I{} bash -c "$CHILD" _ {} "$RUNNER" "$FILTER" "$SHARDS" "$LOGDIR" "$@" || XARGS_RC=$?

if [ "$XARGS_RC" -ne 0 ]; then
    echo "run-catch-sharded: xargs exited $XARGS_RC -- shard dispatch itself failed" >&2
    exit 1
fi

FAILED_SHARDS=()
RAN=0
SKIPPED=0
FAILED_CASES=0
i=0
while [ "$i" -lt "$SHARDS" ]; do
    rc_file="$LOGDIR/$i.rc"
    if [ ! -f "$rc_file" ]; then
        FAILED_SHARDS+=("$i (no exit status recorded)")
    else
        rc="$(cat "$rc_file")"
        if [ "$rc" != "0" ] && [ "$rc" != "4" ]; then
            FAILED_SHARDS+=("$i (exit $rc)")
        fi
    fi
    # Run-level totals: the only <OverallResultsCases> element in the file.
    # A missing or truncated XML contributes nothing, so the sum below comes
    # up short and fails the completeness check.
    xml="$LOGDIR/$i.xml"
    if [ -f "$xml" ]; then
        # Prints "<total> <skips> <failures>" or nothing.
        counts="$(awk '
            /<OverallResultsCases/ {
                for (f = 1; f <= NF; f++) {
                    if (match($f, /^(successes|failures|expectedFailures|skips)="[0-9]+"/)) {
                        v = $f; sub(/^[a-zA-Z]+="/, "", v); sub(/".*/, "", v); s += v
                        if ($f ~ /^skips=/) k += v
                        if ($f ~ /^failures=/) x += v
                    }
                }
                found = 1
            }
            END { if (found) print s + 0, k + 0, x + 0 }
        ' "$xml")"
        if [ -n "$counts" ]; then
            read -r c_total c_skips c_failures <<<"$counts"
            RAN=$((RAN + c_total))
            SKIPPED=$((SKIPPED + c_skips))
            FAILED_CASES=$((FAILED_CASES + c_failures))
        fi
    fi
    i=$((i + 1))
done

STATUS=0
if [ "${#FAILED_SHARDS[@]}" -gt 0 ]; then
    STATUS=1
    echo "run-catch-sharded: ${#FAILED_SHARDS[@]} of $SHARDS shard(s) failed:"
    # Show at most three: a stale tag fails every shard identically.
    shown=0
    for s in "${FAILED_SHARDS[@]}"; do
        shown=$((shown + 1))
        if [ "$shown" -gt 3 ]; then
            echo "--- (remaining failed shards' logs are in $LOGDIR) ---"
            break
        fi
        idx="${s%% *}"
        echo "--- shard $s: $LOGDIR/$idx.log ---"
        # Display only; the failure is already recorded in STATUS.
        tail -25 "$LOGDIR/$idx.log" 2>/dev/null || tail -25 "$LOGDIR/$idx.stdout" 2>/dev/null || true
    done
fi
if [ "$RAN" -ne "$EXPECTED" ]; then
    STATUS=1
    echo "run-catch-sharded: $RAN case(s) ran across $SHARDS shards, but the filter lists $EXPECTED -- incomplete run"
fi
if [ "$FAILED_CASES" -ne 0 ]; then
    STATUS=1
    echo "run-catch-sharded: the XML reports $FAILED_CASES failed case(s)"
fi
if [ "$RAN" -gt 0 ] && [ "$SKIPPED" -eq "$RAN" ]; then
    STATUS=1
    echo "run-catch-sharded: all $RAN case(s) were skipped -- nothing was actually tested"
fi
if [ "$STATUS" -eq 0 ]; then
    echo "run-catch-sharded: $RAN case(s) passed across $SHARDS shards on $JOBS worker(s)"
fi
exit "$STATUS"
