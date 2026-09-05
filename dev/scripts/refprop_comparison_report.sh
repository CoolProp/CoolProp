#!/usr/bin/env bash
#
# Regenerate the HEOS-vs-REFPROP flash-consistency report in one step.
#
# Runs the consistency grid on both backends over every fluid the two have in
# common, then renders the per-point CSVs and the standalone HTML report.
# Intended to be runnable by hand and as a nightly CI job.
#
#   dev/scripts/refprop_comparison_report.sh --out /tmp/refprop-report
#   dev/scripts/refprop_comparison_report.sh --out ./out --fluids Water,R134a   # smoke run
#
# REFPROP is located via --refprop-path or $COOLPROP_REFPROP_ROOT.  Without a
# usable REFPROP this exits 3 WITHOUT producing a report -- a nightly job should
# treat that as "not run" rather than as a passing run with nothing in it.
#
# Exit codes (chosen so CI can distinguish "the tool broke" from "the fluids
# have failures", which they always do -- a zero-failure report is not the goal):
#   0  report generated
#  10  a worker crashed (a backend died in native code); report still generated
#   3  REFPROP unavailable, or Python/CoolProp not importable; nothing generated
#   1  anything else went wrong
set -uo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"
OUT=""
FLUIDS=""
JOBS=""
REFPROP_PATH="${COOLPROP_REFPROP_ROOT:-}"
TITLE="HEOS vs REFPROP Flash Audit"

usage() {
    sed -n '2,/^set -uo/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//;$d'
    exit "${1:-0}"
}

# need_value: `shift 2` with only one argument left FAILS WITHOUT SHIFTING, and with no
# `set -e` that leaves the loop spinning on the same argument forever.  A trailing
# `--fluids` (a templated-away variable in a CI job, say) would burn a core until the
# job timed out, silently.  Check before shifting.
need_value() {
    if [[ $# -lt 2 ]]; then
        echo "error: $1 requires a value" >&2
        exit 1
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out) need_value "$@"; OUT="$2"; shift 2 ;;
        --out=*) OUT="${1#*=}"; shift ;;
        --fluids) need_value "$@"; FLUIDS="$2"; shift 2 ;;
        --fluids=*) FLUIDS="${1#*=}"; shift ;;
        --jobs) need_value "$@"; JOBS="$2"; shift 2 ;;
        --jobs=*) JOBS="${1#*=}"; shift ;;
        --refprop-path) need_value "$@"; REFPROP_PATH="$2"; shift 2 ;;
        --refprop-path=*) REFPROP_PATH="${1#*=}"; shift ;;
        --title) need_value "$@"; TITLE="$2"; shift 2 ;;
        --title=*) TITLE="${1#*=}"; shift ;;
        -h|--help) usage 0 ;;
        *) echo "unknown argument: $1" >&2; usage 1 ;;
    esac
done

if [[ -z "$OUT" ]]; then
    echo "error: --out is required" >&2
    usage 1
fi
mkdir -p "$OUT" || { echo "error: cannot create $OUT" >&2; exit 1; }
OUT="$(cd -- "$OUT" && pwd)"

# ---- preflight: fail loudly and distinctly when the inputs are not there -----
# Check the drivers exist BEFORE running them: python exits 2 on "can't open file",
# which would otherwise be indistinguishable from the driver's own "a worker crashed"
# exit 2 and would be reported as a crash rather than a broken checkout.
for _script in consistency_backend_compare.py consistency_backend_report.py; do
    if [[ ! -f "$SCRIPT_DIR/$_script" ]]; then
        echo "error: $SCRIPT_DIR/$_script is missing" >&2
        exit 1
    fi
done
if ! "$PYTHON" -c "import CoolProp" >/dev/null 2>&1; then
    echo "SKIP: CoolProp is not importable by '$PYTHON'." >&2
    exit 3
fi
if ! "$PYTHON" -c "import pandas, matplotlib" >/dev/null 2>&1; then
    echo "SKIP: pandas and matplotlib are required." >&2
    exit 3
fi

# Probe REFPROP through CoolProp itself rather than looking for files: that is what
# the run will actually depend on, and it catches a present-but-unloadable install.
if ! REFPROP_VERSION=$("$PYTHON" - "$REFPROP_PATH" <<'PY' 2>/dev/null
import os, sys
path = sys.argv[1]
import CoolProp.CoolProp as CP
if path:
    CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, os.path.join(path, ''))
CP.PropsSI('D', 'P', 101325, 'T', 300, 'REFPROP::Water')
v = CP.get_global_param_string('REFPROP_version')
if not v or v == 'n/a':
    raise SystemExit(1)
print(v)
PY
); then
    echo "SKIP: REFPROP is not usable (looked in '${REFPROP_PATH:-<default paths>}')." >&2
    echo "      Set COOLPROP_REFPROP_ROOT or pass --refprop-path." >&2
    exit 3
fi

COOLPROP_VERSION=$("$PYTHON" -c "import CoolProp; print(CoolProp.__version__)" 2>/dev/null)
COOLPROP_VERSION="${COOLPROP_VERSION:-unknown}"
GITREV=$(git -C "$SCRIPT_DIR" rev-parse --short HEAD 2>/dev/null || echo "unknown")

echo "CoolProp $COOLPROP_VERSION @ $GITREV vs REFPROP $REFPROP_VERSION"
echo "Output: $OUT"

# ---- 1. run the grid on both backends ---------------------------------------
COMPARE_ARGS=(--out "$OUT" --backends HEOS,REFPROP)
[[ -n "$FLUIDS" ]] && COMPARE_ARGS+=(--fluids "$FLUIDS")
[[ -n "$JOBS" ]] && COMPARE_ARGS+=(--jobs "$JOBS")
[[ -n "$REFPROP_PATH" ]] && COMPARE_ARGS+=(--refprop-path "$REFPROP_PATH")

# Remove the previous run's aggregate first.  It is the only evidence the shell has
# that the driver got as far as writing results, and nothing else deletes it -- so a
# driver that died on a typo'd --fluids would leave YESTERDAY's json in place, satisfy
# the check below, and re-render a stale report under today's date, exit 0.
rm -f "$OUT/backend_compare.json"

"$PYTHON" "$SCRIPT_DIR/consistency_backend_compare.py" "${COMPARE_ARGS[@]}"
COMPARE_RC=$?
# Driver contract: 0 = clean, 1 = some (fluid, backend) produced no result,
# 10 = a worker crashed.  10 rather than 2 ON PURPOSE: argparse exits 2 for a bad
# option value (--jobs abc), and the driver is invoked with values forwarded from
# this script, so 2 would report tool misuse as a native-code crash.
# 0/1/10 must not stop the report -- the report is where a reader finds out WHICH
# fluid failed, so suppressing it on failure hides the finding.  Anything else is
# the tooling breaking.
if [[ $COMPARE_RC -ne 0 && $COMPARE_RC -ne 1 && $COMPARE_RC -ne 10 ]]; then
    echo "error: the comparison driver failed (exit $COMPARE_RC)" >&2
    exit 1
fi
if [[ ! -f "$OUT/backend_compare.json" ]]; then
    echo "error: the comparison produced no backend_compare.json" >&2
    exit 1
fi

# ---- 2. render the report ----------------------------------------------------
"$PYTHON" "$SCRIPT_DIR/consistency_backend_report.py" \
    --json "$OUT/backend_compare.json" \
    --out "$OUT/report.html" \
    --title "$TITLE" \
    --coolprop "$COOLPROP_VERSION @ $GITREV" \
    --refprop "$REFPROP_VERSION" || {
        echo "error: report generation failed" >&2
        exit 1
    }

echo
echo "Report:      $OUT/report.html"
for _f in backend_compare_pairs.csv backend_compare_points.csv; do
    [[ -f "$OUT/$_f" ]] && echo "             $OUT/$_f"
done

if [[ $COMPARE_RC -eq 10 ]]; then
    echo
    echo "A worker CRASHED -- see the banner at the top of the report." >&2
    exit 10
fi
exit 0
