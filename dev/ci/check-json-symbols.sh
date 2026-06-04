#!/usr/bin/env bash
# Fail if CoolProp's shared library exports any symbol matching the JSON-library
# pattern. The pattern is the 2nd arg (or $JSON_SYMBOL_PATTERN), defaulting to
# 'nlohmann|valijson'. RapidJSON is intentionally NOT in the default pattern
# during the migration (it is still used with default visibility and its
# symbols are expected to be exported); Phase Final tightens the pattern to
# 'nlohmann|valijson|rapidjson' once RapidJSON is removed.
#
# Visibility attributes in a static .a only take effect once linked into a
# shared object, so this MUST inspect a shared product, never the archive.
set -euo pipefail

LIB="${1:-}"
PATTERN="${2:-${JSON_SYMBOL_PATTERN:-nlohmann|valijson}}"
if [[ -z "${LIB}" || ! -f "${LIB}" ]]; then
    echo "usage: $0 <path-to-shared-library (.so/.dylib)> [symbol-regex]" >&2
    exit 2
fi

# Fail fast if symbol extraction itself fails — masking it (|| true) would
# yield an empty list and a false "OK", silently disabling the leak gate.
case "$(uname -s)" in
    Darwin) NM_OPTS=(-gU) ;;
    *)      NM_OPTS=(-D --defined-only) ;;
esac
if ! SYMS_RAW="$(nm "${NM_OPTS[@]}" "${LIB}")"; then
    echo "FAIL: unable to read symbols from ${LIB} (nm failed)" >&2
    exit 1
fi
# Guard against a vacuous pass: a real CoolProp shared library exports many
# defined symbols. An empty list means we were handed the wrong file, a
# stripped/empty library, or an unexpected nm format — i.e. we cannot actually
# validate, which must be treated as a failure (fail-closed), not "OK".
NSYMS="$(printf '%s\n' "${SYMS_RAW}" | grep -c '[^[:space:]]' || true)"
if [[ "${NSYMS}" -eq 0 ]]; then
    echo "FAIL: ${LIB} exports no defined symbols — wrong file or stripped library? Cannot validate." >&2
    exit 1
fi
# Demangle when c++filt is present; mangled names still contain the literal
# namespace text (e.g. "nlohmann"), so the pattern match works either way.
if command -v c++filt >/dev/null 2>&1; then
    SYMS="$(printf '%s\n' "${SYMS_RAW}" | c++filt)"
else
    SYMS="${SYMS_RAW}"
fi

LEAKS="$(printf '%s\n' "${SYMS}" | grep -E "${PATTERN}" || true)"
if [[ -n "${LEAKS}" ]]; then
    echo "FAIL: shared library exports JSON-library symbols matching /${PATTERN}/:" >&2
    printf '%s\n' "${LEAKS}" | head -50 >&2
    exit 1
fi
echo "OK: no symbols matching /${PATTERN}/ exported from ${LIB}"
