#!/usr/bin/env sh
#
# bd-shim.sh - a `bd` that installs itself the first time you actually use it.
#
# dev/ci/bootstrap-beads.sh symlinks this onto PATH as `bd` in ephemeral
# containers, where the real binary is absent.  The first `bd` command then
# pays for the install and the database hydration (dev/ci/beads-install.sh,
# well under a minute via npm); every later command adds about 60 ms of PATH
# resolution before exec'ing the real binary, which itself takes 120-200 ms.
#
# This is why the SessionStart hook is fast: a session that never touches the
# issue tracker never installs anything, and one that does pays the cost at the
# moment it asked for something, rather than everybody paying it up front.
#
# It finds the real binary by walking PATH itself, skipping any entry that
# resolves back to this file, rather than asking `command -v bd` — which on a
# container where this shim IS the `bd` on PATH returns the shim, and exec'ing
# that would loop forever.  A second guard before the exec catches any path
# that still resolves to us.

set -u

# Resolve our own real path through the symlink, so we can find the repo we
# were installed from.  `readlink -f` is coreutils; fall back to walking the
# links by hand where it is missing (BSD/macOS).
bd_shim_realpath() {
    _p="$1"
    if command -v readlink >/dev/null 2>&1 && readlink -f "$_p" 2>/dev/null; then
        return 0
    fi
    # Manual walk, bounded so a symlink loop cannot hang the shim.
    _n=0
    while [ -L "$_p" ] && [ "$_n" -lt 32 ]; do
        _l="$(readlink "$_p")" || break
        case "$_l" in
            /*) _p="$_l" ;;
            *)  _p="$(dirname -- "$_p")/$_l" ;;
        esac
        _n=$((_n + 1))
    done
    printf '%s\n' "$(CDPATH='' cd -- "$(dirname -- "$_p")" && pwd)/$(basename -- "$_p")"
}

bd_shim_self="$(bd_shim_realpath "$0")"
bd_shim_dir="$(CDPATH='' cd -- "$(dirname -- "$bd_shim_self")" && pwd)" || {
    echo "bd: cannot locate the CoolProp checkout this shim belongs to." >&2
    exit 127
}
bd_installer="${bd_shim_dir}/beads-install.sh"

# Absolute path to a real bd binary, or empty.
#
# Walks PATH so it finds bd however it was installed — npm, go install, brew —
# and skips ourselves.  `command -v bd` cannot be used: where this shim is the
# `bd` on PATH, that is precisely what it returns, and exec'ing it would loop.
#
# Intentionally duplicated from beads_binary() in dev/ci/beads-install.sh, so
# that the warm path here never has to source or spawn a helper; the two copies
# must stay in step.  GOPATH/bin is appended because `go install` puts bd there
# whether or not it is on PATH.
# Private npm prefix that beads-install.sh installs into.  Must stay in step
# with BD_NPM_PREFIX there; see the rationale in that file for why the global
# npm tree is avoided.
bd_shim_npm_bin="${XDG_CACHE_HOME:-${HOME:-/tmp}/.cache}/coolprop/beads/node_modules/.bin"

bd_shim_binary() {
    {
        printf '%s\n' "$bd_shim_npm_bin"
        printf '%s' "$PATH" | tr ':' '\n'
        [ -n "${GOBIN:-}" ] && printf '%s\n' "$GOBIN"
        [ -n "${GOPATH:-}" ] && printf '%s\n' "${GOPATH}/bin"
        printf '%s\n' "${HOME:-/nonexistent}/go/bin"
    } | while IFS= read -r _d; do
        [ -n "$_d" ] || continue
        [ -x "${_d}/bd" ] || continue
        _r="$(bd_shim_realpath "${_d}/bd")"
        [ "$_r" = "$bd_shim_self" ] && continue
        printf '%s\n' "${_d}/bd"
        break
    done
}

bd_real="$(bd_shim_binary)"

# Should a missing bd be installed right now, or should we quietly stand down?
#
# Not every `bd` on the system is a deliberate request for the issue tracker.
# All five beads git hooks (.beads/hooks/) guard on `command -v bd`, and that
# guard passes as soon as this shim is on PATH - so without this, an ordinary
# `git commit` would stop to install a tracker the committer never asked for.
# git exports GIT_DIR / GIT_INDEX_FILE when it runs a hook, which is a reliable
# signal for that case.
#
# BEADS_SHIM_NO_INSTALL=1 says the same thing explicitly, for callers like the
# PreCompact hook that want to prime an existing database but must never
# trigger an install.
#
# Standing down is quiet and succeeds: the caller asked for a best-effort sync,
# and "bd is not installed here" is a normal answer to that, the same as before
# this shim existed.
bd_shim_may_install() {
    [ "${BEADS_SHIM_NO_INSTALL:-}" = "1" ] && return 1
    [ -n "${GIT_INDEX_FILE:-}" ] && return 1
    [ -n "${GIT_DIR:-}" ] && return 1
    return 0
}

# The binary and the database go missing independently.  A container can have
# bd installed (from a previous session, or an image that ships it) and still
# have no .beads/embeddeddolt/, because that directory is gitignored and has to
# be rebuilt from the committed JSONL.  Check both, or a warm binary against a
# cold database just fails.
bd_repo="$(CDPATH='' cd -- "${bd_shim_dir}/../.." && pwd)" || bd_repo=""
bd_need_setup=0
[ -z "$bd_real" ] && bd_need_setup=1
if [ -n "$bd_repo" ] &&
   [ -f "${bd_repo}/.beads/issues.jsonl" ] &&
   [ ! -d "${bd_repo}/.beads/embeddeddolt" ]; then
    bd_need_setup=1
fi

if [ "$bd_need_setup" = 1 ]; then
    if ! bd_shim_may_install; then
        # Standing down.  With a binary in hand, still run it and let bd report
        # its own state; with none, succeed quietly, exactly as a system with
        # no bd installed would have behaved before this shim existed.
        [ -n "$bd_real" ] || exit 0
    else
        if [ ! -x "$bd_installer" ]; then
            echo "bd: ${bd_installer} is missing - cannot set bd up." >&2
            exit 127
        fi
        "$bd_installer" || exit 127
        bd_real="$(bd_shim_binary)"
        if [ -z "$bd_real" ]; then
            echo "bd: setup reported success but no bd binary was found." >&2
            exit 127
        fi
    fi
fi

# Guard against ever exec'ing ourselves, however the paths were resolved.
if [ "$(bd_shim_realpath "$bd_real")" = "$bd_shim_self" ]; then
    echo "bd: refusing to exec the shim as itself (resolved to ${bd_real})." >&2
    exit 127
fi

exec "$bd_real" "$@"
