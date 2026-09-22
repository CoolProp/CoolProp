#!/usr/bin/env sh
#
# BEADS-SHIM-IDENTITY-b7f3c1
#
# bd-shim.sh - a `bd` that sets itself up the first time you actually use it.
#
# dev/ci/bootstrap-beads.sh symlinks this onto PATH as `bd` in ephemeral
# containers, where the real binary is absent.  The first `bd` command pays for
# the install and the database hydration (dev/ci/beads-install.sh, about
# fifteen seconds); every later command adds roughly 60 ms of resolution before
# exec'ing the real binary, which itself takes 120-200 ms.
#
# This is why the SessionStart hook is fast: a session that never touches the
# issue tracker installs nothing, and one that does pays at the moment it asked
# rather than making everybody wait up front.
#
# The marker line at the top of this file is load-bearing - see beads-lib.sh.

set -u

_bd_self_raw="$0"

# Find our own real path (we are invoked through a symlink) so we can locate
# the checkout we belong to, and the sibling scripts.
_bd_here_seed="$(dirname -- "$_bd_self_raw")"
if [ -r "${_bd_here_seed}/beads-lib.sh" ]; then
    _bd_lib="${_bd_here_seed}/beads-lib.sh"
else
    # Invoked through a symlink from elsewhere on PATH: resolve it by hand,
    # without the library (which is what we are trying to find).
    _bd_p="$_bd_self_raw"
    _bd_n=0
    while [ -L "$_bd_p" ] && [ "$_bd_n" -lt 32 ]; do
        _bd_l="$(readlink "$_bd_p")" || break
        case "$_bd_l" in
            /*) _bd_p="$_bd_l" ;;
            *)  _bd_p="$(dirname -- "$_bd_p")/$_bd_l" ;;
        esac
        _bd_n=$((_bd_n + 1))
    done
    _bd_lib="$(CDPATH='' cd -- "$(dirname -- "$_bd_p")" 2>/dev/null && pwd)/beads-lib.sh"
fi

if [ ! -r "$_bd_lib" ]; then
    echo "bd: cannot find beads-lib.sh next to the shim - is the checkout intact?" >&2
    exit 127
fi
# shellcheck source=dev/ci/beads-lib.sh
. "$_bd_lib"

bd_shim_dir="$(CDPATH='' cd -- "$(dirname -- "$_bd_lib")" && pwd)" || {
    echo "bd: cannot locate the CoolProp checkout this shim belongs to." >&2
    exit 127
}
bd_installer="${bd_shim_dir}/beads-install.sh"
bd_repo="$(CDPATH='' cd -- "${bd_shim_dir}/../.." && pwd)" || bd_repo=""

# Hard depth guard.  Shim detection below is content-based and should make an
# exec loop impossible, but a loop costs the whole machine, so refuse outright
# to be entered twice rather than relying on detection alone.
if [ -n "${BEADS_SHIM_DEPTH:-}" ]; then
    echo "bd: refusing to re-enter the bd shim (depth ${BEADS_SHIM_DEPTH})." >&2
    exit 127
fi
BEADS_SHIM_DEPTH=1
export BEADS_SHIM_DEPTH

bd_real="$(beads_find_binary)"

# Should a missing bd be set up right now, or should we quietly stand down?
#
# Not every `bd` on the system is a deliberate request for the issue tracker.
# All five hooks in .beads/hooks/ guard on `command -v bd`, which this shim
# satisfies as soon as it is on PATH, so without this an ordinary `git commit`
# - or `git checkout`, or `git push` - would stop to set up a tracker nobody
# asked for.
#
# BD_GIT_HOOK is the right signal and the hooks already export it (see
# .beads/hooks/*, which all set it before invoking bd).  Do NOT key on GIT_DIR:
# measured on git 2.43, GIT_DIR is not exported to hooks at all, and
# GIT_INDEX_FILE is set only for the index-touching ones (pre-commit,
# prepare-commit-msg) - post-checkout, post-merge and pre-push get neither.  It
# stays as a fallback for hooks that are not beads'.
#
# BEADS_SHIM_NO_INSTALL=1 says the same thing explicitly, for the PreCompact
# hook and for the installer's own restore step.
#
# Standing down is quiet and succeeds: the caller asked for a best-effort sync,
# and "bd is not set up here" is a normal answer, exactly as it was before this
# shim existed.
bd_shim_may_setup() {
    [ "${BEADS_SHIM_NO_INSTALL:-}" = "1" ] && return 1
    [ "${BD_GIT_HOOK:-}" = "1" ] && return 1
    [ -n "${GIT_INDEX_FILE:-}" ] && return 1
    return 0
}

# The binary and the database go missing independently.  A container can have
# bd installed - from a previous session, or an image that ships it - and still
# have no database, because .beads/embeddeddolt/ is gitignored and has to be
# rebuilt from the committed JSONL.
bd_need_setup=0
[ -z "$bd_real" ] && bd_need_setup=1
if [ -n "$bd_repo" ] && [ -f "${bd_repo}/.beads/issues.jsonl" ] &&
   ! beads_db_ready "$bd_repo"; then
    bd_need_setup=1
fi

if [ "$bd_need_setup" = 1 ]; then
    if ! bd_shim_may_setup; then
        # Standing down.  With a binary in hand, run it and let bd report its
        # own state; with none, succeed quietly.
        [ -n "$bd_real" ] || exit 0
    else
        if [ ! -x "$bd_installer" ]; then
            echo "bd: ${bd_installer} is missing - cannot set bd up." >&2
            exit 127
        fi
        "$bd_installer" || exit 127
        bd_real="$(beads_find_binary)"
        if [ -z "$bd_real" ]; then
            echo "bd: setup reported success but no bd binary was found." >&2
            exit 127
        fi
    fi
fi

# Last guard before handing over: never exec another shim, however it got here.
if [ -z "$bd_real" ] || beads_is_shim "$bd_real"; then
    echo "bd: refusing to exec '${bd_real:-<none>}' - it is a bd shim, not the real binary." >&2
    exit 127
fi

exec "$bd_real" "$@"
