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

# Find our own real path so we can locate the checkout we belong to, and the
# sibling library.  This has to be done by hand, because the library is what we
# are trying to find.
#
# ALWAYS resolve $0 first; never short-circuit on "is there a beads-lib.sh next
# to the unresolved $0?".  We are normally invoked through a symlink on PATH,
# and that shortcut would source any beads-lib.sh sitting beside the symlink -
# in /usr/local/bin, say - or, when $0 has no slash at all, `./beads-lib.sh`
# from the current working directory.  The library defines the function that
# chooses which binary to exec, so sourcing the wrong one hands over control.
_bd_p="$0"
_bd_n=0
while [ -L "$_bd_p" ] && [ "$_bd_n" -lt 32 ]; do
    _bd_l="$(readlink "$_bd_p")" || break
    case "$_bd_l" in
        /*) _bd_p="$_bd_l" ;;
        *)  _bd_p="$(dirname -- "$_bd_p")/$_bd_l" ;;
    esac
    _bd_n=$((_bd_n + 1))
done
_bd_dir="$(CDPATH='' cd -- "$(dirname -- "$_bd_p")" 2>/dev/null && pwd)" || _bd_dir=""
if [ -z "$_bd_dir" ]; then
    echo "bd: cannot resolve the shim's own location." >&2
    exit 127
fi
_bd_lib="${_bd_dir}/beads-lib.sh"

if [ ! -r "$_bd_lib" ]; then
    echo "bd: cannot find beads-lib.sh next to the shim - is the checkout intact?" >&2
    exit 127
fi
# shellcheck source=dev/ci/beads-lib.sh
. "$_bd_lib"

bd_shim_dir="$_bd_dir"
bd_installer="${bd_shim_dir}/beads-install.sh"
bd_repo="$(CDPATH='' cd -- "${bd_shim_dir}/../.." && pwd)" || bd_repo=""

# Depth counter, as a backstop under the content-based shim detection below.
#
# It must NOT refuse outright on the first re-entry.  A nested `bd` is normal:
# the real binary runs git, git runs .beads/hooks/*, and those run `bd` again.
# Refusing there returns 127, and .beads/hooks/pre-commit and pre-push
# propagate a non-zero hook exit, which aborts the commit or push outright.
#
# So: at depth 1 skip the SETUP only (a nested call must never kick off an
# install), and still resolve and exec the real binary.  Detection at the
# bottom is what actually prevents a loop; this only bounds the damage if
# detection is ever defeated.
BEADS_SHIM_DEPTH=$(( ${BEADS_SHIM_DEPTH:-0} + 1 ))
export BEADS_SHIM_DEPTH
if [ "$BEADS_SHIM_DEPTH" -gt 2 ]; then
    echo "bd: refusing to re-enter the bd shim (depth ${BEADS_SHIM_DEPTH})." >&2
    exit 127
fi

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
    # Nested invocation (see the depth counter above): resolve and run, but
    # never start a setup from inside one.
    [ "${BEADS_SHIM_DEPTH:-1}" -gt 1 ] && return 1
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
