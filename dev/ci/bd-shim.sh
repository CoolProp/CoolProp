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

# Exit status for a HARD setup failure: something is wrong and bd cannot run.
#
# 3, not 127.  All five hooks in .beads/hooks/ special-case exactly two
# statuses - 3 ("database not initialized") and 124 (a timeout) - and propagate
# everything else, and a non-zero status out of pre-commit or pre-push ABORTS
# the commit or the push.  A broken bd setup must never do that.  3 says "bd is
# not usable here", which is true, and which those hooks already know to
# ignore; a human still sees a non-zero status and the reason on stderr.
#
# There is one deliberate exception, further down: when the shim is standing
# down on purpose (a git hook, or BEADS_SHIM_NO_INSTALL) and no binary is
# installed, it exits 0.  Nothing is wrong there - the caller asked for a
# best-effort sync and "bd is not set up here" is the normal answer, exactly as
# it was before this shim existed.
BD_SHIM_UNAVAILABLE=3

# Find our own real path so we can locate the checkout we belong to, and the
# sibling library.  This has to be done by hand, because the library is what we
# are trying to find.
#
# ALWAYS resolve $0 first; never short-circuit on "is there a beads-lib.sh next
# to the unresolved $0?".  We are normally invoked through a symlink on PATH,
# and that shortcut would source any beads-lib.sh sitting beside the symlink -
# in /usr/local/bin, say.  The library defines the function that chooses which
# binary to exec, so sourcing the wrong one hands over control.
#
# Every step below therefore stands down rather than guessing: a path we could
# not resolve is not a path we may source a library from.
_bd_p="$0"

# No slash in $0 means there is nothing to resolve and `dirname` would answer
# ".", i.e. ./beads-lib.sh from whatever directory the caller happened to be
# standing in.  A PATH exec of a #! script always hands the interpreter a path
# with a slash, so this is only reachable when something invoked us oddly.
case "$_bd_p" in
    */*) ;;
    *)
        echo "bd: invoked as '$0' with no path - cannot locate the beads checkout." >&2
        exit "$BD_SHIM_UNAVAILABLE"
        ;;
esac

# Chase the symlink chain.  The bound sits ABOVE Linux's SYMLOOP_MAX of 40 on
# purpose: a bound below it leaves chains the kernel will happily execute but
# that we stop following half way, and the half-resolved path is still a
# symlink - whose DIRECTORY we would then source the library from, which is
# exactly the hole this block exists to close.  The check after the loop is
# what makes that safe; the bound only stops a cycle spinning forever.
_bd_n=0
while [ -L "$_bd_p" ]; do
    if [ "$_bd_n" -ge 64 ]; then
        echo "bd: the symlink chain for '$0' is too long to resolve - standing down." >&2
        exit "$BD_SHIM_UNAVAILABLE"
    fi
    if ! _bd_l="$(readlink "$_bd_p")"; then
        echo "bd: cannot read the symlink '$_bd_p' - standing down." >&2
        exit "$BD_SHIM_UNAVAILABLE"
    fi
    case "$_bd_l" in
        /*) _bd_p="$_bd_l" ;;
        *)  _bd_p="$(dirname -- "$_bd_p")/$_bd_l" ;;
    esac
    _bd_n=$((_bd_n + 1))
done

_bd_dir="$(CDPATH='' cd -- "$(dirname -- "$_bd_p")" 2>/dev/null && pwd)" || _bd_dir=""
if [ -z "$_bd_dir" ]; then
    echo "bd: cannot resolve the shim's own location." >&2
    exit "$BD_SHIM_UNAVAILABLE"
fi
_bd_lib="${_bd_dir}/beads-lib.sh"

if [ ! -r "$_bd_lib" ]; then
    echo "bd: cannot find beads-lib.sh next to the shim - is the checkout intact?" >&2
    exit "$BD_SHIM_UNAVAILABLE"
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
#
# Depth 1 is the ordinary, outermost invocation and may do anything, a setup
# included.  Depth 2 is that expected nesting: resolve and exec, but never
# start a setup from inside it.  Deeper is unexpected, so we stand down - and
# with BD_SHIM_UNAVAILABLE rather than 127, so a hook wrapped around us does
# not abort the git operation.  Detection at the bottom of this file is what
# actually prevents a loop; this only bounds the damage if it is ever defeated.
#
# Sanitize before the arithmetic.  `${x:-0}` covers unset and empty and nothing
# else: a value like "1x" is an arithmetic SYNTAX ERROR, and dash kills the
# script with status 2 on the spot - before the git-hook stand-down below has
# had any chance to run.  Anything that is not a plain number is not ours.
case "${BEADS_SHIM_DEPTH:-}" in
    ''|*[!0-9]*) BEADS_SHIM_DEPTH=0 ;;
esac
BEADS_SHIM_DEPTH=$((BEADS_SHIM_DEPTH + 1))
export BEADS_SHIM_DEPTH
if [ "$BEADS_SHIM_DEPTH" -gt 2 ]; then
    echo "bd: refusing to re-enter the bd shim (depth ${BEADS_SHIM_DEPTH})." >&2
    exit "$BD_SHIM_UNAVAILABLE"
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
    # Depth 2 or more, i.e. a nested invocation (see the depth counter above):
    # resolve and run, but never start a setup from inside one.
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
        # own state; with none, succeed quietly - status 0, not
        # BD_SHIM_UNAVAILABLE, because nothing has gone wrong here (see the
        # note on BD_SHIM_UNAVAILABLE at the top of this file).
        [ -n "$bd_real" ] || exit 0
    else
        if [ ! -x "$bd_installer" ]; then
            echo "bd: ${bd_installer} is missing - cannot set bd up." >&2
            exit "$BD_SHIM_UNAVAILABLE"
        fi
        "$bd_installer" || exit "$BD_SHIM_UNAVAILABLE"
        bd_real="$(beads_find_binary)"
        if [ -z "$bd_real" ]; then
            echo "bd: setup reported success but no bd binary was found." >&2
            exit "$BD_SHIM_UNAVAILABLE"
        fi
    fi
fi

# Last guard before handing over: never exec another shim, however it got here.
if [ -z "$bd_real" ] || beads_is_shim "$bd_real"; then
    echo "bd: refusing to exec '${bd_real:-<none>}' - it is a bd shim, not the real binary." >&2
    exit "$BD_SHIM_UNAVAILABLE"
fi

exec "$bd_real" "$@"
