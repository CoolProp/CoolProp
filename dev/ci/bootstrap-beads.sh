#!/usr/bin/env sh
#
# bootstrap-beads.sh - the SessionStart hook for beads.
#
# It installs NOTHING.  Installing `bd` and hydrating its database takes about
# fifteen seconds - and minutes if it has to fall back to building from source
# - and making every session in every ephemeral container wait for that,
# including the many that never open the issue tracker, is not a cost worth
# paying.
#
# Instead this puts dev/ci/bd-shim.sh on PATH as `bd`, which takes
# milliseconds, and the first actual `bd` command triggers the setup
# (dev/ci/beads-install.sh).  The cost lands on the session that asked for it.
#
# Set BEADS_BOOTSTRAP=0 to opt out; BEADS_BOOTSTRAP=1 additionally starts the
# setup in the background so the first command usually finds it ready.
#
# Strictly non-fatal: every failure warns and lets the session continue, the
# same posture as the git hooks.

set -u

# Derive the checkout from our own location, not from the cwd.  A SessionStart
# hook's cwd is not guaranteed to be the project root (or even inside it), and
# `git rev-parse --show-toplevel` from the wrong place silently yields a
# different repo, or nothing, and bd is quietly unavailable for the session.
BD_HERE="$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)" || exit 0
BD_REPO="$(CDPATH='' cd -- "${BD_HERE}/../.." && pwd)" || exit 0
cd "$BD_REPO" || exit 0

[ -r "${BD_HERE}/beads-lib.sh" ] || exit 0
# shellcheck source=dev/ci/beads-lib.sh
. "${BD_HERE}/beads-lib.sh"

SHIM_SOURCE="${BD_HERE}/bd-shim.sh"

# Kill switch.
if [ "${BEADS_BOOTSTRAP:-}" = "0" ]; then
    exit 0
fi

# A real bd already installed - developer workstation, or a warm container.
#
# Prime it only when the database is also ready.  Priming against a cold
# database is not harmless: `bd prime` creates an empty one as a side effect,
# which then looks like a working tracker while every query returns nothing.
# When the binary is there but the database is not, fall through and install
# the shim so the first real command does the hydration.
BD_REAL="$(beads_find_binary)"
if [ -n "$BD_REAL" ] && beads_db_ready "$BD_REPO"; then
    "$BD_REAL" prime || echo "beads: 'bd prime' failed (exit $?)." >&2
    exit 0
fi

if [ ! -x "$SHIM_SOURCE" ]; then
    echo "beads: ${SHIM_SOURCE} is missing or not executable - bd will be unavailable." >&2
    exit 0
fi

# Put the shim on PATH: first writable directory that is already on PATH,
# preferring the conventional one.  A symlink, so the shim always runs the
# checked-out version rather than a stale copy - and so that only one shim file
# exists no matter how many directories point at it.
beads_install_shim() {
    for _d in /usr/local/bin "${HOME:-/nonexistent}/.local/bin"; do
        case ":${PATH}:" in
            *":${_d}:"*) ;;
            *) continue ;;
        esac
        # Create the $HOME candidate if PATH advertises it but it does not
        # exist yet, which is common on a fresh non-root account.
        [ -d "$_d" ] || mkdir -p "$_d" 2>/dev/null || continue
        [ -w "$_d" ] || continue
        if ln -sfn "$SHIM_SOURCE" "${_d}/bd" 2>/dev/null; then
            printf '%s\n' "${_d}/bd"
            return 0
        fi
    done
    return 1
}

if beads_install_shim >/dev/null; then
    # stdout on purpose: SessionStart output joins the session context, so the
    # agent learns bd exists without anything having to be installed to tell it.
    echo "beads: 'bd' is available; the first command sets it up (~15 s), then runs normally. Start with 'bd prime'."
else
    echo "beads: no writable directory on PATH - could not install the bd shim." >&2
    exit 0
fi

# Optional background warm-up, so the first real command usually finds bd
# ready.  Detached and non-blocking: this hook returns now either way, and
# beads-install.sh takes the same lock the shim would, with a bounded wait, so
# a command issued mid-setup waits for it rather than failing.
if [ "${BEADS_BOOTSTRAP:-}" = "1" ]; then
    if [ -x "${BD_HERE}/beads-install.sh" ]; then
        nohup "${BD_HERE}/beads-install.sh" >/dev/null 2>&1 &
    fi
fi

exit 0
