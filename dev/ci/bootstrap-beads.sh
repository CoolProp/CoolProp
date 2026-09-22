#!/usr/bin/env sh
#
# bootstrap-beads.sh - the SessionStart hook for beads.
#
# It does NOT install anything.  Installing `bd` and hydrating its database
# takes about fifteen seconds - and minutes, if it has to fall back to building
# from source - and making every session in every ephemeral container wait for
# that, including the many that never open the issue tracker, is not a cost
# worth paying.
#
# Instead this puts dev/ci/bd-shim.sh on PATH as `bd`, which takes
# milliseconds, and the first actual `bd` command triggers the install
# (dev/ci/beads-install.sh).  The cost lands on the session that asked for it.
#
# Three situations, in order:
#
#   1. A real bd is already installed (developer workstation, or a warm
#      container) - run `bd prime` as before and get out of the way.
#   2. No bd - symlink the shim onto PATH and say so.  Nothing is installed.
#   3. BEADS_BOOTSTRAP=1 - also start the install in the background, so the
#      first command usually finds it ready.  Still non-blocking; the shim's
#      lock makes the two safe together.
#
# Set BEADS_BOOTSTRAP=0 to opt out entirely.
#
# Strictly non-fatal: every failure warns and lets the session continue, the
# same posture as the git hooks.

set -u

BD_REPO="$(git rev-parse --show-toplevel 2>/dev/null)" || exit 0
[ -n "$BD_REPO" ] || exit 0
cd "$BD_REPO" || exit 0

SHIM_SOURCE="${BD_REPO}/dev/ci/bd-shim.sh"

# Kill switch.
if [ "${BEADS_BOOTSTRAP:-}" = "0" ]; then
    exit 0
fi

# Is the `bd` on PATH our shim rather than a real binary?  Compare against the
# shim source, since the thing on PATH is a symlink to it.
beads_on_path_is_shim() {
    _p="$(command -v bd 2>/dev/null)" || return 1
    [ -n "$_p" ] || return 1
    _r="$(readlink -f "$_p" 2>/dev/null)" || _r="$_p"
    [ "$_r" = "$SHIM_SOURCE" ]
}

# 1. A real bd is already available: behave exactly as before.
if command -v bd >/dev/null 2>&1 && ! beads_on_path_is_shim; then
    bd prime || echo "beads: 'bd prime' failed (exit $?)." >&2
    exit 0
fi

# 2. Put the shim on PATH.  First writable directory that is already on PATH,
#    preferring the conventional one; a symlink so the shim always runs the
#    checked-out version rather than a stale copy.
if [ ! -x "$SHIM_SOURCE" ]; then
    echo "beads: ${SHIM_SOURCE} is missing or not executable - bd will be unavailable." >&2
    exit 0
fi

beads_install_shim() {
    for _d in /usr/local/bin "${HOME:-/nonexistent}/.local/bin"; do
        case ":${PATH}:" in
            *":${_d}:"*) ;;
            *) continue ;;
        esac
        [ -d "$_d" ] && [ -w "$_d" ] || continue
        if ln -sfn "$SHIM_SOURCE" "${_d}/bd" 2>/dev/null; then
            printf '%s\n' "${_d}/bd"
            return 0
        fi
    done
    return 1
}

if _shim_path="$(beads_install_shim)"; then
    # This goes to stdout on purpose: SessionStart output is added to the
    # session context, so the agent learns bd exists without anyone paying for
    # an install to tell it so.
    echo "beads: 'bd' is available; the first command sets it up (~15 s), then runs normally. Start with 'bd prime'."
else
    echo "beads: no writable directory on PATH - could not install the bd shim." >&2
    exit 0
fi

# 3. Optional background warm-up, so the first real command usually finds bd
#    already installed.  Detached and non-blocking: this hook returns now
#    either way, and beads-install.sh takes the same flock the shim would, so
#    a command issued mid-install waits rather than racing.
if [ "${BEADS_BOOTSTRAP:-}" = "1" ]; then
    if [ -x "${BD_REPO}/dev/ci/beads-install.sh" ]; then
        nohup "${BD_REPO}/dev/ci/beads-install.sh" >/dev/null 2>&1 &
    fi
fi

exit 0
