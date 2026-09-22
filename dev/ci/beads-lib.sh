# shellcheck shell=sh
#
# beads-lib.sh - shared helpers for the beads bootstrap, shim and installer.
#
# Sourced, never executed.  It exists so that the shim and the installer cannot
# drift apart: an earlier version kept two copies of the binary resolver "in
# step" by comment alone, and they immediately diverged - the installer's copy
# lost the portable realpath fallback and started returning the shim itself as
# the real binary, which forks until the machine gives up.
#
# Sourcing costs a file read and no fork, so the shim's warm path stays cheap.

# Identity marker.  bd-shim.sh carries this string; nothing else does.  Shim
# detection greps for it rather than comparing paths, because a path comparison
# only recognises THIS shim file: two checkouts, or a copy taken instead of a
# symlink, produce two shim files that each consider the other a real binary
# and exec each other forever.
BEADS_SHIM_MARKER="BEADS-SHIM-IDENTITY-b7f3c1"

# Private npm prefix.  Not `npm install -g`: after an `npm uninstall -g`, the
# next global install of the same package prints "changed 1 package", exits 0
# and installs nothing - not with --force, and not after clearing the empty
# scope directory it leaves behind.  A prefixed install is self-contained,
# repeatable after a wipe, and needs no root.
BEADS_NPM_PREFIX="${XDG_CACHE_HOME:-${HOME:-/tmp}/.cache}/coolprop/beads"
BEADS_NPM_BIN="${BEADS_NPM_PREFIX}/node_modules/.bin"

# Is this candidate our shim rather than a real bd?  Content, not path.
beads_is_shim() {
    [ -f "$1" ] || return 1
    head -n 30 "$1" 2>/dev/null | grep -q "$BEADS_SHIM_MARKER"
}

# Absolute path to a real bd binary, or nothing.
#
# Walks PATH so it finds bd however it was installed - npm, go install, brew -
# plus the private npm prefix and the Go bin directories, which need not be on
# PATH.  `command -v bd` cannot be used: where the shim is the `bd` on PATH,
# that is exactly what it returns.
#
# Requires a regular file: `[ -x ]` alone is true for a DIRECTORY named bd,
# which would be selected and then permanently shadow the real binary, since
# the search stops at the first match.
beads_find_binary() {
    {
        printf '%s\n' "$BEADS_NPM_BIN"
        # The \n in the format string is load-bearing: PATH does not end in a
        # colon, so without it the last PATH entry and the first Go directory
        # run together into one nonsense line and BOTH are lost - which can end
        # a successful install with "no bd binary was found".
        printf '%s\n' "$PATH" | tr ':' '\n'
        [ -n "${GOBIN:-}" ] && printf '%s\n' "$GOBIN"
        [ -n "${GOPATH:-}" ] && printf '%s\n' "${GOPATH}/bin"
        printf '%s\n' "${HOME:-/nonexistent}/go/bin"
    } | while IFS= read -r _bfd; do
        [ -n "$_bfd" ] || continue
        [ -f "${_bfd}/bd" ] && [ -x "${_bfd}/bd" ] || continue
        beads_is_shim "${_bfd}/bd" && continue
        printf '%s\n' "${_bfd}/bd"
        break
    done
}

# Hydration marker, written by the installer only after the database has been
# verified non-empty, and living INSIDE the database directory so that the
# `rm -rf` that tears a bad database down removes it too.  .beads/.gitignore
# already ignores embeddeddolt/, so it never shows up in git status.
#
# Presence of the directory is NOT a usable signal on its own: `bd prime`
# creates an empty database as a side effect, which would then latch "ready"
# forever and leave every `bd ready` / `bd show` reporting nothing, with no
# warning and exit 0.
beads_db_marker() {
    printf '%s\n' "$1/.beads/embeddeddolt/.coolprop-hydrated"
}

beads_db_ready() {
    _bdm="$(beads_db_marker "$1")"
    [ -d "$1/.beads/embeddeddolt" ] && [ -f "$_bdm" ]
}
