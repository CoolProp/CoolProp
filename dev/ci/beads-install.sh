#!/usr/bin/env sh
#
# beads-install.sh - install `bd` and hydrate its database, on demand.
#
# This is the slow part of getting beads working in an ephemeral container:
#
#   * the `bd` binary is absent, and
#   * the embedded database directory .beads/embeddeddolt/ is gitignored and
#     therefore absent too, so it has to be rehydrated from the committed
#     .beads/issues.jsonl, which is the source of truth.
#
# Together those take about fifteen seconds, which is why nothing calls this at
# session start.  dev/ci/bd-shim.sh calls it on the first actual `bd` command,
# so a session that never touches the issue tracker never pays for it.
#
# Idempotent and non-fatal.  A warm container returns almost immediately.  Any
# failure warns on stderr and returns non-zero; the caller decides, and the
# session carries on without bd.
#
# All diagnostics go to stderr; stdout stays clean because the shim may be in a
# pipeline.

set -u

_bd_here="$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)" || exit 1
# shellcheck source=dev/ci/beads-lib.sh
. "${_bd_here}/beads-lib.sh"

BD_REPO="$(CDPATH='' cd -- "${_bd_here}/../.." && pwd)" || exit 1
cd "$BD_REPO" || exit 1

# Pinned.  The two differ because they are different distributions of bd: the
# npm package is versioned 1.3.0 and ships a prebuilt binary, while v1.1.0 is
# the last version verified to build from the `github.com/steveyegge/beads`
# module path (the project has since moved to `gastownhall/beads`).
BD_NPM_VERSION="1.3.0"
BD_GO_VERSION="v1.1.0"

beads_have_binary() {
    [ -n "$(beads_find_binary)" ]
}

# ---------------------------------------------------------------- install ---
#
# npm first, `go install` as a fallback.  npm is not merely faster, it is the
# route that works: the npm package carries a prebuilt, CGO-enabled binary and
# installs in about three seconds from registry.npmjs.org, which the agent
# proxy allows.  Building from source runs into a pincer:
#
#   * a normal (cgo) `go install` needs the ICU development headers that
#     github.com/dolthub/go-icu-regex compiles against, and these containers
#     ship libicu but not libicu-dev, so it dies minutes in with
#     "fatal error: unicode/uregex.h: No such file or directory"; while
#   * CGO_ENABLED=0 builds fine and then refuses to run: "embedded Dolt
#     requires a CGO build, but this bd binary was built with CGO_ENABLED=0" -
#     and .beads/metadata.json selects embedded mode.
#
# So the go path only works where the ICU headers are present.  It is kept for
# exactly that case, and for hosts with no npm.
beads_install_binary() {
    beads_have_binary && return 0

    if command -v npm >/dev/null 2>&1; then
        echo "beads: installing bd ${BD_NPM_VERSION} from npm (first use)..." >&2
        if mkdir -p "$BEADS_NPM_PREFIX" 2>/dev/null &&
           npm install --prefix "$BEADS_NPM_PREFIX" "@beads/bd@${BD_NPM_VERSION}" >&2; then
            # Never trust npm's exit status alone; verify a binary appeared.
            beads_have_binary && return 0
            echo "beads: npm reported success but produced no usable bd." >&2
        else
            echo "beads: 'npm install @beads/bd' failed; falling back to 'go install'." >&2
        fi
    fi

    if ! command -v go >/dev/null 2>&1; then
        echo "beads: neither npm nor go is available - cannot install bd." >&2
        return 1
    fi

    echo "beads: building bd ${BD_GO_VERSION} with 'go install' (several minutes)..." >&2
    # GOTOOLCHAIN=auto lets go fetch the newer toolchain beads' go.mod requires.
    # cgo is left ENABLED on purpose - see the note above.
    if ! GOTOOLCHAIN=auto GOFLAGS=-mod=mod go install \
        "github.com/steveyegge/beads/cmd/bd@${BD_GO_VERSION}" >&2; then
        echo "beads: 'go install bd' failed - bd unavailable.  If it failed on" >&2
        echo "       unicode/uregex.h, install the ICU dev headers (libicu-dev)." >&2
        return 1
    fi
    if ! beads_have_binary; then
        echo "beads: 'go install' reported success but no bd binary was found." >&2
        return 1
    fi
    return 0
}

# ---------------------------------------------------------------- hydrate ---

# Is the database usable?  A health probe, because `bd init` creates the
# directory BEFORE importing and `bd prime` creates an empty database as a side
# effect, so directory presence would latch "hydrated" after a partial import
# and leave every query silently returning nothing.
beads_hydrated() {
    _n="$("$1" count --quiet 2>/dev/null)" || return 1
    [ -n "$_n" ] && [ "$_n" -gt 0 ] 2>/dev/null
}

# `bd init --stealth` still normalizes a few tracked files.  Restore them, but
# only when they were clean beforehand - never clobber a developer's edits.
_bd_did_init=0
_bd_pre_dirty=1   # default "don't restore", so a trap firing early is a no-op

# shellcheck disable=SC2329  # invoked from the EXIT/INT/TERM traps below
beads_restore_tracked() {
    [ "$_bd_did_init" = 1 ] || return 0
    [ "$_bd_pre_dirty" = 0 ] || return 0
    # A path-limited `git checkout -- <path>` still fires the post-checkout
    # hook (measured on git 2.43), and .beads/hooks/post-checkout invokes bd -
    # which, with the shim on PATH and the database just torn down, would
    # re-enter this very script.  The nested run then blocks on the lock this
    # one holds, times out, and every restore fails: the safeguard defeats
    # itself exactly when it is needed.  Tell the shim to stand down.
    BEADS_SHIM_NO_INSTALL=1
    export BEADS_SHIM_NO_INSTALL
    for _f in .beads/config.yaml .beads/.gitignore .gitignore; do
        if git ls-files --error-unmatch "$_f" >/dev/null 2>&1; then
            git checkout -- "$_f" ||
                echo "beads: warning: failed to restore $_f - tree may be left dirty" >&2
        fi
    done
}
trap 'beads_restore_tracked' EXIT
# A trap alone does not terminate the shell on INT/TERM (dash keeps running
# past a caught signal), so exit explicitly.  This re-fires the EXIT trap,
# which is a harmless no-op the second time.
trap 'beads_restore_tracked; exit 130' INT
trap 'beads_restore_tracked; exit 143' TERM

# One mutex around the whole of install AND hydrate.
#
# flock(1) rather than a hand-rolled mkdir+pidfile: the kernel makes the
# acquire atomic, and a holder that dies for any reason - including SIGKILL,
# which no trap can catch - has its lock released when its descriptors close.
#
# -w, not -n.  A non-blocking acquire makes any bd command issued during the
# ~15 s setup window fail outright, which is worst of all for the background
# warm-up (BEADS_BOOTSTRAP=1): the optimisation would break the very command
# it is meant to speed up.  Waiting is what the caller wants; the bound stops
# a wedged holder hanging a session forever.
BEADS_LOCK_WAIT="${BEADS_LOCK_WAIT:-300}"

beads_acquire_lock() {
    if ! command -v flock >/dev/null 2>&1; then
        echo "beads: 'flock' not found - cannot safely coordinate setup." >&2
        return 1
    fi
    mkdir -p .beads 2>/dev/null || true
    # Create the lock file before `exec` touches it: `exec` is a special
    # builtin, so in dash a redirection failure exits the shell outright and
    # the diagnostic below would never print.
    : > .beads/.bootstrap.lock 2>/dev/null || {
        echo "beads: cannot create .beads/.bootstrap.lock - skipping setup." >&2
        return 1
    }
    exec 9>".beads/.bootstrap.lock"
    # fd 9 stays open for the rest of this process, which is what makes the
    # lock self-release on any exit path.
    if ! flock -w "$BEADS_LOCK_WAIT" 9; then
        echo "beads: timed out after ${BEADS_LOCK_WAIT}s waiting for another setup to finish." >&2
        return 1
    fi
    return 0
}

beads_hydrate() {
    _bd="$1"
    [ -f .beads/issues.jsonl ] || return 0
    # Re-check under the lock: another process may have finished while we
    # waited for it.
    if beads_hydrated "$_bd"; then
        # A healthy database with no marker is one this tooling did not create:
        # a developer's existing database, or a container warm from before the
        # marker existed.  ADOPT it - write the marker and leave it alone.
        #
        # It must never fall through to the rm -rf below.  The database can
        # hold issues that were never exported to .beads/issues.jsonl, and
        # re-importing would silently discard them.  With BEADS_BOOTSTRAP=1
        # that would happen in a background job at session start, with output
        # going to /dev/null.
        if [ ! -f "$(beads_db_marker "$BD_REPO")" ]; then
            : > "$(beads_db_marker "$BD_REPO")" 2>/dev/null ||
                echo "beads: could not write the hydration marker; setup will re-check next time." >&2
        fi
        return 0
    fi

    _bd_did_init=1
    # Record whether the tree was already dirty, so the restore only runs on a
    # clean tree.  Checked separately from emptiness: a `git status` that
    # itself errors (index lock, git missing, cwd outside a repo) must not read
    # as "clean".
    if _bd_dirty_out="$(git status --porcelain -- .beads/config.yaml .beads/.gitignore .gitignore 2>/dev/null)"; then
        [ -z "$_bd_dirty_out" ] && _bd_pre_dirty=0
    fi

    echo "beads: hydrating the issue database from .beads/issues.jsonl..." >&2
    # Only reached when `bd count` could not read any issues out of it, so
    # there is nothing here to lose - a healthy database was adopted above and
    # returned before this point.  Clearing is necessary because
    # `bd init --from-jsonl` refuses to run against ANY existing database
    # ("already initialized"), and `bd prime` leaves an empty one behind as a
    # side effect.
    rm -rf .beads/embeddeddolt
    # --stealth keeps beads files out of git (via .git/info/exclude), so init
    #   makes NO commits - critical when this runs from a hook or a shim;
    # --from-jsonl imports the committed JSONL in the same step;
    # --non-interactive / --quiet for a non-TTY container.
    if ! "$_bd" init --stealth --non-interactive --quiet --from-jsonl >/dev/null; then
        echo "beads: 'bd init --from-jsonl' failed - database unavailable." >&2
        rm -rf .beads/embeddeddolt
        return 1
    fi
    if ! beads_hydrated "$_bd"; then
        echo "beads: bd init reported success but the database has no issues - treating as a failed hydration." >&2
        rm -rf .beads/embeddeddolt
        return 1
    fi
    # Only now record success.  The marker lives inside the database directory,
    # so the rm -rf above removes it with everything else.
    : > "$(beads_db_marker "$BD_REPO")" || {
        echo "beads: could not write the hydration marker - setup will rerun next time." >&2
    }
    return 0
}

# ------------------------------------------------------------------- main ---
#
# The lock covers the install too, not just the hydration: two first-use
# commands (or the background warm-up plus an eager first command) would
# otherwise both run `npm install` into the same prefix, and a half-written
# node_modules/.bin/bd can be observed as "installed".
beads_acquire_lock || exit 1

beads_install_binary || exit 1

_bd_path="$(beads_find_binary)"
if [ -z "$_bd_path" ]; then
    echo "beads: no bd binary on PATH after install." >&2
    exit 1
fi

beads_hydrate "$_bd_path" || exit 1

exit 0
