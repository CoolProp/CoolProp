#!/usr/bin/env sh
#
# beads-install.sh - install `bd` and hydrate its database, on demand.
#
# This does the slow part of getting beads working in an ephemeral container:
#
#   * `bd` is not installed, and the upstream `curl | bash` installer pulls a
#     prebuilt binary from GitHub Releases, which the agent proxy blocks (403).
#     `go install` through the allowlisted module proxy works instead, and
#     beads embeds its Dolt engine, so no separate `dolt` binary is needed.
#   * The embedded DB dir (.beads/embeddeddolt/) is gitignored and therefore
#     absent, so it has to be rehydrated from the committed .beads/issues.jsonl,
#     which is the source of truth.
#
# Together those take two to three minutes, which is why nothing calls this at
# session start.  dev/ci/bd-shim.sh calls it on the first actual `bd` command,
# so a session that never touches the issue tracker never pays for it.
#
# Idempotent and non-fatal.  A warm container returns almost immediately (both
# steps no-op).  Any failure warns on stderr and returns non-zero; the caller
# decides what to do about it, and the session carries on without bd.
#
# All diagnostics go to stderr.  Stdout is left clean, because the shim may be
# running inside a pipeline.

set -u

# Pinned: a future release must not silently wedge a cold start.  The two
# differ because they are two different distributions of bd — the npm package
# is versioned 1.3.0 and carries a prebuilt binary, while v1.1.0 is the last
# version verified to build from the `github.com/steveyegge/beads` module path
# (the project has since moved to `gastownhall/beads`).
BD_NPM_VERSION="1.3.0"
BD_GO_VERSION="v1.1.0"

# bd is installed into a private npm prefix rather than the global tree.
#
# `npm install -g` is not reliably repeatable: after `npm uninstall -g`, the
# next global install of the same package reports "changed 1 package", exits 0
# and installs nothing at all - not even with --force, and not after clearing
# the empty scope directory it leaves behind.  A prefixed install writes a
# self-contained tree we own, is repeatable after a wipe, needs no root, and
# cannot disturb anything else on the machine.
#
# Outside the repo on purpose: shared across worktrees, and nothing to
# gitignore.  Kept in step with bd-shim.sh.
BD_NPM_PREFIX="${XDG_CACHE_HOME:-${HOME:-/tmp}/.cache}/coolprop/beads"
BD_NPM_BIN="${BD_NPM_PREFIX}/node_modules/.bin"

# This script lives at <repo>/dev/ci/, so the repo is two levels up.  Derive it
# from our own location rather than the caller's cwd: `bd` can be invoked from
# anywhere, and the JSONL we hydrate from is repo-relative.
_bd_self="$0"
_bd_here="$(CDPATH='' cd -- "$(dirname -- "$_bd_self")" && pwd)" || exit 1
BD_REPO="$(CDPATH='' cd -- "${_bd_here}/../.." && pwd)" || exit 1

cd "$BD_REPO" || exit 1

# Absolute path to a real bd binary, or nothing.
#
# Walks PATH rather than guessing install locations, so it finds bd however it
# got there — npm, go install, brew, a hand-placed binary — and skips this
# repo's own shim, which would otherwise be exec'd as itself forever.
#
# `command -v bd` is deliberately not used: where the shim is installed as
# `bd`, that is exactly what it returns.
#
# Intentionally duplicated in dev/ci/bd-shim.sh.  The shim runs on every single
# bd command and must not fork a subshell to source a shared helper; these two
# copies must stay in step.
beads_binary() {
    _shim="$1"   # realpath of the shim to skip, or empty
    {
        printf '%s\n' "$BD_NPM_BIN"
        printf '%s' "$PATH" | tr ':' '\n'
    } | while IFS= read -r _d; do
        [ -n "$_d" ] || continue
        [ -x "${_d}/bd" ] || continue
        _r="$(readlink -f "${_d}/bd" 2>/dev/null)" || _r="${_d}/bd"
        [ "$_r" = "$_shim" ] && continue
        printf '%s\n' "${_d}/bd"
        break
    done
}

BD_SHIM="$(readlink -f "${_bd_here}/bd-shim.sh" 2>/dev/null)" || BD_SHIM="${_bd_here}/bd-shim.sh"

beads_have_binary() {
    [ -n "$(beads_binary "$BD_SHIM")" ]
}

# ---------------------------------------------------------------- install ---

# npm first, `go install` as a fallback.
#
# npm is not just faster here, it is the route that works.  The npm package
# carries a prebuilt, CGO-enabled binary and installs in about three seconds
# from registry.npmjs.org, which the agent proxy allows.  Building from source
# instead runs into a pincer:
#
#   * a normal (cgo) `go install` needs the ICU development headers that
#     github.com/dolthub/go-icu-regex compiles against, and these containers
#     ship libicu but not libicu-dev, so it dies minutes in with
#     "fatal error: unicode/uregex.h: No such file or directory"; while
#   * CGO_ENABLED=0 builds fine and then refuses to run:
#     "embedded Dolt requires a CGO build, but this bd binary was built with
#     CGO_ENABLED=0" — and .beads/metadata.json selects embedded mode.
#
# So the go path only works where ICU headers are present. It is kept as a
# fallback for exactly that case, and for hosts with no npm.
beads_install_binary() {
    beads_have_binary && return 0

    if command -v npm >/dev/null 2>&1; then
        echo "beads: installing bd ${BD_NPM_VERSION} from npm (first use)..." >&2
        if mkdir -p "$BD_NPM_PREFIX" 2>/dev/null &&
           npm install --prefix "$BD_NPM_PREFIX" "@beads/bd@${BD_NPM_VERSION}" >&2; then
            # Never trust npm's exit status alone; verify the binary exists.
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
    # GOTOOLCHAIN=auto lets go fetch the newer toolchain beads' go.mod requires
    # (through the allowlisted module proxy) when the base toolchain is older.
    # cgo is left ENABLED on purpose - see the note above.
    if ! GOTOOLCHAIN=auto GOFLAGS=-mod=mod go install \
        "github.com/steveyegge/beads/cmd/bd@${BD_GO_VERSION}" >&2; then
        echo "beads: 'go install bd' failed - bd unavailable. If it failed on" >&2
        echo "       unicode/uregex.h, install the ICU dev headers (libicu-dev)." >&2
        return 1
    fi
    # go install drops the binary in GOPATH/bin, which is not necessarily on
    # PATH; add it so the PATH walk above can see it.
    _gp="$(go env GOPATH 2>/dev/null)" || _gp=""
    if [ -n "$_gp" ] && [ -x "${_gp}/bin/bd" ]; then
        case ":${PATH}:" in
            *":${_gp}/bin:"*) ;;
            *) PATH="${_gp}/bin:${PATH}"; export PATH ;;
        esac
    fi
    if ! beads_have_binary; then
        echo "beads: 'go install' reported success but no bd binary was found." >&2
        return 1
    fi
    return 0
}

# ---------------------------------------------------------------- hydrate ---

# "Already hydrated" is a health probe, not a directory check: `bd init`
# creates .beads/embeddeddolt/ *before* importing, so directory presence would
# latch "hydrated" forever after a partial or failed import.  `bd count` also
# exits non-zero when there is no DB at all, so one probe covers both "never
# initialized" and "initialized but empty".
beads_hydrated() {
    _bd="$1"
    _n="$("$_bd" count --quiet 2>/dev/null)" || return 1
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
    for _f in .beads/config.yaml .beads/.gitignore .gitignore; do
        if git ls-files --error-unmatch "$_f" >/dev/null 2>&1; then
            git checkout -- "$_f" ||
                echo "beads: warning: failed to restore $_f - tree may be left dirty" >&2
        fi
    done
}
trap 'beads_restore_tracked' EXIT
# A trap alone does not terminate the shell on INT/TERM (dash keeps running
# past a caught signal), so exit explicitly. This re-fires the EXIT trap, which
# is a harmless no-op the second time.
trap 'beads_restore_tracked; exit 130' INT
trap 'beads_restore_tracked; exit 143' TERM

# Single-shot mutex around the hydration critical section.  Two shims racing in
# one container would otherwise both see "not hydrated" and both rm -rf and
# reinit .beads/embeddeddolt concurrently.
#
# flock(1) rather than a hand-rolled mkdir+pidfile: the kernel makes the
# acquire itself atomic, with no TOCTOU window, and a holder that dies for any
# reason - including SIGKILL, which no trap can catch - has its lock released
# when its descriptors close.  No staleness or steal logic to get subtly wrong.
#
# Fails CLOSED.  The whole point is exclusivity, so "could not tell" has to be
# treated as "someone else holds it".
beads_acquire_lock() {
    if ! command -v flock >/dev/null 2>&1; then
        echo "beads: 'flock' not found - cannot safely coordinate hydration." >&2
        return 1
    fi
    # No 2>/dev/null on this exec: a bare `exec` applies its redirects to the
    # CURRENT shell permanently, so that would silently swallow every later
    # stderr warning for the rest of this script's life, not just a failure of
    # this one open.
    if ! exec 9>".beads/.bootstrap.lock"; then
        echo "beads: could not open .beads/.bootstrap.lock for locking." >&2
        return 1
    fi
    # fd 9 stays open for the rest of this process, which is what makes the
    # lock self-release on any exit path.
    flock -n 9
}

beads_hydrate() {
    _bd="$1"
    [ -f .beads/issues.jsonl ] || return 0
    beads_hydrated "$_bd" && return 0

    if ! beads_acquire_lock; then
        # Either someone else is mid-hydration, or the lock could not be taken
        # (a specific warning was already logged).  Either way the DB is not
        # safe to touch, and not safe to read.
        echo "beads: could not acquire the hydration lock - skipping; retry the command shortly." >&2
        return 1
    fi

    _bd_did_init=1
    # Record whether the tree was already dirty, so the restore above only runs
    # on a clean tree.  Checked separately from emptiness: a `git status` that
    # itself errors (index lock, git missing, cwd outside a repo) must not read
    # as "clean".  stderr is discarded rather than merged, so incidental noise
    # on an otherwise-empty status cannot be mistaken for "dirty".
    if _bd_dirty_out="$(git status --porcelain -- .beads/config.yaml .beads/.gitignore .gitignore 2>/dev/null)"; then
        [ -z "$_bd_dirty_out" ] && _bd_pre_dirty=0
    fi

    echo "beads: hydrating the issue database from .beads/issues.jsonl..." >&2
    # The probe above just said "not hydrated", so nothing of value is here.
    # Clearing matters because `bd prime` lazily creates an empty Dolt DB as an
    # undocumented side effect, and `bd init --from-jsonl` refuses to run
    # against ANY existing DB ("already initialized").
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
    return 0
}

# ------------------------------------------------------------------- main ---

beads_install_binary || exit 1

_bd_path="$(beads_binary "$BD_SHIM")"
if [ -z "$_bd_path" ]; then
    echo "beads: no bd binary on PATH after install." >&2
    exit 1
fi

beads_hydrate "$_bd_path" || exit 1

exit 0
