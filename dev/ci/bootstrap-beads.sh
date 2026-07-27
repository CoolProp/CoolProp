#!/usr/bin/env sh
#
# bootstrap-beads.sh — make `bd` (beads) usable in ephemeral containers.
#
# Claude Code web/CI containers start from a fresh clone, where:
#   * the `bd` binary is not installed, and the upstream `curl | bash`
#     installer downloads a prebuilt binary from GitHub Releases — which the
#     agent proxy blocks (403 on api.github.com / release asset hosts); but
#   * `go install` through the allowlisted module proxy (proxy.golang.org)
#     works, and beads embeds its Dolt engine, so no separate `dolt` binary
#     is needed; and
#   * the embedded DB dir (.beads/embeddeddolt/) is gitignored and therefore
#     absent — it must be rehydrated from the committed .beads/issues.jsonl,
#     which is the source of truth.
#
# This is the ONLY SessionStart hook for beads (see .claude/settings.json).
# Claude Code runs all hooks for one event concurrently — array position does
# not sequence them — so a separate `bd prime` hook would race a `bd` this
# script hasn't finished installing yet. Instead this script runs `bd prime`
# itself as its last step, once bd and the DB are ready.
#
# Restricted to recognized ephemeral containers (CI=true or CLAUDE_CODE_REMOTE
# set): a developer workstation with `go` on PATH must never get an
# unprompted `go install` / toolchain fetch and multi-minute session stall.
# Set BEADS_BOOTSTRAP_FORCE=1 to opt in anywhere else (e.g. to exercise this
# script locally).
#
# It is idempotent (hydration is checked with a `bd count` health probe, not
# directory presence, so a partial/failed import is retried next session
# rather than latching "done" forever) and strictly non-fatal: a warm
# container already has bd on PATH and the DB in place (both steps no-op),
# and any failure logs a warning and lets the session continue without bd —
# the git hooks are already guarded the same way.

set -u

BD_VERSION="v1.1.0"  # pinned: a future release must not silently wedge cold start

cd "$(git rev-parse --show-toplevel 2>/dev/null)" 2>/dev/null || exit 0

# Container gate: only bootstrap in recognized ephemeral environments, unless
# explicitly opted in. Elsewhere, just prime an already-installed bd (if any)
# and get out of the way.
if [ -z "${BEADS_BOOTSTRAP_FORCE:-}" ] \
    && [ "${CI:-}" != "true" ] \
    && [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
    command -v bd >/dev/null 2>&1 && bd prime
    exit 0
fi

# Tracks whether we've started an import this run, and whether the tree was
# clean beforehand — read by the trap below so an interrupted or completed
# hydration never leaves the tree dirty (bd init normalizes a few tracked
# config files even under --stealth; see step 2). Default pre_dirty=1 (i.e.
# "don't restore") so a trap firing before either variable is computed can
# never touch the tree.
_bd_did_init=0
_bd_pre_dirty=1
_bd_lockdir=".beads/.bootstrap.lock"
_bd_lock_held=0

_beads_bootstrap_cleanup() {
    [ "$_bd_lock_held" = 1 ] && rm -rf "$_bd_lockdir" 2>/dev/null
    [ "$_bd_did_init" = 1 ] || return 0
    [ "$_bd_pre_dirty" = 0 ] || return 0
    for _f in .beads/config.yaml .beads/.gitignore .gitignore; do
        if git ls-files --error-unmatch "$_f" >/dev/null 2>&1; then
            git checkout -- "$_f" || echo "beads-bootstrap: warning: failed to restore $_f — tree may be left dirty" >&2
        fi
    done
}
trap '_beads_bootstrap_cleanup' EXIT
# A trap alone does not terminate the shell on INT/TERM (dash keeps running
# past the signal once it's caught) — exit explicitly after cleanup so an
# external kill (e.g. a hook timeout) still takes effect promptly. This also
# re-fires the EXIT trap above, which is a harmless no-op the second time.
trap '_beads_bootstrap_cleanup; exit 130' INT
trap '_beads_bootstrap_cleanup; exit 143' TERM

# bd installs to $GOPATH/bin; make sure that's on PATH for this session.
if command -v go >/dev/null 2>&1; then
    _gobin="$(go env GOPATH 2>/dev/null)/bin"
    # On Windows/Git Bash, `go env GOPATH` prints a Windows path
    # (C:\Users\...\go); a bare PATH-membership test never matches it, and
    # prepending it verbatim would split on the drive-letter colon.
    if command -v cygpath >/dev/null 2>&1; then
        _gobin="$(cygpath -u "$_gobin")"
    fi
    case ":${PATH}:" in
        *":${_gobin}:"*) ;;
        *) PATH="${_gobin}:${PATH}"; export PATH ;;
    esac
    hash -r 2>/dev/null || true
fi

# 1. Install bd if missing.
if ! command -v bd >/dev/null 2>&1; then
    if command -v go >/dev/null 2>&1; then
        echo "beads-bootstrap: installing bd ${BD_VERSION} via 'go install'..." >&2
        # GOTOOLCHAIN=auto lets go fetch the newer toolchain beads' go.mod
        # requires (via the allowlisted module proxy) when the base toolchain
        # is older.
        if ! GOTOOLCHAIN=auto GOFLAGS=-mod=mod go install \
            "github.com/steveyegge/beads/cmd/bd@${BD_VERSION}" >&2; then
            echo "beads-bootstrap: 'go install bd' failed — bd commands unavailable this session." >&2
        fi
        hash -r 2>/dev/null || true
    else
        echo "beads-bootstrap: 'go' not found — cannot install bd." >&2
    fi
fi

if command -v bd >/dev/null 2>&1; then
    # 2. Hydrate the embedded Dolt DB from the committed JSONL if needed.
    #    "already hydrated" is checked with a health probe (bd count), not
    #    directory presence: bd init creates .beads/embeddeddolt/ *before*
    #    importing, so a partial/failed import would otherwise leave a
    #    directory that latches "hydrated" forever. `bd count` also exits
    #    non-zero when no DB exists at all, so a single probe covers both
    #    "never initialized" and "initialized but empty".
    _bd_hydrated() {
        _n="$(bd count --quiet 2>/dev/null)" || return 1
        [ -n "$_n" ] && [ "$_n" -gt 0 ] 2>/dev/null
    }

    # Single-shot mutex around the hydration critical section below: two
    # SessionStart hooks racing in the same container (e.g. two sessions
    # attached to one environment) would otherwise both see "not hydrated"
    # and both rm -rf/reinit .beads/embeddeddolt concurrently. No retry loop
    # — if the lock is held by a live process we just skip hydration this
    # run (self-heals next session via the same `bd count` probe); if the
    # holder's PID is dead (e.g. it was SIGKILLed before it could release —
    # no trap catches that), steal the lock instead of wedging forever.
    _beads_bootstrap_acquire_lock() {
        if mkdir "$_bd_lockdir" 2>/dev/null; then
            echo "$$" >"$_bd_lockdir/pid" 2>/dev/null
            _bd_lock_held=1
            return 0
        fi
        _lock_pid="$(cat "$_bd_lockdir/pid" 2>/dev/null)"
        if [ -n "$_lock_pid" ] && ! kill -0 "$_lock_pid" 2>/dev/null; then
            rm -rf "$_bd_lockdir" 2>/dev/null
            if mkdir "$_bd_lockdir" 2>/dev/null; then
                echo "$$" >"$_bd_lockdir/pid" 2>/dev/null
                _bd_lock_held=1
                return 0
            fi
        fi
        return 1
    }

    if [ -f .beads/issues.jsonl ] && ! _bd_hydrated; then
        if _beads_bootstrap_acquire_lock; then
            _bd_did_init=1
            # Record whether the tree already had changes, so the restore
            # below (via the trap) only runs on a clean fresh container and
            # never touches a developer's in-progress edits to these config
            # files. Checked separately from emptiness: a `git status` that
            # itself errors (index lock contention, git missing, cwd outside
            # a repo) must not be read as "clean". stderr is discarded, not
            # merged into the captured text, so incidental stderr noise on
            # an otherwise-successful, otherwise-empty status can't be
            # mistaken for "dirty".
            if _bd_dirty_out="$(git status --porcelain -- .beads/config.yaml .beads/.gitignore .gitignore 2>/dev/null)"; then
                [ -z "$_bd_dirty_out" ] && _bd_pre_dirty=0
            fi

            echo "beads-bootstrap: hydrating beads DB from .beads/issues.jsonl..." >&2
            # Our own health probe just said "not hydrated", so nothing of
            # value is lost by clearing whatever's here first. This matters
            # because step 3's `bd prime` lazily creates an empty Dolt DB as
            # a side effect whenever none exists (confirmed by testing, not
            # documented behavior) — so a prior failed/skipped attempt in
            # this same container can leave an empty DB behind, and
            # `bd init --from-jsonl` refuses to run against ANY existing DB
            # ("already initialized") without this.
            rm -rf .beads/embeddeddolt
            # --stealth keeps beads files out of git (via .git/info/exclude), so
            #   init makes NO commits — critical in a hook, which must never
            #   auto-commit;
            # --from-jsonl imports the committed .beads/issues.jsonl in the same
            #   step;
            # --non-interactive / --quiet for the non-TTY container.
            if ! bd init --stealth --non-interactive --quiet --from-jsonl >/dev/null; then
                echo "beads-bootstrap: 'bd init --from-jsonl' failed — beads DB unavailable this session." >&2
                rm -rf .beads/embeddeddolt
            elif ! _bd_hydrated; then
                echo "beads-bootstrap: bd init reported success but the DB has no issues — treating as a failed hydration so it retries next session." >&2
                rm -rf .beads/embeddeddolt
            fi
            # The tracked-file restore and the lock release both happen in
            # the EXIT trap above, so they also fire on an interrupted hook,
            # not just a clean finish.
        else
            echo "beads-bootstrap: another bootstrap appears to be in progress (lock held) — skipping hydration this session; it will retry next session." >&2
        fi
    fi

    # 3. Hand off to the beads-managed workflow-context hook.
    bd prime
fi

exit 0
