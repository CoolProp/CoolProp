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
# Opt-in via BEADS_BOOTSTRAP=1, not auto-detected from CI/CLAUDE_CODE_REMOTE:
# the cold-start install+hydrate can take a couple of minutes, a cost that
# should be a deliberate choice for an environment, not a surprise sprung on
# every session in every ephemeral container. Set it once in the persistent
# env config for an environment that actually wants beads tracking; every
# session there then gets bd automatically, cold start included. Without it,
# this script only primes an already-installed bd (near-zero cost) and gets
# out of the way — never installs, never hydrates.
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

# Opt-in gate: only bootstrap when explicitly enabled for this environment.
# Elsewhere, just prime an already-installed bd (if any) and get out of the
# way — no install, no hydration, no cold-start cost.
if [ "${BEADS_BOOTSTRAP:-}" != "1" ]; then
    if command -v bd >/dev/null 2>&1; then
        bd prime || echo "beads-bootstrap: 'bd prime' failed (exit $?)." >&2
    fi
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
_bd_lockfile=".beads/.bootstrap.lock"

_beads_bootstrap_cleanup() {
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
    hash -r || echo "beads-bootstrap: 'hash -r' failed — a stale command lookup cache could mask the freshly-installed bd." >&2
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
        hash -r || echo "beads-bootstrap: 'hash -r' failed — a stale command lookup cache could mask the freshly-installed bd." >&2
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
    # and both rm -rf/reinit .beads/embeddeddolt concurrently. Use flock(1)
    # rather than a hand-rolled mkdir+pid-file scheme: the kernel guarantees
    # the acquire itself is atomic (no TOCTOU window to race), and a holder
    # that dies for any reason — including SIGKILL, which no trap can catch
    # — has its lock released automatically when its file descriptors close,
    # so there's no separate staleness/steal logic to get subtly wrong.
    #
    # Fail CLOSED, not open, if the lock itself can't be obtained (flock
    # missing, or the lock file can't be opened): the entire point of this
    # function is to guarantee exclusivity, so "couldn't tell" must be
    # treated the same as "someone else has it" — skip this session and let
    # the next session's `bd count` probe retry from scratch, rather than
    # silently racing an unprotected rm -rf/bd init. flock ships with
    # util-linux and is present on the Ubuntu-based containers this targets,
    # so this fallback is defensive, not the expected path.
    _beads_bootstrap_acquire_lock() {
        if ! command -v flock >/dev/null 2>&1; then
            echo "beads-bootstrap: 'flock' not found — cannot safely coordinate hydration." >&2
            return 1
        fi
        # No `2>/dev/null` here: a bare `exec` (no command) applies its
        # redirects to the CURRENT shell persistently, not just to this one
        # statement — `2>/dev/null` on it would silently swallow every
        # later stderr warning for the rest of the script's life, not just
        # a failure of this one open. Let a genuine open failure print, and
        # warn explicitly too so "lock open failed" isn't indistinguishable
        # from "lock acquired" in the caller.
        if ! exec 9>"$_bd_lockfile"; then
            echo "beads-bootstrap: could not open $_bd_lockfile for locking." >&2
            return 1
        fi
        # fd 9 stays open (and inherited by bd/go/git children) for the rest
        # of this process's life — that's what makes the flock self-release
        # on any exit. Safe today because bd's embedded Dolt engine runs
        # in-process per invocation rather than forking a long-lived
        # daemon (per the header comment above); revisit if that changes.
        flock -n 9
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
            # The tracked-file restore happens in the EXIT trap above, so it
            # also fires on an interrupted hook, not just a clean finish;
            # the flock itself needs no explicit release — it's tied to fd 9
            # and the kernel drops it when this process's descriptors close,
            # on any exit path including SIGKILL.
        else
            # Either another process holds the lock and is actively
            # hydrating (or tearing down and rebuilding) .beads/embeddeddolt
            # right now, or the lock itself couldn't be obtained (see the
            # specific warning already logged inside
            # _beads_bootstrap_acquire_lock) — both cases mean we can't
            # safely touch the DB this run. Don't call `bd prime` either —
            # it can read mid-rebuild state, or itself write to the DB (see
            # the comment above `rm -rf .beads/embeddeddolt`) — race a
            # concurrent hydration we can't see. Skip priming this session
            # entirely; the next session's `bd count` probe re-evaluates
            # from scratch.
            echo "beads-bootstrap: could not acquire the hydration lock — skipping this session; it will retry next session." >&2
            exit 0
        fi
    fi

    # 3. Hand off to the beads-managed workflow-context hook.
    bd prime || echo "beads-bootstrap: 'bd prime' failed (exit $?)." >&2
fi

exit 0
