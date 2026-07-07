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
# This script is registered as the FIRST SessionStart hook (see
# .claude/settings.json) so it runs before the beads-managed
# `bd prime --hook-json` hook, which then finds an installed bd and a
# populated database.
#
# It is idempotent and strictly non-fatal: a warm container already has bd on
# PATH and the DB in place (both steps no-op), and any failure logs a warning
# and exits 0 so a beads hiccup never blocks the session — the git hooks are
# already guarded the same way.

set -u

BD_VERSION="v1.1.0"  # pinned: a future release must not silently wedge cold start

# bd installs to $GOPATH/bin; make sure that's on PATH for this session.
if command -v go >/dev/null 2>&1; then
    _gobin="$(go env GOPATH 2>/dev/null)/bin"
    case ":${PATH}:" in
        *":${_gobin}:"*) ;;
        *) PATH="${_gobin}:${PATH}"; export PATH ;;
    esac
fi

# 1. Install bd if missing.
if ! command -v bd >/dev/null 2>&1; then
    if command -v go >/dev/null 2>&1; then
        echo "beads-bootstrap: installing bd ${BD_VERSION} via 'go install'..." >&2
        # GOTOOLCHAIN=auto lets go fetch the newer toolchain beads' go.mod
        # requires (via the allowlisted module proxy) when the base toolchain
        # is older.
        if ! GOTOOLCHAIN=auto GOFLAGS=-mod=mod go install \
            "github.com/steveyegge/beads/cmd/bd@${BD_VERSION}" >&2 2>&1; then
            echo "beads-bootstrap: 'go install bd' failed — bd commands unavailable this session." >&2
        fi
    else
        echo "beads-bootstrap: 'go' not found — cannot install bd." >&2
    fi
fi
command -v bd >/dev/null 2>&1 || exit 0

# 2. Hydrate the embedded Dolt DB from the committed JSONL if it isn't there.
#    Presence of the gitignored DB dir is the reliable "already hydrated" test
#    (bd's own commands exit 0 even when no DB exists, so exit codes can't be
#    used here).
if [ ! -d .beads/embeddeddolt ]; then
    # Record whether the tree already had changes, so the post-init restore
    # below only runs on a clean fresh container and never touches a
    # developer's in-progress edits to these config files.
    _bd_pre_dirty="$(git status --porcelain -- .beads/config.yaml .beads/.gitignore .gitignore 2>/dev/null)"
    echo "beads-bootstrap: hydrating beads DB from .beads/issues.jsonl..." >&2
    # --stealth keeps beads files out of git (via .git/info/exclude), so init
    #   makes NO commits and rewrites NO tracked files — critical in a hook,
    #   which must never dirty the tree or auto-commit;
    # --from-jsonl imports the committed .beads/issues.jsonl in the same step;
    # --non-interactive / --quiet for the non-TTY container.
    if ! bd init --stealth --non-interactive --quiet --from-jsonl >/dev/null 2>&1; then
        echo "beads-bootstrap: 'bd init --from-jsonl' failed — beads DB unavailable this session." >&2
    fi
    # --stealth avoids commits and the CLAUDE.md/AGENTS.md/settings rewrites,
    # but bd init still normalizes a few tracked config files
    # (.beads/config.yaml, .beads/.gitignore, .gitignore). Restore them so a
    # fresh container's tree stays clean; the gitignored embedded DB is kept
    # and remains usable (config.yaml holds only commented defaults). Guarded
    # so it only runs when the pre-hook tree was clean (fresh container), never
    # clobbering real in-progress edits.
    if [ -z "$_bd_pre_dirty" ]; then
        git checkout -- .beads/config.yaml .beads/.gitignore .gitignore 2>/dev/null || true
    fi
fi

exit 0
