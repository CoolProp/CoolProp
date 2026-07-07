# Project Instructions for AI Agents

This file provides instructions and context for AI coding agents working on this project.

<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:ca08a54f -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

## Session Completion

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   bd dolt push
   git push
   git status  # MUST show "up to date with origin"
   ```
5. **Clean up** - Clear stashes, prune remote branches
6. **Verify** - All changes committed AND pushed
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
<!-- END BEADS INTEGRATION -->

### `bd` in ephemeral (Claude Code web/CI) containers

Fresh containers have no `bd` binary and no Dolt DB (the DB dir
`.beads/embeddeddolt/` is gitignored). `dev/ci/bootstrap-beads.sh` — wired
as the first `SessionStart` hook in `.claude/settings.json`, ahead of the
beads-managed `bd prime` — installs `bd` via `go install` (the `curl|bash`
installer's GitHub-release download is blocked by the agent proxy; the Go
module proxy is allowed) and rehydrates the DB from the committed
`.beads/issues.jsonl`. It is idempotent and non-fatal. First cold start
takes ~2 min (Go toolchain fetch + build); warm containers skip both steps.
To persist newly filed issues, append their `bd export --include-memories`
lines to `.beads/issues.jsonl` and commit — do **not** overwrite the file
wholesale (bd 1.1.0 re-serializes the `dependencies` field on unrelated
issues, and a plain `bd export` drops the persisted memories). `bd dolt
push` in the session-completion steps above is a no-op here — there is no
Dolt remote in web containers; committing `issues.jsonl` is the sync.


## Build & Test

```bash
# Configure (one-time).  Always set -DCMAKE_BUILD_TYPE=Release: the
# default is an UNOPTIMIZED build, which makes the [slow] SVDSBTL surface
# builds (a dense SVD at production resolution, NT=200/NR=800/rank=20)
# crawl for minutes per fresh fluid.  Release turns that into seconds.
cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release

# Build the Catch2 runner
cmake --build build_catch --target CatchTestRunner -j8

# Run the test suite — see "Test filter discipline" below for tag scope
./build_catch/CatchTestRunner [SBTL]            # SBTL adapter layer
./build_catch/CatchTestRunner [SVDSBTL]         # backend-level tests
./build_catch/CatchTestRunner "[!slow]"         # everything fast
```

## Pre-Push Gate — REQUIRED before every `git push`

Run `./dev/ci/preflight.sh` before any `git push` (or install the
pre-push hook once via `ln -s ../../dev/ci/pre-push.sample
.git/hooks/pre-push`).  The script mirrors what CI runs:

- clang-format dry-run vs `origin/master` (version pinned from `.pre-commit-config.yaml`)
- build CatchTestRunner
- Catch2 tests with auto-selected tag scope based on changed paths
- cppcheck (`--enable=warning`) on changed files
- clang-tidy diff-only (requires LLVM 18+ on PATH)
- semgrep `p/security-audit` + local `.semgrep/` rules (uvx-resolved)

If preflight passes, CI should pass with high probability.  See
`dev/ci/README.md#preflightsh--local-pre-push-gate` for details.

**`git commit --no-verify` only skips pre-commit hooks (clang-format,
bd auto-export).  It does NOT skip the pre-push gate.**  If you must
push without preflight, use `git push --no-verify` and document why.

## Conventions & Patterns

### Test filter discipline

When changes touch files under `src/SBTL/`, `include/CoolProp/sbtl/`,
`src/Backends/SVDSBTL/`, or `src/Region/`, run the **umbrella**
`[SBTL]` tag locally — NOT just `[SVDSBTL]`.  The SBTL adapter layer
(serializer round-trip, multi-fluid PH preset tests) lives under
`[SBTL]` only; narrowing to `[SVDSBTL]` misses tests that bite in CI.
`./dev/ci/preflight.sh` auto-selects the right umbrella tag from the
changed-file paths; running preflight is the safe default.

### Pre-PR adversarial review — MANDATORY

Before any `gh pr create` invocation, **you MUST** invoke an adversarial
review subagent against the diff.  The pre-push shell hook can't
mechanically gate this (subagents are a Claude Code construct, not a CLI),
so the gate is procedural — but it is NOT optional.  Every recent PR where
CodeRabbit found a blocking issue (null-deref, FD-out-of-range,
noexcept-on-throwing, NaN-absorbed-by-std::max, and `|| true` masking that
silently disabled a CI gate) would have been caught by the reviewer first
at zero CI latency.  Skip this step and the same class of findings keeps
recurring.

**Agent availability:** prefer `subagent_type: "superpowers:code-reviewer"`;
if that agent type isn't registered in the current environment, fall back to
`subagent_type: "general-purpose"` with the same prompt.  Do NOT skip the
review because the named agent is missing.

**Pre-`gh pr create` checklist (MUST complete in order):**

1. `./dev/ci/preflight.sh` passes (or you can explain each `--skip`).
2. The review subagent returns with no blocking findings, OR you've
   addressed/justified each one.  Review the **final** diff vs the branch's
   actual base (for a stacked branch that's the parent branch, NOT always
   `origin/master`).
3. `git push` (the pre-push hook re-runs preflight as a safety net).
4. THEN `gh pr create`.
5. **Re-review the delta.**  Any commit pushed AFTER step 2 (review-feedback
   fixes, follow-on changes, CI fixes) is unreviewed — re-run the review on
   those commits before treating the PR as done.  New files added late (e.g.
   a benchmark program) are the easiest to ship unreviewed.

Canonical review invocation:

```
Agent({
  subagent_type: "superpowers:code-reviewer",   // fall back to "general-purpose" if unavailable
  description: "Pre-PR review of <branch>",
  prompt: "Adversarial review of the diff between <branch> and its actual base.
           Project conventions are in CLAUDE.md.  Check for:
             - null-deref risk on shared_ptr inputs
             - noexcept on functions that can throw (resize/allocate/etc.)
             - FD stencils that step outside their valid range
             - non-finite values silently absorbed by std::max/std::min
             - bit-exact compare where ULP-class noise exists
             - fopen without permission restriction
             - preconditions checked in the factory but not the public constructor
             - integer narrowing / unsigned-overflow slipping past a range check
             - CI/shell gate robustness: `|| true` or `2>/dev/null` that masks a
               failure, swallowed exit codes, skip-on-error in a gate, `grep -c`
               without `|| echo 0` under `set -e`, unguarded `find -o` precedence
             - implementation that diverges from its own spec/plan without
               updating the doc
             - bot-comment-class issues that recur across recent PRs
           For EVERY check, gate, guard, or validation step in the diff, answer
           explicitly: 'what input or failure makes this pass when it should
           fail?'  A safeguard that can silently no-op (fail-open) is a BLOCKING
           finding even if it does NOT fail CI.
           Report blocking findings.  If you notice an issue and judge it an
           acceptable trade-off, STILL report it as a flagged trade-off — never
           silently dismiss it.  Skip pure style nits."
})
```

### REFPROP is available locally

REFPROP is installed on Ian's primary dev machine — `[refprop]`-tagged
Catch2 tests *run* locally, they don't silently SKIP.  Don't narrow
test filters assuming REFPROP-only tests are CI-gated; they're catchable
locally and should be in the pre-push sweep.

### `.beads/issues.jsonl` should not be in source PRs

The `bd` pre-commit hook auto-exports + stages `.beads/issues.jsonl`.
For source-code PRs, restore it explicitly before committing:

```bash
git restore --staged .beads/issues.jsonl
git checkout .beads/issues.jsonl
git commit --no-verify -m "..."   # --no-verify so the hook doesn't re-add it
```

The `--no-verify` here is intentional and scoped to ONE commit.  Always
run `./dev/ci/preflight.sh` separately afterwards since `--no-verify`
also skips clang-format checking.
