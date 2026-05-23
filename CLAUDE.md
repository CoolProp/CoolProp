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


## Build & Test

```bash
# Configure (one-time)
cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON

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

### Pre-PR adversarial review

Before opening a PR (`gh pr create`), run the `superpowers:code-reviewer`
subagent against the diff with the project conventions in this file as
the standard.  It catches the same null-deref / edge-stencil / dead-arg
class of findings that CodeRabbit flags post-push, at zero CI-cycle
latency.  Example invocation:

```
Agent({
  subagent_type: "superpowers:code-reviewer",
  description: "Pre-PR review of <branch>",
  prompt: "Adversarial review of the diff between <branch> and origin/master.
           Check for: null-deref risk on shared_ptr inputs, FD stencils that
           step outside their valid range, fopen without permission restriction,
           tests that use bit-exact compare where ULP-class noise exists,
           bot-comment-class issues that recur across the SVDSBTL PRs.
           Report blocking findings only."
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
