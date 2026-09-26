# Project Instructions for AI Agents

This file provides instructions and context for AI coding agents working on this project.

## Attribution — REQUIRED

Work you produce here must be attributed. Every commit whose content you
generated or materially shaped carries a `Co-authored-by:` trailer naming the
model, and the PR body states what you did and what the human verified. Full
policy:
[`AGENTS.md`](AGENTS.md#ai-assistance-and-attribution).

```
Co-authored-by: Claude Opus 5 <noreply@anthropic.com>
```

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

A fresh container has no `bd` binary and no database.  The embedded Dolt dir
`.beads/embeddeddolt/` is gitignored, so it has to be rehydrated from the
committed `.beads/issues.jsonl`, which is the source of truth.  Installing
`bd` and hydrating takes about fifteen seconds (and minutes, if it ever has to
fall back to building from source).

**That cost is paid on first use, not at session start.**  The `SessionStart`
hook (`dev/ci/bootstrap-beads.sh`) does no installing: it puts
`dev/ci/bd-shim.sh` on PATH as `bd`, which takes milliseconds, and returns.
(If a real `bd` is already reachable as `bd` AND the database is already
hydrated, it primes that instead and installs no shim.  It never symlinks over
a real `bd` binary; it takes the name only when it is free, already one of our
own shims, or a dangling link - that last so a checkout that was moved or
renamed does not leave a broken `bd` nobody can repair.  It can still *shadow*
a real `bd` that sits later on PATH, which costs the 60 ms of resolution and
nothing else, since the shim then execs that very binary.)

The first actual `bd` command runs `dev/ci/beads-install.sh`: an npm install of
`@beads/bd` into a private prefix under `~/.cache/coolprop/beads` (about 3 s;
the upstream `curl | bash` installer's GitHub-Releases download is
proxy-blocked, and building from source needs ICU headers these containers
lack, while a `CGO_ENABLED=0` build cannot open an embedded Dolt database at
all), then `bd init --from-jsonl` (about 12 s), and then it execs the real
binary.  Every later command adds roughly 60 ms of resolution.  A session that
never opens the tracker installs nothing.

So **just run `bd prime`** when you need the tracker; the first command does
the setup and says so.  Nothing needs enabling.

- `BEADS_BOOTSTRAP=1` additionally starts the setup in the background at
  session start, so the first command usually finds it ready.  Still
  non-blocking: a command issued mid-setup waits on the installer's lock
  rather than failing.
- `BEADS_BOOTSTRAP=0` disables the hook entirely.
- `BEADS_SHIM_NO_INSTALL=1` primes an already-installed `bd` but never
  triggers an install.  Used by the `PreCompact` hook, and by the installer's
  own restore step.
- `BEADS_LOCK_WAIT` (default 300) is how many seconds the installer waits for
  another setup to finish before giving up.  Anything that is not a whole
  number of seconds is ignored with a warning.
- `BEADS_SHIM_DEPTH` is set and incremented by the shim itself and is not for
  callers.  Depth 1 is the ordinary outermost call and may do anything, a
  setup included; at depth 2 (the expected nesting, bd -> git -> hook -> bd)
  the shim still resolves and execs `bd` but starts no setup; above depth 2 it
  stands down.  It is only a backstop: what actually prevents an exec loop is
  that the shim recognises another shim by content before handing over.
- When the shim cannot run `bd` at all it exits **3**, never 127.  All five
  hooks in `.beads/hooks/` neutralise exactly two statuses, 3 ("database not
  initialized") and 124 (timeout), and propagate everything else, so a 127 out
  of `pre-commit` or `pre-push` would abort the commit or the push.
- The shim stands down inside git hooks, keyed on `BD_GIT_HOOK`, which all
  five hooks in `.beads/hooks/` already export.  They guard on
  `command -v bd`, which the shim satisfies, so without this an ordinary
  `git commit` (or `git checkout`, or `git push`) would stop to set up a
  tracker nobody asked for.  Do **not** key on `GIT_DIR`: measured on git
  2.43 it is not exported to hooks at all, and `GIT_INDEX_FILE` is set only
  for the index-touching hooks, so `post-checkout`, `post-merge` and
  `pre-push` would slip straight through.

The install goes into a private npm prefix rather than `npm install -g`,
because the global tree is not reliably repeatable: after `npm uninstall -g`,
the next global install of the same package prints "changed 1 package", exits
0 and installs nothing, not with `--force`, and not after clearing the empty
scope directory it leaves behind.  A prefixed install is self-contained,
repeatable, and needs no root.

Hydration is idempotent and non-fatal.  The installer checks it with a
`bd count` health probe and, on success, writes a marker inside
`.beads/embeddeddolt/`; the shim gates on that marker.  Neither uses directory
presence, because `bd init` creates the directory before importing and
`bd prime` creates an empty database as a side effect, and either would latch
"ready" and leave every query silently returning nothing.  A database that is
healthy but has no marker (a developer's own, or one from before the marker
existed) is **adopted**, never rebuilt: it can hold issues that were never
exported to `.beads/issues.jsonl`, and re-importing would discard them.  A
populated database that the probe cannot read at all is **refused** rather
than rebuilt, with a message telling you to look at it, because "`bd count`
failed" and "there are no issues in here" are not the same statement.  Only an
absent, empty or provably issue-free database is cleared and re-imported.  A
single lock (`flock -w`, bounded by `BEADS_LOCK_WAIT`) covers both the install
and the hydration, so two concurrent first commands wait for one another
instead of racing or failing.  Any failure warns and lets the session continue
without bd.

A hydration killed by SIGKILL or a hard container teardown can leave
`.beads/config.yaml`, `.beads/.gitignore` and `.gitignore` modified (no trap
survives SIGKILL), and the next run reads that as a developer's edit and will
not auto-restore it; clear it with
`git checkout -- .beads/config.yaml .beads/.gitignore .gitignore`.

To persist newly filed issues, append their `bd export --include-memories`
lines to `.beads/issues.jsonl` and commit.  Do **not** overwrite the file
wholesale (bd 1.1.0 re-serializes the `dependencies` field on unrelated
issues, and a plain `bd export` drops the persisted memories).  `bd dolt push`
in the session-completion steps above is a no-op here: there is no Dolt remote
in web containers, so committing `issues.jsonl` is the sync.


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
./build_catch/CatchTestRunner "~[slow]"         # everything except [slow]
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
`src/Backends/SVDSBTL/`, `src/Region/`, **or `dev/fluids/` and
`dev/mixtures/`**, run the **umbrella** `[SBTL]` tag locally — NOT just
`[SVDSBTL]`, and NOT just `~[slow]`.  The SVD tables are sampled from the
fluid data, so changing a fluid silently invalidates its cached table; the
tests that would catch it are tagged `[slow]`.  The SBTL adapter layer
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
