# CoolProp CI scripts and contributor tooling

This directory holds CI helper scripts and the central code-quality workflow
doc for contributors. Each section below maps to a CI workflow under
`.github/workflows/`.

If you're new, the **fast path** is:

```bash
cmake -G Ninja -B build -S .                    # one-time configure
pip install pre-commit && pre-commit install    # one-time hook setup
cmake --build build --target format-check       # check formatting any time

# Local pre-push gate (CoolProp-6r6) — mirrors what CI runs, catches
# clang-format / cppcheck / semgrep / test failures before they hit a PR.
# Install once by symlinking the hook (it's NOT auto-installed since
# .git/hooks/ isn't tracked):
ln -s ../../dev/ci/pre-push.sample .git/hooks/pre-push   # one-time
./dev/ci/preflight.sh                                    # run any time
```

Then commit normally — the pre-commit hook will block formatting violations
on staged C/C++ files, and the pre-push hook will run the full preflight
(clang-format + build + tests + cppcheck + clang-tidy + semgrep) before
any `git push` succeeds.

---

## preflight.sh — local pre-push gate

`dev/ci/preflight.sh` is a single script that runs the same checks CI
runs against the diff between HEAD and the upstream branch.  Designed
to be invoked from a pre-push git hook so a passing preflight strongly
predicts a green CI.

| Check | What it does | Skip flag |
|---|---|---|
| clang-format | uvx clang-format (version pinned from `.pre-commit-config.yaml`) dry-run on changed `.cpp` / `.h` files | `--skip=clang-format` |
| build | cmake builds `CatchTestRunner` in `build_catch/` (auto-configures on first run) | `--skip=build` |
| tests | Catch2 runner with auto-selected tag scope — `[SBTL]`, `[SVDSBTL]`, etc. picked from the changed paths | `--skip=tests` |
| cppcheck | `--enable=warning` (real-bug-class) on changed files, `--language=c++ --std=c++17` to handle headers | `--skip=cppcheck` |
| clang-tidy | diff-only via existing `run-clang-tidy-staged.sh`, requires `build_catch/compile_commands.json` | `--skip=clang-tidy` |
| semgrep | `p/security-audit` + local `.semgrep/` rules (uvx-resolved, Python 3.12 pinned) | `--skip=semgrep` |
| incomp-sanity | `dev/incompressible_liquids/test_json_sanity.py` on the committed incompressible JSON — rejects optimizer starting guesses, all-zero fits, non-finite/boolean coefficients and cleared vital properties. Only runs when `dev/incompressible_liquids/` is in the diff; falls back to calling the `test_*` functions directly when pytest is unavailable | `--skip=incomp-sanity` |

Invocation:

```bash
./dev/ci/preflight.sh                          # check vs origin/master
./dev/ci/preflight.sh --base=HEAD~1            # check vs an earlier ref
./dev/ci/preflight.sh --skip=cppcheck,semgrep  # subset
```

Tools missing locally (`semgrep`, `clang-tidy`) are *gracefully skipped*
rather than blocking — but the skip count is reported in the summary so
agents see what's actually being checked.

### Custom semgrep rules

Rules under `.semgrep/` are loaded in addition to the public registry.
Current local rules:

- `cpp-fopen-without-restricted-permissions` — flags `std::fopen("...", "w*")`
  patterns without an obvious permission restriction.  Suppress with
  `// nosemgrep: cpp-fopen-without-restricted-permissions` on the line
  when the surrounding code restricts via `std::filesystem::permissions`
  after the close.  See PR #2947 for the originating incident.

Add new rules as the CI surfaces new CodeQL-class findings that didn't
get caught locally.

---

## bootstrap-beads.sh — bd in ephemeral containers

`dev/ci/bootstrap-beads.sh` is the `SessionStart` hook (see
`.claude/settings.json`) that makes `bd` (the beads issue tracker) usable in
Claude Code web/CI containers, which start from a fresh clone with no `bd`
binary and no Dolt database. It installs `bd` via `go install` (the
`curl|bash` installer's GitHub-release download is blocked by the agent
proxy; the Go module proxy is allowed), rehydrates the embedded DB from the
committed `.beads/issues.jsonl`, and then runs `bd prime` itself — see
CLAUDE.md's "`bd` in ephemeral (Claude Code web/CI) containers" section for
the full story (gating, idempotency, and why it's the only `SessionStart`
hook for beads).

It's opt-in: the install+hydrate steps only run when `BEADS_BOOTSTRAP=1` is
set in an environment's persistent config (the cold-start cost is a
per-environment choice, not auto-detected from `CI`/`CLAUDE_CODE_REMOTE`).
Without it, the hook just primes an already-installed `bd`, if any.

---

## clang-format

`.clang-format` at repo root is the source of truth for formatting rules.
Three runnable paths, all pinned to clang-format 18.1.x:

| Path | Scope | Modifies files? | Use when |
|---|---|---|---|
| `cmake --build build --target format-check` | whole tree (src/, include/) | no (dry-run) | spot-check before committing |
| `cmake --build build --target format` | whole tree | yes (`-i`) | apply repo-wide formatting |
| `pre-commit run` | staged files only | no (dry-run) | before each commit (auto via hook) |
| `./dev/ci/clang-format.sh HEAD <base>` | files changed vs `<base>` | yes (`-i`) | format a PR-sized diff against a ref |

CI uses `dev/ci/clang-format.sh` against the PR base SHA, so any of the
above will mirror what CI does for you.

### pre-commit framework

The `.pre-commit-config.yaml` at repo root pins clang-format via the
`pre-commit/mirrors-clang-format` upstream. One-time setup:

```bash
pip install pre-commit
pre-commit install
```

After that, every `git commit` runs clang-format on staged C/C++ files.
Manual runs:

```bash
pre-commit run                  # staged files
pre-commit run --all-files      # whole repo (slow)
```

To skip the hook for a single commit (e.g. WIP):

```bash
git commit --no-verify
```

To uninstall fully, see the header comment in `.pre-commit-config.yaml`.

### Excluded paths

`.pre-commit-config.yaml` and the CMake `format`/`format-check` targets
both skip `externals/` (vendored third-party). The CMake targets also
skip auto-generated `include/*_JSON*.h`, `include/gitrevision.h`, and
`include/miniz.h` (vendored). If you regenerate the JSON-as-string-
literal headers, format the *generator template*, not the output.

---

## compile_commands.json

`clang-tidy`, `include-what-you-use`, and most editor LSP integrations need a
JSON compilation database. CoolProp configures `CMAKE_EXPORT_COMPILE_COMMANDS`
ON unconditionally in `CMakeLists.txt`, so any cmake configure produces it:

```bash
cmake -G Ninja -B build -S .            # Ninja or Makefile generator
ls build/compile_commands.json          # 90+ entries covering src/
```

### macOS contributors

The Homebrew `llvm@18` and `llvm` packages ship `clang-tidy`, but they don't
know about the Xcode SDK that Apple's `/usr/bin/c++` links against. If
`clang-tidy` reports `'iterator' file not found` or `__builtin_clzg`-style
errors against system headers, point it at the Xcode sysroot:

```bash
SDK=$(xcrun --show-sdk-path)
clang-tidy -p build --extra-arg=--sysroot=$SDK src/CPstrings.cpp
```

clang-tidy 19+ is recommended on macOS — earlier versions don't recognize new
libc++ builtins (`__builtin_clzg`, `__builtin_ctzg`) that Apple's libc++
headers use.

CI runs on Ubuntu where this issue doesn't apply.

### Running clang-tidy locally

Once `build/compile_commands.json` exists, the manual-stage pre-commit
hook runs clang-tidy with the strict repo `.clang-tidy` config:

```bash
pre-commit run --hook-stage manual --all-files clang-tidy
# or per file:
pre-commit run --hook-stage manual --files src/CPstrings.cpp clang-tidy
```

It's manual-stage because clang-tidy on the CoolProp codebase is too slow
to run on every commit (single-file runs can be tens of seconds with the
full check set). Invoke it before pushing instead.

The hook auto-skips with a warning when clang-tidy isn't installed or no
`build/` directory exists. Override the build dir via
`COOLPROP_BUILD_DIR=<path>` if you have multiple build trees.

The CI counterpart (the `clang-tidy` job in
`.github/workflows/dev_checks.yml`) runs `clang-tidy-diff.py` on
PR-touched lines only; the local hook runs clang-tidy on whole files,
so local output is a strict superset of what CI surfaces.

---

## Other CI tooling (warning-only)

All of these now live as parallel jobs in a single
`.github/workflows/dev_checks.yml` workflow (CoolProp-rog consolidated
the previously-separate `dev_*.yml` files). They produce artifacts on
every PR but never fail the build — they exist to surface signal, not
to gate. The exception is `asan`, which does fail.

- **cppcheck** (`cppcheck` job in `dev_checks.yml`) — uploads a colorized
  cppcheck report. Run locally with `cppcheck --std=c++17 ./src` after
  `apt install cppcheck`.
- **clang-tidy diff** (`clang-tidy` job in `dev_checks.yml`) — runs
  `clang-tidy-diff.py` on PR-touched lines and uploads
  `clang-tidy-diff.log`. Empirical noise survey on representative `src/`
  files showed the strict `.clang-tidy` config produces a high cascade
  of `misc-include-cleaner` / `misc-const-correctness` /
  `readability-isolate-declaration` warnings before any meaningful
  bug-finding signal — informational by design, no plans to gate.
- **CodeQL** (`codeql` job in `dev_checks.yml`) — runs the
  `security-and-quality` query suite on every PR. Findings appear in the
  repo's Security tab.
- **Coverity** (`.github/workflows/dev_coverity.yml`) — schedule-only
  (twice weekly) due to free-tier quota. Lives in its own workflow file
  because its trigger model is fundamentally different from the
  consolidated dev_checks jobs. Uploads `coverity-defects.json`
  (machine-readable) for AI-agent consumption; see also the Coverity
  Scan web UI for the curated view.
- **IWYU** (`iwyu` job in `dev_checks.yml`) — runs include-what-you-use
  via the `COOLPROP_IWYU` CMake opt-in and uploads `iwyu.log`.
- **AddressSanitizer** (`asan` job in `dev_checks.yml`) — the Catch test
  suite under ASan on every PR, **excluding `[slow]`**. This one *does* fail
  builds, since memory bugs are real bugs. See
  [ASan scope and budget](#asan-scope-and-budget) for why `[slow]` is
  excluded and what to do when the job times out.

## ASan scope and budget

The `asan` job is the longest-running job in CI, so its scope is a
deliberate trade-off rather than an accident. Two things to know before
changing it.

### It runs `~[slow]`, not everything

`./CatchTestRunner "~[slow]"` excludes the 58 `[slow]`-tagged cases, 49 of
which are in the SVDSBTL suite. Those build SVD surfaces at **production
resolution** (`SurfaceSpec` defaults `NT=200`, `NR=800`, `rank=20`), and the
dense BDCSVD over a 200x800 grid is the bulk of the job's wall time.

Measured on one machine, non-ASan, with a cold surface cache (no
`*.svd.bin.z` present, as in CI):

| Scope | Cases | Wall time |
| --- | --- | --- |
| `[slow]` | 58 (6 skipped, no REFPROP locally) | 887 s |
| `~[slow]` | 410 | >16 min |

So `[slow]` is roughly 40-45% of the local total. That is a *lower* bound on
the CI saving: the 6 cases that skip locally for lack of REFPROP include the
5 `SVDSBTL&REFPROP` production-resolution builds, which do run in CI. Note
these are non-ASan figures; ASan scales everything up but the ratio is what
matters here.

Note also that excluding `[slow]` does **not** remove production-resolution
surface builds entirely: 17 non-`[slow]` `TEST_CASE`s construct SVDSBTL
surfaces with no `grid` option, so they build at 200x800 too. Surfaces are cached per
`(fluid, source, input_pair, options)` key, so that cost is paid once rather
than per test -- which is why the saving comes from the excluded tests'
repeat work and their REFPROP-sourced surfaces, not from eliminating dense
SVD altogether.

**Why not just shrink the grid under ASan?** The resolution is load-bearing
for those tests. They deliberately pass no `grid` option so that they run at
production resolution, then assert agreement with a truth source (HEOS /
REFPROP / IF97) at about `1e-3` relative. `SurfaceSpec` records that the
defaults deliver ~`7e-5` fractional error for Water, so there is only ~14x
headroom. Dropping to `NT=40, NR=80, rank=10` — 50x fewer grid points and
half the SVD rank — would breach that and force the tolerances open. The
result would be a real accuracy gate quietly becoming a vacuous one *in the
ASan build only*, which is the worst place to lose signal.

**Why the ASan coverage cost is small.** Because production resolution still
runs. As noted above, 17 non-`[slow]` `TEST_CASE`s build surfaces at the
200x800/rank-20 defaults, so full-size dense BDCSVD, its temporaries and its
allocation pattern are all still exercised under ASan. That is what makes the
exclusion cheap for the classes ASan detects -- heap overflow, use-after-free,
leaks, ODR.

**Do not restate this as "the small-grid tests keep that path covered."** They
do not run. All 8 of the `NT=40, NR=80, rank=10` grid specs in
`src/Tests/CoolProp-Tests-SVDSBTL.cpp` sit inside `[slow]` `TEST_CASE`s, so
`~[slow]` excludes every one of them; the only non-`[slow]` references to a
40x80 grid are schema/round-trip tests in `CoolProp-Tests-SVDSBTLOptions.cpp`,
which build no surface. An earlier version of this section had it backwards --
it counted the 8 option-string literals and assumed their tags without
checking them.

What is actually given up is the excluded tests' repeat work and their
REFPROP-sourced surfaces (the 6 `SVDSBTL&REFPROP` truth-comparison cases), not
the large-matrix code path. Note the consequence: because 200x800 still runs, a
size-dependent bug reachable only at production extent is still reachable here,
so there is no "only visible at full resolution" residual risk to trade away.
If the trade ever looks wrong,
the answer is more runner minutes, not a resolution cut (see above for why
lowering the grid breaks assertions on both sides of the filter).

If you add a test that is the *only* coverage of a memory-safety-relevant
path, do not tag it `[slow]`, or ASan will not see it.

One consequence worth knowing before anyone proposes lowering the grid for
the *surviving* tests instead: some of them also assert accuracy against
HEOS (`SVDSBTL DT two-phase dome` at `1e-3`, `DT fast_evaluate` at `1e-2`),
so a global resolution cut is not safe even with `[slow]` excluded.

### A timeout is a step failure, not a cancellation

The test step is wrapped in `timeout` with a 70-minute budget, inside the
job's 90-minute `timeout-minutes`. That ordering is intentional:

- A job that exceeds `timeout-minutes` is reported by GitHub with conclusion
  **`cancelled`**, which is indistinguishable from a concurrency-supersede or
  from `cancel_merged_pr_runs.yml` cancelling a merged PR's runs. That is
  exactly the ambiguity the job cap was introduced to remove (CoolProp-xuu),
  so leaving the cap as the only bound reintroduces it.
- Bounding the *step* makes an overrun exit `124` (or `137` after
  `--kill-after`), which fails the step with an explicit `::error::`
  annotation saying it was a timeout rather than a memory error.

So when this job goes red, read the annotation first:

| Symptom | Meaning |
| --- | --- |
| `::error::ASan run exceeded its budget` (exit 124) | Timeout. Investigate what got slower; do not just raise the budget. |
| `::error::ASan run was killed (SIGKILL)` (exit 137) | Ambiguous: either the budget elapsed and SIGINT was ignored, **or** the kernel OOM-killed it. ASan roughly doubles memory, so OOM is realistic on a 2-vCPU runner — check for an ASan OOM report before assuming a timeout. |
| ASan report in the log (`ERROR: AddressSanitizer: ...`) | A real memory bug. Fix it. |
| Conclusion `cancelled`, no annotation | Most likely a superseded push or a merged PR (`cancel_merged_pr_runs.yml`). But **check which step was running first**: only the test step is wrapped in `timeout`, so if REFPROP setup or the compile overran the job cap you also get an unannotated `cancelled`. |

### Measured wall time (the first clean `~[slow]` run)

PR #3299, run `30364245551`, job `90291277987` — the first ASan job to finish
rather than hit the cap. Conclusion `success`:

| Step | Wall time |
| --- | --- |
| REFPROP build | 36 s |
| cmake config | 28 s |
| compile (`-j$(nproc)`) | 8.8 min |
| **test run (`~[slow]`)** | **50.7 min (3041 s)** |
| job total | 60.6 min |

That settles the sizing, and both current values turn out to be about right:

- `budget=4200` is **1.38x** the measured 3041 s, leaving 19 min unused —
  essentially the 1.3x target, so leave it. Runners are shared and vary by
  tens of percent between runs; a tighter budget would start flagging normal
  variance as a timeout, which is the failure mode this gate exists to avoid.
- `timeout-minutes: 90` is 1.49x the 60.6 min job total. Also fine to leave.

Both runs of the same scope, before and after `verbosity=1` was dropped:

| Run | `verbosity=1` | Job total |
| --- | --- | --- |
| `30364245551` (`0d2af4f9`) | set | 60.6 min |
| `30371990208` (`3a0f63a5`) | not set | 60.2 min |

So the sizing above holds for the shipped config. Read the 0.4 min gap as
nothing more than that: one observed difference between two runs that also
differ in one setting. It is not a variance bound in either direction — two
observations cannot establish one — and it is not enough to infer a
performance effect. If you ever need the real spread, collect repeated runs
at a fixed config; do not tighten the budget off a pair of numbers.

For contrast, the full-scope (`[slow]` included) predecessor on PR #3294 did
not finish at all: it ran 09:33:03 -> 11:03:20, hit the 90 min cap, and its
log ends `Terminate orphan process: pid (4681) (CatchTestRunner)`.

### Why `verbosity=1` is not set

The job used to run with `ASAN_OPTIONS=verbosity=1:...`. That was removed
because it buys nothing and produces an enormous amount of log noise.

**It is not a performance fix, despite appearances.** Removing it did not
measurably speed the job up: 60.6 min with it, 60.2 min without. If you are
looking for the wall-time lever, it is the `~[slow]` scope, not this. The
volume figures below are why the noise is worth removing on its own terms —
they are not evidence of a time saving, and the extrapolation that suggested
one turned out not to hold (the 4365 lines/s rate is clearly not sustained
across the whole run, or the difference would have shown up).

At `verbosity>=1` ASan reports its own container-annotation bookkeeping, which
for this Eigen- and `std::vector`-heavy suite is a sustained stream of

```text
==4681==poisoning: 0x7f7bed013040 800
==4681==unpoisoning: 0x7f7bed0141c0 400
```

In PR #3294's ASan log that ran at **4365 lines/s** (measured across a 0.68 s
window; ~38 bytes per line). Extrapolated over a 70-minute test window that
is on the order of 18 M lines and ~690 MB of stderr — treat that total as an
upper bound, since the rate was sampled in one window rather than averaged
over the whole run. Every one of those lines is a formatted `write()` into a
pipe that the Actions runner reads, timestamps and uploads.

It is not a diagnostic setting. ASan's error reports, leak reports and
ODR-violation reports print unconditionally, at any verbosity; `verbosity`
only adds ASan's internal chatter. Confirmed directly: a trivial
`std::vector` program emits **0** lines at the default verbosity and **236**
at `verbosity=1`.

`detect_odr_violation=1` is kept — that one is doing real work (see the ODR
note in `dev_checks.yml` for why it is level 1 and not 2).

### The exit code does propagate test failures (two ASan options exist)

There are two ASan switches in `CMakeLists.txt` and only one is used by CI:

- **`COOLPROP_ASAN`** (what `dev_checks.yml` passes, with
  `-DCMAKE_BUILD_TYPE=Asan`) adds an `Asan` build type at `-O2 -g -DNDEBUG`.
  `CatchTestRunner` links `Catch2::Catch2WithMain`, so a failed assertion
  makes the process exit non-zero. Specifically it exits **42** --
  Catch2's `TestFailureExitCode`, a fixed sentinel, *not* the number of
  failures (`catch_session.cpp`: `if (totals.assertions.failed) { return
  TestFailureExitCode; }`). Verified by injecting one bad `CHECK`: 5 failed
  assertions still yields 42. Catch2's other codes are 1 unspecified error,
  2 no tests run, 3 unmatched test spec, 4 all skipped, 5 invalid spec.
  This matters for the `124`/`137` handling below: since 42 is fixed, a test
  failure can never be mistaken for a `timeout` exit code.
- **`COOLPROP_CLANG_ADDRESS_SANITIZER`** is a separate, older buildbot-era
  switch. It appends `src/Tests/catch_always_return_success.cxx`, whose
  `main` runs the session and then returns `EXIT_SUCCESS` **unconditionally**
  — deliberately, so that buildbot saw ASan's exit code rather than the Catch
  failure count.

So the `exit $rc` gate in the workflow is meaningful, but only because CI
uses the first option. If anyone ever switches the job to
`COOLPROP_CLANG_ADDRESS_SANITIZER`, ordinary test failures would stop failing
the job while ASan's own aborts still would — a silent narrowing of the gate.
Don't.

### Reproducing locally

```bash
cmake -B build_asan -S . -DCMAKE_BUILD_TYPE=Asan -DCOOLPROP_ASAN=ON \
      -DCOOLPROP_CATCH_MODULE=ON -DCMAKE_CXX_COMPILER=clang++
cmake --build build_asan --target CatchTestRunner -j$(nproc)
ASAN_OPTIONS=detect_odr_violation=1 \
ASAN_SYMBOLIZER_PATH=$(command -v llvm-symbolizer) \
  timeout --signal=INT --kill-after=60 4200 \
  ./build_asan/CatchTestRunner "~[slow]" --benchmark-samples 1
```

That mirrors what CI runs. Drop the `timeout` wrapper and
`--benchmark-samples 1` if you only want a functional check and do not care
about reproducing the budget or the benchmark sampling.

Note `detect_stack_use_after_return` is deliberately **not** enabled, and
`detect_odr_violation` is lowered to `1`; both have specific justifications
in the job's comments in `dev_checks.yml`.

## `git blame` and the one-shot reformat

The repo's `.git-blame-ignore-revs` lists SHAs of pure-formatting / mechanical
commits (the whole-repo clang-format pass in PR #2803, plus future
`clang-tidy --fix` passes). To make local `git blame` skip them, opt in
once per clone:

```bash
git config blame.ignoreRevsFile .git-blame-ignore-revs
```

GitHub's blame view honours this file automatically — no setup needed in the
web UI.
