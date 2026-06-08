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
- **AddressSanitizer** (`asan` job in `dev_checks.yml`) — full Catch test
  suite under ASan on every PR. This one *does* fail builds, since memory
  bugs are real bugs.

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
