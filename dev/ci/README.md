# CoolProp CI scripts and contributor tooling

This directory holds CI helper scripts and the central code-quality workflow
doc for contributors. Each section below maps to a CI workflow under
`.github/workflows/`.

If you're new, the **fast path** is:

```bash
cmake -G Ninja -B build -S .                    # one-time configure
pip install pre-commit && pre-commit install    # one-time hook setup
cmake --build build --target format-check       # check formatting any time
```

Then commit normally — the pre-commit hook will block formatting violations
on staged C/C++ files.

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

---

## Other CI tooling (warning-only)

These workflows produce artifacts on every PR but never fail the build —
they exist to surface signal, not to gate.

- **cppcheck** (`.github/workflows/dev_cppcheck.yml`) — uploads a colorized
  cppcheck report. Run locally with `cppcheck --std=c++17 ./src` after
  `apt install cppcheck`.
- **clang-tidy diff** (`.github/workflows/dev_clangtidy.yml`) — runs
  `clang-tidy-diff.py` on PR-touched lines and uploads
  `clang-tidy-diff.log`. Empirical noise survey on representative `src/`
  files showed the strict `.clang-tidy` config produces a high cascade
  of `misc-include-cleaner` / `misc-const-correctness` /
  `readability-isolate-declaration` warnings before any meaningful
  bug-finding signal — informational by design, no plans to gate.
- **CodeQL** (`.github/workflows/dev_codeql.yml`) — runs the
  `security-and-quality` query suite on every PR. Findings appear in the
  repo's Security tab.
- **Coverity** (`.github/workflows/dev_coverity.yml`) — schedule-only
  (twice weekly) due to free-tier quota. The workflow uploads
  `coverity-defects.json` (machine-readable) for AI-agent consumption;
  see also the Coverity Scan web UI for the curated view.
- **IWYU** (`.github/workflows/dev_iwyu.yml`) — runs include-what-you-use
  via the `COOLPROP_IWYU` CMake opt-in and uploads `iwyu.log`.
- **AddressSanitizer** (`.github/workflows/dev_asan.yml`) — full Catch test
  suite under ASan on every PR. This one *does* fail builds, since memory
  bugs are real bugs.

## Future tooling

These will land as the `CoolProp-2uw` epic progresses; cross-references
will be added here as each lands:

- `.git-blame-ignore-revs` for the one-shot reformat (Tier 4.1, `CoolProp-2uw.12`)
