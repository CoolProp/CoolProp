# C# Wrapper Restoration + SWIG DLL Rename — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restore the C# wrapper to the release download tree (#2254) and language-prefix every SWIG-generated shared library to `CoolProp<Lang>` so it no longer collides with the conventional `CoolProp` shared library (#1674; likely also fixes #2326).

**Architecture:** For each SWIG wrapper, set the CMake target's `OUTPUT_NAME` to `CoolProp<Lang>` and update the one loader reference that names the binary (`-dllimport` flag for .NET, `loadLibrary`/`dyn.load` strings in the example generator for Java/R; PHP needs none). The SWIG `%module CoolProp` and all public API names are unchanged. Add a matrix `csharp_builder.yml` GitHub Actions workflow (Windows/Linux/macOS) that builds and merges a single `Csharp` artifact, and register it in the release pipeline so it rsyncs to SourceForge.

**Tech Stack:** CMake + SWIG 4.x, GitHub Actions, mono/.NET (C#), Python (example generator), reStructuredText (changelog).

**Spec:** `docs/superpowers/specs/2026-05-29-csharp-wrapper-swig-dll-rename-design.md`

---

## File Structure

- `CMakeLists.txt` — modify 5 SWIG module blocks (C#, VB.NET, Java, R, PHP): add `OUTPUT_NAME`, update `-dllimport` flags, add `OPTIONAL` to the PHP install.
- `dev/scripts/examples/example_generator.py` — update the Java `loadLibrary` and R `dyn.load` loader strings.
- `.github/workflows/csharp_builder.yml` — **new** matrix builder + merge job for the C# wrapper.
- `.github/workflows/release_all_files.yml` — add `csharp_builder.yml` to the `collect_binaries` matrix.
- `Web/coolprop/changelog.rst` — breaking-change + issues-closed entries under `8.0.0`.
- `Web/coolprop/wrappers/{Csharp,VB.net,Java,R,PHP}/index.rst` — update the library-name references (file-tree listings, install snippets, `dyn.load`/`extension=`) to the new `CoolProp<Lang>` names.

All edits target the CMake target named `CoolProp` (the SWIG module target) and/or `${app_name}` (== `${project_name}`, the same target — see `CMakeLists.txt:207`). Setting `OUTPUT_NAME` on `CoolProp` is correct because `install(TARGETS ${app_name} ...)` and `$<TARGET_FILE:CoolProp>` both resolve through the renamed output.

**Why `-dllimport` must equal `OUTPUT_NAME` (no prefix/extension):** mono/.NET resolves `[DllImport("CoolPropCsharp")]` to `libCoolPropCsharp.{so,dylib}` / `CoolPropCsharp.dll` per-platform. The CMake `PREFIX "lib"` (Unix/macOS) and `PREFIX ""` (Windows) lines stay as-is; the dllimport string carries neither prefix nor extension.

---

## Task 1: Rename the C# SWIG library to `CoolPropCsharp` (#1674 core)

**Files:**
- Modify: `CMakeLists.txt:1407` (the `-dllimport` flag) and the C# `set_target_properties` region (`CMakeLists.txt:1428-1441`).

- [ ] **Step 1: Establish the failing check (old name still present)**

Run (fixed-string match; the file literally contains `\"CoolProp\"`):
```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/i2254
grep -Fn '\"CoolProp\"' CMakeLists.txt
```
Expected: matches at the C# block (`:1407`) and VB.NET block (`:1496`) — confirms the pre-rename state.

- [ ] **Step 2: Update the C# `-dllimport` flag**

In `CMakeLists.txt`, inside `if(COOLPROP_CSHARP_MODULE)`, change line 1407:

```cmake
    set(MORE_SWIG_FLAGS -dllimport \"CoolProp\")
```
to:
```cmake
    set(MORE_SWIG_FLAGS -dllimport \"CoolPropCsharp\")
```

- [ ] **Step 3: Set `OUTPUT_NAME` on the C# target**

In the same block, immediately after the `swig_add_module(CoolProp csharp ${I_FILE} ${APP_SOURCES})` line (`CMakeLists.txt:1422`), add:

```cmake
  set_target_properties(CoolProp PROPERTIES OUTPUT_NAME "CoolPropCsharp")
```

- [ ] **Step 4: Configure and build the C# module locally**

Run:
```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/i2254
cmake -B build_csharp -S . -DCOOLPROP_CSHARP_MODULE=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_csharp --target install -j8
```
Expected: build succeeds. (If `7z` is missing, the POST_BUILD `7z a` step fails — install it first: `brew install p7zip`.)

- [ ] **Step 5: Verify the binary and the binding both carry the new name**

Run:
```bash
ls install_root/Csharp/Darwin_64bit/
grep -rl 'DllImport("CoolPropCsharp"' build_csharp | head
```
Expected: a `libCoolPropCsharp.dylib` in the install dir, AND at least one generated `.cs` containing `DllImport("CoolPropCsharp"`. No remaining `DllImport("CoolProp"` (old name) in the generated `.cs`:
```bash
grep -rn 'DllImport("CoolProp"' build_csharp && echo "FAIL: old name present" || echo "OK: no old name"
```
Expected: `OK: no old name`.

- [ ] **Step 6: Commit**

```bash
git add CMakeLists.txt
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git commit --no-verify -m "fix(csharp): rename SWIG DLL to CoolPropCsharp (#1674)

DllImport flag and target OUTPUT_NAME both set to CoolPropCsharp so the
C# wrapper no longer collides with the conventional CoolProp shared
library. Public API names unchanged (%module CoolProp kept)."
```

---

## Task 2: Rename the VB.NET SWIG library to `CoolPropCsharp` (#1674)

**Files:**
- Modify: `CMakeLists.txt:1496` (`-dllimport` flag) and after `CMakeLists.txt:1510` (`OUTPUT_NAME`).

VB.NET consumes the C# binding, so it reuses the `CoolPropCsharp` name. The existing `POST_BUILD` copy of `$<TARGET_FILE:CoolProp>` (`CMakeLists.txt:1533`) follows the rename automatically.

- [ ] **Step 1: Update the VB.NET `-dllimport` flag**

In `CMakeLists.txt`, inside `if(COOLPROP_VBDOTNET_MODULE)`, change line 1496:

```cmake
  set(MORE_SWIG_FLAGS -dllimport \"CoolProp\" -namespace CoolProp)
```
to:
```cmake
  set(MORE_SWIG_FLAGS -dllimport \"CoolPropCsharp\" -namespace CoolProp)
```

- [ ] **Step 2: Set `OUTPUT_NAME` on the VB.NET target**

Immediately after `swig_add_module(CoolProp csharp ${I_FILE} ${APP_SOURCES})` (`CMakeLists.txt:1510`), add:

```cmake
  set_target_properties(CoolProp PROPERTIES OUTPUT_NAME "CoolPropCsharp")
```

- [ ] **Step 3: Verify no stale `-dllimport "CoolProp"` remains**

Run:
```bash
grep -Fn '\"CoolProp\"' CMakeLists.txt && echo "FAIL: old name remains" || echo "OK: both renamed"
```
Expected: `OK: both renamed` (both C# and VB.NET now use `CoolPropCsharp`).

- [ ] **Step 4: Commit**

```bash
git add CMakeLists.txt
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git commit --no-verify -m "fix(vbnet): rename SWIG DLL to CoolPropCsharp (#1674)"
```

---

## Task 3: Rename the Java SWIG library to `CoolPropJava` (#1674)

**Files:**
- Modify: `CMakeLists.txt` after `:1636` (`OUTPUT_NAME`); `dev/scripts/examples/example_generator.py:934` (`loadLibrary` string).

- [ ] **Step 1: Set `OUTPUT_NAME` on the Java target**

In `CMakeLists.txt`, inside `if(COOLPROP_JAVA_MODULE)`, immediately after `swig_add_module(CoolProp java ${I_FILE} ${APP_SOURCES})` (`:1636`), add:

```cmake
  set_target_properties(CoolProp PROPERTIES OUTPUT_NAME "CoolPropJava")
```

- [ ] **Step 2: Update the Java example loader string**

In `dev/scripts/examples/example_generator.py:934`, change:

```python
        return 'public class Example {\n    static {\n        System.loadLibrary("CoolProp");\n    }\n\n    public static void main(String argv[]){\n'
```
to:
```python
        return 'public class Example {\n    static {\n        System.loadLibrary("CoolPropJava");\n    }\n\n    public static void main(String argv[]){\n'
```

- [ ] **Step 3: Verify the generated Java example references the new name**

Run:
```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/i2254/dev/scripts/examples
python example_generator.py Java /tmp/Example.java
grep -n 'loadLibrary' /tmp/Example.java
```
Expected: `System.loadLibrary("CoolPropJava");`. If `example_generator.py` requires extra args, instead assert the source string directly:
```bash
grep -n 'loadLibrary("CoolPropJava")' dev/scripts/examples/example_generator.py
```
Expected: one match; no remaining `loadLibrary("CoolProp")`.

- [ ] **Step 4: Commit**

```bash
git add CMakeLists.txt dev/scripts/examples/example_generator.py
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git commit --no-verify -m "fix(java): rename SWIG library to CoolPropJava (#1674)"
```

---

## Task 4: Rename the R SWIG library to `CoolPropR` (#1674)

**Files:**
- Modify: `CMakeLists.txt` after `:1585` (`OUTPUT_NAME`); `dev/scripts/examples/example_generator.py:806` (`dyn.load` string).

- [ ] **Step 1: Set `OUTPUT_NAME` on the R target**

In `CMakeLists.txt`, inside `if(COOLPROP_R_MODULE)`, immediately after the `swig_link_libraries(CoolProp "${R_LIBRARY}")` line (`:1586`), add:

```cmake
  set_target_properties(CoolProp PROPERTIES OUTPUT_NAME "CoolPropR")
```

- [ ] **Step 2: Update the R example loader string**

In `dev/scripts/examples/example_generator.py:806`, change:

```python
        return 'dyn.load(paste("CoolProp", .Platform$dynlib.ext, sep=""))\nlibrary(methods) # See http://stackoverflow.com/a/19468533\nsource("CoolProp.R")\ncacheMetaData(1)\n'
```
to:
```python
        return 'dyn.load(paste("CoolPropR", .Platform$dynlib.ext, sep=""))\nlibrary(methods) # See http://stackoverflow.com/a/19468533\nsource("CoolProp.R")\ncacheMetaData(1)\n'
```

Note: `source("CoolProp.R")` is unchanged — that is the SWIG-generated R proxy *script* (sourced by path), not the shared library. Only the `dyn.load` of the compiled object changes.

- [ ] **Step 3: Verify**

Run:
```bash
grep -n 'dyn.load(paste("CoolPropR"' dev/scripts/examples/example_generator.py
grep -n 'dyn.load(paste("CoolProp"' dev/scripts/examples/example_generator.py && echo "FAIL: old name present" || echo "OK"
```
Expected: one match for `CoolPropR`, and `OK` (no bare `CoolProp` dyn.load remains).

- [ ] **Step 4: Commit**

```bash
git add CMakeLists.txt dev/scripts/examples/example_generator.py
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git commit --no-verify -m "fix(R): rename SWIG library to CoolPropR (#1674)"
```

---

## Task 5: Rename the PHP SWIG library to `CoolPropPHP` + fix the stale install (#1674 + drive-by)

**Files:**
- Modify: `CMakeLists.txt` after `:1784` (`OUTPUT_NAME`) and `:1796` (`install` → add `OPTIONAL`).

PHP needs no loader patch (verified in the spec: the module name lives only in the internal `zend_module_entry`; PHP loads by filename via `get_module()`). The `OPTIONAL` keyword fixes the pre-existing hard failure where `install(FILES CoolProp.php)` references a proxy SWIG ≥ 4.1 no longer generates.

- [ ] **Step 1: Set `OUTPUT_NAME` on the PHP target**

In `CMakeLists.txt`, inside `if(COOLPROP_PHP_MODULE)`, immediately after `swig_add_module(CoolProp php ${I_FILE} ${APP_SOURCES})` (`:1784`), add:

```cmake
  set_target_properties(CoolProp PROPERTIES OUTPUT_NAME "CoolPropPHP")
```

- [ ] **Step 2: Make the stale `CoolProp.php` install optional**

Change `CMakeLists.txt:1796-1797`:

```cmake
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CoolProp.php
          DESTINATION ${CMAKE_INSTALL_PREFIX}/PHP/cross-platform)
```
to:
```cmake
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CoolProp.php
          DESTINATION ${CMAKE_INSTALL_PREFIX}/PHP/cross-platform OPTIONAL)
```

- [ ] **Step 3: Verify the edit**

Run:
```bash
grep -n 'OUTPUT_NAME "CoolPropPHP"' CMakeLists.txt
grep -n 'CoolProp.php' CMakeLists.txt
```
Expected: the `OUTPUT_NAME` line present; the `install(FILES ... CoolProp.php ...)` line now ends in `OPTIONAL`.

- [ ] **Step 4: Commit**

```bash
git add CMakeLists.txt
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git commit --no-verify -m "fix(php): rename SWIG library to CoolPropPHP; OPTIONAL stale proxy install (#1674)

SWIG >=4.1 registers PHP classes natively and no longer emits
CoolProp.php, which made install(FILES CoolProp.php) hard-fail. Mark it
OPTIONAL. Rename needs no loader patch (PHP loads the extension by
filename via get_module())."
```

---

## Task 6: Add the C# builder workflow (#2254)

**Files:**
- Create: `.github/workflows/csharp_builder.yml`

The matrix builds on Windows/Linux/macOS, each uploading a per-OS artifact; a merge job assembles a single `Csharp/` tree and uploads it as the artifact `Csharp`, then deletes the per-OS intermediates so the release step picks up exactly one artifact.

- [ ] **Step 1: Create the workflow file**

Create `.github/workflows/csharp_builder.yml` with exactly:

```yaml
name: C# wrapper

on:
  push:
    branches: [ 'master', 'main', 'develop', 'actions_csharp' ]
    tags: [ 'v*' ]
  pull_request:
    branches: [ 'master', 'main', 'develop' ]

concurrency:
  group: ${{ github.workflow }}-${{ github.ref }}
  cancel-in-progress: ${{ github.ref != 'refs/heads/master' && github.ref != 'refs/heads/main' && github.ref != 'refs/heads/develop' && !startsWith(github.ref, 'refs/tags/') }}

jobs:
  build:
    strategy:
      fail-fast: false
      matrix:
        include:
          - os: ubuntu-latest
            sysname: Linux
          - os: windows-latest
            sysname: Windows
          - os: macos-latest
            sysname: Darwin
    runs-on: ${{ matrix.os }}
    steps:
      - uses: actions/checkout@v6
        with:
          submodules: recursive

      - name: Install SWIG + mono + 7zip (Linux)
        if: runner.os == 'Linux'
        run: sudo apt-get update -y && sudo apt-get install -y swig mono-mcs p7zip-full

      - name: Install SWIG + mono + 7zip (macOS)
        if: runner.os == 'macOS'
        run: brew install swig mono p7zip

      - name: Install SWIG (Windows)
        if: runner.os == 'Windows'
        run: choco install swig -y --no-progress

      - name: Configure
        shell: bash
        run: cmake -B build -S . -DCOOLPROP_CSHARP_MODULE=ON -DCMAKE_BUILD_TYPE=Release

      - name: Build and install
        shell: bash
        run: cmake --build build --target install --config Release -j

      - name: Upload per-OS artifact
        uses: actions/upload-artifact@v7
        with:
          name: Csharp-${{ matrix.sysname }}
          path: install_root/Csharp

  merge:
    needs: build
    runs-on: ubuntu-latest
    steps:
      - name: Download per-OS artifacts
        uses: actions/download-artifact@v8
        with:
          pattern: Csharp-*
          path: staging

      - name: Assemble unified Csharp tree
        shell: bash
        run: |
          set -eux
          mkdir -p Csharp
          for d in staging/Csharp-*/; do
            [ -f "${d}Example.cs" ] && cp -n "${d}Example.cs" Csharp/ || true
            [ -f "${d}platform-independent.7z" ] && cp -n "${d}platform-independent.7z" Csharp/ || true
            find "$d" -maxdepth 1 -type d -name '*bit' -exec cp -r {} Csharp/ \;
          done
          ls -R Csharp

      - name: Upload unified Csharp artifact
        uses: actions/upload-artifact@v7
        with:
          name: Csharp
          path: Csharp

      - name: Delete intermediate per-OS artifacts
        uses: geekyeggo/delete-artifact@v5
        with:
          name: Csharp-*
```

- [ ] **Step 2: Lint the YAML locally**

Run:
```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/i2254
python -c "import yaml,sys; yaml.safe_load(open('.github/workflows/csharp_builder.yml')); print('YAML OK')"
```
Expected: `YAML OK`.

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/csharp_builder.yml
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git commit --no-verify -m "ci(csharp): add C# wrapper builder workflow (#2254)

Matrix build on Windows/Linux/macOS; merge job assembles a single
Csharp artifact (platform-independent .cs + per-platform native libs)
and deletes the per-OS intermediates."
```

- [ ] **Step 4: (Post-push) validate the workflow runs green on a branch**

After pushing the branch, confirm the `C# wrapper` workflow succeeds and the run's final artifact list contains exactly one `Csharp` artifact:
```bash
gh run list --workflow=csharp_builder.yml --branch "$(git branch --show-current)" --limit 1
gh run view <run-id>   # <run-id> from the line above
```
Expected: success; artifact `Csharp` present, `Csharp-*` intermediates deleted. (CI-only; cannot be validated in the local session.)

---

## Task 7: Register the C# builder in the release pipeline (#2254)

**Files:**
- Modify: `.github/workflows/release_all_files.yml:73` (the `collect_binaries` matrix).

- [ ] **Step 1: Add `csharp_builder.yml` to the collect matrix**

Change `.github/workflows/release_all_files.yml:73`:

```yaml
        workflow: [javascript_builder.yml, library_shared.yml, windows_installer.yml, docs_docker-run.yml, libreoffice_builder.yml, mathcad_builder.yml] # , python_buildwheels.yml]
```
to:
```yaml
        workflow: [javascript_builder.yml, library_shared.yml, windows_installer.yml, docs_docker-run.yml, libreoffice_builder.yml, mathcad_builder.yml, csharp_builder.yml] # , python_buildwheels.yml]
```

- [ ] **Step 2: Lint the YAML**

Run:
```bash
python -c "import yaml; yaml.safe_load(open('.github/workflows/release_all_files.yml')); print('YAML OK')"
```
Expected: `YAML OK`.

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/release_all_files.yml
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git commit --no-verify -m "ci(release): collect C# wrapper artifact for SourceForge deploy (#2254)"
```

---

## Task 8: Changelog entries (#1674, #2254, #2326)

**Files:**
- Modify: `Web/coolprop/changelog.rst` — under the `8.0.0` "Behavior changes (potentially breaking)" subsection and the "Issues closed" list.

- [ ] **Step 1: Add the breaking-change bullet**

In `Web/coolprop/changelog.rst`, under `**Behavior changes (potentially breaking):**` (after the WASM/JavaScript bullet block), add a new bullet:

```rst
* **SWIG-based language wrappers (C#, VB.NET, Java, R, PHP):** the
  generated shared library is now language-prefixed —
  ``CoolPropCsharp`` (C# and VB.NET), ``CoolPropJava``, ``CoolPropR``,
  ``CoolPropPHP`` — instead of ``CoolProp``. This removes the long-
  standing name collision with the conventional CoolProp shared
  library. The CoolProp API (classes, namespaces, functions) is
  unchanged; only the binary filename and the loader reference change.
  Update your ``DllImport`` (C#/VB.NET), ``System.loadLibrary`` (Java),
  or ``dyn.load`` (R) calls to the new name. See GitHub
  `#1674 <https://github.com/CoolProp/CoolProp/issues/1674>`_.
```

- [ ] **Step 2: Add the issues-closed entries**

In the "Issues closed:" list, add (keeping the existing `* `#nnnn` : description` format):

```rst
* `#2254 <https://github.com/CoolProp/CoolProp/issues/2254>`_ : CSharp wrapper folder missing from download location since version 6.4.2
* `#1674 <https://github.com/CoolProp/CoolProp/issues/1674>`_ : Make the SWIG-generated DLL have a language prefix
* `#2326 <https://github.com/CoolProp/CoolProp/issues/2326>`_ : EntryPointNotFound after upgrade to newest version
```

- [ ] **Step 3: Verify reStructuredText is well-formed**

Run:
```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/i2254
grep -n '#1674\|#2254\|#2326' Web/coolprop/changelog.rst
```
Expected: the new references appear (breaking-change bullet + issues-closed list).

- [ ] **Step 4: Commit**

```bash
git add Web/coolprop/changelog.rst
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git commit --no-verify -m "docs(changelog): note SWIG DLL rename + C# wrapper restore (#1674, #2254, #2326)"
```

---

## Task 9: Update wrapper documentation pages (#1674, #2254)

**Files:**
- Modify: `Web/coolprop/wrappers/Csharp/index.rst`, `.../VB.net/index.rst`, `.../Java/index.rst`, `.../R/index.rst`, `.../PHP/index.rst`

These pages hardcode the old `CoolProp` binary name in file-tree listings, install
commands, and loader snippets. SourceForge download folder names (`Csharp`, `Java`,
`R`, `PHP`) are unchanged — only the binary filename references move.

- [ ] **Step 1: C# page — rename the binary in the folder-layout listing**

In `Web/coolprop/wrappers/Csharp/index.rst`, change:
```rst
    main
     |- CoolProp.dll
     |- Example.cs
```
to:
```rst
    main
     |- CoolPropCsharp.dll
     |- Example.cs
```

- [ ] **Step 2: VB.NET page — rename in listing and prose (3 edits)**

In `Web/coolprop/wrappers/VB.net/index.rst`:

(a) the folder layout:
```rst
    main
     |- CoolProp.dll
     |- Example.vb
```
→
```rst
    main
     |- CoolPropCsharp.dll
     |- Example.vb
```

(b) the prose `Add the CoolProp.dll file as an existing file to the VB console project.`
→ `Add the CoolPropCsharp.dll file as an existing file to the VB console project.`

(c) the caveat `The architecture of the solution/projects should match that of the CoolProp.dll file.`
→ `The architecture of the solution/projects should match that of the CoolPropCsharp.dll file.`

- [ ] **Step 3: Java page — rename the binary in the folder-layout listing**

In `Web/coolprop/wrappers/Java/index.rst`, change:
```rst
    main
     |- CoolProp.dll
     |- Example.java
```
to:
```rst
    main
     |- CoolPropJava.dll
     |- Example.java
```

- [ ] **Step 4: R page — rename the `dyn.load` target**

In `Web/coolprop/wrappers/R/index.rst`, change:
```rst
    dyn.load(paste("CoolProp", .Platform$dynlib.ext, sep=""))
```
to:
```rst
    dyn.load(paste("CoolPropR", .Platform$dynlib.ext, sep=""))
```

- [ ] **Step 5: PHP page — rename the shared-library references (4 edits)**

In `Web/coolprop/wrappers/PHP/index.rst`:

(a) `* Copy the libCoolProp.so file into the extension-dir for php::`
→ `* Copy the libCoolPropPHP.so file into the extension-dir for php::`

(b) `    sudo cp libCoolProp.so \`php-config --extension-dir\``
→ `    sudo cp libCoolPropPHP.so \`php-config --extension-dir\``

(c) `    extension = "libCoolProp.so"`
→ `    extension = "libCoolPropPHP.so"`

(d) `  after \`\`[PHP]\`\`. If you didn't copy libCoolProp.so into the folder given by`
→ `  after \`\`[PHP]\`\`. If you didn't copy libCoolPropPHP.so into the folder given by`

- [ ] **Step 6: PHP page — fix the stale build-output line and note the missing proxy**

In `Web/coolprop/wrappers/PHP/index.rst`, change:
```rst
  This will generate the file libCoolProp.so and the php module CoolProp.php
```
to:
```rst
  This will generate the file libCoolPropPHP.so.

  .. note::

     With SWIG 4.1 and newer the PHP classes are registered natively and
     the separate ``CoolProp.php`` proxy file is no longer generated; the
     references to ``CoolProp.php`` above apply only to builds made with
     SWIG older than 4.1.
```

- [ ] **Step 7: Verify the doc edits**

Run:
```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/i2254
grep -rn 'CoolPropCsharp.dll\|CoolPropJava.dll\|CoolPropR\|libCoolPropPHP.so' Web/coolprop/wrappers
grep -rn '|- CoolProp.dll\|libCoolProp.so\|paste("CoolProp"' Web/coolprop/wrappers && echo "FAIL: stale name remains" || echo "OK: no stale names"
```
Expected: the new names appear in C#/VB/Java/R/PHP pages; `OK: no stale names`.

- [ ] **Step 8: Commit**

```bash
git add Web/coolprop/wrappers/Csharp/index.rst Web/coolprop/wrappers/VB.net/index.rst Web/coolprop/wrappers/Java/index.rst Web/coolprop/wrappers/R/index.rst Web/coolprop/wrappers/PHP/index.rst
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git commit --no-verify -m "docs(wrappers): update SWIG library names to CoolProp<Lang> (#1674, #2254)"
```

---

## Final verification (before PR)

- [ ] **Step 1: Confirm no stale `CoolProp` SWIG names remain in build config**

```bash
grep -Fn '\"CoolProp\"' CMakeLists.txt && echo "FAIL: dllimport stale" || echo "OK: dllimport renamed"
grep -En 'loadLibrary\("CoolProp"\)|dyn.load\(paste\("CoolProp"' dev/scripts/examples/example_generator.py && echo "FAIL: loaders stale" || echo "OK: loaders renamed"
```
Expected: both `OK`.

- [ ] **Step 2: Re-run the local C# build end-to-end**

```bash
rm -rf build_csharp install_root/Csharp
cmake -B build_csharp -S . -DCOOLPROP_CSHARP_MODULE=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_csharp --target install -j8
ls install_root/Csharp/Darwin_64bit/
```
Expected: `libCoolPropCsharp.dylib` present.

- [ ] **Step 3: Pre-PR adversarial review (per CLAUDE.md — MANDATORY)**

Invoke the `superpowers:code-reviewer` subagent against the diff vs `origin/master` before `gh pr create`, per the project checklist in `CLAUDE.md`.

- [ ] **Step 4: Run the pre-push gate**

```bash
./dev/ci/preflight.sh
```
Note: preflight auto-selects test scope from changed paths; this change touches only build/CI/docs files, so the C++ Catch2 scope may be empty — that is expected. Explain any `--skip` used.
```
