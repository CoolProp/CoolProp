# Design: restore C# wrapper downloads (#2254) + language-prefix SWIG libraries (#1674)

Date: 2026-05-29
Issues: [#2254](https://github.com/CoolProp/CoolProp/issues/2254), [#1674](https://github.com/CoolProp/CoolProp/issues/1674); likely also resolves [#2326](https://github.com/CoolProp/CoolProp/issues/2326)
Target: next major release

## Problem

1. **#2254** — Since v6.4.2 the `Csharp/` wrapper folder (and other SWIG language
   wrappers) disappeared from the SourceForge download tree. Root cause: the build
   was migrated to GitHub Actions, but no C# builder workflow was ever ported. The
   release pipeline assembles the download tree from a fixed list of builder
   workflows in `.github/workflows/release_all_files.yml` (the `collect_binaries`
   matrix, currently `release_all_files.yml:73`); no `csharp_builder.yml` exists and
   nothing references `COOLPROP_CSHARP_MODULE`, so no C# artifact is ever produced or
   deployed. The CMake build target itself (`CMakeLists.txt:1396`) is intact.

2. **#1674** — Every SWIG wrapper is declared `swig_add_module(CoolProp <lang> ...)`,
   so each emits a binary named `CoolProp`, colliding with the conventional shared
   library produced by `library_shared.yml`. This name collision is a recurring
   source of user confusion. Affected: C#, VB.NET, Java, R, PHP. Octave is exempt
   (`.oct` suffix). The collision is the suspected cause of #2326
   (`EntryPointNotFound` on 64-bit C#).

A major release is the right time to make the breaking rename in #1674.

## Decisions (locked during brainstorming)

- **#1674 scope**: rename all 5 SWIG wrappers (C#, VB.NET, Java, R, PHP), not just C#.
- **Naming convention**: `CoolProp<Lang>` — `CoolPropCsharp`, `CoolPropJava`,
  `CoolPropR`, `CoolPropPHP`. VB.NET reuses the C# binding → `CoolPropCsharp`.
- **C# CI build matrix**: Windows (x64), Linux, macOS. Legacy Win32 dropped.
- **CI builders for Java/R/PHP**: out of scope. The rename lands in the build system
  for all 5 so manual/other builds get correct names, but only C# gets a CI builder.

## Guiding principle for the rename

Change only the **shared-library filename** (`OUTPUT_NAME`) and the **loader
reference** that points at it. Keep `%module CoolProp` in `src/CoolProp.i`
untouched, so all user-facing class/namespace/API names (`CoolProp.PropsSI`, etc.)
are unchanged. The CMake target remains named `CoolProp`; only `OUTPUT_NAME` changes.

## Part A — #1674: rename the SWIG shared libraries

Per wrapper in `CMakeLists.txt` (and `dev/scripts/examples/example_generator.py`):

| Wrapper | Module guard | Output name | Loader reference to update |
|---|---|---|---|
| C# | `COOLPROP_CSHARP_MODULE` (`:1396`) | `CoolPropCsharp` | `-dllimport "CoolProp"` → `-dllimport "CoolPropCsharp"` (`:1407`) |
| VB.NET | `COOLPROP_VBDOTNET_MODULE` (`:1486`) | `CoolPropCsharp` | `-dllimport "CoolProp"` → `-dllimport "CoolPropCsharp"` (`:1496`) |
| Java | `COOLPROP_JAVA_MODULE` (`:1606`) | `CoolPropJava` | `example_generator.py:934` `System.loadLibrary("CoolProp")` → `loadLibrary("CoolPropJava")` |
| R | `COOLPROP_R_MODULE` (`:1555`) | `CoolPropR` | `example_generator.py:806` `dyn.load(paste("CoolProp", ...))` → `"CoolPropR"` |
| PHP | `COOLPROP_PHP_MODULE` (`:1751`) | `CoolPropPHP` | none required — see PHP finding below |

Mechanics:
- Add `set_target_properties(CoolProp PROPERTIES OUTPUT_NAME "CoolProp<Lang>")` to
  each module block. The existing `PREFIX ""` / `PREFIX "lib"` lines stay.
- VB.NET: the existing `POST_BUILD` copy of `$<TARGET_FILE:CoolProp>`
  (`CMakeLists.txt:1533`) follows the rename automatically — no path edit needed.
- Install destinations (`Csharp/`, `Java/`, `R/`, `PHP/`, `VB.NET/`) are unchanged.

**PHP finding (tested with SWIG 4.3.0, 2026-05-29)**: the rename is safe with no
loader patch. Verified by running `swig -c++ -php7` on a minimal `%module CoolProp`
interface:
- The module name `"CoolProp"` appears only in the internal `zend_module_entry`,
  `SWIG_name`, and `INIT_CLASS_ENTRY` (`CoolProp_wrap.cxx`). That is the API-level
  extension name (used by `extension_loaded("CoolProp")` and `get_module()`), which
  we keep unchanged.
- There is **no** hardcoded `.so` filename and **no** `dl()` call in the generated
  wrapper. PHP loads a C extension by filename (`extension=CoolPropPHP.so` or
  `dl("CoolPropPHP.so")`), then calls `get_module()`; filename and internal module
  name are decoupled, so renaming the `.so` to `CoolPropPHP` requires no patch.

**Pre-existing PHP bug — fixed as a drive-by in this branch**: SWIG 4.1+ rewrote the
PHP module to register classes as native PHP 8 Zend classes (verified: with an
`AbstractState` class in the interface, SWIG 4.3 emits `PHP_METHOD(AbstractState,...)`
directly in `php_CoolProp.h`/`_wrap.cxx` and produces **no** `CoolProp.php` proxy).
But `CMakeLists.txt:1796` still does
`install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CoolProp.php ...)`, so
`cmake --build --target install` for the PHP module hard-fails on any SWIG ≥ 4.1.

Fix: add the `OPTIONAL` keyword to that `install(FILES ...)`. This stops the failure
on modern SWIG (file absent → skipped) while still installing the proxy on a pre-4.1
SWIG that does generate it. Minimal and version-agnostic.

Note on PHP wrapper viability: the generated extension itself is functional — the
`.so` compiles and registers `PropsSI`/`AbstractState` natively for PHP 8 (SWIG ≥ 4.1
dropped PHP 5/7). The wrapper is shippable once the install step is fixed; it is
simply not built in CI (no PHP builder workflow, which remains out of scope).

## Part B — #2254: C# builder workflow

New `.github/workflows/csharp_builder.yml`, modeled on `libreoffice_builder.yml`
(same `on:` triggers and `concurrency` block).

- **Build matrix**: `windows-latest`, `ubuntu-latest`, `macos-latest` (all x64).
- **Per-OS job**:
  - Install SWIG and a C# toolchain (mono / .NET) sufficient for `FindCsharp.cmake`.
  - `cmake -B build -S . -DCOOLPROP_CSHARP_MODULE=ON`
  - `cmake --build build --target install`
  - Upload a per-OS artifact (`Csharp-windows`, `Csharp-linux`, `Csharp-macos`) so
    `actions/upload-artifact` v4 unique-name rules are satisfied. Each contains that
    OS's `install_root/Csharp/` (platform-independent `.cs` + `platform-independent.7z`
    + the `Csharp/<System>_64bit/` native lib).
- **Merge job** (in the same workflow, `needs:` the matrix): download the 3 per-OS
  artifacts, assemble a single `Csharp/` tree — platform-independent files once
  (deduped), each platform's native-lib subfolder alongside — and upload a single
  artifact named `Csharp`.

**Release wiring**: add `csharp_builder.yml` to the `collect_binaries` matrix in
`release_all_files.yml:73`. `release_get_artifact.yml` downloads the latest
successful run's single `Csharp` artifact, which then flows through `deploy_files`
(rsync over SSH) into the SourceForge tree — restoring the `Csharp/` folder.

## Part C — Changelog

Add to `Web/coolprop/changelog.rst` under the next major version:
- **Breaking change**: SWIG-generated shared libraries renamed `CoolProp` →
  `CoolProp<Lang>` (C#, VB.NET, Java, R, PHP). Consumers must update their
  `DllImport` / `System.loadLibrary` / `dyn.load` references accordingly. The
  CoolProp API names (classes, namespaces, functions) are unchanged.
- **Fixed**: C# (and other SWIG) wrapper downloads restored to the release tree
  (#2254). Resolves the DLL-name collision with the conventional CoolProp shared
  library (#1674) and the associated `EntryPointNotFound` on 64-bit C# (#2326).

## Testing & verification

- **C#**: CMake already carries Example.cs test scaffolding (`CMakeLists.txt`
  ~`:1460+`). CI runs it on each OS, directly validating that the `-dllimport` name
  matches the renamed DLL — this is the concrete check for both #1674 and #2326.
- **Java / R**: existing example tests (`R_test`, the Java example run at
  `CMakeLists.txt:1683`) validate the loader-name updates.
- **PHP**: rename verified safe by SWIG-generation inspection (see PHP finding); no
  runtime test needed for the rename itself.
- **Local**: the full SWIG toolchain set is not all installable locally; rely on the
  CI matrix for cross-platform confirmation.

## Drive-by fix (folded in)

- Add `OPTIONAL` to `install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CoolProp.php ...)` at
  `CMakeLists.txt:1796` so the PHP module's `--target install` no longer hard-fails
  on SWIG ≥ 4.1 (which no longer generates the proxy). See the PHP finding above.

## Out of scope

- Building Java / R / PHP wrappers in CI (no builder workflows added).
- Legacy Win32 C# binaries.
- Any change to the public CoolProp API surface.
