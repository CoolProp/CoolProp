# Option A — link-time JSON-symbol export control — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the fragile compile-time `#pragma GCC visibility push(hidden)` + `CP_JSON_LOCAL` + the rj2i pre-include stopgap with **link-time export control**: a per-shared-product linker hide-list for nlohmann/valijson symbols. Mechanism-proof against the libc-poisoning failure (CoolProp-rj2i) and immune to dependency bumps.

**Architecture:** A reusable CMake helper `coolprop_hide_json_symbols(<target>)` attaches, per shared product, an ELF `--version-script` (`{ local: *nlohmann*; *valijson*; };`) or a Mach-O `-unexported_symbols_list` (`*nlohmann*`/`*valijson*`). Windows (MSVC) is a no-op (opt-in `.def` exports). With link-time hiding in place, the compile-time visibility machinery is deleted: the pragma + stopgap pre-includes in `detail/json.h`, and the `CP_JSON_LOCAL` macro + its uses across the 4 loader headers.

**Tech Stack:** CMake (`target_link_options`, ELF version scripts / Mach-O `-unexported_symbols_list`), the symbol-leak gate (`dev/ci/check-json-symbols.sh`).

**Spec/design:** `docs/superpowers/specs/2026-06-07-json-symbol-visibility-strategy-reassessment-design.md` (§2-5, Option A).
**Issue:** `CoolProp-xa8w.6`. **Branch:** `ihb/json-optionA-link-export-control` (stacked on the Phase 5 branch).

**Verification reality:** this is an ELF/glibc-class change; the dev machine is macOS/Mach-O. **Local proof = the Mach-O symbol gate on `libCoolProp.dylib` stays green AFTER the pragma is deleted** (which proves the link-time hide-list works, since compile-time hiding is gone). The ELF version script + the Octave `.oct` link are verified in **CI** (the docs build + the symbol gate). Iterate against CI.

---

## Task 1: Symbol-list files + the CMake helper

**Files:**
- Create: `dev/linker/coolprop_hide_json.map` (ELF version script)
- Create: `dev/linker/coolprop_hide_json.exp` (Mach-O unexported list)
- Create: `cmake/CoolPropJSONVisibility.cmake` (the helper)
- Modify: `CMakeLists.txt` (include the helper module once, near the top after `CMAKE_MODULE_PATH` is set)

- [ ] **Step 1: ELF version script `dev/linker/coolprop_hide_json.map`**
```
{
  local:
    *nlohmann*;
    *valijson*;
};
```
(Anonymous version with only a `local:` list — GNU ld keeps all unmatched symbols at their default (exported) binding, and forces matching ones local/hidden. Matches MANGLED names, which contain the literal `nlohmann`/`valijson`.)

- [ ] **Step 2: Mach-O unexported list `dev/linker/coolprop_hide_json.exp`**
```
*nlohmann*
*valijson*
```
(`-unexported_symbols_list` glob patterns over mangled symbol names; `*nlohmann*` matches `__ZN8nlohmann…`.)

- [ ] **Step 3: `cmake/CoolPropJSONVisibility.cmake`**
```cmake
# Hide nlohmann/valijson symbols from a shared product's dynamic export table at
# LINK time (replaces the compile-time visibility pragma; see CoolProp-xa8w.6).
# ELF: --version-script hide-list; Mach-O: -unexported_symbols_list. MSVC: no-op
# (exports are opt-in via src/CoolPropLib.def).
function(coolprop_hide_json_symbols target)
    if(NOT TARGET ${target})
        return()
    endif()
    if(APPLE)
        set(_exp "${CMAKE_CURRENT_SOURCE_DIR}/dev/linker/coolprop_hide_json.exp")
        target_link_options(${target} PRIVATE "LINKER:-unexported_symbols_list,${_exp}")
        set_property(TARGET ${target} APPEND PROPERTY LINK_DEPENDS "${_exp}")
    elseif(UNIX)
        set(_map "${CMAKE_CURRENT_SOURCE_DIR}/dev/linker/coolprop_hide_json.map")
        target_link_options(${target} PRIVATE "LINKER:--version-script=${_map}")
        set_property(TARGET ${target} APPEND PROPERTY LINK_DEPENDS "${_map}")
    endif()
    # MSVC/WIN32: nothing — only .def-listed symbols are exported.
endfunction()
```
(`CMAKE_CURRENT_SOURCE_DIR` resolves to the repo root for the top-level CMakeLists and to the wrapper dir for `wrappers/Python` — so use an absolute path via `${CMAKE_SOURCE_DIR}` instead if called from a subdirectory. **Use `${CMAKE_SOURCE_DIR}` in the helper** so it works from any CMakeLists. Adjust the two `set(_exp …)`/`set(_map …)` lines to `${CMAKE_SOURCE_DIR}/dev/linker/…`.)

- [ ] **Step 4: include the module in `CMakeLists.txt`**

After the `list(APPEND CMAKE_MODULE_PATH …)` / project setup near the top, add:
```cmake
include(CoolPropJSONVisibility)
```
(Confirm `CMAKE_MODULE_PATH` includes `${CMAKE_SOURCE_DIR}/cmake`; if the helper isn't found, add `list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")` first. The Python wrapper already does `list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")` then includes CPM — it can `include(CoolPropJSONVisibility)` after adding the ROOT cmake dir to its module path.)

- [ ] **Step 5: configure-only smoke**

Run: `cmake -B build_shared -S . -DCOOLPROP_SHARED_LIBRARY=ON -DCMAKE_BUILD_TYPE=Release` — expect clean configure (the helper parses; no target yet wired — that's Task 2). Commit after Task 2 (these files are inert until applied).

---

## Task 2: Apply the helper to every shared product

**Files:** `CMakeLists.txt` (libCoolProp + Octave/Java/C#/R/PHP swig targets), `wrappers/Python/CMakeLists.txt` (Cython + nanobind modules)

For EACH shared product, add `coolprop_hide_json_symbols(<target>)` right after the target is created / its libraries are linked. The SWIG targets are all named `CoolProp`; the main lib is `${LIB_NAME}`; the Python modules are `CoolProp_module`.

- [ ] **Step 1: libCoolProp** — after `add_library(${LIB_NAME} SHARED …)` (~CMakeLists.txt:651) and its `target_link_libraries`, add `coolprop_hide_json_symbols(${LIB_NAME})`.
- [ ] **Step 2: Octave** — after `swig_add_library(CoolProp LANGUAGE octave …)` + `swig_link_libraries` (~1409), add `coolprop_hide_json_symbols(CoolProp)`. (This is THE product whose link broke; it must be covered.)
- [ ] **Step 3: Java** (`swig_add_module(CoolProp java …)` ~1680), **C#** (~1455/1544), **R** (~1620), **PHP** (~1829): add `coolprop_hide_json_symbols(CoolProp)` after each module's creation/link (inside each `if(COOLPROP_*_MODULE)` block, so only configured ones apply).
- [ ] **Step 4: Python** (`wrappers/Python/CMakeLists.txt`): add `coolprop_hide_json_symbols(CoolProp_module)` after `Python_add_library(CoolProp_module …)` (~184) and after `nanobind_add_module(CoolProp_module …)` (~218). Ensure the helper is `include()`d there (add `list(APPEND CMAKE_MODULE_PATH "${ROOT_DIR}/cmake")` + `include(CoolPropJSONVisibility)` near its top).
- [ ] **Step 5: build + Mach-O gate (the local proof that hiding still works WITH the pragma still present)**

Run: `cmake -B build_shared -S . -DCOOLPROP_SHARED_LIBRARY=ON -DCMAKE_BUILD_TYPE=Release && cmake --build build_shared -j8 && ./dev/ci/check-json-symbols.sh "$(find build_shared -name 'libCoolProp.dylib' | head -1)"`
Expected: `OK` (0 nlohmann/valijson exported). At this point BOTH the pragma AND the version script hide them — the next task removes the pragma and the gate must STILL be green (proving the link-time hide is doing the work).

- [ ] **Step 6: Commit** (Tasks 1+2 together — the files + their application)
```bash
git add dev/linker/ cmake/CoolPropJSONVisibility.cmake CMakeLists.txt wrappers/Python/CMakeLists.txt
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit -m "build(json): add link-time export control for nlohmann/valijson symbols (CoolProp-xa8w.6)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Delete the compile-time visibility machinery

**Files:** `include/CoolProp/detail/json.h` (pragma + stopgap), the 4 loader headers (`CP_JSON_LOCAL`).

- [ ] **Step 1: `detail/json.h` — delete the stopgap pre-include block AND the pragma**

Remove the `// STOPGAP (CoolProp-rj2i) …` comment block and ALL the pre-includes it added (`<cassert> <cmath> <cstddef> <cstdint> <cstdio> <cstdlib> <cstring> <exception> <iostream> <limits> <memory> <new> <stdexcept> <string> <string_view> <system_error> <typeinfo> <vector>`), AND the `#if defined(__GNUC__)… pragma GCC visibility push(hidden) …#endif` / `…pop…#endif` wrappers (lines ~17-25 and ~49-51 and ~59-61). Keep the actual includes:
```cpp
#include <nlohmann/json.hpp>
#include <valijson/adapters/nlohmann_json_adapter.hpp>
#include <valijson/schema.hpp>
#include <valijson/schema_parser.hpp>
#include <valijson/validator.hpp>
```
(now at default visibility — hidden at link by the version script). Re-add only the standard headers the file's OWN code actually needs (`<string>`, `<string_view>`, `<vector>`, `<cstdint>` if used by the `cpjson` helpers below — check what `detail/json.h`'s own code uses and keep those; drop the rest). Update the file-top comment block (the "Symbol-leak strategy" note) to describe the link-time approach.

- [ ] **Step 2: Delete `CP_JSON_LOCAL` in the 4 loader headers**

In each of `src/Backends/Helmholtz/Fluids/FluidLibrary.h`, `src/Backends/Incompressible/IncompressibleLibrary.h`, `src/Backends/PCSAFT/PCSAFTLibrary.h`, `src/Backends/Cubics/UNIFACLibrary.h`: delete the `#if defined(__GNUC__)… #define CP_JSON_LOCAL …#endif` macro definition block, remove the `CP_JSON_LOCAL ` prefix from every method declaration that has it, and delete the `#undef CP_JSON_LOCAL` line. Also update the comment in `src/Backends/Helmholtz/Fluids/FluidLibraryFactories.h` (~line 15) that references the `CP_JSON_LOCAL` treatment.

- [ ] **Step 3: Build + the decisive Mach-O gate**

Run: `cmake --build build_shared -j8 && ./dev/ci/check-json-symbols.sh "$(find build_shared -name 'libCoolProp.dylib' | head -1)"`
Expected: clean build AND `OK` (0 nlohmann/valijson exported). **This is the key local proof:** with the pragma and `CP_JSON_LOCAL` GONE, the only thing hiding the symbols is the link-time version script — a green gate confirms it works. If the gate now FAILS (symbols leaked), the version script / unexported list isn't matching — fix the pattern (check `nm -gU libCoolProp.dylib | c++filt | grep -E 'nlohmann|valijson'` to see what leaked and adjust the `.exp`/`.map` globs).
Also build the Catch runner (`cmake --build build_catch --target CatchTestRunner -j8`) — confirm deleting the pragma/CP_JSON_LOCAL doesn't break compilation.

- [ ] **Step 4: Commit**
```bash
git add include/CoolProp/detail/json.h src/Backends/Helmholtz/Fluids/FluidLibrary.h src/Backends/Helmholtz/Fluids/FluidLibraryFactories.h src/Backends/Incompressible/IncompressibleLibrary.h src/Backends/PCSAFT/PCSAFTLibrary.h src/Backends/Cubics/UNIFACLibrary.h
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit -m "refactor(json): delete the visibility pragma + CP_JSON_LOCAL — link-time export control replaces them (CoolProp-xa8w.6)

Removes the rj2i pre-include stopgap, the GCC visibility pragma around the
nlohmann/valijson includes (which poisoned __assert_fail in non-NDEBUG loader
builds), and the CP_JSON_LOCAL per-method attributes. Symbols are now hidden at
link time per shared product (Task 1-2), which cannot poison system externs.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Extend the symbol gate to a wrapper product

**Files:** `dev/ci/preflight.sh` (the json-symbols step), optionally `dev/ci/check-json-symbols.sh` (unchanged — it already takes any product path).

The gate today inspects only `libCoolProp` — which is why the `.oct` leak/poison was invisible. Add the **Octave `.oct`** (the product that broke) to the gate, gracefully skipping when Octave/SWIG aren't installed.

- [ ] **Step 1: preflight — build + gate the Octave module when available**

In the `json-symbols` step of `dev/ci/preflight.sh`, after the existing `libCoolProp` check, add a sub-check: if `swig` and `octave-config` (or `octave`) are on PATH, configure+build a `build_octave` with `-DCOOLPROP_OCTAVE_MODULE=ON` and run `./dev/ci/check-json-symbols.sh` on the resulting `CoolProp.oct`. Skip (not fail) with a clear message when the Octave toolchain is absent (so it runs in CI, skips on dev machines without Octave). Fail-closed if the build is attempted and fails, or if the `.oct` exports nlohmann/valijson.

- [ ] **Step 2: Verify the step skips cleanly locally** (no Octave on the dev machine → skip message, not fail) and document that CI (which has Octave) runs it.

- [ ] **Step 3: Commit**
```bash
git add dev/ci/preflight.sh
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit -m "ci(json): extend the symbol-leak gate to the Octave .oct wrapper (CoolProp-xa8w.6)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final Verification (before PR)

- [ ] **Step 1: Mach-O symbol gate green with the pragma gone** (done in Task 3) — re-confirm `OK` on a fresh `build_shared`.
- [ ] **Step 2: No compile-time visibility machinery left**
  Run: `grep -rn "CP_JSON_LOCAL\|pragma GCC visibility" include src | grep -v "FluidLibrary.h:.*comment"` → expect nothing (all pragma/CP_JSON_LOCAL gone; the only `visibility` left should be the linker files / none).
- [ ] **Step 3: Full fast suite** `./build_catch/CatchTestRunner "~[slow]"` → 0 failures (deleting the pragma is behaviour-neutral).
- [ ] **Step 4: Preflight + adversarial review (MANDATORY).** Focus: (a) the version-script/`.exp` globs actually hide nlohmann/valijson AND keep the public C API + `CoolProp::*` exported (the gate's nm check + the Catch tests linking against the lib prove the public API survives); (b) the helper is applied to EVERY shared product (none missed — a missed product re-leaks/re-exposes); (c) `target_link_options` LINKER: prefix is correct for both ELF and Mach-O; (d) no fail-open in the extended gate (Octave skip is legitimate-absent only, not masking a build failure); (e) deleting the pragma did not leave `detail/json.h` needing a header it dropped.
- [ ] **Step 5: Push + PR (stacked on the Phase 5 branch).** `git push -u origin ihb/json-optionA-link-export-control`; `gh pr create --base ihb/json-phase5-configuration-off-rapidjson …` (the parent in the stack). In the PR body, be explicit that the ELF link + the Octave `.oct` are verified by **CI** (the docs build + the symbol gate), not locally. After CI runs, **re-check the docs "Documentation builds (HTML)" job goes green** — that is the end-to-end proof Option A fixed the poisoning durably (and the stopgap can later be considered redundant). Re-review any post-review commits.
- [ ] **Step 6: bd.** `bd update CoolProp-xa8w.6 --notes "PR #…"`; note that `w0cm` (the poisoning guard) is now moot — the pragma is gone — and can close once this merges; note `rj2i` is permanently fixed (the stopgap is deleted here, superseded by link-time control).

---

## Self-Review (plan vs design)

**Design coverage:** reassessment §2 (Option A mechanism) → Tasks 1-2 (helper + per-product); §3/§5 (delete pragma + CP_JSON_LOCAL + stopgap) → Task 3; "extend the gate to wrappers" → Task 4; "w0cm moot / rj2i permanent" → Final Verification §6.

**Placeholder scan:** the per-product line numbers are approximate (the implementer confirms each target's exact creation site); the symbol-list globs and the helper are given in full. No "TODO" left in shipped code.

**Risk/verification consistency:** every step that can be checked on macOS (configure, build, Mach-O symbol gate, Catch suite) is a local gate; the ELF version script + Octave `.oct` are explicitly deferred to CI, stated in Task intro, Task 3 Step 3, and Final Verification Step 5.
