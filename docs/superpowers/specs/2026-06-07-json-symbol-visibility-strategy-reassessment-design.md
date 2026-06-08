# JSON symbol-visibility strategy — reassessment & durable fix — Design

**Date:** 2026-06-07
**Epic:** `CoolProp-xa8w` (RapidJSON → nlohmann/json migration)
**Issues:** `CoolProp-rj2i` (the bug — docs/Octave link failure), `CoolProp-w0cm` (poisoning CI guard), plus a new Option-A task to be filed.
**Supersedes:** the hidden-visibility approach in `migration design §4 (layer 3)`.
**Root-cause analysis (authoritative, empirically verified in a gcc:12 container):**
`docs/superpowers/specs/2026-06-07-json-hidden-visibility-link-failure.md`.

---

## 1. What went wrong

The migration kept nlohmann/valijson symbols out of CoolProp's exported ABI by
wrapping their includes in `#pragma GCC visibility push(hidden)` (in
`include/CoolProp/detail/json.h`). Phase 2-4 (PR #3094) switched the
Incompressible/PC-SAFT/UNIFAC loaders to include `detail/json.h`, which spread
that pragma into the loader translation units. Consequence (root-cause doc §2):

- nlohmann (via `detail/macro_scope.hpp`) and valijson (via `utf8_utils.hpp`)
  include `<cassert>` **inside** the hidden region.
- `<cassert>` re-declares `__assert_fail` **unguarded on every include** (by
  design, to support `NDEBUG` toggling). GCC pins a symbol's ELF visibility at
  its **first declaration** in the TU — so `__assert_fail` became
  `GLOBAL HIDDEN UND`.
- A hidden undefined symbol cannot bind to a shared library; ELF takes the most
  restrictive visibility across input objects, so one poisoned loader object
  made `__assert_fail` hidden for the whole `CoolProp.oct` link.
- Only the **docs CI** broke: it is the one build with no `-DNDEBUG`, so
  `assert()`/`__assert_fail` are live. Any non-Release shared build on glibc
  would hit it.

The pragma is the wrong tool: it operates at compile time and **poisons system
externs**, it is **fail-open across nlohmann/valijson version bumps** (a future
header that first-declares a new extern inside the region re-breaks non-NDEBUG
links), and it was **incomplete** — `superancillary.h` leaked nlohmann directly
(now fixed by PR #3104) and only `libCoolProp` is gate-checked, never the
wrappers (the `.oct` that actually broke).

## 2. Strategy shift: link-time export control (the report's Option A)

Replace compile-time visibility with **link-time export control**, applied per
shared product:

- **ELF (Linux):** `-Wl,--version-script=<map>` with `{ local: *nlohmann*; *valijson*; };`
- **Mach-O (macOS):** `-Wl,-unexported_symbols_list,<list>` matching `*nlohmann*` / `*valijson*`
- **Windows (MSVC):** no-op — exports are opt-in via `src/CoolPropLib.def`; nothing leaks.

This operates on the output's dynamic symbol table at link time — **orthogonal
to compile-time visibility**, so the libc-poisoning failure mode cannot exist.
It is a *hide-list* (hide only nlohmann/valijson; keep the public API exported),
so it needs no enumeration of the public API. It covers valijson and any
header-only JSON leak (incl. superancillary) uniformly, and lets the pragma, the
`detail/json.h` pre-include list, and all `CP_JSON_LOCAL` annotations be deleted.

**Build-structure facts** (mapped): there is no central shared-target helper —
~6 products build independently (`libCoolProp` SHARED, Octave `.oct`, Java, C#,
R via SWIG; Python via `Python_add_library`/`nanobind_add_module`), each
compiling the CoolProp sources directly. So the export control is applied via a
**reusable CMake helper `coolprop_hide_json_symbols(<target>)`** called once per
product (ELF/Mach-O branch inside). (nanobind already compiles
`-fvisibility=hidden`, which — per the verified table — hides symbols *without*
poisoning externs; the helper is still applied for uniformity/belt-and-suspenders.)

## 3. The two-part fix

### Part 1 — Tactical stopgap (unblocks docs CI now; this PR)

Complete the partial pre-include list in `detail/json.h`: pre-include the system
headers nlohmann/valijson reach (critically `<cassert>` — the unguarded
re-declaration that poisons `__assert_fail`) **before** the pragma push, pinning
those externs at default visibility. Near-zero-risk, ~a dozen `#include` lines.

`<cassert>` is the essential, principled fix (it is uniquely unguarded among
standard headers; most others are pinned default by CoolProp's broad earlier
includes). A modest set of additional libc/runtime headers
(`<cstdlib> <cstring> <cstdio> <cmath> <memory> <new> <exception>`) is added as
defense-in-depth for any build whose include order differs.

This stopgap is explicitly temporary — **Part 2 deletes the entire pre-include
block and the pragma.** It carries a `// TODO(CoolProp-<optionA>)` marker.

### Part 2 — Durable Option A (separate follow-up PR, own plan)

1. Add `cmake/CoolPropJSONVisibility.cmake` (or inline helper) defining
   `coolprop_hide_json_symbols(<target>)` + the two symbol-list files
   (`dev/linker/coolprop_json.version`, `dev/linker/coolprop_json_unexport.txt`).
2. Call it for every shared product (libCoolProp, Octave, Java, C#, R, Python).
3. **Delete** the `#pragma GCC visibility push/pop` + the pre-include stopgap in
   `detail/json.h`, and the `CP_JSON_LOCAL` macro + its 12 usages across
   `FluidLibrary.h`, `IncompressibleLibrary.h`, `PCSAFTLibrary.h`,
   `UNIFACLibrary.h`.
4. **Extend the symbol gate** to the wrapper products (Octave at minimum) — the
   gate today only inspects `libCoolProp`, which is why the `.oct` leak/poison
   was invisible. (Build the Octave `.oct` in the gate job and run
   `check-json-symbols.sh` on it too.)
5. Verify: `check-json-symbols.sh` stays green on every product; the docs
   (non-NDEBUG) Octave build links; full suite + symbol gate green.

## 4. Epic reconciliation

- **`rj2i`** (bug) → closed by Part 1 (stopgap) landing; Part 2 makes the fix
  permanent.
- **`w0cm`** (CI guard for hidden-undefined libc externs) → **becomes moot once
  the pragma is deleted in Part 2** (the poisoning mechanism is gone). Close it
  when Part 2 lands. (Its intent is subsumed by the extended symbol gate +
  "no compile-time visibility regions around system headers" rule.)
- **Phase Final (`xa8w.3`)** "tighten gates / decide getter-layer fate": the
  wrapper-gate extension and the `CP_JSON_LOCAL` deletion belong here or in the
  Option-A PR; reference this design from `xa8w.3`.
- **No conflict** with in-flight `rxxx` (#3100, install-exclude) or the merged
  superancillary de-leak (#3104).

## 5. Constraints / risks

- **Local repro is limited:** the bug is ELF/glibc (Linux); this dev machine is
  macOS/Mach-O. The stopgap is verifiable by reasoning + a Mach-O build; the
  ELF link and the docs Octave build are verified in **CI**. Option A's
  per-wrapper link must be iterated against CI (the docs workflow is the real
  test for the `.oct`).
- **Version-script glob matches mangled names** (`*nlohmann*` matches
  `_ZN8nlohmann...`); confirm in CI that the public C API + `CoolProp::*` stay
  exported (the symbol gate's existing nm check + a smoke import per wrapper).
- The stopgap remains fail-open across dependency bumps — acceptable because it
  is short-lived and Part 2 removes the mechanism entirely.

## 6. Scope of THIS spec's first PR

Part 1 (the stopgap) only: edit `include/CoolProp/detail/json.h` to pre-include
the system headers before the pragma. Part 2 (Option A) gets its own plan + PR.
