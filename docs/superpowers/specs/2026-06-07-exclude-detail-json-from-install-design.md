# Exclude `detail/json.h` from installed headers — Design

**Date:** 2026-06-07
**Status:** Approved design, ready for implementation planning
**Epic:** `CoolProp-xa8w` (RapidJSON → nlohmann/json migration)
**Issue:** `CoolProp-rxxx` (narrowed to the install-exclude; the superancillary de-leak is split into a new sibling issue)

---

## 1. Goal

Stop shipping `include/CoolProp/detail/json.h` in the installed header set. That
header `#include`s `nlohmann/json.hpp` + valijson (build-only CPM dependencies),
so shipping it puts a header pulling those deps into the public install tree even
though no consumer needs it. This is one of the two installed headers that pull
nlohmann; `superancillary/superancillary.h` is the other and is handled
separately.

## 2. Key facts (verified)

- **No installed header `#include`s `detail/json.h`.** A repo-wide grep of
  `include/CoolProp/**` finds only a *comment* reference in `detail/rapidjson.h`
  (line ~303). The actual `#include`rs of `detail/json.h` are all non-installed:
  `src/Backends/Helmholtz/Fluids/FluidLibrary.{h,cpp}`, `FluidLibraryFactories.h`,
  `src/Backends/PCSAFT/PCSAFTLibrary.{h,cpp}`, `src/Backends/Incompressible/…`,
  `src/Backends/Cubics/CubicsLibrary.cpp`, and tests.
- **`detail/rapidjson.h` does NOT include `detail/json.h`.** Line 14 includes
  `<rapidjson/rapidjson.h>` (the library). The two CoolProp wrappers each define
  the `cpjson::schema_validation_code` enum independently, both guarded by
  `CPJSON_SCHEMA_VALIDATION_CODE_DEFINED` so a TU including both does not get a
  redefinition error. **There is no include coupling between them**, so no enum
  decoupling is required.
- **`detail/json.h` is currently installed.** The header-install rules ship
  `include/CoolProp/**/*.h`; the existing `PATTERN "*_JSON*.h" EXCLUDE` does NOT
  match `json.h`, and nothing else excludes it. So it lands in the install tree
  today as an orphan (shipped, never `#include`d by any shipped header).

Because nothing installed reaches `detail/json.h`, removing it from the shipped
set cannot break any consumer's include graph. This is pure hygiene.

## 3. The change (CMake-only; no source edits)

Add a precise exclude for `detail/json.h` to every header-install rule:

- **`CMakeLists.txt`** — three `install(DIRECTORY ${...}/include/CoolProp
  DESTINATION <dest>/include FILES_MATCHING PATTERN "*.h")` clauses (static at
  ~L631, shared at ~L668, cmake-install at ~L857). The *first* clause in each
  trio ships the `CoolProp/` tree (incl. `detail/`); that is the one needing the
  exclude. The *second* clause (`include/` with `PATTERN "CoolProp" EXCLUDE`)
  already skips the `CoolProp` directory, so it needs no change.
- **`wrappers/Python/CMakeLists.txt`** — the `install(DIRECTORY
  "${ROOT_DIR}/include/" DESTINATION CoolProp/include … PATTERN "*_JSON.h"
  EXCLUDE PATTERN "*_JSON_z.h" EXCLUDE)` clause ships all of `include/`; add the
  exclude here too.

Matcher: `REGEX "detail/json\\.h$" EXCLUDE` — anchored so it targets exactly
`detail/json.h` and does **not** match `detail/rapidjson.h` (which must keep
shipping until Phase Final) or any `*_JSON.h` data header. (`PATTERN "json.h"`
would also work since only one file is named `json.h`, but the anchored REGEX is
unambiguous and self-documenting.)

## 4. Regression guard — `dev/ci/check-installed-headers.sh`

A lightweight, **fail-closed** script that proves the exclude actually holds (not
just that the CMakeLists text looks right), wired into `dev/ci/preflight.sh`:

- Takes a configured+built shared build dir (reuses `build_shared`, which
  preflight already builds for the symbol-leak gate) and runs `cmake --install
  <dir> --prefix <tmpdir>`.
- **Asserts:** `detail/json.h` is **absent** from the installed tree, AND
  `detail/rapidjson.h` is **present** (a positive control proving the exclude is
  specific, not over-broad, and that the install actually produced headers).
- **Fail-closed:** if `cmake --install` fails, or the install tree contains zero
  headers (wrong dir / empty), the check FAILS rather than vacuously passing —
  same discipline as `check-json-symbols.sh` (no `|| true` that masks a failure;
  no empty-result false pass).
- Preflight wiring: a new step near the symbol-gate step (so it reuses the
  already-built `build_shared`); honors a `--skip=install-headers` flag and skips
  gracefully when the shared build is unavailable (mirrors the json-symbols
  step's skip logic).

**Scope of the assertion:** intentionally narrow — *"`detail/json.h` is not
installed."* The broader invariant *"no shipped header includes nlohmann/
valijson"* cannot pass yet because `superancillary.h` still ships
`nlohmann/json.hpp`; the script is written so that broader assertion can be added
later (when superancillary is de-leaked) without restructuring. A `# TODO`
documents that follow-on, so the narrowing is explicit, not silent.

## 5. Testing & verification

- **Behavior:** install-time-only change; compilation is unaffected. A normal
  `[!slow]` Catch run is a sanity check, not the point.
- **The guard is the test:** `dev/ci/check-installed-headers.sh build_shared`
  passes (json.h absent, rapidjson.h present).
- **Manual cross-check** in the plan: `cmake --install` to a temp prefix in the
  shared config and `find <prefix> -name json.h -path '*detail*'` returns nothing.
- **Preflight:** `./dev/ci/preflight.sh` runs the new step; clang-format/cppcheck/
  clang-tidy are unaffected (no source change). Pre-PR adversarial review per
  CLAUDE.md (focus: the fail-closed guarantees of the new script — no path where
  it passes when `detail/json.h` is actually shipped).

## 6. bd bookkeeping

- Narrow `CoolProp-rxxx` to this install-exclude (update its title/description).
- File a new sibling issue under `CoolProp-xa8w` for the superancillary de-leak
  (de-leak `nlohmann::json` from the installed templated `superancillary.h` via
  explicit instantiation of the construction path — `get_superanc()` constructs
  inline in `CoolPropFluid.h`, so the string/JSON ctor must move out-of-line into
  `libCoolProp`). Make it block Phase Final (`CoolProp-xa8w.3`) alongside `rxxx`.

## 7. Out of scope

- `superancillary.h` nlohmann de-leak (separate issue/effort).
- `detail/rapidjson.h` install-exclude and the dual-enum guard removal
  (Phase Final, `CoolProp-xa8w.3`).
- The broader "no shipped header pulls a build-only JSON dep" consumer-compile
  gate (deferred; coordinated with the nanobind interface).
