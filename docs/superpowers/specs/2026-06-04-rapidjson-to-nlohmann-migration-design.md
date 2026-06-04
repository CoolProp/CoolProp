# RapidJSON → nlohmann/json Migration — Design

**Date:** 2026-06-04
**Status:** Approved design, ready for implementation planning
**Scope:** Replace RapidJSON entirely with nlohmann/json across the CoolProp C++
library, with zero JSON-library symbols exposed in the public ABI.

---

## 1. Motivation

RapidJSON has not had a meaningful release in roughly a decade. Because it is a
header-only library whose symbols are embedded into every consumer, version
skew between CoolProp's copy and another loaded package's copy has produced ODR
violations / symbol clashes in the field. Consolidating onto a single,
actively-maintained JSON library (nlohmann/json) removes the stale dependency —
but only if we *also* prevent the new library's symbols from leaking into the
ABI, or we will have simply traded one clash source for another.

This migration is a breaking change. It does not have to ship in v8, but since
the v8/v9 line already breaks API/ABI in other ways, this is a natural time to
land it.

## 2. Goals & Constraints

**Goals**

- Remove RapidJSON from the codebase and build entirely.
- Standardize on nlohmann/json internally, using its idiomatic API at call
  sites (tight coupling is acceptable and desired).

**Hard constraints**

1. **Zero JSON symbols in the public ABI.** No `nlohmann::` (or `rapidjson::`)
   type appears in any *installed* header, and no such symbol is exported from
   the built shared library. This is the ODR/version-clash fix that motivates
   the effort, and it is enforced mechanically (Section 7).
2. **No maintenance-heavy wrapper layer.** Symbol-hiding comes from structure
   plus build flags, not a hand-written facade that must be babysat.
3. **CPM-managed dependencies.** nlohmann/json and the schema validator are
   fetched via CPM in `cmake/dependencies.cmake`, exactly as RapidJSON is today.
   Nothing is hand-vendored under `externals/`.

## 3. Current State (as surveyed 2026-06-04)

- RapidJSON is fetched via CPM (`cmake/dependencies.cmake:33-38`,
  `DOWNLOAD_ONLY`), include dir added at `CMakeLists.txt:336`. It is **not** a
  git submodule or an `externals/` checkout.
- A `cpjson` helper namespace lives in
  `include/CoolProp/detail/rapidjson.h` (~15 getter/setter/serialize/validate
  helpers) — the natural migration seam.
- **~15 files** include RapidJSON headers. The per-backend library headers
  (`src/Backends/.../FluidLibrary.h`, `PCSAFTLibrary.h`,
  `IncompressibleLibrary.h`, `UNIFACLibrary.h`) live under `src/` and are
  **not installed**.
- **RapidJSON types leak into installed headers today:** `Configuration.h`,
  `SchemaValidation.h`, `fluids/Helmholtz.h`, `fluids/Ancillaries.h`,
  `fluids/PCSAFTFluid.h`, and `detail/rapidjson.h` itself. The install rules
  (`CMakeLists.txt:626-631`) ship the **entire `include/**/*.h` tree** (only
  generated `*_JSON.h`, `gitrevision.h`, `cpversion.h` are excluded), so
  `detail/` is currently public.
- Schema validation uses `rapidjson::SchemaDocument` / `SchemaValidator` in
  `src/SchemaValidation.cpp`, invoked by the PC-SAFT and cubic loaders to
  validate **user-supplied** fluid JSON. Schema source files:
  `dev/cubics/cubic_fluids_schema.json`,
  `dev/pcsaft/pcsaft_fluids_schema.json`,
  `dev/mixtures/mixture_departure_functions_schema.json`.
- Major ingestion paths: Helmholtz (`all_fluids.json`, ~20 MB, embedded as
  `all_fluids_JSON_z.h`, parsed on first library use), Incompressible, PC-SAFT,
  Cubics/UNIFAC, plus SVDSBTL surface-option validation and Configuration
  round-trip.
- A single-header nlohmann placeholder exists at `externals/nlohmann-json/`
  (added to dodge the heavy full checkout) with its include path at
  `CMakeLists.txt:331`. It is otherwise unused.

## 4. Architecture — Keeping Symbols Out of the Binary

Three cheap, mostly-mechanical layers. No standing wrapper library.

1. **Do not expose JSON types in installed headers.** Every class or function
   whose signature mentions a JSON node is handled one of two ways:
   - **Internal parsing classes** (the JSON-node-taking constructors such as
     `SaturationAncillaryFunction(json&)`, `SurfaceTensionCorrelation(json&)`,
     the `to_json(...)` methods, `PCSAFTFluid(...)`): their declarations move
     into `src/` (not installed) or receive an install `EXCLUDE`. These were
     installed by accident; they are implementation detail. nlohmann may appear
     in their signatures because the header never ships.
   - **Genuinely public surface** (`Configuration`, `SchemaValidation`):
     converted to **string in / string out** — JSON text is parsed internally
     and never crosses the boundary as a typed node.

   Net effect: no JSON type — rapidjson or nlohmann — appears in any installed
   header.

2. **Custom nlohmann namespace via macro.** Define the
   `NLOHMANN_JSON_NAMESPACE_*` macros (as a compile definition on the CoolProp
   targets, set through the CPM package options — *not* by editing a vendored
   header) so CoolProp's nlohmann lives in a private inline namespace (e.g.
   `nlohmann_coolprop`). Even an accidental instantiation cannot ODR-collide
   with another library's `nlohmann::json`. This is the belt to the
   "don't expose" suspenders.

3. **Localized hidden visibility around the JSON includes.** Wrap the
   nlohmann/Valijson `#include`s in `#pragma GCC visibility push(hidden)` /
   `pop` (equivalently, compile only the JSON-heavy translation units with
   hidden visibility) so that *only* the JSON-library symbols become local while
   the rest of CoolProp keeps its current default-visibility export behavior.
   Internally-instantiated JSON template symbols stay local rather than
   exported.

   **Why localized, not global `-fvisibility=hidden`.** CoolProp does **not**
   curate its exports today: there is no linker version script and no
   visibility annotations. (The `exports.txt` / `headers.txt` outputs at
   `CMakeLists.txt:683-701` are MSVC-only `dumpbin` *dumps* for inspection —
   they do not filter exports.) On Unix/macOS every symbol with external
   linkage is exported by default visibility, which is exactly why RapidJSON's
   weak/inline symbols leak. Flipping `-fvisibility=hidden` *globally* would
   therefore hide the **entire un-annotated public API** (the C API in
   `CoolPropLib.h` and the whole C++ `AbstractState` surface), requiring a
   `visibility("default")` export macro on every public symbol — a large,
   separate ABI project. Localized hiding gets the no-leak result without that
   overhaul. A proper project-wide export map (global `-fvisibility=hidden` +
   annotated public API, or a version script / `-exported_symbols_list`) is
   **deferred** — see Section 9.

   **Platform note.** This is a Unix/macOS concern. MSVC is hidden-by-default
   and only exports `dllexport`ed symbols, so nlohmann is already not exported
   on Windows; the namespace macro (layer 2) is the defense-in-depth there.

   **Exception-translation becomes load-bearing.** With JSON typeinfo hidden,
   an nlohmann exception that propagated across a shared-library boundary could
   not be caught by type. Every nlohmann/Valijson exception must therefore be
   caught and translated to CoolProp's own `ValueError` (or similar) *before it
   leaves a `.cpp`*. This is required for correctness, not just API hygiene, and
   is covered by an explicit test (Section 7).

The existing `cpjson` helpers are reimplemented as thin nlohmann one-liners and
kept **only** where they aid readability. They are convenience, not a required
abstraction, and they live in a non-installed `src/` header.

## 5. Schema Validation

Two tiers, splitting the work between build time and runtime.

- **Build-time (the real gate).** A Python step validates every `dev/*.json`
  fluid file against its schema using the mature `jsonschema` package, wired
  into the data-generation tooling and CI. Malformed *built-in* data is caught
  here, before embedding. Embedded blobs are then trusted at runtime and parsed
  without re-validation.

- **Runtime (user-supplied fluids only).** `add_fluids_from_JSON_string` for
  PC-SAFT and cubics keeps real schema validation via **Valijson**
  (`tristanpenman/valijson`): header-only, validates an `nlohmann::json`
  instance directly through its bundled nlohmann adapter, fetched via CPM and
  symbol-hidden the same way nlohmann is. Draft-7 is sufficient for CoolProp's
  schemas (required keys, types, array shapes, enums).

- The `rapidjson::SchemaValidator` path in `SchemaValidation.cpp` is replaced by
  Valijson. The public `SchemaValidation.h` API becomes **string-based**
  (`validate_json_against_schema(std::string instance, std::string schema)`), so
  no JSON node crosses the ABI. Schema source files are unchanged; their
  embedded `*_JSON.h` forms remain.

**Validator choice rationale.** pboettch/json-schema-validator (the current
"not super robust" candidate) is not header-only. jsoncons is the most
spec-complete (through Draft 2020-12) but uses its own parallel JSON DOM, which
would reintroduce the "two JSON libraries in one binary" problem this migration
exists to kill. Valijson validates nlohmann instances directly and is
header-only, making it the cleanest fit; Draft-7 capability is adequate here.

## 6. Data Loading & Parse Performance

nlohmann parses more slowly than RapidJSON (commonly 3–8× on large documents),
and `all_fluids.json` (~20 MB) is parsed on first Helmholtz-library use.

Approach: **straight swap then benchmark.** Each loader's parse path
(`Document::Parse` → `nlohmann::json::parse`) is converted directly. A focused
micro-benchmark times first-load of `all_fluids.json` before/after and is
reported in the relevant migration PR. If the regression is small, accept it.
Only if it is material do we reach for mitigations (nlohmann parse tuning, or a
pre-parsed CBOR/MessagePack embedded blob that nlohmann reads natively) — those
are deferred, not designed in up front. No pre-optimization.

`cpjson::to_canonical_json` (deterministic key-sorted serialization) is
reimplemented on nlohmann's ordered/sorted serialization primitives, preserving
byte-stable output.

## 7. Testing & Verification

- **Symbol-leak gate (headline invariant).** A post-build check runs
  `nm -D` / `objdump -T` and **fails the build** if any exported symbol matches
  `nlohmann` or `rapidjson`. The check must target a **shared product** (the
  built shared library, and ideally the Python extension `.so`) — *not* the
  static `.a`, because visibility attributes in a `.a` only take effect once it
  is linked into a shared object, so `nm` on the archive would show JSON symbols
  regardless and give a false reading. Active in CI from Phase 0 onward — the
  no-leak constraint becomes an enforced invariant, not a hope.
- **Behavior parity.** Existing Catch2 suites — the `[SBTL]` umbrella,
  fluid-loading, configuration round-trip, PC-SAFT/cubic user-fluid add — must
  stay green at every phase. Serializer round-trip and canonical-JSON tests
  guard byte-stable output.
- **Schema-error parity.** Add a test that feeds a deliberately-malformed user
  fluid and asserts a useful Valijson error message, so error quality does not
  silently regress.
- **No escaping JSON exceptions.** Add a test asserting that malformed input at
  the public boundary surfaces a CoolProp exception type (e.g. `ValueError`),
  never a raw `nlohmann`/Valijson exception — guarding the catch-and-translate
  contract that localized hidden visibility makes load-bearing (Section 4).
- **Per-phase preflight.** `./dev/ci/preflight.sh` plus the
  `superpowers:code-reviewer` subagent before every PR, per CLAUDE.md. When
  changes touch `src/SBTL/`, `include/CoolProp/sbtl/`,
  `src/Backends/SVDSBTL/`, or `src/Region/`, run the `[SBTL]` umbrella tag.

## 8. Rollout Plan

Staged. Each PR keeps the tree green and bisectable.

**Phase 0 — Seam & scaffolding (additive, nothing migrated):**
- Add nlohmann/json and Valijson via `CPMAddPackage` in
  `cmake/dependencies.cmake`, each pinned to a tag.
- Apply the `NLOHMANN_JSON_NAMESPACE_*` custom-namespace compile definition and
  the localized `#pragma GCC visibility push(hidden)` wrapping around the
  JSON-library includes (Section 4) — not a global `-fvisibility=hidden`.
- Delete `externals/nlohmann-json/` and remove its include path at
  `CMakeLists.txt:331`.
- Wire the symbol-leak CI check.
- Add the Python build-time schema validator to the data-generation tooling.
- Stand up the new `cpjson`-on-nlohmann helper header under `src/`.
- RapidJSON remains present and in use; build stays green.

**Phases 1…N — One loader per PR (dependency order):**
Helmholtz → Incompressible → PC-SAFT → Cubics/UNIFAC → SVDSBTL →
Configuration. Each PR converts one backend's `.cpp`, de-publishes/relocates
its JSON-typed headers (into `src/` or via install `EXCLUDE`), converts any
public surface to strings, and runs the relevant backend / `[SBTL]` tags.

**Phase Final — Delete RapidJSON:**
- Remove the rapidjson `CPMAddPackage` block from `cmake/dependencies.cmake`
  and its include dir from `CMakeLists.txt:336`.
- Delete `include/CoolProp/detail/rapidjson.h` and the `rapidjson_include.h`
  deprecation shim.
- Confirm the symbol-leak check is green and grep the tree for any surviving
  `rapidjson::`.

## 9. Out of Scope / Deferred

- Parse-performance mitigations (binary embedded format, parse tuning) — only
  if benchmarking shows a material regression.
- Migrating wrapper-side JSON (e.g. the Rust/Tauri GUI's `serde_json`) — those
  do not use the C++ RapidJSON and are unaffected.
- Adopting JSON Schema drafts beyond Draft-7.
- **Project-wide export map.** A curated public-export surface for CoolProp
  (global `-fvisibility=hidden` + `visibility("default")` annotations on the C
  and C++ public API, or a linker version script / macOS
  `-exported_symbols_list`). This would supersede the localized hiding in
  Section 4 and harden the ABI generally, but it is a large standalone ABI
  project independent of the JSON migration and is not required to meet the
  no-leak constraint here.
