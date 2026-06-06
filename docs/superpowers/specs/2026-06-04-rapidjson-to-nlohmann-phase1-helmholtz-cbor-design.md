# RapidJSON → nlohmann/json — Phase 1 (Helmholtz loader + CBOR) — Design

**Date:** 2026-06-04
**Status:** Approved design, ready for implementation planning
**Parent design:** `docs/superpowers/specs/2026-06-04-rapidjson-to-nlohmann-migration-design.md`
**Phase 0 (merged):** PR #3084 — scaffolding (CPM nlohmann+valijson, hidden-visibility `cpjson` wrapper, symbol-leak gate, Python schema validation).

---

## 1. Goal

Migrate the **Helmholtz fluid loader** off RapidJSON and onto nlohmann/json, and
embed the `all_fluids` database as **uncompressed CBOR** (dropping the
zlib/miniz step from this load path). This is the first per-loader phase; it
also establishes the conversion pattern every later loader phase will follow.

Success = RapidJSON fully removed from the Helmholtz load path; no nlohmann type
in any installed header; every fluid loads identically (behavior parity);
first-load latency at or below today's; the embedded blob is provably
byte-equivalent to the source JSON.

## 2. Why the loader swap forces the public-header change

`FluidLibrary` constructs fluid-model objects directly from JSON nodes, e.g.
`SaturationAncillaryFunction(ancillaries[...])` and
`SurfaceTensionCorrelation(surface_tension)` (`FluidLibrary.h` `parse_ancillaries`
/ `parse_surface_tension`). Today those nodes are `rapidjson::Value&`. The moment
the loader parses with nlohmann, the nodes become `nlohmann::json`, so those
constructor signatures must change — and `Ancillaries.h` / `Helmholtz.h` are
**installed** headers (CMake ships all `include/CoolProp/**/*.h` except
`*_JSON*.h`). Taking `nlohmann::json` there would re-leak exactly the symbols
Phase 0 hid. So Phase 1 inescapably bundles the loader swap with de-publishing
the JSON-typed public members.

`FluidLibrary.h`/`.cpp` themselves live under `src/Backends/Helmholtz/Fluids/`
(NOT installed), so the ~23 `parse_*` methods may take `nlohmann::json` directly
with no de-publishing concern.

## 3. Decisions

> **Implementation note (shipped, supersedes the `cbor2` references below):** the
> encoder was changed during implementation from `cbor2` to a vendored stdlib-only
> CBOR codec, `dev/cbor_min.py`. `cbor2` is a compiled (Rust) package, and getting
> it installed into the exact interpreter CMake invokes across CoolProp's ~12
> heterogeneous build workflows (manylinux/emscripten/MSVC/conda) proved
> intractable. `cbor_min` needs no third-party dependency on the build path. A CI
> job (`dev/check_cbor_min_vs_cbor2.py`, wired into the `test_catch2` hard gate)
> asserts `cbor_min` is **byte-identical** to `cbor2`, so the wording below stands
> with `cbor2` → `cbor_min` substituted on the build path.

| Area | Decision |
|------|----------|
| **`parse_*` conversion** | Convert the ~23 `parse_*` methods to `const nlohmann::json&`, using the **nlohmann `cpjson` getters** (`cpjson::get_double(j,"x")`, `get_long_double_array(j,"n")`, …) — same names as the rapidjson ones, so it's a near-mechanical include-swap + member-access/iterator tweaks. Use native nlohmann only where clearly cleaner. |
| **JSON-taking constructors** | De-publish as **free-function factories** in non-installed detail headers, e.g. `cpjson::make_saturation_ancillary(const nlohmann::json&) -> SaturationAncillaryFunction`. The classes keep their non-JSON public constructors and stay in the installed headers. nlohmann appears only in the non-installed factory headers the loader includes. |
| **`to_json` methods** | **Drop** all ~13 `to_json(rapidjson::Value&, rapidjson::Document&)` and the single debug-only caller (`HelmholtzEOSMixtureBackend.cpp` `_calc_all_critical_points`, `if (get_debug_level() > 10)` branch, ~lines 4290-4299). Their output was never a stable/public format. |
| **Embedded format** | `all_fluids` becomes **uncompressed CBOR**. zlib/miniz dropped from this path. |
| **Encoder** | `dev/generate_headers.py` emits CBOR via `cbor2`, with a build-time round-trip self-check (decode own output, compare to source) before writing the header. |
| **Loader read** | `nlohmann::json::from_cbor(embedded_bytes)` replaces `uncompress()` + `rapidjson::Document::Parse`. |
| **Conversion idiom (long-term)** | Hybrid: getters now through all loader phases; an explicit Phase-Final task decides getters-stay vs inline-to-native (with `validate_schema`/`parse` kept regardless). |

## 4. Data Flow

**Build time:** `dev/fluids/*.json` → `generate_headers.py::combine_json` combines
into one array → `cbor2.dumps(array)` → `dev/all_fluids.cbor` → hex byte array →
`include/all_fluids_CBOR.h` (incbin template, same mechanism as the current
`all_fluids_JSON_z.h`). The generator **round-trips its own output**
(`cbor2.loads(blob) == array`) and aborts on mismatch before writing.

**Runtime (`FluidLibrary::load()`):**
```
incbin bytes (all_fluids_CBOR.h)
  → nlohmann::json doc = cpjson::from_cbor(bytes)   // new helper, exception-translated
  → add_many(doc)
  → parse_*(doc[...])                                // nlohmann walks, cpjson getters
```
No decompress step. The `parse_EOS` `cpjson::json2string(EOS_json["SUPERANCILLARY"])`
rapidjson→string bridge is replaced by passing the `nlohmann::json` subobject (or
`.dump()`) directly to the superancillary system (already nlohmann-based).

## 5. Components & Boundaries

- **`include/CoolProp/detail/json.h`** (existing, non-installed wrapper): add
  `cpjson::from_cbor(const std::vector<std::uint8_t>&) -> nlohmann::json` (and/or
  a span/pointer overload for the incbin symbol), exception-translated to
  `CoolProp::ValueError` like `cpjson::parse`.
- **New non-installed factory header(s)** (under `src/Backends/Helmholtz/Fluids/`
  or `include/CoolProp/detail/` with an install `EXCLUDE`): the
  `cpjson::make_*` factories that take `nlohmann::json` and return the
  ancillary/correlation objects. The loader includes these; nothing installed does.
- **`Ancillaries.h` / `Helmholtz.h`** (installed): lose the JSON-taking
  constructors and all `to_json` methods; stop including `detail/rapidjson.h`.
  They keep the class definitions and non-JSON constructors.
- **`FluidLibrary.h`/`.cpp`** (not installed): `parse_*` take `nlohmann::json`;
  include `detail/json.h` + the factory header instead of `detail/rapidjson.h`.
- **`dev/generate_headers.py`**: CBOR emission + self-check; drop the
  `all_fluids.json.z` zlib path (the `zvalues` entry and the `zlib.compress` line).
- **`include/all_fluids_CBOR.h`** (generated, must NOT ship): replaces
  `all_fluids_JSON_z.h`. **Note:** the existing install rules exclude
  `*_JSON*.h`, which `all_fluids_CBOR.h` does NOT match — so Phase 1 must add a
  `PATTERN "*_CBOR*.h" EXCLUDE` (or equivalent) to the three `install(DIRECTORY
  include/...)` rules in `CMakeLists.txt`, or the generated CBOR header will be
  shipped. (A `[cbor]`-or-install check should guard this.)

## 6. Sequencing (one PR, two commits)

Isolate the large mechanical change from the format change; each independently
verifiable; no shipped slow intermediate (both land together).

- **Commit A — parser swap, format unchanged.** Loader still inflates zlib→JSON,
  but parses with **nlohmann** (`nlohmann::json::parse`) instead of rapidjson; all
  `parse_*` converted; constructors de-published to factories; `to_json` + debug
  branch dropped; Helmholtz/Ancillaries stop including `detail/rapidjson.h`.
  RapidJSON is gone from the Helmholtz path. **Verify:** full fluid-load test
  suite green (behavior parity); first-load benchmark (the in-context ~+70 ms
  naive point).
- **Commit B — format swap.** `generate_headers.py` emits CBOR; loader uses
  `cpjson::from_cbor`; `all_fluids_JSON_z.h`/`.json.z` and the zlib path removed.
  **Verify:** `[cbor]` byte-equivalence test; first-load benchmark (the ~35 ms
  point); blob-size check.

## 7. Testing & Verification

- **Behavior parity (the primary safety net for the 23-method conversion):** the
  existing Helmholtz fluid-loading suite must stay green at commit A — every
  fluid loads with identical numbers. Tag scope per CLAUDE.md.
- **`[cbor]` byte-equivalence gate:** a Catch2 test decodes the embedded CBOR
  blob (`gall_fluids_CBOR*` symbol, available in tests via `COOLPROP_NO_INCBIN`)
  into `nlohmann::json` and asserts value-equality against
  `nlohmann::json::parse(dev/all_fluids.json)`. The source path is provided via a
  **CMake-configured path macro** (robust against CWD), not a relative guess.
  This catches encoder/decoder drift and source-vs-blob staleness.
- **Generator self-check:** build-time round-trip assertion in
  `generate_headers.py` (decode own CBOR, compare to the source dict).
- **Symbol-leak gate (now meaningful):** Phase 1 is the first time nlohmann is
  compiled INTO `libCoolProp` (Phase 0 only used it in tests). The gate confirms
  no `nlohmann`/`valijson` symbols are exported from the shared library.
- **In-context first-load benchmark** at both commits, reported in the PR
  (decompress+parse at A; from_cbor at B), validating the spec §6 projection
  (~37 ms today → ~107 ms at A → ~35 ms at B).
- **Per-commit preflight** + the strengthened adversarial review (CLAUDE.md
  "Pre-PR adversarial review", including the fail-open gate question and
  re-review-the-delta rule).

## 8. Out of Scope / Deferred (tracked)

- **Phase Final:** delete `detail/rapidjson.h` + RapidJSON (CPM block, include
  dir), add the `detail/` install EXCLUDE, tighten the symbol-gate pattern to
  include `rapidjson`.
- **Phase-Final getter-layer decision:** keep the nlohmann `cpjson` getter
  conveniences vs inline them to native nlohmann (keeping `validate_schema`/
  `parse`). New tracked task; decided with the full call-site surface visible.
- Other loaders (Incompressible, PC-SAFT, Cubics/UNIFAC, SVDSBTL, Configuration)
  — their own later phases; they follow the pattern this phase establishes.
- `FluidLibrary::get_JSONstring` (Helmholtz): appears to have no live caller;
  whether to delete it is a separate cleanup, not required here unless it blocks
  the rapidjson removal from this file (in which case it converts with the rest).
