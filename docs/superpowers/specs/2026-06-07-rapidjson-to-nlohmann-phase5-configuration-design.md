# RapidJSON → nlohmann — Phase 5 (Configuration / SchemaValidation / SVDSBTL off rapidjson) — Design

**Date:** 2026-06-07
**Status:** Approved design, ready for implementation planning
**Epic:** `CoolProp-xa8w` (RapidJSON → nlohmann/json migration) · child `CoolProp-xa8w.2`
**Predecessors (merged):** Phase 0-4 (loaders, value-type boundaries, superancillary de-leak).

---

## 1. Goal

Move the **last rapidjson consumers** off RapidJSON: `Configuration`,
`SchemaValidation`, the SVDSBTL options path, `MixtureParameters`, and a REFPROP
test parse. Make the public `Configuration` and `SchemaValidation` installed
headers **string-based** (no rapidjson — or nlohmann — type in any installed
signature). After this phase, the only file that still references rapidjson is
`detail/rapidjson.h` itself (plus the deprecated `rapidjson_include.h` shim and
the enum-guard comment in `detail/json.h`) — all deleted in Phase Final
(`CoolProp-xa8w.3`).

The **payoff**: with `Configuration.h` no longer pulling `detail/rapidjson.h`
into the loader TUs, the two Phase-2-4 `validate_schema` function-pointer casts
(PCSAFT, Cubics) become unambiguous and are **deleted** — the proof the
migration has fully landed.

## 2. Key facts (verified)

- **The rapidjson-typed public functions have no real C++ callers.**
  - `Configuration::get_config_as_json(rapidjson::Document&)` /
    `set_config_json(rapidjson::Document&)` (Configuration.h:346/365): the only
    references are `%ignore` lines in `src/CoolProp.i:13-14`. → **removed outright.**
  - `SchemaValidation::validate_against_schema(const rapidjson::Value&, …)` /
    `to_canonical_json(const rapidjson::Value&)` (SchemaValidation.h:26/32/48):
    two internal callers — `SVDSBTLBackend.cpp:180-181` and the `[SchemaValidation]`
    test — both migrate to the **string-based** API. → **removed outright.**
- **The string-based public API already exists** and is the SWIG-exposed,
  supported surface: `get_config_as_json_string()`, `set_config_as_json_string(std::string)`,
  `validate_json_against_schema(std::string, std::string)`, `to_canonical_json_str(std::string)`.
- **`Configuration.h` is the transitive rapidjson source** for the loader TUs:
  `FluidLibrary.h` → `Configuration.h` → `detail/rapidjson.h`, which is why
  PCSAFTLibrary.cpp / CubicsLibrary.cpp see both `validate_schema` overloads and
  needed the cast. Dropping `detail/rapidjson.h` from `Configuration.h` removes
  that path.
- **`to_canonical_json` simplifies under nlohmann:** `nlohmann::json` is
  `std::map`-backed, so object keys serialize **sorted at every nesting level** —
  the canonical form is just `cpjson::parse(s).dump()` (no manual recursive sort).
- **SchemaValidation's validator** moves rapidjson `SchemaDocument`/`SchemaValidator`
  → **Valijson** (the same engine Phase 2-4 adopted for the loaders), reusing
  `cpjson::validate_schema` semantics.

## 3. Components

### Configuration (`include/CoolProp/Configuration.h`, `src/Configuration.cpp`)
- **Header:** delete `get_config_as_json(rapidjson::Document&)` and
  `set_config_json(rapidjson::Document&)`; drop `#include "CoolProp/detail/rapidjson.h"`.
  Retype `ConfigurationItem::add_to_json`/`set_from_json` to nlohmann **and move
  their bodies out of the installed header into the `.cpp`** (so no
  `nlohmann::json` type appears in the installed header — these become declared
  in the header only if needed, or fully `.cpp`-internal). Net: the installed
  `Configuration.h` carries **no JSON-library type** — only the `std::string`
  config API.
- **.cpp:** convert serialize (`add_to_json`/`get_config_as_json`) and
  deserialize (`set_from_json`/`set_config_as_json`) to nlohmann; keep the
  `*_as_json_string` public functions (now nlohmann-backed). Include
  `detail/json.h` (non-installed) in the `.cpp`.
- **`src/CoolProp.i`:** remove the two `%ignore` lines for the deleted rapidjson
  overloads.

### SchemaValidation (`include/CoolProp/SchemaValidation.h`, `src/SchemaValidation.cpp`)
- **Header:** delete the three rapidjson-typed functions; keep the string API
  (`validate_json_against_schema`, `to_canonical_json_str`); drop the rapidjson
  include. Installed header becomes JSON-library-free.
- **.cpp:** `validate_json_against_schema` → Valijson (via `cpjson::validate_schema`
  or a direct Valijson call); `to_canonical_json_str` → `cpjson::parse(s).dump()`
  (sorted-key canonical). Include `detail/json.h`.

### SVDSBTL (`src/Backends/SVDSBTL/SVDSBTLBackend.cpp`)
- Replace `validate_against_schema(opts /*rapidjson*/, schema)` +
  `to_canonical_json(opts)` with the **string** API on the options string from
  the factory. `resolve_grid(const rapidjson::Document&)` → takes
  `const nlohmann::json&` and reads `grid.NT/NR/rank` via nlohmann. Drop
  `#include "CoolProp/detail/rapidjson.h"`; include `detail/json.h`.

### MixtureParameters (`src/Backends/Helmholtz/MixtureParameters.cpp`)
- Convert the predefined-mixtures JSON parsing/iteration from rapidjson to
  nlohmann (`cpjson::parse` + range-for + `cpjson::get_*`), idiom-mapped as in
  Phase 2-4.

### REFPROP (`src/Backends/REFPROP/REFPROPMixtureBackend.cpp:2920`)
- Convert the single rapidjson parse of a fluid-param JSON string to nlohmann.

### The cast deletion (the payoff)
- `src/Backends/PCSAFT/PCSAFTLibrary.cpp` and
  `src/Backends/Cubics/CubicsLibrary.cpp`: delete the
  `static_cast<… (*)(std::string_view, …)>(&cpjson::validate_schema)` indirection
  and call `cpjson::validate_schema(...)` directly (now the only overload visible),
  removing the explanatory comments about the ambiguity.

### Tests
- `src/Tests/CoolProp-Tests-SchemaValidation.cpp`: drop its `rapidjson::Document
  parse()` helper; call `validate_json_against_schema(instance_str, schema_str)`
  and `to_canonical_json_str(str)`. Assertions unchanged in intent (valid passes,
  invalid throws, canonical is sorted/deterministic).

## 4. Sequencing (one PR, logical commits)

1. **SchemaValidation** → Valijson + string-only public API (+ test).
2. **Configuration** → nlohmann internal + string-only public API (+ `CoolProp.i`).
3. **SVDSBTL + MixtureParameters + REFPROP** → nlohmann.
4. **Delete the two `validate_schema` casts** (verify they compile away cleanly).

Each commit keeps the tree green.

## 5. Testing & verification

- **Behavior parity:** `[SchemaValidation]`, `[SVDSBTL]`/`[SVDSBTL][options]`,
  `[SBTL]` umbrella (per CLAUDE.md, SBTL changes run the umbrella tag), the
  mixture tests (predefined mixtures load identically), `[PCSAFT]`/`[cubic]`
  (the cast change is behaviour-neutral). Config get/set round-trip
  (`get_config_as_json_string`/`set_config_as_json_string`) unchanged.
- **Canonical-JSON equivalence:** confirm `to_canonical_json_str` produces
  sorted-key output matching the previous rapidjson canonical form (the test
  asserts determinism + sorting).
- **Installed-header check:** `Configuration.h` and `SchemaValidation.h` contain
  no `rapidjson`/`nlohmann` token.
- **Exhaustive grep:** after this phase,
  `git grep -l "rapidjson" include src` returns only `detail/rapidjson.h`,
  `rapidjson_include.h`, the `detail/json.h` enum-guard comment, and
  `PCSAFTFluidFactory.h`'s comment — nothing else.
- **Symbol gate** stays green; full `~[slow]` suite green.
- Per-commit `./dev/ci/preflight.sh` + the mandatory adversarial review.

## 6. Out of scope (Phase Final, `xa8w.3`)

Delete `detail/rapidjson.h` + `rapidjson_include.h`, drop the
`CPJSON_SCHEMA_VALIDATION_CODE_DEFINED` dual-enum guard (make the enum
unconditional in `detail/json.h`), tighten the symbol gate pattern to include
`rapidjson`, and the deferred getter-layer decision. Phase 5 leaves the tree in
the state where Phase Final is a near-mechanical deletion.
