# RapidJSON → nlohmann/json — Phase 2-4 (Incompressible + PC-SAFT + Cubics/UNIFAC loaders) — Design

**Date:** 2026-06-06
**Status:** Approved design, ready for implementation planning
**Epic:** `CoolProp-xa8w` (RapidJSON → nlohmann/json migration) · this is child `CoolProp-xa8w.1`
**Parent design:** `docs/superpowers/specs/2026-06-04-rapidjson-to-nlohmann-migration-design.md`
**Predecessors (merged):** Phase 0 scaffolding (PR #3084), Phase 1 Helmholtz loader + CBOR (PR #3086).

---

## 1. Goal

Convert the three remaining JSON-string-fed backend loaders — **Incompressible**,
**PC-SAFT**, and **Cubics/UNIFAC** — off RapidJSON onto nlohmann/json, in one PR.
After this lands, RapidJSON survives only in SVDSBTL / Configuration /
SchemaValidation (Phase 5) ahead of the Final deletion. No nlohmann (or
rapidjson) type may appear in any installed header, and no nlohmann/valijson
symbol may be exported from `libCoolProp`.

These loaders embed **plain, uncompressed `*_JSON.h` string headers** (no zlib),
all small (incompressibles ~175 KB, cubic/pcsaft smaller). So unlike Phase 1's
all_fluids, **there is no CBOR conversion** — just swap the parser. The data
stays human-readable in the embedded header; nlohmann parses ~175 KB in well
under a millisecond.

## 2. Why the conversion is mechanical

Every `cpjson::get_*` helper has an **identical signature** in `detail/json.h`
(nlohmann) and `detail/rapidjson.h`. So each loader is: swap the include, retype
parse-entry parameters, translate the small set of direct rapidjson API calls,
and the `cpjson::get_*` call sites are unchanged. The idiom map (same as Phase 1):

| rapidjson | nlohmann |
|-----------|----------|
| `rapidjson::Value& x` (param) | `const nlohmann::json& x` |
| `v.HasMember("k")` | `v.contains("k")` |
| `v["k"]` (read) | `v.at("k")` |
| `v["k"].GetString()/GetInt()/GetDouble()` | `cpjson::get_string/get_integer/get_double(v,"k")` (or `v.at("k").get<T>()`) |
| `v.IsObject()/IsArray()` | `v.is_object()/is_array()` |
| `for (auto itr=v.Begin(); itr!=v.End(); ++itr) { Value& c=*itr; … }` | `for (const auto& c : v) { … }` |
| `rapidjson::Document dd; dd.Parse<0>(s.data(), s.size()); if (dd.HasParseError()) …` | `nlohmann::json dd = cpjson::parse(s);` (throws `ValueError` on bad JSON) |
| `cpjson::json2string(v)` | `v.dump()` |

**Optional-member discipline:** rapidjson `v["k"]` (UB on missing) → nlohmann
`v.at("k")` (throws). Every optional key must keep its `contains("k")` guard, as
in the rapidjson original — verified per method, with the behavior-parity suite
as the net.

## 3. Schema validation → Valijson (the casts delete themselves)

PC-SAFT (`PCSAFTLibrary.cpp:36-38`) and Cubics (`CubicsLibrary.cpp:154-156`)
validate user-supplied fluids against their embedded schemas. Today they call
`cpjson::validate_schema` through an explicit **function-pointer cast** to
disambiguate the rapidjson and nlohmann overloads (a Phase-1 bridge, because the
loader TU saw both). Once these loaders stop including `detail/rapidjson.h`, only
the nlohmann (Valijson) `cpjson::validate_schema` overload is visible — so the
cast is **deleted**, and validation runs on Valijson with no other change (same
enum, same signature, same `schema_validation_code` return).

## 4. Delete the Phase-1 temporary bridges

Phase 1 left removable bridges in these files (grep "Removed once" /
"migrates off rapidjson"); they delete cleanly once the loaders are nlohmann-native:

| File:line | Bridge |
|-----------|--------|
| `CubicsLibrary.cpp:50-57` | `alpha0` rapidjson→`json2string`→`cpjson::parse`→nlohmann round-trip feeding `JSONFluidLibrary::parse_alpha0` |
| `UNIFACLibrary.cpp:63-70` | same `alpha0` round-trip |
| `CubicsLibrary.cpp:154-156` | `validate_schema` function-pointer cast |
| `PCSAFTLibrary.cpp:36-37` | `validate_schema` function-pointer cast |

After conversion, the cubic/UNIFAC loaders pass their (now nlohmann) `alpha0`
node directly to `JSONFluidLibrary::parse_alpha0(const nlohmann::json&)`.

## 5. The de-publishing boundary — POD constructor, zero nlohmann in installed headers

`include/CoolProp/fluids/PCSAFTFluid.h` (installed) has a public constructor
`PCSAFTFluid(rapidjson::Value::ValueIterator)` and `#include
"CoolProp/detail/rapidjson.h"`. Its fields are **protected** (getters only).

**Decision: POD-constructor boundary (not `json_fwd`+friend).** Add to the
installed header a plain constructor taking already-parsed values —
`PCSAFTFluid(std::string name, std::string CAS, CoolPropDbl molemass,
std::vector<std::string> aliases, PCSAFTValues params)` — all POD/std types
(`PCSAFTValues` is already a public struct in this header). The nlohmann parsing
lives **entirely** in a new non-installed `src/` factory
`cpjson::make_pcsaft_fluid(const nlohmann::json&)` (under
`src/Backends/PCSAFT/`), which extracts the fields and calls the POD ctor. The
installed header **drops the rapidjson include and gains no nlohmann at all —
not even `<nlohmann/json_fwd.hpp>`**.

**Why POD and not `json_fwd`+friend:** the `json_fwd` pattern (used for Phase 1's
`SaturationAncillaryFunction`) leaves `#include <nlohmann/json_fwd.hpp>` in the
installed header. CoolProp's install ships only `include/CoolProp/**` — it does
**not** ship nlohmann's headers (nlohmann is a build-only CPM dep) — so any
downstream consumer whose TU reaches that header would need `json_fwd.hpp` on its
own include path to compile. The POD boundary makes the installed header
self-contained: **downstream consumers need nothing from nlohmann.** (Not the
full library, not even the forward-decl header.)

**Retrofit `Ancillaries.h` to match.** Phase 1's
`include/CoolProp/fluids/Ancillaries.h` currently `#include
<nlohmann/json_fwd.hpp>` and friends `cpjson::make_saturation_ancillary(const
nlohmann::json&)`. It is reachable from `include/CoolProp/CoolPropFluid.h`, so it
imposes the same latent `json_fwd` dependency on downstream consumers. This phase
**converts `SaturationAncillaryFunction` to the same POD-constructor boundary**:
a POD/Eigen-typed constructor in `Ancillaries.h` (no nlohmann include, no friend),
with `cpjson::make_saturation_ancillary` (in the existing non-installed
`FluidLibraryFactories.h`) parsing and calling it. Net: the installed
fluid-model headers `Ancillaries.h` and `PCSAFTFluid.h` are nlohmann-free.

Incompressible and Cubics/UNIFAC need **no** de-publishing: `IncompressibleFluid`
constructs internally via setters (no JSON ctor in any installed header), and
Cubics/UNIFAC have no per-fluid public JSON constructors (string-based interface).

## 6. Components & Boundaries

- **`src/Backends/Incompressible/IncompressibleLibrary.{h,cpp}`** — swap include
  to `detail/json.h`; `parse_coefficients`/`parse_value`/`parse_ifrac`/`add_many`
  take `const nlohmann::json&`; `load_incompressible_library()` uses
  `cpjson::parse`.
- **`src/Backends/PCSAFT/PCSAFTLibrary.{h,cpp}` + `PCSAFTFluid.cpp`** — swap
  include; `add_many`/binary-pairs loops → nlohmann; `add_fluids_from_JSON_string`
  uses `cpjson::parse` + Valijson `validate_schema` (drop the cast); construct
  fluids via `cpjson::make_pcsaft_fluid`.
- **`include/CoolProp/fluids/PCSAFTFluid.h`** (installed) — drop rapidjson
  include; replace the JSON ctor with the POD ctor; **no nlohmann**.
- **New `src/Backends/PCSAFT/PCSAFTFluidFactory.h`** (non-installed) —
  `cpjson::make_pcsaft_fluid(const nlohmann::json&)` (the old ctor body, idiom-mapped).
- **`src/Backends/Cubics/CubicsLibrary.cpp` + `UNIFACLibrary.{h,cpp}`** — swap
  includes; `add_many`/`populate`/`jsonize` → nlohmann (`cpjson::parse`); drop the
  2 alpha0 bridges + the cast; Valijson validation.
- **`include/CoolProp/fluids/Ancillaries.h` + `Ancillaries.cpp` +
  `src/Backends/Helmholtz/Fluids/FluidLibraryFactories.h`** — retrofit
  `SaturationAncillaryFunction` to the POD-constructor boundary (remove `json_fwd`
  + friend; factory calls the POD ctor).

## 7. Sequencing (one PR, four commits)

Each commit keeps the tree green and bisectable:
1. **Incompressible loader** — simplest; no schema/bridges/factory. Validates the pattern.
2. **Ancillaries.h POD retrofit** — small, self-contained; removes the Phase-1
   `json_fwd` downstream dependency. (Independent of the loaders; sequenced early
   so the de-publishing pattern is proven before PCSAFTFluid reuses it.)
3. **PC-SAFT loader** — parser swap + `PCSAFTFluid` POD ctor + `make_pcsaft_fluid`
   factory + Valijson validation + drop its cast.
4. **Cubics/UNIFAC loader** — parser swap + Valijson validation + delete the 2
   alpha0 bridges + the cast.

## 8. Testing & Verification

- **Behavior parity (the safety net):** the existing `[INCOMP]`, `[cubic*]`,
  `[PCSAFT]`, `[UNIFAC]` Catch2 suites must stay green at every commit — every
  fluid loads with identical numbers.
- **User-fluid validation test (new):** the schema-validation path for
  user-supplied fluids moves engine (rapidjson → Valijson), so add a Catch2 test
  feeding a deliberately-malformed PC-SAFT and a malformed cubic payload to the
  public `add_fluids_as_JSON` entry points and asserting a clean failure with a
  useful (Valijson) error message — guarding error quality across the engine swap.
- **Symbol-leak gate:** the shared-library gate stays green — no `nlohmann`/
  `valijson` symbol exported. The new `make_pcsaft_fluid` factory must not export
  (inline in Release; mark `CP_JSON_LOCAL` if any out-of-line method of a
  JSON-taking class would otherwise emit a symbol).
- **Installed-header check:** grep confirms `PCSAFTFluid.h` and `Ancillaries.h`
  contain no `nlohmann`/`rapidjson` token after the retrofit.
- **Grep:** the three loaders' files are rapidjson-free; the 4 named Phase-1
  bridges are gone.
- **Per-commit:** `./dev/ci/preflight.sh` + the strengthened adversarial review
  (CLAUDE.md), including the fail-open gate question.

## 9. Out of Scope / Deferred (tracked)

- **Consumer-compile end-to-end gate** — a dependent library compiling against
  the install prefix with **no** nlohmann/rapidjson on its include path. This can
  only pass once the *entire* reachable installed surface is self-contained:
  `Ancillaries.h` + `PCSAFTFluid.h` (this phase), **`superancillary.h`** (native
  nlohmann, `SuperAncillary(const nlohmann::json&)` public ctor — `CoolProp-rxxx`),
  and **`detail/rapidjson.h` deleted + `detail/` install-excluded** (Phase Final
  `CoolProp-xa8w.3`). Added as an acceptance item on `xa8w.3`. **Coordinate with
  the in-flight nanobind interface session** — that interface is a genuine
  downstream consumer of CoolProp's public headers and can serve as (or inform)
  this gate, rather than building a duplicate synthetic consumer project.
- Phase 5 (SVDSBTL + Configuration → string-based public surface) and Phase Final
  (delete RapidJSON, drop the dual-enum guard, getter-layer decision) — `xa8w.2`,
  `xa8w.3`.
