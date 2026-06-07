# RapidJSON → nlohmann — Phase 5 (Configuration / SchemaValidation / SVDSBTL) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the last rapidjson consumers (Configuration, SchemaValidation, SVDSBTL options, MixtureParameters, a REFPROP test parse) onto nlohmann/Valijson; make `Configuration.h`/`SchemaValidation.h` JSON-library-free; delete the two Phase-2-4 `validate_schema` casts.

**Architecture:** The rapidjson-typed public functions are removed outright (no C++ callers). SchemaValidation re-backs its string API with Valijson (`cpjson::validate_schema`) + nlohmann canonicalization (`cpjson::parse(s).dump()` — nlohmann's `std::map`-backed objects serialize keys sorted at every level). Configuration's per-item serialize/deserialize moves to non-installed `.cpp` free functions using `ConfigurationItem`'s accessors. The loaders' `validate_schema` casts become unambiguous once `Configuration.h` stops pulling `detail/rapidjson.h`.

**Tech Stack:** C++17, nlohmann/json + Valijson via `include/CoolProp/detail/json.h`'s `cpjson::`, Catch2 v3, CMake (Release).

**Spec:** `docs/superpowers/specs/2026-06-07-rapidjson-to-nlohmann-phase5-configuration-design.md`
**Epic:** `CoolProp-xa8w.2`. **Branch:** `ihb/json-phase5-configuration-off-rapidjson`.

---

## Idiom map (for the mechanical conversions — same as Phase 2-4)

| rapidjson | nlohmann |
|-----------|----------|
| `rapidjson::Document dd; dd.Parse<0>(s.c_str()); … HasParseError()` | `nlohmann::json dd = cpjson::parse(s);` (throws `ValueError`) |
| `v.IsObject()/IsArray()/IsBool()/IsInt()/IsDouble()/IsString()` | `v.is_object()/is_array()/is_boolean()/is_number_integer()/is_number()/is_string()` |
| `v.GetBool()/GetInt()/GetDouble()/GetString()` | `cpjson::get_*` or `v.get<bool>()/get<int>()/get<double>()/get<std::string>()` |
| `for (auto it=v.MemberBegin(); it!=v.MemberEnd(); ++it) { it->name.GetString(); it->value … }` | `for (auto& [name, value] : v.items()) { … }` |
| `for (auto it=v.Begin(); it!=v.End(); ++it) { *it }` | `for (const auto& el : v) { … }` |
| `val.AddMember(name, v, alloc)` | `obj[name] = v;` |
| `cpjson::to_string(doc)` / writer | `doc.dump()` |

---

## Task 1: SchemaValidation → Valijson + string-only API

**Files:**
- Modify: `include/CoolProp/SchemaValidation.h` (remove rapidjson include + the 3 rapidjson-typed decls)
- Rewrite: `src/SchemaValidation.cpp`
- Modify: `src/Tests/CoolProp-Tests-SchemaValidation.cpp`

- [ ] **Step 1: Header — string-only, JSON-library-free**

In `SchemaValidation.h`: delete `#include "CoolProp/detail/rapidjson.h"` (line 7); delete the declarations of `validate_against_schema(const rapidjson::Value&, const rapidjson::Document&)` (26), `validate_against_schema(const rapidjson::Value&, const std::string&)` (32), and `to_canonical_json(const rapidjson::Value&)` (48). Keep only `validate_json_against_schema(const std::string&, const std::string&)` and `to_canonical_json_str(const std::string&)`. Remove any now-unused `#if !defined(SWIG)` guards that only wrapped the deleted decls.

- [ ] **Step 2: Rewrite `src/SchemaValidation.cpp`**

Replace the entire file body with the collapsed nlohmann/Valijson version:
```cpp
#include "CoolProp/SchemaValidation.h"

#include <string>

#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/json.h"

namespace CoolProp {

void validate_json_against_schema(const std::string& instance_json, const std::string& schema_json) {
    std::string errstr;
    // cpjson::validate_schema(schemaJson, inputJson, errstr) — Valijson-backed.
    cpjson::schema_validation_code code = cpjson::validate_schema(schema_json, instance_json, errstr);
    if (code != cpjson::SCHEMA_VALIDATION_OK) {
        throw ValueError(std::string("schema validation failed: ") + errstr);
    }
}

std::string to_canonical_json_str(const std::string& json_str) {
    // nlohmann::json is std::map-backed, so object keys serialize sorted at every
    // level — this is the canonical (deterministic, key-sorted) form. Arrays keep order.
    return cpjson::parse(json_str).dump();
}

}  // namespace CoolProp
```
(Deletes `copy_sorted`, `pointer_to_string`, `parse_json`, and both rapidjson-typed overloads. NOTE the argument order: `cpjson::validate_schema(schema, instance, errstr)` — schema first.)

- [ ] **Step 3: Update the `[SchemaValidation]` test**

In `src/Tests/CoolProp-Tests-SchemaValidation.cpp`: delete the local `rapidjson::Document parse(const std::string&)` helper (lines ~52-56) and the rapidjson include. Replace each `validate_against_schema(parse(inst), std::string(kSchema))` with `validate_json_against_schema(inst, std::string(kSchema))`, and each `to_canonical_json(parse(s))` with `to_canonical_json_str(s)`. The `REQUIRE_NOTHROW`/`REQUIRE_THROWS` structure is unchanged. For canonical-output assertions: if the test compares an EXACT canonical string, update the expected value to nlohmann's `dump()` output (number formatting may differ from rapidjson) — run the test to capture the actual canonical string and confirm it is sorted/deterministic; if the test only asserts "two equivalent inputs produce equal output" (determinism), no expected-string change is needed.

- [ ] **Step 4: Build + test**

Run: `cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[SchemaValidation]" -s`
Expected: clean build; all sections pass (valid passes, invalid throws, canonical sorted/deterministic).

- [ ] **Step 5: Commit**
```bash
git add include/CoolProp/SchemaValidation.h src/SchemaValidation.cpp src/Tests/CoolProp-Tests-SchemaValidation.cpp
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit -m "refactor(json): SchemaValidation off rapidjson — Valijson + nlohmann, string-only API (CoolProp-xa8w.2)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Configuration → nlohmann internal, string-only public API

**Files:**
- Modify: `include/CoolProp/Configuration.h` (remove rapidjson include + the public `*_as_json(Document&)` decls + the `ConfigurationItem::add_to_json`/`set_from_json` methods)
- Modify: `src/Configuration.cpp` (nlohmann internals; de-published per-item free functions)
- Modify: `src/CoolProp.i` (remove the two `%ignore`s)

- [ ] **Step 1: Confirm `ConfigurationItem`'s public accessors**

Read `include/CoolProp/Configuration.h` around the `ConfigurationItem` class (the cast operators / getters used by `get_config_bool/integer/double/string` — e.g. `static_cast<double>(item)`, `get_key()`, the `type`). Confirm there is a public way to read the item's **type** and its typed value. If a public `ConfigurationDataTypes get_type() const { return type; }` accessor is absent, ADD one (small, non-JSON) — the de-published free functions need it.

- [ ] **Step 2: Header — remove the rapidjson surface**

In `Configuration.h`: delete `#include "CoolProp/detail/rapidjson.h"` (line 10); delete the public `void get_config_as_json(rapidjson::Document&)` (346) and `void set_config_json(rapidjson::Document&)` (365) declarations; delete the `ConfigurationItem::add_to_json(rapidjson::Value&, rapidjson::Document&)` (124-154) and `set_from_json(rapidjson::Value&)` (155-190) methods entirely (and their `#if !defined(SWIG)` wrapper if it now wraps nothing). Keep `get_config_as_json_string()` / `set_config_as_json_string(const std::string&)`. The installed `Configuration.h` must end up with **no `rapidjson`/`nlohmann` token**.

- [ ] **Step 3: `Configuration.cpp` — de-published per-item free functions (nlohmann)**

Add `#include "CoolProp/detail/json.h"` to the `.cpp`. Add file-local free functions reproducing the old methods' logic, idiom-mapped, using `ConfigurationItem`'s accessors:
```cpp
namespace {
void item_to_json(const ConfigurationItem& item, nlohmann::json& obj) {
    const std::string name = config_key_to_string(item.get_key());
    switch (item.get_type()) {
        case CONFIGURATION_BOOL_TYPE:    obj[name] = static_cast<bool>(item); break;
        case CONFIGURATION_INTEGER_TYPE: obj[name] = static_cast<int>(item); break;
        case CONFIGURATION_DOUBLE_TYPE:  obj[name] = static_cast<double>(item); break;
        case CONFIGURATION_STRING_TYPE:  obj[name] = static_cast<std::string>(item); break;
        case CONFIGURATION_ENDOFLIST_TYPE:
        case CONFIGURATION_NOT_DEFINED_TYPE: throw ValueError();
    }
}
void item_from_json(ConfigurationItem& item, const nlohmann::json& val) {
    switch (item.get_type()) {
        case CONFIGURATION_BOOL_TYPE:
            if (!val.is_boolean()) throw ValueError(format("Input is not boolean"));
            item.set_bool(val.get<bool>()); break;
        case CONFIGURATION_INTEGER_TYPE:
            if (!val.is_number_integer()) throw ValueError(format("Input is not integer"));
            item.set_integer(val.get<int>()); break;
        case CONFIGURATION_DOUBLE_TYPE:
            if (!val.is_number()) throw ValueError(format("Input [%s] is not double (or castable)", val.dump().c_str()));
            item.set_double(val.get<double>()); break;
        case CONFIGURATION_STRING_TYPE:
            if (!val.is_string()) throw ValueError(format("Input is not string"));
            item.set_string(val.get<std::string>()); break;
        case CONFIGURATION_ENDOFLIST_TYPE:
        case CONFIGURATION_NOT_DEFINED_TYPE: throw ValueError();
    }
}
}  // namespace
```
IMPORTANT: use whatever the ACTUAL `ConfigurationItem` typed getters/setters are named (the cast operators `static_cast<T>(item)` are confirmed used elsewhere; the setters may be `set_bool`/`set_integer`/… or assignment — read the class and match). `val.is_number()` covers the old `IsDouble() || IsInt()`.

- [ ] **Step 4: `Configuration.cpp` — string functions on nlohmann; delete the rapidjson ones**

Replace `get_config_as_json(rapidjson::Document&)` + `get_config_as_json_string()` (105-117) with:
```cpp
std::string get_config_as_json_string() {
    nlohmann::json doc = nlohmann::json::object();
    for (auto& kv : _get_config()->get_items()) {
        item_to_json(kv.second, doc);
    }
    return doc.dump();
}
```
Replace `set_config_as_json(rapidjson::Value&)` + `set_config_as_json_string(const std::string&)` (118-153) with:
```cpp
void set_config_as_json_string(const std::string& s) {
    nlohmann::json doc = cpjson::parse(s);
    // First check that all keys are valid
    for (auto& [name, value] : doc.items()) {
        try {
            _get_config()->get_item(config_string_to_key(name));
        } catch (std::exception& e) {
            throw ValueError(format("Unable to parse json file with error: %s", e.what()));
        }
    }
    // Now set the values
    for (auto& [name, value] : doc.items()) {
        ConfigurationItem& item = _get_config()->get_item(config_string_to_key(name));
        try {
            item_from_json(item, value);
        } catch (std::exception& e) {
            throw ValueError(format("Unable to parse json file with error: %s", e.what()));
        }
    }
}
```
(The public `get_config_as_json(Document&)`/`set_config_as_json(Value&)` are gone; the string functions are the only public JSON config API.)

- [ ] **Step 5: `src/CoolProp.i` — remove the dead `%ignore`s**

Delete lines 13-14:
```
%ignore CoolProp::set_config_json(rapidjson::Document &);
%ignore CoolProp::get_config_as_json(rapidjson::Document &);
```

- [ ] **Step 6: Build + config round-trip test**

Run: `cmake --build build_catch --target CatchTestRunner -j8`
Then exercise the config JSON round-trip (find the covering test; config is tested via backends + any `[configuration]`/`[config]` tag):
`./build_catch/CatchTestRunner "[configuration],[config]"` (if no such tag, run `"~[slow]"` later in Final Verification — the change is exercised broadly). Confirm `get_config_as_json_string`/`set_config_as_json_string` round-trip a config unchanged.
Then: `grep -n "rapidjson\|nlohmann" include/CoolProp/Configuration.h` → no matches.

- [ ] **Step 7: Commit**
```bash
git add include/CoolProp/Configuration.h src/Configuration.cpp src/CoolProp.i
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit -m "refactor(json): Configuration off rapidjson — string-only public API, nlohmann internal (CoolProp-xa8w.2)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: SVDSBTL + MixtureParameters + REFPROP → nlohmann

**Files:**
- Modify: `src/Backends/SVDSBTL/SVDSBTLBackend.cpp` (include; `resolve_grid`; the validate/canonical calls)
- Modify: `src/Backends/Helmholtz/MixtureParameters.cpp` (predefined-mixtures parse)
- Modify: `src/Backends/REFPROP/REFPROPMixtureBackend.cpp:~2920` (one parse)

- [ ] **Step 1: SVDSBTL**

In `SVDSBTLBackend.cpp`: replace `#include "CoolProp/detail/rapidjson.h"` (38) with `#include "CoolProp/detail/json.h"`. The options path currently parses the validated options string into a `rapidjson::Document opts`, then calls `validate_against_schema(opts, kSVDSBTLOptionsSchemaJson)` + `to_canonical_json(opts)` (lines ~180-181). Switch to the **string** API on the options string `opts_str`:
```cpp
CoolProp::validate_json_against_schema(opts_str, std::string(kSVDSBTLOptionsSchemaJson));
std::string canonical = CoolProp::to_canonical_json_str(opts_str);
```
Convert `resolve_grid(const rapidjson::Document& opts)` → `resolve_grid(const nlohmann::json& opts)`: `opts["grid"].IsObject()` → `opts.contains("grid") && opts.at("grid").is_object()`; `grid["NT"].GetInt()` → `cpjson::get_integer(grid, "NT")` (and NR/rank). Parse the (canonical or validated) options string with `cpjson::parse` where the code currently built the rapidjson Document, and pass the nlohmann node to `resolve_grid`. Read the surrounding factory/options code to wire `opts_str` through correctly.

- [ ] **Step 2: MixtureParameters**

In `MixtureParameters.cpp`: replace the rapidjson include with `detail/json.h`; convert the predefined-mixtures `rapidjson::Document`/`Value` parse + iteration (lines ~25, 34-39) to `cpjson::parse` + range-for + `cpjson::get_*`, per the idiom map. Read the functions and apply the same mechanical conversion used for the Phase 2-4 loaders (preserve every `HasMember`→`contains` guard).

- [ ] **Step 3: REFPROP**

In `REFPROPMixtureBackend.cpp` near line 2920: convert the single `rapidjson` parse of a fluid-param JSON string to `cpjson::parse(...)` + nlohmann access. (This is inside a test/diagnostic path — keep the assertion behavior identical.)

- [ ] **Step 4: Build + tests**

Run: `cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[SBTL]"`  (the SBTL **umbrella** per CLAUDE.md, not just `[SVDSBTL]`), then `./build_catch/CatchTestRunner "[SVDSBTL][options]"` and the mixture tests (e.g. `"[mixture]"` smoke). Expected: pass; SVDSBTL options validate/round-trip identically; predefined mixtures load identically.
Then: `grep -rn "rapidjson" src/Backends/SVDSBTL/ src/Backends/Helmholtz/MixtureParameters.cpp src/Backends/REFPROP/REFPROPMixtureBackend.cpp` → no matches.

- [ ] **Step 5: Commit**
```bash
git add src/Backends/SVDSBTL/SVDSBTLBackend.cpp src/Backends/Helmholtz/MixtureParameters.cpp src/Backends/REFPROP/REFPROPMixtureBackend.cpp
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit -m "refactor(json): SVDSBTL options + MixtureParameters + REFPROP parse off rapidjson (CoolProp-xa8w.2)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Delete the two `validate_schema` casts (the payoff)

**Files:**
- Modify: `src/Backends/PCSAFT/PCSAFTLibrary.cpp` (the cast in `add_fluids_from_JSON_string`)
- Modify: `src/Backends/Cubics/CubicsLibrary.cpp` (the cast in `add_fluids_from_JSON_string`)

- [ ] **Step 1: PCSAFT — call `validate_schema` directly**

In `PCSAFTLibrary.cpp`, replace the cast block:
```cpp
    auto validate_schema_nlohmann =
      static_cast<cpjson::schema_validation_code (*)(std::string_view, std::string_view, std::string&)>(&cpjson::validate_schema);
    cpjson::schema_validation_code val_code = validate_schema_nlohmann(pcsaft_fluids_schema_JSON, JSON, errstr);
```
with:
```cpp
    cpjson::schema_validation_code val_code = cpjson::validate_schema(pcsaft_fluids_schema_JSON, JSON, errstr);
```
and delete the explanatory comment about the ambiguity (the FluidLibrary.h→Configuration.h→rapidjson chain that caused it is now gone).

- [ ] **Step 2: Cubics — same**

In `CubicsLibrary.cpp`, replace its `validate_schema_nlohmann` cast block with the direct `cpjson::validate_schema(cubic_fluids_schema_JSON, JSON, errstr);` call; delete the comment.

- [ ] **Step 3: Build + verify the casts compiled away cleanly**

Run: `cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[PCSAFT],[cubic],[json_validation]"`
Expected: clean build (no ambiguity error — proving the rapidjson overload is no longer visible in these TUs), tests pass.

- [ ] **Step 4: Commit**
```bash
git add src/Backends/PCSAFT/PCSAFTLibrary.cpp src/Backends/Cubics/CubicsLibrary.cpp
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit -m "refactor(json): drop the validate_schema overload-disambiguation casts (CoolProp-xa8w.2)

The casts existed only because Configuration.h pulled detail/rapidjson.h into
these TUs (via FluidLibrary.h), making both validate_schema overloads visible.
Phase 5 removed that path, so the nlohmann/Valijson overload is unambiguous.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final Verification (before PR)

- [ ] **Step 1: The exhaustive rapidjson grep — only the detail headers remain**

Run: `git grep -l "rapidjson" -- 'include/**' 'src/**'`
Expected: ONLY `include/CoolProp/detail/rapidjson.h`, `include/rapidjson_include.h`, `include/CoolProp/detail/json.h` (the enum-guard comment), and `src/Backends/PCSAFT/PCSAFTFluidFactory.h` (a comment). NO Configuration/SchemaValidation/SVDSBTL/MixtureParameters/REFPROP/loader-cast references. If anything else appears, it was missed — convert it.

- [ ] **Step 2: Installed headers JSON-library-free**

Run: `grep -n "rapidjson\|nlohmann" include/CoolProp/Configuration.h include/CoolProp/SchemaValidation.h`
Expected: no matches.

- [ ] **Step 3: Full fast suite + symbol gate**

Run: `./build_catch/CatchTestRunner "~[slow]"` (0 failures) and `./dev/ci/check-json-symbols.sh "$(find build_shared -name 'libCoolProp.*' | head -1)"` after a shared build (still 0 nlohmann/valijson exported).

- [ ] **Step 4: Preflight + adversarial review (MANDATORY)**

`./dev/ci/preflight.sh`, then dispatch the review vs `master`. Focus: (a) `item_from_json` type checks match the old rapidjson ones (`is_number()` = old `IsDouble()||IsInt()`; bool/int/string guards preserved); (b) the canonical-JSON output is genuinely sorted/deterministic and the SVDSBTL cache-key/opthash change (rapidjson→nlohmann number formatting) is acknowledged (cache invalidation only, not correctness); (c) the two casts compiled away because the overload is truly gone (not silently still-ambiguous); (d) `ConfigurationItem` accessors used match the real interface. Then push + `gh pr create`, then re-review the delta.

- [ ] **Step 5: bd**

`bd update CoolProp-xa8w.2 --notes "PR #<n>…"`. Note in `xa8w.3` (Phase Final) that the tree is now one mechanical deletion away (only `detail/rapidjson.h` + the shim + the enum guard remain).

---

## Self-Review (plan vs spec)

**Spec coverage:** §3 Configuration → Task 2; §3 SchemaValidation → Task 1; §3 SVDSBTL/MixtureParameters/REFPROP → Task 3; §3 cast deletion → Task 4; §5 verification (exhaustive grep, installed-header check, suite, symbol gate) → Final Verification; §1 string-only public API → Tasks 1-2 (rapidjson decls removed, string API kept).

**Placeholder scan:** none — SchemaValidation.cpp is given in full; the Configuration free functions + string fns are given in full (with the explicit caveat to match the real `ConfigurationItem` setter names); SVDSBTL/MixtureParameters/REFPROP use the idiom map (the implementer reads each site and applies it, as in Phase 2-4).

**Type consistency:** `cpjson::validate_schema(schema, instance, errstr)` arg order is used consistently (Task 1, and the loaders in Task 4); `item_to_json`/`item_from_json` signatures match their call sites in `get/set_config_as_json_string`; `to_canonical_json_str`/`validate_json_against_schema` are the string API referenced by SVDSBTL (Task 3) and the test (Task 1).
