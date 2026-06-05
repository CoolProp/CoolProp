# RapidJSON → nlohmann/json — Phase 1 (Helmholtz loader + CBOR) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Migrate the Helmholtz fluid loader off RapidJSON onto nlohmann/json, then embed `all_fluids` as uncompressed CBOR (dropping zlib/miniz from this load path), with no nlohmann type in any installed header and provable behavior + byte equivalence.

**Architecture:** Two commits in one PR. **Commit A**: convert `FluidLibrary`'s `parse_*` methods and `add_many`/`add_one` to `nlohmann::json` (format unchanged — still zlib+JSON, but parsed by nlohmann), de-publish the JSON-taking constructors as free-function factories in a non-installed header, and drop the `to_json` methods. **Commit B**: emit the fluid blob as CBOR via `cbor2` and load it with `from_cbor`, dropping zlib.

**Tech Stack:** C++17 (CoolProp), nlohmann/json, CMake, Python `cbor2`, Catch2, incbin/miniz (B removes miniz from this path).

**Spec:** `docs/superpowers/specs/2026-06-04-rapidjson-to-nlohmann-phase1-helmholtz-cbor-design.md`

**Prerequisites:** configure + build once:
```bash
cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_catch --target CatchTestRunner -j8
```
Behavior-parity tag scope for the Helmholtz loader: run `./build_catch/CatchTestRunner "[fluid]"` plus a broad `"[!slow]"` sweep (see CLAUDE.md test-filter discipline; `./dev/ci/preflight.sh` auto-selects).

---

## File Structure

| File | Responsibility | Action |
|------|----------------|--------|
| `src/Backends/Helmholtz/Fluids/FluidLibraryFactories.h` | Non-installed: `cpjson::make_*` factories that build ancillary/surface-tension objects from `nlohmann::json` | Create |
| `include/CoolProp/fluids/Ancillaries.h` | Public (installed): drop JSON ctors, stop including `detail/rapidjson.h` | Modify |
| `src/Backends/Helmholtz/Fluids/Ancillaries.cpp` | `SaturationAncillaryFunction` JSON-ctor body → moved/retyped into the factory | Modify |
| `include/CoolProp/fluids/Helmholtz.h` | Public (installed): drop the ~13 `to_json` methods, stop including `detail/rapidjson.h` | Modify |
| `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp:4290-4299` | Drop the `debug_level>10` `to_json` debug branch | Modify |
| `src/Backends/Helmholtz/Fluids/FluidLibrary.h` | Convert ~23 `parse_*` + `add_many`/`add_one` to `nlohmann::json`; swap include | Modify |
| `src/Backends/Helmholtz/Fluids/FluidLibrary.cpp` | `load()`: nlohmann parse (A), then `from_cbor` (B) | Modify |
| `include/CoolProp/detail/json.h` | Add `cpjson::from_cbor` (B) | Modify |
| `dev/generate_headers.py` | Emit CBOR + round-trip self-check; drop the zlib path (B) | Modify |
| `include/all_fluids_CBOR.h` | Generated CBOR incbin header (B); replaces `all_fluids_JSON_z.h` | Generated |
| `CMakeLists.txt` | Add `PATTERN "*_CBOR*.h" EXCLUDE` to the 3 header-install rules (B); CBOR source-path macro for the test (B) | Modify |
| `src/Tests/CoolProp-Tests-CBOR.cpp` | `[cbor]` byte-equivalence test (B) | Create |

---

## Conversion Reference — rapidjson → nlohmann idiom map

Use this table for every `parse_*` body. **The `cpjson::get_*` helpers are unchanged** (same names/signatures in `detail/json.h`), so lines that already call them need only their *input node type* changed. Translate only the direct rapidjson API:

| rapidjson | nlohmann |
|-----------|----------|
| `rapidjson::Value& x` (param) | `const nlohmann::json& x` (or non-const where mutated — loaders don't mutate) |
| `v.HasMember("k")` | `v.contains("k")` |
| `v["k"]` (read) | `v.at("k")` |
| `v["k"].GetString()` | `v.at("k").get<std::string>()` (or `cpjson::get_string(v,"k")`) |
| `v["k"].GetDouble()` | `cpjson::get_double(v,"k")` |
| `v.IsArray()` | `v.is_array()` |
| `el.GetString()` (in a loop) | `el.get<std::string>()` |
| `for (auto itr = v.Begin(); itr != v.End(); ++itr) { rapidjson::Value& c = *itr; ... }` | `for (const auto& c : v) { ... }` |
| `rapidjson::Value::ValueIterator` / `ConstValueIterator` | range-`for` over `const auto&` |
| `cpjson::get_long_double_array(v["n"])` | **unchanged** — `cpjson::get_long_double_array(v.at("n"))` (use `.at`) |
| `cpjson::json2string(v["SUPERANCILLARY"])` | `v.at("SUPERANCILLARY").dump()` |

Member access: prefer `v.at("k")` (throws `nlohmann::type_error`/`out_of_range` on miss) but the surrounding code already guards with `contains`; where a guard exists, `v.at("k")` is safe. The `cpjson` getters keep the friendly `CoolProp::ValueError` messages, so prefer them for scalars.

---

# COMMIT A — Parser swap (format unchanged)

## Task A1: Create the factory header (SurfaceTensionCorrelation)

**Files:**
- Create: `src/Backends/Helmholtz/Fluids/FluidLibraryFactories.h`

- [ ] **Step 1: Create the non-installed factory header**

Create `src/Backends/Helmholtz/Fluids/FluidLibraryFactories.h`:

```cpp
#ifndef FLUIDLIBRARY_FACTORIES_H
#define FLUIDLIBRARY_FACTORIES_H

// NON-INSTALLED (lives under src/). Builds fluid-model objects from
// nlohmann::json nodes, keeping nlohmann out of the installed public headers
// (Ancillaries.h / Helmholtz.h). The loader includes this; nothing public does.

#include "CoolProp/detail/json.h"
#include "CoolProp/fluids/Ancillaries.h"

namespace cpjson {

/// Build a SurfaceTensionCorrelation from its JSON node.
inline CoolProp::SurfaceTensionCorrelation make_surface_tension_correlation(const nlohmann::json& j) {
    CoolProp::SurfaceTensionCorrelation out;
    out.a = cpjson::get_long_double_array(j.at("a"));
    out.n = cpjson::get_long_double_array(j.at("n"));
    out.Tc = cpjson::get_double(j, "Tc");
    out.BibTeX = cpjson::get_string(j, "BibTeX");
    out.N = out.n.size();
    out.s = out.n;
    return out;
}

}  // namespace cpjson

#endif  // FLUIDLIBRARY_FACTORIES_H
```

> The members (`a`, `n`, `Tc`, `BibTeX`, `N`, `s`) are already public on `SurfaceTensionCorrelation` (see `Ancillaries.h`). The default constructor stays.

- [ ] **Step 2: Add the SaturationAncillaryFunction factory**

Read the existing `SaturationAncillaryFunction(rapidjson::Value&)` constructor body in `src/Backends/Helmholtz/Fluids/Ancillaries.cpp`. Add a factory to `FluidLibraryFactories.h` (inside `namespace cpjson`) named `make_saturation_ancillary(const nlohmann::json& j) -> CoolProp::SaturationAncillaryFunction` whose body is that constructor's body with the idiom map applied (the `cpjson::get_*` calls are unchanged; `json_code["x"]` → `j.at("x")`, `HasMember` → `contains`, etc.). If the class has no suitable public setter path, give it a private constructor taking the parsed fields and `friend` the factory, OR add a `from_json`-free path — choose the minimal change that compiles; report which you used as a concern.

- [ ] **Step 3: Build to confirm the header compiles standalone**

Add a temporary `.cpp` or rely on Task A4's include. For now:
```bash
cmake --build build_catch --target CatchTestRunner -j8 2>&1 | tail -5
```
Expected: builds (the header is not yet included anywhere, so this just confirms no syntax error once included in A4). Defer the real build to A4.

- [ ] **Step 4: Commit (interim — squashed into Commit A later)**

```bash
git add src/Backends/Helmholtz/Fluids/FluidLibraryFactories.h
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit --no-verify -m "feat(json): add nlohmann factories for fluid-model objects (Phase 1A)"
```

## Task A2: De-publish the JSON constructors from Ancillaries.h

**Files:**
- Modify: `include/CoolProp/fluids/Ancillaries.h`
- Modify: `src/Backends/Helmholtz/Fluids/Ancillaries.cpp`

- [ ] **Step 1: Remove the SurfaceTensionCorrelation JSON ctor**

In `include/CoolProp/fluids/Ancillaries.h`, delete the `SurfaceTensionCorrelation(rapidjson::Value& json_code) { ... }` constructor (lines ~36-49). Keep the default constructor and all members.

- [ ] **Step 2: Remove the SaturationAncillaryFunction JSON ctor declaration + definition**

Delete the `SaturationAncillaryFunction(rapidjson::Value& json_code);` declaration in `Ancillaries.h` (~line 121) and its definition in `Ancillaries.cpp` (now superseded by the factory from A1 Step 2).

- [ ] **Step 3: Drop the rapidjson include**

In `include/CoolProp/fluids/Ancillaries.h`, remove `#include "CoolProp/detail/rapidjson.h"` (line 7). If any remaining code in the header needs `CoolPropDbl`/helpers, include `CoolProp/CPnumerics.h` or the minimal header that provides it instead — verify by compiling.

- [ ] **Step 4: Build (will fail at call sites until A4)**

```bash
cmake --build build_catch --target CatchTestRunner -j8 2>&1 | tail -15
```
Expected: errors only at the `parse_ancillaries`/`parse_surface_tension` call sites in `FluidLibrary.h` (no more `SaturationAncillaryFunction(...)`/`SurfaceTensionCorrelation(...)` ctor). That confirms the de-publish is complete; A4 fixes the call sites.

- [ ] **Step 5: Commit (interim)**

```bash
git add include/CoolProp/fluids/Ancillaries.h src/Backends/Helmholtz/Fluids/Ancillaries.cpp
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit --no-verify -m "refactor(json): de-publish Ancillaries JSON constructors (Phase 1A)"
```

## Task A3: Drop the to_json methods + debug branch

**Files:**
- Modify: `include/CoolProp/fluids/Helmholtz.h`
- Modify: `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp`

- [ ] **Step 1: Remove the debug caller**

In `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp`, delete the `if (get_debug_level() > 10) { ... }` block in `_calc_all_critical_points` (~lines 4290-4299) — the one that builds a `rapidjson::Document` and calls `mat[0][1]->phi.to_json(val, doc)`.

- [ ] **Step 2: Remove every to_json declaration/definition in Helmholtz.h**

In `include/CoolProp/fluids/Helmholtz.h`, delete each `void to_json(rapidjson::Value& el, rapidjson::Document& doc);` declaration and any inline definitions (there are ~13; grep `to_json` in the file). Remove `#include "CoolProp/detail/rapidjson.h"` (line 8) once no rapidjson type remains.

- [ ] **Step 3: Grep-confirm no remaining to_json / rapidjson in these public headers**

```bash
grep -n "to_json\|rapidjson" include/CoolProp/fluids/Helmholtz.h include/CoolProp/fluids/Ancillaries.h
```
Expected: no matches.

- [ ] **Step 4: Build**

```bash
cmake --build build_catch --target CatchTestRunner -j8 2>&1 | tail -15
```
Expected: still only the FluidLibrary call-site errors from A2 (no new errors from the to_json removal — it had a single, now-deleted caller).

- [ ] **Step 5: Commit (interim)**

```bash
git add include/CoolProp/fluids/Helmholtz.h src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit --no-verify -m "refactor(json): drop debug-only to_json serializers from Helmholtz (Phase 1A)"
```

## Task A4: Convert FluidLibrary parse_* + add_many/add_one to nlohmann

This is the bulk mechanical conversion. Apply the Conversion Reference table to every `parse_*` body and the `add_*` signatures.

**Files:**
- Modify: `src/Backends/Helmholtz/Fluids/FluidLibrary.h`

- [ ] **Step 1: Swap the include + add the factory include**

In `FluidLibrary.h`, replace `#include "CoolProp/detail/rapidjson.h"` (line 7) with:
```cpp
#include "CoolProp/detail/json.h"
#include "FluidLibraryFactories.h"
```

- [ ] **Step 2: Convert the two factory call sites (worked example)**

`parse_surface_tension` (currently `FluidLibrary.h:1185-1187`):
```cpp
void parse_surface_tension(const nlohmann::json& surface_tension, CoolPropFluid& fluid) {
    fluid.ancillaries.surface_tension = cpjson::make_surface_tension_correlation(surface_tension);
}
```
`parse_ancillaries` (currently `1130-1182`): change the parameter to `const nlohmann::json& ancillaries`, replace every `ancillaries.HasMember("k")` with `ancillaries.contains("k")` and every `SaturationAncillaryFunction(ancillaries["k"])` with `cpjson::make_saturation_ancillary(ancillaries.at("k"))`. The structure/branches are otherwise unchanged. Example for the first lines:
```cpp
void parse_ancillaries(const nlohmann::json& ancillaries, CoolPropFluid& fluid) {
    if (!ancillaries.contains("rhoL") || !ancillaries.contains("rhoV")) {
        throw ValueError("Ancillary curves for either rhoL or rhoV are missing");
    }
    fluid.ancillaries.rhoL = cpjson::make_saturation_ancillary(ancillaries.at("rhoL"));
    fluid.ancillaries.rhoV = cpjson::make_saturation_ancillary(ancillaries.at("rhoV"));
    // ... remaining branches: HasMember -> contains, ["k"] -> .at("k") ...
}
```

- [ ] **Step 3: Convert the remaining parse_* methods**

Apply the Conversion Reference table to each of these (all in `FluidLibrary.h`), changing the parameter to `const nlohmann::json&` and translating direct rapidjson API; leave `cpjson::get_*` calls intact (switch their `v["k"]` arg to `v.at("k")`):
`parse_alphar`, `parse_alpha0`, `parse_environmental`, `parse_EOS`, `parse_EOS_listing`, `parse_dilute_viscosity`, `parse_initial_density_viscosity`, `parse_higher_order_viscosity`, `parse_ECS_conductivity`, `parse_ECS_viscosity`, `parse_Chung_viscosity`, `parse_rhosr_viscosity`, `parse_viscosity`, `parse_dilute_conductivity`, `parse_residual_conductivity`, `parse_critical_conductivity`, `parse_thermal_conductivity`, `parse_transport`, `parse_melting_line`, `parse_states`.
In `parse_EOS`, replace `cpjson::json2string(EOS_json["SUPERANCILLARY"])` with `EOS_json.at("SUPERANCILLARY").dump()`.

- [ ] **Step 4: Convert add_many / add_one signatures**

In `FluidLibrary.h`:
```cpp
void add_many(const nlohmann::json& listing);
void add_one(const nlohmann::json& fluid_json);
```
Convert their bodies (iterate `listing` with range-`for`; `fluid_json.at("...")`). Update the `get_JSONstring` body if it uses rapidjson `Document`/allocator — if it has no live caller, you may leave it temporarily ONLY IF it no longer references rapidjson; otherwise convert it to build an `nlohmann::json` and `.dump()`. Grep first: `grep -rn "get_JSONstring" src include wrappers | grep -v CubicsLibrary`.

- [ ] **Step 5: Build**

```bash
cmake --build build_catch --target CatchTestRunner -j8 2>&1 | tail -20
```
Expected: compiles except `FluidLibrary.cpp::load()` (still hands a `rapidjson::Document` to `add_many`). A5 fixes load().

- [ ] **Step 6: Grep-confirm FluidLibrary.h is rapidjson-free**

```bash
grep -n "rapidjson" src/Backends/Helmholtz/Fluids/FluidLibrary.h
```
Expected: no matches.

- [ ] **Step 7: Commit (interim)**

```bash
git add src/Backends/Helmholtz/Fluids/FluidLibrary.h
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit --no-verify -m "refactor(json): convert FluidLibrary parse_* to nlohmann (Phase 1A)"
```

## Task A5: Switch load() to nlohmann parse (format still zlib+JSON) + parity

**Files:**
- Modify: `src/Backends/Helmholtz/Fluids/FluidLibrary.cpp`

- [ ] **Step 1: Parse the decompressed string with nlohmann**

In `FluidLibrary.cpp::load()`, replace the rapidjson parse block (lines 61-72) with:
```cpp
    try {
        nlohmann::json dd = cpjson::parse(buf);
        library.add_many(dd);
    } catch (std::exception& e) {
        std::cout << e.what() << '\n';
    }
```
Add `#include "CoolProp/detail/json.h"` at the top; the zlib `uncompress` block (lines 47-53) and the `all_fluids_JSON_z.h` include stay (format unchanged in Commit A).

- [ ] **Step 2: Build**

```bash
cmake --build build_catch --target CatchTestRunner -j8 2>&1 | tail -10
```
Expected: clean build.

- [ ] **Step 3: Behavior parity — full fluid suite must stay green**

```bash
./build_catch/CatchTestRunner "[fluid]" 2>&1 | tail -5
./build_catch/CatchTestRunner "[!slow]" 2>&1 | tail -5
```
Expected: PASS — every fluid loads with identical numbers. This is the safety net for the whole conversion. If any fluid-property test fails, a parse_* conversion changed behavior — debug against the rapidjson original before proceeding.

- [ ] **Step 4: First-load benchmark (in-context)**

Time a cold first-load (process start → first `PropsSI("D","T",300,"P",1e5,"Water")`), compare to master. Expected: ~+70 ms vs master (the naive nlohmann-parse cost; spec §6). Record the number in the PR.

- [ ] **Step 5: Symbol-leak gate (now load-bearing — nlohmann is in libCoolProp)**

```bash
cmake -B build_shared -S . -DCOOLPROP_SHARED_LIBRARY=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_shared -j8
SHARED_LIB="$(find build_shared \( -name 'libCoolProp.so' -o -name 'libCoolProp.dylib' \) | head -1)"
./dev/ci/check-json-symbols.sh "$SHARED_LIB"
```
Expected: `OK: no symbols matching /nlohmann|valijson/ exported …`. (If it FAILS, a public header still includes nlohmann — find and de-publish it.)

- [ ] **Step 6: Commit (Commit A)**

```bash
git add src/Backends/Helmholtz/Fluids/FluidLibrary.cpp
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit --no-verify -m "refactor(json): Helmholtz loader parses via nlohmann (Phase 1A, format unchanged)"
```

---

# COMMIT B — Format swap (JSON → CBOR, drop zlib)

## Task B1: Add cpjson::from_cbor helper

**Files:**
- Modify: `include/CoolProp/detail/json.h`
- Modify: `src/Tests/CoolProp-Tests-JSONHelpers.cpp`

- [ ] **Step 1: Write the failing test**

Append to `src/Tests/CoolProp-Tests-JSONHelpers.cpp` before the final `#endif`:
```cpp
TEST_CASE("cpjson::from_cbor round-trips a document", "[json]") {
    nlohmann::json j = cpjson::parse(R"({"a": 1.5, "b": [1, 2, 3]})");
    std::vector<std::uint8_t> blob = nlohmann::json::to_cbor(j);
    nlohmann::json back = cpjson::from_cbor(blob.data(), blob.size());
    REQUIRE(back == j);
}

TEST_CASE("cpjson::from_cbor throws ValueError on garbage", "[json]") {
    std::vector<std::uint8_t> bad = {0xFF, 0xFF, 0xFF};
    REQUIRE_THROWS_AS(cpjson::from_cbor(bad.data(), bad.size()), CoolProp::ValueError);
}
```

- [ ] **Step 2: Run to verify failure**

```bash
cmake --build build_catch --target CatchTestRunner -j8
```
Expected: compile error — `from_cbor` not a member of `cpjson`.

- [ ] **Step 3: Implement from_cbor**

In `include/CoolProp/detail/json.h`, inside `namespace cpjson` after `parse`, add:
```cpp
/// Decode a CBOR byte buffer into an nlohmann::json document.
/// Throws CoolProp::ValueError (never a raw nlohmann exception) on failure.
inline nlohmann::json from_cbor(const std::uint8_t* data, std::size_t size) {
    try {
        return nlohmann::json::from_cbor(data, data + size);
    } catch (const std::exception& e) {
        throw CoolProp::ValueError(std::string("Unable to decode CBOR: ") + e.what());
    }
}
```

- [ ] **Step 4: Run to verify pass**

```bash
cmake --build build_catch --target CatchTestRunner -j8
./build_catch/CatchTestRunner "[json]" 2>&1 | tail -3
```
Expected: all `[json]` tests PASS.

- [ ] **Step 5: Commit**

```bash
git add include/CoolProp/detail/json.h src/Tests/CoolProp-Tests-JSONHelpers.cpp
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit --no-verify -m "feat(json): add cpjson::from_cbor helper (Phase 1B)"
```

## Task B2: Emit CBOR from generate_headers.py (with self-check)

**Files:**
- Modify: `dev/generate_headers.py`

- [ ] **Step 1: Switch combine_json to write CBOR + self-check**

In `dev/generate_headers.py::combine_json`, change the fluids `DependencyManager` destination from `all_fluids.json.z` to `all_fluids.cbor` (line 361), and replace the write block (line 387-388):
```python
        import cbor2
        blob = cbor2.dumps(master)
        # Round-trip self-check: a bad encode must never be committed.
        if cbor2.loads(blob) != master:
            raise ValueError("CBOR round-trip self-check failed for all_fluids")
        with (Path(root_dir) / 'dev' / 'all_fluids.cbor').open('wb') as fp:
            fp.write(blob)
```
Keep the `all_fluids.json` / `all_fluids_verbose.json` writes (the source JSON stays for the byte-equivalence test and human inspection).

- [ ] **Step 2: Update the zvalues entry + the incbin comment**

Change `zvalues` (line 62-64) to:
```python
zvalues = [
    ('all_fluids.cbor', 'all_fluids_CBOR.h', 'all_fluids_CBOR'),
]
```
Update the hardcoded comment in `incbin_template` (line 72) from `all_fluids.json.z` to `all_fluids.cbor`.

- [ ] **Step 3: Regenerate + verify the CBOR header is produced**

```bash
python3 -c "import cbor2" || uvx --from cbor2 python -c "import cbor2"   # confirm cbor2 available
python3 dev/generate_headers.py   # or the project's invocation; check how CMake calls it
ls -lh include/all_fluids_CBOR.h dev/all_fluids.cbor
```
Expected: both produced; `all_fluids.cbor` ~4.4 MB; the header is a hex incbin array. Confirm the self-check passed (no exception).

> Note: `cbor2` must be available to the build's Python. If CMake drives `generate_headers.py`, ensure `cbor2` is on that interpreter (add to the dev requirements / the CMake Python invocation). Report how the build invokes the generator.

- [ ] **Step 4: Commit**

```bash
git add dev/generate_headers.py dev/all_fluids.cbor include/all_fluids_CBOR.h
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit --no-verify -m "build(json): generate all_fluids as CBOR with round-trip self-check (Phase 1B)"
```

## Task B3: Load from CBOR, drop zlib, exclude the CBOR header from install

**Files:**
- Modify: `src/Backends/Helmholtz/Fluids/FluidLibrary.cpp`
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Add the install EXCLUDE for the CBOR header**

In `CMakeLists.txt`, in each of the three `install(DIRECTORY .../include ... FILES_MATCHING PATTERN "*.h" ...)` rules that already have `PATTERN "*_JSON*.h" EXCLUDE` (the static, shared, and Linux/usr blocks — grep `_JSON*.h` to find all three), add alongside it:
```cmake
            PATTERN "*_CBOR*.h" EXCLUDE
```
Otherwise the generated `all_fluids_CBOR.h` ships in the installed tree.

- [ ] **Step 2: Switch load() to from_cbor; remove zlib**

In `FluidLibrary.cpp`, replace the include block (lines 9-24) so it incbin's `all_fluids_CBOR.h` (symbol `all_fluids_CBOR`) instead of `all_fluids_JSON_z.h`, and rewrite `load()` (lines 46-73) to:
```cpp
void load() {
    if (getenv("COOLPROP_DISABLE_SUPERANCILLARIES_ENTIRELY")) {
        std::cout << "CoolProp: superancillaries have been disabled because the COOLPROP_DISABLE_SUPERANCILLARIES_ENTIRELY environment variable has been defined" << '\n';
    }
    try {
        nlohmann::json dd = cpjson::from_cbor(gall_fluids_CBORData, gall_fluids_CBORSize);
        library.add_many(dd);
    } catch (std::exception& e) {
        std::cout << e.what() << '\n';
    }
}
```
Remove `#include "miniz.h"` and the `uncompress`/`outbuffer` logic. (Confirm the incbin symbol names: the generator emits `g<SYMBOL>Data` / `g<SYMBOL>Size`, i.e. `gall_fluids_CBORData` / `gall_fluids_CBORSize`.)

- [ ] **Step 3: Build + parity + byte-size**

```bash
cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_catch --target CatchTestRunner -j8
./build_catch/CatchTestRunner "[fluid]" 2>&1 | tail -5
./build_catch/CatchTestRunner "[!slow]" 2>&1 | tail -5
```
Expected: clean build; all fluid tests PASS (parity preserved through the format change).

- [ ] **Step 4: First-load benchmark**

Re-time cold first-load. Expected: ~35 ms (spec §6) — at or below master's ~37 ms, and well below Commit A's ~107 ms. Record in the PR.

- [ ] **Step 5: Commit**

```bash
git add src/Backends/Helmholtz/Fluids/FluidLibrary.cpp CMakeLists.txt
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit --no-verify -m "perf(json): load fluid data from embedded CBOR, drop zlib (Phase 1B)"
```

## Task B4: Byte-equivalence [cbor] test

**Files:**
- Create: `src/Tests/CoolProp-Tests-CBOR.cpp`
- Modify: `CMakeLists.txt` (register test + source-path macro)

- [ ] **Step 1: Provide the source-JSON path via a CMake macro**

In `CMakeLists.txt`, near the `CatchTestRunner` target definition, add:
```cmake
  target_compile_definitions(CatchTestRunner PRIVATE COOLPROP_ALL_FLUIDS_JSON_PATH="${CMAKE_CURRENT_SOURCE_DIR}/dev/all_fluids.json")
```

- [ ] **Step 2: Write the byte-equivalence test**

Create `src/Tests/CoolProp-Tests-CBOR.cpp`:
```cpp
// Byte-equivalence gate (RapidJSON->nlohmann Phase 1): the embedded CBOR blob
// must decode to exactly the source dev/all_fluids.json (value-equality). This
// catches encoder/decoder drift (cbor2/nlohmann version bumps) and staleness
// (source JSON edited without regenerating the blob).
#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>
#    include <fstream>
#    include <sstream>
#    include <string>

#    include "CoolProp/detail/json.h"

extern "C" {
extern const unsigned char gall_fluids_CBORData[];
extern const unsigned int gall_fluids_CBORSize;
}

TEST_CASE("embedded CBOR decodes to the source all_fluids.json", "[cbor]") {
    nlohmann::json from_blob = cpjson::from_cbor(gall_fluids_CBORData, gall_fluids_CBORSize);

    std::ifstream f(COOLPROP_ALL_FLUIDS_JSON_PATH, std::ios::binary);
    REQUIRE(f.good());
    std::ostringstream ss;
    ss << f.rdbuf();
    nlohmann::json from_src = cpjson::parse(ss.str());

    REQUIRE(from_blob == from_src);
}

#endif  // ENABLE_CATCH
```

- [ ] **Step 3: Register the test file**

In `CMakeLists.txt`, in the `COOLPROP_CATCH_MODULE` `APP_SOURCES` list, add:
```cmake
  list(APPEND APP_SOURCES
       "${CMAKE_CURRENT_SOURCE_DIR}/src/Tests/CoolProp-Tests-CBOR.cpp")
```

- [ ] **Step 4: Reconfigure, build, run**

```bash
cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_catch --target CatchTestRunner -j8
./build_catch/CatchTestRunner "[cbor]" 2>&1 | tail -5
```
Expected: PASS — embedded CBOR == source JSON. (Note: under `COOLPROP_NO_INCBIN`, the symbols come from the included header; confirm the `gall_fluids_CBOR*` symbols link in the test runner. If incbin is disabled for tests, the header form must define the array — it does.)

- [ ] **Step 5: Commit**

```bash
git add src/Tests/CoolProp-Tests-CBOR.cpp CMakeLists.txt
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit --no-verify -m "test(json): [cbor] byte-equivalence gate vs source all_fluids.json (Phase 1B)"
```

## Task B5: Phase-1 verification sweep

- [ ] **Step 1: Full fast suite + cbor + json**

```bash
cmake --build build_catch --target CatchTestRunner -j8
./build_catch/CatchTestRunner "[!slow]" 2>&1 | tail -5
./build_catch/CatchTestRunner "[cbor][json]" 2>&1 | tail -5
```
Expected: PASS.

- [ ] **Step 2: Symbol gate on the shared lib**

```bash
cmake --build build_shared -j8
SHARED_LIB="$(find build_shared \( -name 'libCoolProp.so' -o -name 'libCoolProp.dylib' \) | head -1)"
./dev/ci/check-json-symbols.sh "$SHARED_LIB"
```
Expected: `OK: no symbols matching /nlohmann|valijson/ exported …`.

- [ ] **Step 3: Grep — rapidjson gone from the Helmholtz path**

```bash
grep -rn "rapidjson" src/Backends/Helmholtz/Fluids/FluidLibrary.* include/CoolProp/fluids/Helmholtz.h include/CoolProp/fluids/Ancillaries.h
```
Expected: no matches. (RapidJSON remains elsewhere — other loaders — until their phases.)

- [ ] **Step 4: Pre-PR adversarial review (per CLAUDE.md, MANDATORY)**

Invoke the review subagent (`superpowers:code-reviewer`, or `general-purpose` fallback) on the diff vs the branch base. Apply the fail-open-gate question and re-review-the-delta rule. Address/justify each blocking finding.

- [ ] **Step 5: Preflight + push + PR**

```bash
./dev/ci/preflight.sh
git push -u origin ihb/rapidjson-nlohmann-phase1-helmholtz
```
Then `gh pr create` against master. Note the two first-load benchmark numbers (A and B) in the PR body.

---

## Self-Review (completed during authoring)

- **Spec coverage:** parse_* conversion → A4; de-publish constructors via factories → A1/A2; drop to_json + debug branch → A3; CBOR emit + self-check → B2; from_cbor load + drop zlib → B1/B3; install EXCLUDE for CBOR header → B3 Step 1 (the §5 gotcha); `[cbor]` byte-equivalence + CMake path macro → B4; symbol gate load-bearing → A5/B5; A/B sequencing → commit structure; benchmarks → A5/B3; Phase-Final deferrals → spec §8 (not in this plan). All spec sections map to tasks.
- **Placeholder scan:** the two genuinely-read-at-execution spots are the `SaturationAncillaryFunction` ctor body (A1 Step 2 — its exact body is in Ancillaries.cpp; the transform is fully specified by the idiom map) and the bulk `parse_*` bodies (A4 Step 3 — specified by the idiom map + two fully-worked examples + the parity-test gate). This is the correct way to specify a large mechanical conversion; the parity suite (A5 Step 3) is the objective correctness check.
- **Type consistency:** `cpjson::make_surface_tension_correlation`, `cpjson::make_saturation_ancillary`, `cpjson::from_cbor`, `cpjson::parse`; incbin symbols `gall_fluids_CBORData`/`gall_fluids_CBORSize`; `const nlohmann::json&` params throughout. Consistent across tasks.

---

## Roadmap: subsequent phases (unchanged, each its own plan)

Phase 2 Incompressible · Phase 3 PC-SAFT · Phase 4 Cubics/UNIFAC · Phase 5 SVDSBTL+Configuration · Phase Final (delete rapidjson + `detail/rapidjson.h`, `detail/` install EXCLUDE, tighten symbol-gate pattern, **decide nlohmann getter-layer keep-vs-inline**).
