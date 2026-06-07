# RapidJSON → nlohmann/json — Phase 2-4 (Incompressible + PC-SAFT + Cubics/UNIFAC loaders) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert the Incompressible, PC-SAFT, and Cubics/UNIFAC backend loaders off RapidJSON onto nlohmann/json, with a POD-constructor de-publishing boundary that leaves zero nlohmann/rapidjson tokens in any installed header.

**Architecture:** Each loader is an include-swap (`detail/rapidjson.h` → `detail/json.h`) plus a mechanical idiom translation (the `cpjson::get_*` helpers are signature-identical across both wrappers). Schema validation moves from rapidjson `SchemaValidator` to Valijson automatically once the rapidjson overload is no longer visible. Two installed headers that currently expose JSON construction (`PCSAFTFluid.h`, `Ancillaries.h`) are converted to a POD-constructor boundary: a plain-typed constructor in the installed header, with all nlohmann parsing done by a non-installed `src/` factory. The four Phase-1 bridges (2 `alpha0` round-trips + 2 `validate_schema` function-pointer casts) are deleted.

**Tech Stack:** C++17, nlohmann/json (via the `cpjson` wrapper in `include/CoolProp/detail/json.h`), Valijson, Catch2 v3, CMake (Release), CoolProp `cpjson::` helper layer.

**Spec:** `docs/superpowers/specs/2026-06-06-rapidjson-to-nlohmann-phase2-4-loaders-design.md`
**Epic:** `CoolProp-xa8w.1`. **Branch:** `ihb/rapidjson-nlohmann-phase2-4-loaders`.

---

## Idiom map (applies to every mechanical conversion in this plan)

| rapidjson | nlohmann |
|-----------|----------|
| `rapidjson::Value& x` (param/local) | `const nlohmann::json& x` |
| `x.HasMember("k")` | `x.contains("k")` |
| `x["k"]` (read) | `x.at("k")` |
| `x["k"].GetInt()/GetDouble()/GetString()` | `x.at("k").get<int>()` / `.get<double>()` / `.get<std::string>()` (or the matching `cpjson::get_*`) |
| `x["k"].IsNumber()/IsObject()/IsArray()` | `x.at("k").is_number()/is_object()/is_array()` |
| `for (auto itr=L.Begin(); itr!=L.End(); ++itr) { …(*itr)… }` | `for (const auto& el : L) { …el… }` |
| `rapidjson::Document dd; dd.Parse<0>(s.data(), s.size()); if (dd.HasParseError()) throw …;` | `nlohmann::json dd = cpjson::parse(s);` (throws `CoolProp::ValueError` on bad JSON) |
| `cpjson::json2string(v)` | `v.dump()` (or keep `cpjson::json2string(v)` — identical in `detail/json.h`) |
| `cpjson::get_string(*itr,"k")` etc. | `cpjson::get_string(el,"k")` — **call sites unchanged**, only the operand type changes |

**Optional-member rule:** rapidjson `x["k"]` is UB on a missing key; nlohmann `x.at("k")` throws. Wherever the original guards a key with `HasMember` (or `IsNumber`), keep an equivalent `contains("k")` (and `is_number()`) guard. Behavior parity is the safety net.

---

## Task 1: Incompressible loader → nlohmann

**Files:**
- Modify: `src/Backends/Incompressible/IncompressibleLibrary.h:9,151-153,165-166`
- Modify: `src/Backends/Incompressible/IncompressibleLibrary.cpp:5-6 (includes), 348-426, 554-566`

This loader has no schema, no user-supplied fluids, and no JSON-typed installed header — the simplest case. It proves the pattern.

- [ ] **Step 1: Swap the include in the header**

In `IncompressibleLibrary.h`, replace line 9:
```cpp
#include "CoolProp/detail/rapidjson.h"
```
with:
```cpp
#include "CoolProp/detail/json.h"
```

- [ ] **Step 2: Retype the parse_* / add_* declarations in the header**

In `IncompressibleLibrary.h`, change the three `parse_*` declarations (151-153) and the two `add_*` declarations (165-166) from `rapidjson::Value&` to `const nlohmann::json&`:
```cpp
    IncompressibleData parse_coefficients(const nlohmann::json& obj, const std::string& id, bool vital);
    double parse_value(const nlohmann::json& obj, const std::string& id, bool vital, double def);
    composition_types parse_ifrac(const nlohmann::json& obj, const std::string& id);
```
```cpp
    /// Add all the fluid entries in the nlohmann::json array passed in
    void add_many(const nlohmann::json& listing);
    void add_one(const nlohmann::json& fluid_json);
```
Also update the doc comment at lines 132-133 (`rapidjson::Value`/`rapidjson array` → `nlohmann::json`).

- [ ] **Step 3: Convert the include and parse_coefficients in the .cpp**

In `IncompressibleLibrary.cpp`, ensure the include is `detail/json.h` (the `#include "CoolProp/detail/rapidjson.h"` near the top, ~line 6, becomes `detail/json.h`; if `all_incompressibles_JSON.h` is the only other JSON include, leave it).

Replace `parse_coefficients` (348-393) — note the nested-member access `obj[id]["coeffs"]` → `obj.at(id).at("coeffs")` and the `HasMember` guards → `contains`:
```cpp
IncompressibleData JSONIncompressibleLibrary::parse_coefficients(const nlohmann::json& obj, const std::string& id, bool vital) {
    IncompressibleData fluidData;
    if (obj.contains(id)) {
        const nlohmann::json& entry = obj.at(id);
        if (entry.contains("type")) {
            if (entry.contains("coeffs")) {
                std::string type = cpjson::get_string(entry, "type");
                if (!type.compare("polynomial")) {
                    fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL;
                    fluidData.coeffs = vec_to_eigen(cpjson::get_double_array2D(entry.at("coeffs")));
                    return fluidData;
                } else if (!type.compare("exponential")) {
                    fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL;
                    fluidData.coeffs = vec_to_eigen(cpjson::get_double_array(entry.at("coeffs")));
                    return fluidData;
                } else if (!type.compare("logexponential")) {
                    fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL;
                    fluidData.coeffs = vec_to_eigen(cpjson::get_double_array(entry.at("coeffs")));
                    return fluidData;
                } else if (!type.compare("exppolynomial")) {
                    fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL;
                    fluidData.coeffs = vec_to_eigen(cpjson::get_double_array2D(entry.at("coeffs")));
                    return fluidData;
                } else if (!type.compare("polyoffset")) {
                    fluidData.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYOFFSET;
                    fluidData.coeffs = vec_to_eigen(cpjson::get_double_array(entry.at("coeffs")));
                    return fluidData;
                } else if (vital) {
                    throw ValueError(format("The type [%s] is not understood for [%s] of incompressible fluids. Please check your JSON file.",
                                            type.c_str(), id.c_str()));
                }
            } else {
                throw ValueError(format("Your file does not have an entry for \"coeffs\" in [%s], which is vital for this function.", id.c_str()));
            }
        } else {
            throw ValueError(format("Your file does not have an entry for \"type\" in [%s], which is vital for this function.", id.c_str()));
        }
    } else {
        if (vital) {
            throw ValueError(format("Your file does not have information for [%s], which is vital for an incompressible fluid.", id.c_str()));
        }
    }
    return fluidData;
}
```

- [ ] **Step 4: Convert parse_value, parse_ifrac, add_many in the .cpp**

`parse_value` (396-406):
```cpp
double JSONIncompressibleLibrary::parse_value(const nlohmann::json& obj, const std::string& id, bool vital, double def = 0.0) {
    if (obj.contains(id)) {
        return cpjson::get_double(obj, id);
    } else {
        if (vital) {
            throw ValueError(format("Your file does not have information for [%s], which is vital for an incompressible fluid.", id.c_str()));
        } else {
            return def;
        }
    }
}
```
`parse_ifrac` (409): only the parameter type changes — `rapidjson::Value& obj` → `const nlohmann::json& obj`; body unchanged.
`add_many` (422-426):
```cpp
void JSONIncompressibleLibrary::add_many(const nlohmann::json& listing) {
    for (const auto& fluid_json : listing) {
        add_one(fluid_json);
    }
}
```
`add_one` (428): change the signature `rapidjson::Value& fluid_json` → `const nlohmann::json& fluid_json`; the body uses only `cpjson::get_string(fluid_json, …)` and the `parse_*` helpers, so it is otherwise unchanged.

- [ ] **Step 5: Convert load_incompressible_library in the .cpp**

Replace (554-566):
```cpp
void load_incompressible_library() {
    // This json formatted string comes from the all_incompressibles_JSON.h header
    nlohmann::json dd = cpjson::parse(all_incompressibles_JSON);
    try {
        library.add_many(dd);
    } catch (std::exception& e) {
        std::cout << e.what() << '\n';
    }
}
```
(`cpjson::parse` throws `ValueError` on malformed JSON, replacing the explicit `HasParseError` check.)

- [ ] **Step 6: Build**

Run: `cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release && cmake --build build_catch --target CatchTestRunner -j8`
Expected: clean build.

- [ ] **Step 7: Run the incompressible suite (behavior parity)**

Run: `./build_catch/CatchTestRunner "[INCOMP]"`
Expected: all assertions pass (every incompressible fluid loads with identical numbers). The `[1374]` Melinder-source-data and `[1578]`/`[2209]` tests are the key parity checks.

- [ ] **Step 8: Verify the loader is rapidjson-free**

Run: `grep -rn "rapidjson" src/Backends/Incompressible/`
Expected: no matches.

- [ ] **Step 9: Commit**

```bash
git add src/Backends/Incompressible/IncompressibleLibrary.h src/Backends/Incompressible/IncompressibleLibrary.cpp
git commit -m "refactor(json): incompressible loader off rapidjson onto nlohmann (CoolProp-xa8w.1)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Ancillaries.h POD-constructor retrofit (remove the json_fwd downstream dependency)

**Files:**
- Modify: `include/CoolProp/fluids/Ancillaries.h:11-23 (remove json_fwd + friend), 89-126 (enum public + Values struct + ctor decl)`
- Modify: `src/Backends/Helmholtz/Fluids/Ancillaries.cpp` (add the out-of-line ctor)
- Modify: `src/Backends/Helmholtz/Fluids/FluidLibraryFactories.h:32-70 (rewrite make_saturation_ancillary to call the POD ctor)`

`SaturationAncillaryFunction` currently friends `cpjson::make_saturation_ancillary(const nlohmann::json&)` and `#include <nlohmann/json_fwd.hpp>`. That header is reachable downstream via `CoolPropFluid.h`, imposing a `json_fwd.hpp` build dependency on consumers. Replace the friend boundary with a POD `Values` bundle + a constructor, so the installed header has no nlohmann token.

- [ ] **Step 1: Remove the json_fwd include and the cpjson friend forward-declaration**

In `Ancillaries.h`, delete lines 11-23 (the `// Forward declaration …` comment block, `#include <nlohmann/json_fwd.hpp>`, the `namespace CoolProp { class SaturationAncillaryFunction; }` forward decl, and the `namespace cpjson { … make_saturation_ancillary … }` block).

- [ ] **Step 2: Make the type enum public and add the POD Values bundle + constructor declaration**

In `SaturationAncillaryFunction`, move the `ancillaryfunctiontypes` enum (108-114) into the `public:` section, and add a POD `Values` bundle plus a constructor that consumes it. Replace the `private:`…`public:` boundary around the enum and the existing default ctor (89-126) so it reads:

```cpp
class SaturationAncillaryFunction
{
   public:
    enum ancillaryfunctiontypes
    {
        TYPE_NOT_SET = 0,
        TYPE_NOT_EXPONENTIAL,     ///< It is a non-exponential type of equation
        TYPE_EXPONENTIAL,         ///< It is an exponential type equation, with or without the T_c/T term
        TYPE_RATIONAL_POLYNOMIAL  ///< It is a rational polynomial equation
    };

    /// Plain-typed bundle of already-parsed ancillary values. This is the
    /// de-publishing boundary: the non-installed factory
    /// cpjson::make_saturation_ancillary (FluidLibraryFactories.h) does all the
    /// nlohmann parsing and hands one of these across, keeping every JSON type
    /// out of this installed header.
    struct Values
    {
        ancillaryfunctiontypes type = TYPE_NOT_SET;
        Eigen::MatrixXd num_coeffs, den_coeffs;
        std::vector<double> n, t;
        CoolPropDbl max_abs_error = _HUGE;
        bool using_tau_r = false;
        CoolPropDbl reducing_value = _HUGE, T_r = _HUGE;
        CoolPropDbl Tmin = _HUGE, Tmax = _HUGE;
    };

   private:
    Eigen::MatrixXd num_coeffs,   ///< Coefficients for numerator in rational polynomial
      den_coeffs;                 ///< Coefficients for denominator in rational polynomial
    std::vector<double> n, t, s;  // For TYPE_NOT_EXPONENTIAL & TYPE_EXPONENTIAL
    union
    {
        CoolPropDbl max_abs_error;  ///< For TYPE_RATIONAL_POLYNOMIAL
        struct
        {
            bool using_tau_r;
            CoolPropDbl reducing_value, T_r;
            std::size_t N;
        };
    };
    CoolPropDbl Tmax, Tmin;
    ancillaryfunctiontypes type;  ///< The type of ancillary curve being used

   public:
    SaturationAncillaryFunction()
      : type(TYPE_NOT_SET), Tmin(_HUGE), Tmax(_HUGE){};

    /// Build from a plain-typed Values bundle (defined out-of-line in
    /// Ancillaries.cpp). All JSON parsing happens in the non-installed factory.
    explicit SaturationAncillaryFunction(const Values& v);
```
Then continue with the rest of the existing class body (`bool enabled()` onward at line 128) unchanged. Delete the old `friend …` declaration (126).

- [ ] **Step 3: Define the constructor out-of-line in Ancillaries.cpp**

Add to `src/Backends/Helmholtz/Fluids/Ancillaries.cpp` (in `namespace CoolProp`), reproducing the field assignments the factory used to make, including the computed `N` and `s`:
```cpp
SaturationAncillaryFunction::SaturationAncillaryFunction(const Values& v) {
    type = v.type;
    Tmin = v.Tmin;
    Tmax = v.Tmax;
    if (type == TYPE_RATIONAL_POLYNOMIAL) {
        num_coeffs = v.num_coeffs;
        den_coeffs = v.den_coeffs;
        max_abs_error = v.max_abs_error;
    } else {
        n = v.n;
        t = v.t;
        N = n.size();
        s = n;
        reducing_value = v.reducing_value;
        using_tau_r = v.using_tau_r;
        T_r = v.T_r;
    }
}
```
(The `union` is written through the matching active member per branch — exactly as the old friend factory did.)

- [ ] **Step 4: Rewrite make_saturation_ancillary to build a Values and call the ctor**

In `FluidLibraryFactories.h`, replace `make_saturation_ancillary` (32-70) so it populates a `SaturationAncillaryFunction::Values` and constructs from it (no private-member access, so the friend is gone):
```cpp
/// Build a SaturationAncillaryFunction from its JSON description. Parses into a
/// plain-typed SaturationAncillaryFunction::Values bundle and constructs from
/// it, so no nlohmann type crosses the installed-header boundary.
inline CoolProp::SaturationAncillaryFunction make_saturation_ancillary(const nlohmann::json& j) {
    CoolProp::SaturationAncillaryFunction::Values vals;
    std::string type = cpjson::get_string(j, "type");
    if (!type.compare("rational_polynomial")) {
        vals.type = CoolProp::SaturationAncillaryFunction::TYPE_RATIONAL_POLYNOMIAL;
        vals.num_coeffs = CoolProp::vec_to_eigen(cpjson::get_double_array(j.at("A")));
        vals.den_coeffs = CoolProp::vec_to_eigen(cpjson::get_double_array(j.at("B")));
        vals.max_abs_error = cpjson::get_double(j, "max_abs_error");
        try {
            vals.Tmin = cpjson::get_double(j, "Tmin");
            vals.Tmax = cpjson::get_double(j, "Tmax");
        } catch (...) {
            vals.Tmin = _HUGE;
            vals.Tmax = _HUGE;
        }
    } else {
        if (!type.compare("rhoLnoexp"))
            vals.type = CoolProp::SaturationAncillaryFunction::TYPE_NOT_EXPONENTIAL;
        else
            vals.type = CoolProp::SaturationAncillaryFunction::TYPE_EXPONENTIAL;
        vals.n = cpjson::get_double_array(j.at("n"));
        vals.t = cpjson::get_double_array(j.at("t"));
        if (vals.n.size() != vals.t.size()) {
            throw CoolProp::ValueError("Ancillary 'n' and 't' arrays must have equal length");
        }
        vals.Tmin = cpjson::get_double(j, "Tmin");
        vals.Tmax = cpjson::get_double(j, "Tmax");
        vals.reducing_value = cpjson::get_double(j, "reducing_value");
        vals.using_tau_r = cpjson::get_bool(j, "using_tau_r");
        vals.T_r = cpjson::get_double(j, "T_r");
    }
    return CoolProp::SaturationAncillaryFunction(vals);
}
```
Update the file-top comment block (4-8) to drop the "friends these free functions … populate the otherwise-private … fields directly" wording (the `make_surface_tension_correlation` factory already needs no friend since `SurfaceTensionCorrelation` has public members; `make_saturation_ancillary` now also needs none).

- [ ] **Step 5: Build**

Run: `cmake --build build_catch --target CatchTestRunner -j8`
Expected: clean build.

- [ ] **Step 6: Run the ancillary-touching suites (behavior parity)**

Run: `./build_catch/CatchTestRunner "[ancillaries],[fast]" ; ./build_catch/CatchTestRunner "[ancillary]"`
Expected: pass. If neither tag matches in this checkout, run a saturation-pressure smoke instead: `./build_catch/CatchTestRunner "[helmholtz]"` and confirm no regressions. The decisive parity check is that ancillary-derived saturation values are unchanged.

- [ ] **Step 7: Verify the installed header is nlohmann-free**

Run: `grep -n "nlohmann\|json_fwd\|rapidjson" include/CoolProp/fluids/Ancillaries.h`
Expected: no matches.

- [ ] **Step 8: Commit**

```bash
git add include/CoolProp/fluids/Ancillaries.h src/Backends/Helmholtz/Fluids/Ancillaries.cpp src/Backends/Helmholtz/Fluids/FluidLibraryFactories.h
git commit -m "refactor(json): Ancillaries.h POD-ctor boundary — drop json_fwd downstream dep (CoolProp-xa8w.1)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: PC-SAFT loader → nlohmann + PCSAFTFluid POD ctor + factory + Valijson

**Files:**
- Modify: `include/CoolProp/fluids/PCSAFTFluid.h:8 (drop rapidjson include), 35-37 (replace JSON ctor with POD ctor)`
- Create: `src/Backends/PCSAFT/PCSAFTFluidFactory.h` (non-installed `cpjson::make_pcsaft_fluid`)
- Modify: `src/Backends/PCSAFT/PCSAFTFluid.cpp` (drop the rapidjson ctor; keep calc_water_sigma)
- Modify: `src/Backends/PCSAFT/PCSAFTLibrary.h:10,28,38 (includes + member signatures)`
- Modify: `src/Backends/PCSAFT/PCSAFTLibrary.cpp:7 (include), 30-56 (validation→Valijson), 111-213 (add_many), 340-403 (binary pairs)`

- [ ] **Step 1: Replace the JSON ctor with a POD ctor in the installed header**

In `PCSAFTFluid.h`, delete line 8 (`#include "CoolProp/detail/rapidjson.h"`). Replace the JSON constructor (36):
```cpp
    PCSAFTFluid(rapidjson::Value::ValueIterator itr);
```
with a POD constructor:
```cpp
    PCSAFTFluid(std::string name, std::string CAS, CoolPropDbl molemass,
                std::vector<std::string> aliases, PCSAFTValues params)
      : name(std::move(name)), CAS(std::move(CAS)), molemass(molemass),
        aliases(std::move(aliases)), params(std::move(params)) {}
```
(`PCSAFTValues` is already declared above in this header. No JSON type appears.)

- [ ] **Step 2: Create the non-installed factory PCSAFTFluidFactory.h**

Create `src/Backends/PCSAFT/PCSAFTFluidFactory.h` carrying the old ctor body, idiom-mapped, building the POD pieces and calling the POD ctor:
```cpp
#ifndef PCSAFTFLUIDFACTORY_H
#define PCSAFTFLUIDFACTORY_H

// NON-INSTALLED factory that builds a PCSAFTFluid from nlohmann::json. Lives
// under src/ so the nlohmann::json type never appears in the installed
// PCSAFTFluid.h. Mirrors the (now-removed) rapidjson constructor.

#include "CoolProp/detail/json.h"
#include "CoolProp/fluids/PCSAFTFluid.h"

namespace cpjson {

inline CoolProp::PCSAFTFluid make_pcsaft_fluid(const nlohmann::json& fluid) {
    CoolProp::PCSAFTValues params;
    params.m = cpjson::get_double(fluid, "m");
    params.sigma = cpjson::get_double(fluid, "sigma");
    params.u = cpjson::get_double(fluid, "u");

    params.uAB = (fluid.contains("uAB") && fluid.at("uAB").is_number()) ? cpjson::get_double(fluid, "uAB") : 0.;
    params.volA = (fluid.contains("volA") && fluid.at("volA").is_number()) ? cpjson::get_double(fluid, "volA") : 0.;
    params.assocScheme = fluid.contains("assocScheme") ? cpjson::get_string_array(fluid, "assocScheme") : std::vector<std::string>{};
    params.dipm = (fluid.contains("dipm") && fluid.at("dipm").is_number()) ? cpjson::get_double(fluid, "dipm") : 0.;
    params.dipnum = (fluid.contains("dipnum") && fluid.at("dipnum").is_number()) ? cpjson::get_double(fluid, "dipnum") : 0.;
    params.z = (fluid.contains("charge") && fluid.at("charge").is_number()) ? cpjson::get_double(fluid, "charge") : 0.;

    return CoolProp::PCSAFTFluid(cpjson::get_string(fluid, "name"), cpjson::get_string(fluid, "CAS"),
                                 cpjson::get_double(fluid, "molemass"), cpjson::get_string_array(fluid, "aliases"),
                                 std::move(params));
}

}  // namespace cpjson

#endif  // PCSAFTFLUIDFACTORY_H
```

- [ ] **Step 3: Remove the rapidjson ctor from PCSAFTFluid.cpp**

In `PCSAFTFluid.cpp`, delete the `#include "CoolProp/detail/rapidjson.h"` (6) and the entire `PCSAFTFluid::PCSAFTFluid(rapidjson::Value::ValueIterator itr)` definition (10-55). Keep `calc_water_sigma` (57-65) and the includes it needs (`<cmath>`, `PCSAFTFluid.h`).

- [ ] **Step 4: Convert the PCSAFTLibrary header**

In `PCSAFTLibrary.h`: replace line 10 `#include "CoolProp/detail/rapidjson.h"` with `#include "CoolProp/detail/json.h"`. Retype:
```cpp
    void load_from_JSON(const nlohmann::json& doc);   // was rapidjson::Document&
```
```cpp
    int add_many(const nlohmann::json& listing);       // was rapidjson::Value&
```

- [ ] **Step 5: Convert validation to Valijson + parse in PCSAFTLibrary.cpp**

In `PCSAFTLibrary.cpp`: replace line 7 `#include "CoolProp/detail/rapidjson.h"` with `#include "CoolProp/detail/json.h"`, and add `#include "PCSAFTFluidFactory.h"`. Rewrite `add_fluids_from_JSON_string` (30-56), dropping the function-pointer cast (the rapidjson overload is no longer visible, so the call binds to Valijson):
```cpp
void add_fluids_from_JSON_string(PCSAFTLibraryClass& dest, const std::string_view& JSON) {
    std::string errstr;
    cpjson::schema_validation_code val_code = cpjson::validate_schema(pcsaft_fluids_schema_JSON, JSON, errstr);
    if (val_code == cpjson::SCHEMA_VALIDATION_OK) {
        nlohmann::json dd = cpjson::parse(JSON);
        try {
            dest.add_many(dd);
        } catch (std::exception& e) {
            std::cout << e.what() << '\n';
        }
    } else {
        if (get_debug_level() > 0) {
            throw ValueError(format("Unable to load PC-SAFT library with error: %s", errstr.c_str()));
        }
    }
}
```

- [ ] **Step 6: Convert add_many in PCSAFTLibrary.cpp**

Change the signature (111) and the loop/construction (114-116); the rest of the body (the already-present checks, removal, index bookkeeping) is unchanged because it uses `fluid.getName()`/`getCAS()`/`getAliases()`:
```cpp
int PCSAFTLibraryClass::add_many(const nlohmann::json& listing) {
    int counter = 0;
    std::string fluid_name;
    for (const auto& fluid_json : listing) {
        try {
            PCSAFTFluid fluid = cpjson::make_pcsaft_fluid(fluid_json);
            fluid_name = fluid.getName();
            // … (lines 119-207 unchanged) …
        } catch (const std::exception& e) {
            throw ValueError(format("Unable to load fluid [%s] due to error: %s", fluid_name.c_str(), e.what()));
        }
    }
    return counter;
}
```

- [ ] **Step 7: Convert the binary-pairs loaders in PCSAFTLibrary.cpp**

`load_from_JSON` (340-394): change `rapidjson::Document& doc` → `const nlohmann::json& doc`, the loop to `for (const auto& el : doc)`, each `*itr` → `el`, and the two `itr->HasMember("kij"/"kijT")` → `el.contains("kij"/"kijT")`. The `cpjson::get_string(el, …)` / `cpjson::get_double(el, …)` calls are otherwise unchanged.
`load_from_string` (396-403):
```cpp
void PCSAFTLibraryClass::load_from_string(const std::string_view& str) {
    nlohmann::json doc = cpjson::parse(str);
    load_from_JSON(doc);
}
```

- [ ] **Step 8: Build**

Run: `cmake --build build_catch --target CatchTestRunner -j8`
Expected: clean build.

- [ ] **Step 9: Run the PCSAFT suite (behavior parity)**

Run: `./build_catch/CatchTestRunner "[PCSAFT]"`
Expected: pass (the `[Qmass][PCSAFT][mixture]` test and any PCSAFT flash/property checks load fluids identically).

- [ ] **Step 10: Verify rapidjson-free + installed header clean**

Run: `grep -rn "rapidjson" src/Backends/PCSAFT/ ; grep -n "nlohmann\|rapidjson" include/CoolProp/fluids/PCSAFTFluid.h`
Expected: no matches in either.

- [ ] **Step 11: Commit**

```bash
git add include/CoolProp/fluids/PCSAFTFluid.h src/Backends/PCSAFT/PCSAFTFluidFactory.h src/Backends/PCSAFT/PCSAFTFluid.cpp src/Backends/PCSAFT/PCSAFTLibrary.h src/Backends/PCSAFT/PCSAFTLibrary.cpp
git commit -m "refactor(json): PC-SAFT loader off rapidjson; PCSAFTFluid POD ctor + Valijson (CoolProp-xa8w.1)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Cubics/UNIFAC loaders → nlohmann + delete the alpha0 bridges + Valijson

**Files:**
- Modify: `src/Backends/Cubics/CubicsLibrary.cpp:6 (include), 29-94 (add_many + alpha0 bridge), 111-116 (get_JSONstring), 144-172 (validation→Valijson)`
- Modify: `src/Backends/Cubics/UNIFACLibrary.h:6,91,94 (include + populate/jsonize sigs)`
- Modify: `src/Backends/Cubics/UNIFACLibrary.cpp:7-90 (jsonize/populate + alpha0 bridge)`

- [ ] **Step 1: Swap the Cubics include**

In `CubicsLibrary.cpp`, replace line 6 `#include "CoolProp/detail/rapidjson.h"` with `#include "CoolProp/detail/json.h"`.

- [ ] **Step 2: Convert CubicsLibraryClass::add_many and delete its alpha0 bridge**

Replace `add_many` (29-94). Note: `alpha0` now passes the nlohmann node directly to `parse_alpha0` (bridge deleted), and `json2string(*itr)` → `el.dump()`:
```cpp
    int add_many(const nlohmann::json& listing) {
        int counter = 0;
        for (const auto& el : listing) {
            CubicsValues val;
            val.Tc = cpjson::get_double(el, "Tc");
            val.pc = cpjson::get_double(el, "pc");
            val.acentric = cpjson::get_double(el, "acentric");
            val.molemass = cpjson::get_double(el, "molemass");
            val.name = cpjson::get_string(el, "name");
            val.aliases = cpjson::get_string_array(el, "aliases");
            val.CAS = cpjson::get_string(el, "CAS");
            if (el.contains("rhomolarc") && el.at("rhomolarc").is_number()) {
                val.rhomolarc = cpjson::get_double(el, "rhomolarc");
            }
            if (el.contains("alpha") && el.at("alpha").is_object()) {
                const nlohmann::json& alpha = el.at("alpha");
                val.alpha_type = cpjson::get_string(alpha, "type");
                val.alpha_coeffs = cpjson::get_double_array(alpha, "c");
            } else {
                val.alpha_type = "default";
            }
            if (el.contains("alpha0") && el.at("alpha0").is_array()) {
                val.alpha0 = JSONFluidLibrary::parse_alpha0(el.at("alpha0"));
            }
            std::pair<std::map<std::string, CubicsValues>::iterator, bool> ret;
            ret = fluid_map.emplace(upper(val.name), val);
            if (ret.second == false && get_config_bool(OVERWRITE_FLUIDS)) {
                fluid_map.erase(ret.first);
                ret = fluid_map.emplace(upper(val.name), val);
                if (get_debug_level() > 0) {
                    std::cout << "added the cubic fluid: " + val.name << '\n';
                }
                assert(ret.second == true);
            }
            for (auto it = val.aliases.begin(); it != val.aliases.end(); ++it) {
                if (aliases_map.find(*it) == aliases_map.end()) {
                    aliases_map.emplace(*it, upper(val.name));
                }
            }
            std::pair<std::map<std::string, std::string>::iterator, bool> addJson;
            addJson = JSONstring_map.emplace(upper(val.name), el.dump());
            if (addJson.second == false && get_config_bool(OVERWRITE_FLUIDS)) {
                JSONstring_map.erase(addJson.first);
                addJson = JSONstring_map.emplace(upper(val.name), el.dump());
                if (get_debug_level() > 0) {
                    std::cout << "added the cubic fluid: " + val.name << '\n';
                }
                assert(addJson.second == true);
            }
            counter++;
        }
        return counter;
    };
```

- [ ] **Step 3: Convert get_JSONstring (wrap-in-array) in CubicsLibrary.cpp**

Replace the rapidjson array-wrapping (111-116) with nlohmann:
```cpp
        // Then, wrap the single fluid object in an array, as callers expect
        nlohmann::json doc = cpjson::parse(it->second);
        nlohmann::json doc2 = nlohmann::json::array();
        doc2.push_back(doc);
        return doc2.dump();
```
(Deletes the `cpjson::JSON_string_to_rapidjson` / `doc.GetAllocator()` calls.)

- [ ] **Step 4: Convert Cubics validation to Valijson + parse**

Rewrite `add_fluids_from_JSON_string` (148-172), dropping the function-pointer cast:
```cpp
void add_fluids_from_JSON_string(CubicsLibraryClass& dest, const std::string_view& JSON) {
    std::string errstr;
    cpjson::schema_validation_code val_code = cpjson::validate_schema(cubic_fluids_schema_JSON, JSON, errstr);
    if (val_code == cpjson::SCHEMA_VALIDATION_OK) {
        nlohmann::json dd = cpjson::parse(JSON);
        try {
            dest.add_many(dd);
        } catch (std::exception& /* e */) {
            throw ValueError(format("Unable to load cubics library with error: %s", errstr.c_str()));
        }
    } else {
        throw ValueError(format("Unable to validate cubics library against schema with error: %s", errstr.c_str()));
    }
}
```

- [ ] **Step 5: Convert the UNIFAC header**

In `UNIFACLibrary.h`: replace line 6 `#include "CoolProp/detail/rapidjson.h"` with `#include "CoolProp/detail/json.h"`. Retype:
```cpp
    void jsonize(std::string& s, nlohmann::json& doc);
```
```cpp
    /// Populate internal data structures based on nlohmann::json arrays
    void populate(const nlohmann::json& group_data, const nlohmann::json& interaction_data, const nlohmann::json& decomp_data);
```

- [ ] **Step 6: Convert UNIFACLibrary.cpp jsonize + populate, delete the alpha0 bridge**

`jsonize` (7-14):
```cpp
void UNIFACParameterLibrary::jsonize(std::string& s, nlohmann::json& d) {
    d = cpjson::parse(s);
}
```
`populate(nlohmann::json…)` (15-82): change the three params to `const nlohmann::json&`; convert each `for (rapidjson::Value::ValueIterator itr = X.Begin(); …)` to `for (const auto& el : X)`; each `(*itr)["k"].GetInt()/GetDouble()/GetString()` → `el.at("k").get<int>()/get<double>()/get<std::string>()`; `(*itr).HasMember("userid")` → `el.contains("userid")`; the `alpha` block uses `el.at("alpha")` with `is_object()`; **delete the alpha0 bridge** (63-70) and replace with the direct call:
```cpp
        if (el.contains("alpha0") && el.at("alpha0").is_array()) {
            c.alpha0 = CoolProp::JSONFluidLibrary::parse_alpha0(el.at("alpha0"));
        }
```
and the nested `groups` loop:
```cpp
        const nlohmann::json& groups = el.at("groups");
        for (const auto& g : groups) {
            int count = g.at("count").get<int>();
            int sgi = g.at("sgi").get<int>();
            if (has_group(sgi)) {
                ComponentGroup cg(count, get_group(sgi));
                c.groups.push_back(cg);
            }
        }
```
`populate(std::string…)` (83-92):
```cpp
void UNIFACParameterLibrary::populate(std::string& group_data, std::string& interaction_data, std::string& decomp_data) {
    nlohmann::json group_JSON, interaction_JSON, decomp_JSON;
    jsonize(group_data, group_JSON);
    jsonize(interaction_data, interaction_JSON);
    jsonize(decomp_data, decomp_JSON);
    populate(group_JSON, interaction_JSON, decomp_JSON);
    m_populated = true;
}
```

- [ ] **Step 7: Build**

Run: `cmake --build build_catch --target CatchTestRunner -j8`
Expected: clean build.

- [ ] **Step 8: Run the Cubics + UNIFAC suites (behavior parity)**

Run: `./build_catch/CatchTestRunner "[cubic],[cubic_superanc],[cubic_DmolarT],[UNIFAC]"`
Expected: pass (the Poling UNIFAC example `[UNIFAC]`, cubic superancillary `[2739]`, and cubic DmolarT `[2673]` checks load identically).

- [ ] **Step 9: Verify rapidjson-free + bridges gone**

Run: `grep -rn "rapidjson\|Removed once\|migrates off rapidjson\|validate_schema_rapidjson" src/Backends/Cubics/`
Expected: no matches.

- [ ] **Step 10: Commit**

```bash
git add src/Backends/Cubics/CubicsLibrary.cpp src/Backends/Cubics/UNIFACLibrary.h src/Backends/Cubics/UNIFACLibrary.cpp
git commit -m "refactor(json): Cubics/UNIFAC loaders off rapidjson; delete alpha0 bridges + Valijson (CoolProp-xa8w.1)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: User-fluid validation test (guard the schema-engine swap)

**Files:**
- Modify: `src/Tests/CoolProp-Tests.cpp` (add a new `TEST_CASE` near the other `[PCSAFT]`/`[cubic]` cases)

The user-supplied-fluid validation path moved from rapidjson `SchemaValidator` to Valijson. This test feeds a structurally-invalid payload to the public `add_fluids_as_JSON` entry points and asserts a clean throw — proving the Valijson path rejects bad input rather than crashing or silently accepting.

- [ ] **Step 1: Write the test**

Add to `src/Tests/CoolProp-Tests.cpp`:
```cpp
TEST_CASE("User-fluid schema validation rejects malformed payloads (Valijson)", "[json_validation]") {
    SECTION("PC-SAFT: missing required fields throws") {
        // valid JSON, but a fluid object missing required PC-SAFT parameters (m, sigma, u, …)
        std::string bad = R"([{"name": "Bogus", "CAS": "000-00-0"}])";
        CHECK_THROWS(CoolProp::PCSAFTLibrary::add_fluids_as_JSON(bad));
    }
    SECTION("PC-SAFT: structurally invalid JSON throws") {
        std::string notjson = R"(this is not json)";
        CHECK_THROWS(CoolProp::PCSAFTLibrary::add_fluids_as_JSON(notjson));
    }
    SECTION("Cubic: missing required fields throws") {
        std::string bad = R"([{"name": "Bogus"}])";
        CHECK_THROWS(CoolProp::CubicLibrary::add_fluids_as_JSON(bad));
    }
}
```
Ensure the test file already includes the loader headers; if not, add `#include "Backends/PCSAFT/PCSAFTLibrary.h"` and `#include "Backends/Cubics/CubicsLibrary.h"` near the other backend includes. (The cubic loader throws on schema-validation failure unconditionally; the PCSAFT loader throws on bad input when `get_debug_level() > 0`, but `add_fluids_as_JSON` parse failure throws regardless via `cpjson::parse` — the missing-fields case relies on the schema gate. If the PCSAFT missing-fields `SECTION` does not throw at default debug level, wrap it with `set_debug_level(1)` before and restore after.)

- [ ] **Step 2: Build**

Run: `cmake --build build_catch --target CatchTestRunner -j8`
Expected: clean build.

- [ ] **Step 3: Run the new test**

Run: `./build_catch/CatchTestRunner "[json_validation]" -s`
Expected: all three sections pass (each `CHECK_THROWS` fires). If the PCSAFT missing-fields section does not throw, apply the `set_debug_level(1)` guard noted in Step 1 and re-run.

- [ ] **Step 4: Commit**

```bash
git add src/Tests/CoolProp-Tests.cpp
git commit -m "test(json): user-fluid schema validation rejects malformed PCSAFT/cubic payloads (CoolProp-xa8w.1)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final Verification (before PR)

- [ ] **Step 1: Full fast suite**

Run: `./build_catch/CatchTestRunner "~[slow]"`
Expected: no regressions across all backends.

- [ ] **Step 2: Symbol-leak gate (shared library)**

Build the shared library and run the gate (it greps `nm -D` for exported `nlohmann`/`valijson`):
Run: `./dev/ci/check-json-symbols.sh`
Expected: PASS — zero exported nlohmann/valijson symbols. (If `make_pcsaft_fluid` or a factory leaks, confirm it is `inline` and only used within `libCoolProp`; mark `CP_JSON_LOCAL` on any out-of-line method that would otherwise emit a symbol.)

- [ ] **Step 3: Installed-header cleanliness**

Run: `grep -rn "nlohmann\|rapidjson\|json_fwd" include/CoolProp/fluids/PCSAFTFluid.h include/CoolProp/fluids/Ancillaries.h`
Expected: no matches.

- [ ] **Step 4: Loaders rapidjson-free**

Run: `grep -rln "rapidjson" src/Backends/Incompressible/ src/Backends/PCSAFT/ src/Backends/Cubics/`
Expected: no matches.

- [ ] **Step 5: Preflight**

Run: `./dev/ci/preflight.sh`
Expected: PASS (clang-format, build, the SBTL umbrella tag is not implicated here but the changed-path auto-selection will pick the right Catch scope; cppcheck/clang-tidy/semgrep on changed files clean).

- [ ] **Step 6: Adversarial review (MANDATORY before `gh pr create`)**

Per CLAUDE.md, dispatch the review subagent against the diff vs `master` (this branch's base). Address or justify every blocking finding, paying attention to: the optional-member `contains`-guard parity (a dropped guard turns a missing key into a thrown exception — a behavior change), the union write in the `SaturationAncillaryFunction(Values)` ctor (correct active member per branch), and any factory that could export an nlohmann symbol.

- [ ] **Step 7: Update the bd issue + finish the branch**

Use superpowers:finishing-a-development-branch to push and open the PR. After PR creation, re-review any commits added after Step 6 (per CLAUDE.md "re-review the delta").

---

## Self-Review (plan vs spec)

**Spec coverage:** §2 mechanical conversion → Tasks 1/3/4 + idiom map; §3 Valijson cast-deletion → Tasks 3.5, 4.4; §4 bridge deletion → Tasks 4.2, 4.6 (+ casts in 3.5/4.4); §5 POD boundary (PCSAFTFluid) → Task 3.1-3.2; §5 Ancillaries retrofit → Task 2; §6 components → all tasks; §7 four-commit sequencing → Tasks 1-4 each commit (Task 5 adds the test commit); §8 testing/symbol-gate/grep → Tasks' parity steps + Final Verification; §9 deferred consumer-compile gate → not implemented here (correctly out of scope, tracked on `xa8w.3`).

**Type consistency:** `cpjson::make_pcsaft_fluid(const nlohmann::json&) → PCSAFTFluid` (Task 3.2) matches its call site (Task 3.6) and the POD ctor (Task 3.1). `SaturationAncillaryFunction::Values` fields (Task 2.2) match the ctor assignments (Task 2.3) and the factory population (Task 2.4). `add_many(const nlohmann::json&)`, `load_from_JSON(const nlohmann::json&)`, `populate(const nlohmann::json&…)`, `jsonize(std::string&, nlohmann::json&)` signatures are consistent between their header declarations and .cpp definitions.

**Placeholder scan:** none — every code step shows complete code; the one conditional (PCSAFT debug-level guard in Task 5) is spelled out with its fallback.
