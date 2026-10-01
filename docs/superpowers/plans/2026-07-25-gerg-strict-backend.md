# GERG-2004 / GERG-2008 Strict Backend Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `GERG2004` and `GERG2008` as strict CoolProp backend families that reproduce teqp's GERG models to within 1e-10 relative.

**Architecture:** A single `GERGMixtureBackend` derived from `HelmholtzEOSMixtureBackend` builds its `CoolPropFluid` objects, `GERG2008ReducingFunction`, and `ExcessTerm` from self-contained coefficient tables in a backend-private header, bypassing CoolProp's global fluid and binary-pair JSON libraries entirely.

**Tech Stack:** C++17, CMake, Catch2, Python 3 (offline data generation only, via `teqp>=0.18`).

**Spec:** `docs/superpowers/specs/2026-07-25-gerg-strict-backend-design.md`
**Issue:** `bd CoolProp-p8ub`

## Global Constraints

- C++17 (`CMAKE_CXX_STANDARD 17`, CMakeLists.txt:104).
- `GERGData.h` lives under `src/Backends/GERG/` and is **never** installed or included from `include/CoolProp/`. It is private to the backend.
- Gas constant: `R = 8.314472 J/mol/K` for every GERG component. Ideal-gas ratio `R*/R = 8.314510 / 8.314472`.
- Reference implementation and source of every coefficient: `/Users/ianbell/Code/teqp/include/teqp/models/GERG/GERG.hpp` (teqp 0.23.1, local checkout). Every table transcribed from it carries a comment naming the teqp source lines and the GERG monograph table.
- teqp is NIST-authored and the coefficients are published data from the GERG-2004 monograph and Kunz & Wagner (2012). Each table block gets an attribution comment; no license header is copied.
- SI only. Densities `mol/m^3`, temperatures `K`, pressures `Pa`, molar mass `kg/mol`.
- Catch2 tag for every test in this plan: `[GERG]`. Model-specific sub-tags `[GERG2004]` / `[GERG2008]` where useful.
- Build and test commands, used verbatim throughout:
  ```bash
  cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
  cmake --build build_catch --target CatchTestRunner -j8
  ./build_catch/CatchTestRunner "[GERG]"
  ```
- `.beads/issues.jsonl` must not appear in these commits. Before each commit:
  ```bash
  git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
  ```
  and commit with `--no-verify`. Run `./dev/ci/preflight.sh` before pushing.
- Before `gh pr create`, the mandatory adversarial review in CLAUDE.md applies.

---

## File Structure

| File | Responsibility |
|---|---|
| `src/Backends/GERG/GERGData.h` | All coefficient tables and lookup functions. Namespaces `CoolProp::GERG::GERG2004` / `GERG2008`. No CoolProp types beyond `ValueError`. |
| `src/Backends/GERG/GERGBackend.h` | `GERGModel` enum, `GERGMixtureBackend` declaration. |
| `src/Backends/GERG/GERGBackend.cpp` | Fluid construction, mixture wiring, strictness, generator registration. |
| `src/Backends/GERG/GERGReferenceValues.h` | Generated test fixtures (teqp reference values). Test-only. |
| `src/Tests/CoolProp-Tests-GERG.cpp` | All `[GERG]` tests. |
| `dev/gerg/generate_reference_values.py` | Offline: emits `GERGReferenceValues.h` from teqp. |
| `dev/gerg/fit_ancillaries.py` | Offline: emits the ancillary coefficient table. |
| `include/CoolProp/DataStructures.h` | Two new `backend_families`, two new `backends` enumerators. |
| `src/DataStructures.cpp` | Registry entries for the above. |
| `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.h` | `set_mixture_parameters()` becomes `virtual`. |
| `CMakeLists.txt` | `GERG` in `COOLPROP_ENABLED_BACKENDS`; test source appended. |

---

## Task 1: Build wiring and backend registration

Establishes the factory path end to end before any thermodynamics exists. A `GERG2008` factory call must reach *our* code and fail there, not fail in `extract_backend_families`.

**Files:**
- Modify: `include/CoolProp/DataStructures.h:493-526`
- Modify: `src/DataStructures.cpp:803-821`
- Create: `src/Backends/GERG/GERGBackend.h`
- Create: `src/Backends/GERG/GERGBackend.cpp`
- Modify: `CMakeLists.txt:330-337` (backend list), `CMakeLists.txt:2296` (test source)
- Test: `src/Tests/CoolProp-Tests-GERG.cpp`

**Interfaces:**
- Consumes: nothing.
- Produces: `CoolProp::GERGModel` enum (`GERG_2004`, `GERG_2008`); `CoolProp::GERGMixtureBackend(GERGModel model, const std::vector<std::string>& names)`; enumerators `GERG2004_BACKEND_FAMILY`, `GERG2008_BACKEND_FAMILY`, `GERG2004_BACKEND`, `GERG2008_BACKEND`.

- [ ] **Step 1: Write the failing test**

Create `src/Tests/CoolProp-Tests-GERG.cpp`. Every existing test file in
`src/Tests/` wraps its whole body in `#if defined(ENABLE_CATCH)` — match that,
and close the guard at the end of the file. The repo root is on the include
path (`CMakeLists.txt:370`), so backend headers can be included as
`"src/Backends/GERG/GERGData.h"`; check what neighbouring test files do and
follow them.

```cpp
#if defined(ENABLE_CATCH)
#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "CoolProp/AbstractState.h"
#include "CoolProp/DataStructures.h"
#include "CoolProp/Exceptions.h"

using namespace CoolProp;

TEST_CASE("GERG backend families are registered", "[GERG]") {
    backend_families f1 = INVALID_BACKEND_FAMILY, f2 = INVALID_BACKEND_FAMILY;

    extract_backend_families("GERG2004", f1, f2);
    CHECK(f1 == GERG2004_BACKEND_FAMILY);

    extract_backend_families("GERG2008", f1, f2);
    CHECK(f1 == GERG2008_BACKEND_FAMILY);
}

TEST_CASE("GERG factory reaches the GERG backend", "[GERG]") {
    // Task 1 ships a skeleton whose constructor throws NotImplementedError.
    // The point of this test is that we get THAT exception, proving dispatch
    // works, rather than a ValueError about an unknown backend.
    CHECK_THROWS_AS(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}), NotImplementedError);
    CHECK_THROWS_AS(AbstractState::factory("GERG2004", std::vector<std::string>{"Methane"}), NotImplementedError);
}

#endif /* ENABLE_CATCH */
```

**This test is scaffolding with a scheduled death.** Task 8 implements the
constructor, at which point asserting `NotImplementedError` becomes false and
this test must be deleted — Task 8 Step 1 says so explicitly. Do not leave it
behind; a test asserting that a feature is unimplemented, kept past the point
where it is implemented, is worse than no test.

- [ ] **Step 2: Run test to verify it fails**

```bash
cmake --build build_catch --target CatchTestRunner -j8
```
Expected: compile failure — `GERG2004_BACKEND_FAMILY` not declared.

- [ ] **Step 3: Add the enumerators**

In `include/CoolProp/DataStructures.h`, append to `enum backend_families` (after `SVDSBTL_BACKEND_FAMILY`, keeping existing values stable):

```cpp
    SVDSBTL_BACKEND_FAMILY,
    GERG2004_BACKEND_FAMILY,
    GERG2008_BACKEND_FAMILY
};
```

and to `enum backends` (after `SVDSBTL_BACKEND`):

```cpp
    SVDSBTL_BACKEND,
    GERG2004_BACKEND,
    GERG2008_BACKEND
};
```

Appending only — never reorder. These values are ABI for the wrappers.

- [ ] **Step 4: Add the registry entries**

In `src/DataStructures.cpp`, append to `backend_family_list`:

```cpp
  {SVDSBTL_BACKEND_FAMILY, "SVDSBTL"}, {GERG2004_BACKEND_FAMILY, "GERG2004"}, {GERG2008_BACKEND_FAMILY, "GERG2008"}};
```

and to `backend_list`:

```cpp
                                                {SVDSBTL_BACKEND, "SVDSBTLBackend", SVDSBTL_BACKEND_FAMILY},
                                                {GERG2004_BACKEND, "GERG2004Backend", GERG2004_BACKEND_FAMILY},
                                                {GERG2008_BACKEND, "GERG2008Backend", GERG2008_BACKEND_FAMILY}};
```

- [ ] **Step 5: Create the skeleton backend header**

`src/Backends/GERG/GERGBackend.h`:

```cpp
#ifndef GERGBACKEND_H_
#define GERGBACKEND_H_

#include <string>
#include <vector>

#include "../Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "CoolProp/DataStructures.h"

namespace CoolProp {

/// Which GERG model year this backend instance represents.
enum class GERGModel
{
    GERG_2004,
    GERG_2008
};

/**
 * \brief Strict GERG-2004 / GERG-2008 backend.
 *
 * Admits only the components published with the selected model (18 for
 * GERG-2004, 21 for GERG-2008) and uses only parameters published with it.
 * See docs/superpowers/specs/2026-07-25-gerg-strict-backend-design.md.
 */
class GERGMixtureBackend : public HelmholtzEOSMixtureBackend
{
   public:
    GERGMixtureBackend(GERGModel model, const std::vector<std::string>& names);
    ~GERGMixtureBackend() override = default;

    std::string backend_name() override {
        return get_backend_string(m_model == GERGModel::GERG_2004 ? GERG2004_BACKEND : GERG2008_BACKEND);
    }

    GERGModel model() const {
        return m_model;
    }

   private:
    GERGModel m_model;
};

} /* namespace CoolProp */
#endif /* GERGBACKEND_H_ */
```

- [ ] **Step 6: Create the skeleton backend source with generator registration**

`src/Backends/GERG/GERGBackend.cpp`:

```cpp
#include "GERGBackend.h"

#include "CoolProp/Exceptions.h"

namespace CoolProp {

GERGMixtureBackend::GERGMixtureBackend(GERGModel model, const std::vector<std::string>& names) : m_model(model) {
    (void)names;
    throw NotImplementedError("GERG backend is not yet implemented");
}

class GERG2004Generator : public AbstractStateGenerator
{
   public:
    AbstractState* get_AbstractState(const std::vector<std::string>& fluid_names) override {
        return new GERGMixtureBackend(GERGModel::GERG_2004, fluid_names);
    }
};
static GeneratorInitializer<GERG2004Generator> gerg2004_gen(GERG2004_BACKEND_FAMILY);

class GERG2008Generator : public AbstractStateGenerator
{
   public:
    AbstractState* get_AbstractState(const std::vector<std::string>& fluid_names) override {
        return new GERGMixtureBackend(GERGModel::GERG_2008, fluid_names);
    }
};
static GeneratorInitializer<GERG2008Generator> gerg2008_gen(GERG2008_BACKEND_FAMILY);

} /* namespace CoolProp */
```

Check the exact virtual signature of `AbstractStateGenerator::get_AbstractState` in `include/CoolProp/AbstractState.h` near line 1685 before compiling — if the options-string overload is the pure virtual one, implement that signature and throw `NotImplementedError` when `options_json` is non-empty, matching what TTSE does in `src/AbstractState.cpp`.

- [ ] **Step 7: Wire CMake**

In `CMakeLists.txt`, add `GERG` to `COOLPROP_ENABLED_BACKENDS` (line 330):

```cmake
set(COOLPROP_ENABLED_BACKENDS
    Cubics
    IF97
    Helmholtz
    REFPROP
    Incompressible
    Tabular
    PCSAFT
    GERG)
```

and append the test source alongside the others near line 2296:

```cmake
  list(APPEND APP_SOURCES
       "${CMAKE_CURRENT_SOURCE_DIR}/src/Tests/CoolProp-Tests-GERG.cpp")
```

- [ ] **Step 8: Run tests to verify they pass**

```bash
cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_catch --target CatchTestRunner -j8
./build_catch/CatchTestRunner "[GERG]"
```
Expected: PASS, 2 test cases.

- [ ] **Step 9: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add include/CoolProp/DataStructures.h src/DataStructures.cpp src/Backends/GERG src/Tests/CoolProp-Tests-GERG.cpp CMakeLists.txt
git commit --no-verify -m "feat(gerg): register GERG2004/GERG2008 backend families (CoolProp-p8ub)"
```

---

## Task 2: Component tables and name resolution

**Files:**
- Create: `src/Backends/GERG/GERGData.h`
- Test: `src/Tests/CoolProp-Tests-GERG.cpp`

**Interfaces:**
- Consumes: `GERGModel` from Task 1.
- Produces:
  ```cpp
  namespace CoolProp { namespace GERG {
    struct PureInfo { double rhoc_molm3, Tc_K, M_kgmol; };
    const std::vector<std::string>& component_names(GERGModel model);
    PureInfo get_pure_info(GERGModel model, const std::string& gerg_name);
    std::string resolve_component(GERGModel model, const std::string& user_name);
  }}
  ```
  `resolve_component` maps a CoolProp name/alias/CAS to a GERG component name, throwing `ValueError` if unresolvable or outside the model's set.

- [ ] **Step 1: Write the failing test**

Append to `src/Tests/CoolProp-Tests-GERG.cpp`:

```cpp
#include "../Backends/GERG/GERGData.h"

using namespace CoolProp::GERG;

TEST_CASE("GERG component sets have the published sizes", "[GERG]") {
    CHECK(component_names(GERGModel::GERG_2004).size() == 18);
    CHECK(component_names(GERGModel::GERG_2008).size() == 21);
}

TEST_CASE("GERG pure info matches the monograph", "[GERG]") {
    // Table A3.5, GERG-2004 monograph. Tabulated in mol/dm^3 and kg/kmol;
    // our accessor returns mol/m^3 and kg/mol.
    auto methane = get_pure_info(GERGModel::GERG_2004, "methane");
    CHECK_THAT(methane.rhoc_molm3, Catch::Matchers::WithinRel(10.139342719e3, 1e-14));
    CHECK_THAT(methane.Tc_K, Catch::Matchers::WithinRel(190.564, 1e-14));
    CHECK_THAT(methane.M_kgmol, Catch::Matchers::WithinRel(16.042460e-3, 1e-14));

    // GERG-2008 changes carbon monoxide and isopentane relative to GERG-2004.
    CHECK_THAT(get_pure_info(GERGModel::GERG_2004, "carbonmonoxide").Tc_K, Catch::Matchers::WithinRel(132.800, 1e-14));
    CHECK_THAT(get_pure_info(GERGModel::GERG_2008, "carbonmonoxide").Tc_K, Catch::Matchers::WithinRel(132.860, 1e-14));

    // Added in GERG-2008 only.
    CHECK_THROWS_AS(get_pure_info(GERGModel::GERG_2004, "n-decane"), CoolProp::ValueError);
    CHECK_NOTHROW(get_pure_info(GERGModel::GERG_2008, "n-decane"));
}

TEST_CASE("GERG component resolution accepts CoolProp names and aliases", "[GERG]") {
    CHECK(resolve_component(GERGModel::GERG_2008, "Methane") == "methane");
    CHECK(resolve_component(GERGModel::GERG_2008, "METHANE") == "methane");
    CHECK(resolve_component(GERGModel::GERG_2008, "CO2") == "carbondioxide");
    CHECK(resolve_component(GERGModel::GERG_2008, "124-38-9") == "carbondioxide");
    CHECK(resolve_component(GERGModel::GERG_2008, "n-Butane") == "n-butane");

    // Known to CoolProp, outside the GERG set.
    CHECK_THROWS_AS(resolve_component(GERGModel::GERG_2008, "R134a"), CoolProp::ValueError);
    // In GERG-2008 but not GERG-2004.
    CHECK_THROWS_AS(resolve_component(GERGModel::GERG_2004, "HydrogenSulfide"), CoolProp::ValueError);
    CHECK_NOTHROW(resolve_component(GERGModel::GERG_2008, "HydrogenSulfide"));
    // Not a fluid at all.
    CHECK_THROWS_AS(resolve_component(GERGModel::GERG_2008, "NOT A FLUID"), CoolProp::ValueError);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cmake --build build_catch --target CatchTestRunner -j8
```
Expected: compile failure — `GERGData.h` does not exist.

- [ ] **Step 3: Create GERGData.h with the component tables**

Transcribe from `~/Code/teqp/include/teqp/models/GERG/GERG.hpp`:
- `GERG2004::component_names` — line 431
- `GERG2004::get_pure_info` — lines 434-465
- `GERG2008::component_names` — line 973
- `GERG2008::get_pure_info` — lines 975-997

Structure (the `CAS_to_gerg` table is new — it does not exist in teqp):

```cpp
#ifndef GERGDATA_H_
#define GERGDATA_H_

// Coefficient tables for the GERG-2004 and GERG-2008 equations of state.
//
// Transcribed from teqp (https://github.com/usnistgov/teqp),
// include/teqp/models/GERG/GERG.hpp, which is the reference implementation
// this backend is validated against.  Original data: Kunz, Klimeck, Wagner,
// Jaeschke, "The GERG-2004 Wide-Range Equation of State for Natural Gases and
// Other Mixtures", GERG TM15 (2007), and Kunz & Wagner, J. Chem. Eng. Data 57
// (2012) 3032-3091.
//
// This header is PRIVATE to the GERG backend.  It is not installed and must
// not be included from include/CoolProp/.

#include <map>
#include <string>
#include <vector>

#include "CoolProp/Exceptions.h"
#include "GERGBackend.h"

namespace CoolProp {
namespace GERG {

struct PureInfo
{
    double rhoc_molm3;  ///< Reducing density, mol/m^3
    double Tc_K;        ///< Reducing temperature, K
    double M_kgmol;     ///< Molar mass, kg/mol
};

namespace detail {

/// Table A3.5, GERG-2004 monograph.  Tabulated in mol/dm^3, K, kg/kmol;
/// converted to mol/m^3, K, kg/mol on the way out of get_pure_info.
inline const std::map<std::string, PureInfo>& pure_info_2004() {
    static const std::map<std::string, PureInfo> data = {
      {"methane", {10.139342719, 190.564000000, 16.042460}},
      // ... transcribe all 18 entries from teqp GERG.hpp:437-457
    };
    return data;
}

/// Entries that GERG-2008 changes or adds.  Everything else falls through
/// to pure_info_2004().
inline const std::map<std::string, PureInfo>& pure_info_2008_overrides() {
    static const std::map<std::string, PureInfo> data = {
      {"carbonmonoxide", {10.85, 132.86, 28.010100}},  // changed from GERG-2004
      {"isopentane", {3.271, 460.35, 72.148780}},      // changed from GERG-2004
      {"n-nonane", {1.81, 594.55, 128.2551}},          // new in GERG-2008
      {"n-decane", {1.64, 617.7, 142.28168}},          // new in GERG-2008
      {"hydrogensulfide", {10.19, 373.1, 34.08088}}    // new in GERG-2008
    };
    return data;
}

/// CAS number -> GERG component name.  Used by resolve_component so that
/// CoolProp aliases (CO2, R744, 124-38-9, ...) reach the right component.
inline const std::map<std::string, std::string>& cas_to_gerg() {
    static const std::map<std::string, std::string> data = {
      {"74-82-8", "methane"},        {"7727-37-9", "nitrogen"},        {"124-38-9", "carbondioxide"},
      {"74-84-0", "ethane"},         {"74-98-6", "propane"},           {"106-97-8", "n-butane"},
      {"75-28-5", "isobutane"},      {"109-66-0", "n-pentane"},        {"78-78-4", "isopentane"},
      {"110-54-3", "n-hexane"},      {"142-82-5", "n-heptane"},        {"111-65-9", "n-octane"},
      {"1333-74-0p", "hydrogen"},    {"7782-44-7", "oxygen"},          {"630-08-0", "carbonmonoxide"},
      {"7732-18-5", "water"},        {"7440-59-7", "helium"},          {"7440-37-1", "argon"},
      {"7783-06-4", "hydrogensulfide"}, {"111-84-2", "n-nonane"},      {"124-18-5", "n-decane"}};
    return data;
}

}  // namespace detail

inline const std::vector<std::string>& component_names(GERGModel model) {
    static const std::vector<std::string> names_2004 = {
      "methane", "nitrogen", "carbondioxide", "ethane", "propane", "n-butane", "isobutane", "n-pentane", "isopentane",
      "n-hexane", "n-heptane", "n-octane", "hydrogen", "oxygen", "carbonmonoxide", "water", "helium", "argon"};
    static const std::vector<std::string> names_2008 = [] {
        std::vector<std::string> v = names_2004;
        v.insert(v.end(), {"hydrogensulfide", "n-nonane", "n-decane"});
        return v;
    }();
    return (model == GERGModel::GERG_2004) ? names_2004 : names_2008;
}

inline PureInfo get_pure_info(GERGModel model, const std::string& gerg_name) {
    PureInfo data{};
    bool found = false;
    if (model == GERGModel::GERG_2008) {
        const auto& ov = detail::pure_info_2008_overrides();
        auto it = ov.find(gerg_name);
        if (it != ov.end()) {
            data = it->second;
            found = true;
        }
    }
    if (!found) {
        // GERG-2004 does not contain the three fluids added in GERG-2008.
        const auto& names = component_names(model);
        if (std::find(names.begin(), names.end(), gerg_name) == names.end()) {
            throw ValueError(format("[%s] is not a component of this GERG model", gerg_name.c_str()));
        }
        const auto& base = detail::pure_info_2004();
        auto it = base.find(gerg_name);
        if (it == base.end()) {
            throw ValueError(format("Unable to load GERG pure info for [%s]", gerg_name.c_str()));
        }
        data = it->second;
    }
    data.rhoc_molm3 *= 1000;  // mol/dm^3 -> mol/m^3
    data.M_kgmol /= 1000;     // kg/kmol -> kg/mol
    return data;
}

std::string resolve_component(GERGModel model, const std::string& user_name);

}  // namespace GERG
}  // namespace CoolProp
#endif /* GERGDATA_H_ */
```

The `1333-74-0p` CAS for hydrogen is CoolProp's normal-hydrogen CAS — confirm against `dev/fluids/Hydrogen.json` before trusting it, and likewise verify every other CAS in the table with:

```bash
python3 -c "
import CoolProp.CoolProp as CP
for n in ['Methane','Nitrogen','CarbonDioxide','Ethane','Propane','n-Butane','IsoButane',
          'n-Pentane','Isopentane','n-Hexane','n-Heptane','n-Octane','Hydrogen','Oxygen',
          'CarbonMonoxide','Water','Helium','Argon','HydrogenSulfide','n-Nonane','n-Decane']:
    print(n, CP.get_fluid_param_string(n,'CAS'))
"
```

Fix any mismatch in the table rather than in the test.

- [ ] **Step 4: Implement resolve_component**

Add to `src/Backends/GERG/GERGBackend.cpp`:

```cpp
#include "GERGData.h"
#include "CoolProp/CoolProp.h"

namespace CoolProp {
namespace GERG {

std::string resolve_component(GERGModel model, const std::string& user_name) {
    // Resolve through CoolProp's normal alias/CAS machinery first, so users
    // can spell components the way they do everywhere else in CoolProp.
    std::string cas;
    try {
        cas = get_fluid_param_string(user_name, "CAS");
    } catch (const std::exception&) {
        throw ValueError(format("[%s] is not a fluid CoolProp recognises, so it cannot be a GERG component", user_name.c_str()));
    }
    const auto& table = detail::cas_to_gerg();
    auto it = table.find(cas);
    if (it == table.end()) {
        throw ValueError(format("[%s] (CAS %s) is not a component of this GERG model", user_name.c_str(), cas.c_str()));
    }
    const std::string& gerg_name = it->second;
    const auto& names = component_names(model);
    if (std::find(names.begin(), names.end(), gerg_name) == names.end()) {
        throw ValueError(format("[%s] is a GERG-2008 component but not a GERG-2004 component", user_name.c_str()));
    }
    return gerg_name;
}

}  // namespace GERG
}  // namespace CoolProp
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[GERG]"
```
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/Backends/GERG src/Tests/CoolProp-Tests-GERG.cpp
git commit --no-verify -m "feat(gerg): component tables and CoolProp-name resolution (CoolProp-p8ub)"
```

---

## Task 3: Pure-fluid residual coefficient tables

**Files:**
- Modify: `src/Backends/GERG/GERGData.h`
- Test: `src/Tests/CoolProp-Tests-GERG.cpp`

**Interfaces:**
- Produces:
  ```cpp
  struct PureCoeffs { std::vector<double> n, t, d, c, l; };
  PureCoeffs get_pure_coeffs(GERGModel model, const std::string& gerg_name);
  ```

The GERG pure residual is `alpha^r = sum_i n_i delta^{d_i} tau^{t_i} exp(-c_i delta^{l_i})`, with `c_i = 1` exactly when `l_i > 0`. This is the form `ResidualHelmholtzGeneralizedExponential::add_Power` already implements (`include/CoolProp/fluids/Helmholtz.h:377`), which sets `c = 1.0` when `l > 0` — so `c` need not be passed through, but keep it in the struct for one-to-one comparison against teqp.

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("GERG pure residual coefficients are internally consistent", "[GERG]") {
    for (auto model : {GERGModel::GERG_2004, GERGModel::GERG_2008}) {
        for (const auto& name : component_names(model)) {
            CAPTURE(name);
            auto pc = get_pure_coeffs(model, name);
            REQUIRE(pc.n.size() > 0);
            CHECK(pc.t.size() == pc.n.size());
            CHECK(pc.d.size() == pc.n.size());
            CHECK(pc.c.size() == pc.n.size());
            CHECK(pc.l.size() == pc.n.size());
            // c_i is 1 exactly when l_i > 0
            for (std::size_t i = 0; i < pc.n.size(); ++i) {
                CHECK(pc.c[i] == ((pc.l[i] > 0) ? 1.0 : 0.0));
            }
        }
    }
}

TEST_CASE("GERG pure residual coefficients match the monograph", "[GERG]") {
    // Methane, Table A3.2 (24 terms), first and last n.
    auto ch4 = get_pure_coeffs(GERGModel::GERG_2004, "methane");
    REQUIRE(ch4.n.size() == 24);
    CHECK_THAT(ch4.n[0], Catch::Matchers::WithinRel(0.57335704239162, 1e-14));

    // The generalised 12-term set shared by most fluids.
    auto c3h8 = get_pure_coeffs(GERGModel::GERG_2004, "propane");
    REQUIRE(c3h8.n.size() == 12);
    CHECK(c3h8.t == std::vector<double>{0.250, 1.125, 1.500, 1.375, 0.250, 0.875, 0.625, 1.750, 3.625, 3.625, 14.500, 12.000});
    CHECK(c3h8.d == std::vector<double>{1, 1, 1, 2, 3, 7, 2, 5, 1, 4, 3, 4});
    CHECK(c3h8.l == std::vector<double>{0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 3, 3});

    // GERG-2008 changed carbon monoxide and isopentane.
    CHECK(get_pure_coeffs(GERGModel::GERG_2004, "carbonmonoxide").n
          != get_pure_coeffs(GERGModel::GERG_2008, "carbonmonoxide").n);
    CHECK_THAT(get_pure_coeffs(GERGModel::GERG_2008, "carbonmonoxide").n[0], Catch::Matchers::WithinRel(0.90554, 1e-14));

    // Unchanged fluids fall through identically.
    CHECK(get_pure_coeffs(GERGModel::GERG_2004, "methane").n == get_pure_coeffs(GERGModel::GERG_2008, "methane").n);
}
```

Before writing the two spot values above, read them out of teqp so the test asserts the real numbers rather than transcription noise:

```bash
sed -n '511,591p' ~/Code/teqp/include/teqp/models/GERG/GERG.hpp
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cmake --build build_catch --target CatchTestRunner -j8
```
Expected: compile failure — `get_pure_coeffs` undeclared.

- [ ] **Step 3: Transcribe the tables**

Source: `~/Code/teqp/include/teqp/models/GERG/GERG.hpp` lines 511-591 (`GERG2004::get_pure_coeffs`) and 1105-1131 (`GERG2008::get_pure_coeffs`).

Note the structure teqp uses and preserve it: most fluids share one generalized
12-term `t/d/c/l` exponent set and differ only in `n`, while a minority carry
their own longer sets. **Read teqp to determine which fluids fall in which
group — do not trust a list written here.** (An earlier draft of this plan
asserted that oxygen, carbon monoxide, and argon have their own sets; they do
not, they share the 12-term set. teqp is the authority, not this document.)
Keep the same split so the two files can be diffed by eye.

Add to `GERGData.h`:

```cpp
struct PureCoeffs
{
    std::vector<double> n, t, d, c, l;
};

PureCoeffs get_pure_coeffs(GERGModel model, const std::string& gerg_name);
```

with the definition in `GERGBackend.cpp` (the tables are large; keeping them out of line keeps the header compile time sane).

- [ ] **Step 4: Verify the transcription mechanically**

Do not eyeball 23 tables. Dump both sides and diff:

```bash
mkdir -p /tmp/gergcheck
python3 - <<'PY' > /tmp/gergcheck/teqp_pure.txt
import re, pathlib
src = pathlib.Path.home()/'Code/teqp/include/teqp/models/GERG/GERG.hpp'
text = src.read_text()
for m in re.finditer(r'\{"([a-z\-]+)",\s*\{([-0-9.e+, ]+)\}\}', text):
    nums = [float(x) for x in m.group(2).split(',') if x.strip()]
    print(m.group(1), len(nums), ' '.join('%.14g' % v for v in nums))
PY
wc -l /tmp/gergcheck/teqp_pure.txt
```

Then add a temporary `[.gergdump]` (hidden-tag) Catch2 test that prints the same rows from `get_pure_coeffs` and `diff` the two files. Delete the dump test before committing; the permanent guard is the consistency test from Step 1.

- [ ] **Step 5: Run tests to verify they pass**

```bash
cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[GERG]"
```
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/Backends/GERG src/Tests/CoolProp-Tests-GERG.cpp
git commit --no-verify -m "feat(gerg): pure-fluid residual coefficient tables (CoolProp-p8ub)"
```

---

## Task 4: Ideal-gas coefficients and integration-constant recomputation

The subtlest task in the plan. Get this wrong and `p`, `c_v`, and `w` still match teqp while `h` and `s` silently do not.

**Files:**
- Modify: `src/Backends/GERG/GERGData.h`, `src/Backends/GERG/GERGBackend.cpp`
- Test: `src/Tests/CoolProp-Tests-GERG.cpp`

**Interfaces:**
- Produces:
  ```cpp
  struct AlphaigCoeffs { std::vector<double> n0, theta0; };  // both size 8, indices 1..7 used
  AlphaigCoeffs get_alphaig_coeffs(GERGModel model, const std::string& gerg_name);
  /// Returns {n0_1, n0_2} solving h0 = s0 = 0 for the ideal gas at (T0, p0).
  std::pair<double, double> recalc_integration_constants(const AlphaigCoeffs& c, double T0, double Tci, double rho0, double rhoci,
                                                         double Rstar_R);
  constexpr double R_GERG = 8.314472;   // J/mol/K
  constexpr double RSTAR_GERG = 8.314510;  // J/mol/K
  ```

`get_alphaig_coeffs` returns coefficients **with `n0[1]` and `n0[2]` already replaced** by the recomputed values, exactly as teqp does in `GERG200XAlphaig::get_coeffs` (GERG.hpp:371-383). The published values are discarded.

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("GERG ideal-gas integration constants zero h and s at the reference state", "[GERG]") {
    // teqp GERG.hpp:376-379 — h = s = 0 for the IDEAL GAS at 298.15 K,
    // 101325 Pa, with rho0 = p0/(R*T0).
    const double T0 = 298.15, p0 = 101325.0;
    const double rho0 = p0 / (R_GERG * T0);
    const double Rstar_R = RSTAR_GERG / R_GERG;

    for (auto model : {GERGModel::GERG_2004, GERGModel::GERG_2008}) {
        for (const auto& name : component_names(model)) {
            CAPTURE(name);
            auto info = get_pure_info(model, name);
            auto c = get_alphaig_coeffs(model, name);

            // alphaig and its tau-derivative at the reference state, per
            // teqp GERG200XAlphaig::alphaig_pure (GERG.hpp:395-410).
            auto alphaig00 = [&](double T, double rho) {
                double s = c.n0[1] + c.n0[2] * info.Tc_K / T + c.n0[3] * std::log(info.Tc_K / T);
                if (c.theta0[4] != 0) s += c.n0[4] * std::log(std::abs(std::sinh(c.theta0[4] * info.Tc_K / T)));
                if (c.theta0[6] != 0) s += c.n0[6] * std::log(std::abs(std::sinh(c.theta0[6] * info.Tc_K / T)));
                if (c.theta0[5] != 0) s -= c.n0[5] * std::log(std::abs(std::cosh(c.theta0[5] * info.Tc_K / T)));
                if (c.theta0[7] != 0) s -= c.n0[7] * std::log(std::abs(std::cosh(c.theta0[7] * info.Tc_K / T)));
                return std::log(rho / info.rhoc_molm3) + Rstar_R * s;
            };
            // Aig10 = tau * d(alphaig)/d(tau); evaluate by central difference in T.
            const double dT = 1e-5 * T0;
            double dalpha_dT = (alphaig00(T0 + dT, rho0) - alphaig00(T0 - dT, rho0)) / (2 * dT);
            double Aig10 = -T0 * dalpha_dT;

            double h_over_RT = 1 + Aig10;                  // ideal gas: h/(RT) = 1 + Aig10
            double s_over_R = Aig10 - alphaig00(T0, rho0);  // s/R = Aig10 - Aig00

            CHECK_THAT(h_over_RT, Catch::Matchers::WithinAbs(0.0, 1e-8));
            CHECK_THAT(s_over_R, Catch::Matchers::WithinAbs(0.0, 1e-8));
        }
    }
}

TEST_CASE("GERG ideal-gas theta values match the monograph", "[GERG]") {
    // Table A3.1: methane theta values. n0[1], n0[2] are recomputed, so only
    // n0[3..7] and theta0[4..7] are comparable with the published table.
    auto c = get_alphaig_coeffs(GERGModel::GERG_2004, "methane");
    REQUIRE(c.n0.size() == 8);
    REQUIRE(c.theta0.size() == 8);
    // NOTE: an earlier draft of this plan had 3.00160 and 4.30632556 here.
    // Both were wrong. These are teqp's actual values — read them from teqp,
    // do not trust this document.
    CHECK_THAT(c.n0[3], Catch::Matchers::WithinRel(3.000880, 1e-12));
    CHECK_THAT(c.theta0[4], Catch::Matchers::WithinRel(4.306474465, 1e-12));
}
```

Read the true methane values out of teqp before committing them into the test:

```bash
sed -n '468,509p' ~/Code/teqp/include/teqp/models/GERG/GERG.hpp
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cmake --build build_catch --target CatchTestRunner -j8
```
Expected: compile failure — `get_alphaig_coeffs` undeclared.

- [ ] **Step 3: Transcribe the ideal-gas tables**

Source: teqp GERG.hpp lines 468-509 (`GERG2004::get_alphaig_coeffs`) and 1133-1163 (`GERG2008::get_alphaig_coeffs`). teqp stores `{n0, theta0}` as a pair of vectors, both length 8, with index 0 unused.

- [ ] **Step 4: Port recalc_integration_constants**

Direct port of teqp GERG.hpp:47-75. teqp solves a 2x2 system with Eigen's `colPivHouseholderQr`. CoolProp has Eigen available, but a 2x2 system does not need it — solve by hand and keep the code readable:

```cpp
std::pair<double, double> recalc_integration_constants(const AlphaigCoeffs& c, double T0, double Tci, double rho0, double rhoci,
                                                       double Rstar_R) {
    const double th = Tci / T0;
    auto sinh_term = [&](std::size_t i) { return (c.n0[i] != 0) ? c.n0[i] * std::log(std::abs(std::sinh(c.theta0[i] * th))) : 0.0; };
    auto cosh_term = [&](std::size_t i) { return (c.n0[i] != 0) ? c.n0[i] * std::log(std::abs(std::cosh(c.theta0[i] * th))) : 0.0; };
    auto sinh_dterm = [&](std::size_t i) {
        return (c.n0[i] != 0) ? c.n0[i] * c.theta0[i] * th / std::tanh(c.theta0[i] * th) : 0.0;
    };

    // Coefficients of {n0_1, n0_2, constant} for Aig00 and Aig10.
    const double a00_0 = Rstar_R;
    const double a00_1 = Rstar_R * th;
    const double a00_2 = std::log(rho0 / rhoci) + Rstar_R * (c.n0[3] * std::log(th) + sinh_term(4) + sinh_term(6) - cosh_term(5) - cosh_term(7));

    const double a10_0 = 0.0;
    const double a10_1 = Rstar_R * th;
    const double a10_2 = Rstar_R * (c.n0[3] + sinh_dterm(4) + sinh_dterm(6) - c.n0[5] * c.theta0[5] * th * std::tanh(c.theta0[5] * th)
                                    - c.n0[7] * c.theta0[7] * th * std::tanh(c.theta0[7] * th));

    // Row 0: h0 = 0  =>  Aig10 = -1
    // Row 1: s0 = 0  =>  Aig10 - Aig00 = 0
    const double A00 = a10_0, A01 = a10_1, b0 = -1.0 - a10_2;
    const double A10 = a10_0 - a00_0, A11 = a10_1 - a00_1, b1 = -a10_2 + a00_2;

    const double det = A00 * A11 - A01 * A10;
    if (std::abs(det) < 1e-300) {
        throw ValueError("GERG ideal-gas integration constants: singular 2x2 system");
    }
    const double n1 = (b0 * A11 - A01 * b1) / det;
    const double n2 = (A00 * b1 - b0 * A10) / det;
    return {n1, n2};
}
```

Note that teqp's `Aig10` expression applies the `!= 0` guard to the sinh terms but **not** to the two cosh terms (GERG.hpp:63-65). Preserve that asymmetry exactly — where `n0[5]` or `n0[7]` is zero the term contributes zero anyway, so the behaviour is identical, but keeping the shape identical makes the two implementations diffable.

`get_alphaig_coeffs` calls this and overwrites `n0[1]`, `n0[2]` before returning, using `T0 = 298.15`, `p0 = 101325`, `rho0 = p0/(R_GERG*T0)`, `Tci`/`rhoci` from `get_pure_info`, and `Rstar_R = RSTAR_GERG/R_GERG`.

- [ ] **Step 5: Run tests to verify they pass**

```bash
cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[GERG]"
```
Expected: PASS. If `h_over_RT` is nonzero for every fluid by the same factor, the `R*/R` ratio is misplaced. If it is wrong for a subset, the transcription of those fluids' `theta0` is wrong.

- [ ] **Step 6: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/Backends/GERG src/Tests/CoolProp-Tests-GERG.cpp
git commit --no-verify -m "feat(gerg): ideal-gas coefficients and integration-constant recomputation (CoolProp-p8ub)"
```

---

## Task 5: Reducing-parameter tables

The largest transcription in the plan: 153 pairs for GERG-2004 plus GERG-2008's changes.

**Files:**
- Modify: `src/Backends/GERG/GERGData.h`, `src/Backends/GERG/GERGBackend.cpp`
- Test: `src/Tests/CoolProp-Tests-GERG.cpp`

**Interfaces:**
- Produces:
  ```cpp
  struct BetasGammas { double betaV, gammaV, betaT, gammaT; };
  BetasGammas get_betasgammas(GERGModel model, const std::string& f1, const std::string& f2);
  ```
  Field order matches teqp's struct (`betaV, gammaV, betaT, gammaT`) so brace-initialised rows can be copied verbatim. When the stored pair is `(f2, f1)`, the accessor returns `betaV` and `betaT` **reciprocated**, matching teqp GERG.hpp:757-763.

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("GERG reducing parameters exist for every binary pair", "[GERG]") {
    for (auto model : {GERGModel::GERG_2004, GERGModel::GERG_2008}) {
        const auto& names = component_names(model);
        for (std::size_t i = 0; i < names.size(); ++i) {
            for (std::size_t j = i + 1; j < names.size(); ++j) {
                CAPTURE(names[i], names[j]);
                CHECK_NOTHROW(get_betasgammas(model, names[i], names[j]));
                CHECK_NOTHROW(get_betasgammas(model, names[j], names[i]));
            }
        }
    }
}

TEST_CASE("GERG reducing parameters invert correctly when the pair is reversed", "[GERG]") {
    auto fwd = get_betasgammas(GERGModel::GERG_2008, "methane", "nitrogen");
    auto rev = get_betasgammas(GERGModel::GERG_2008, "nitrogen", "methane");
    CHECK_THAT(rev.betaT, Catch::Matchers::WithinRel(1.0 / fwd.betaT, 1e-14));
    CHECK_THAT(rev.betaV, Catch::Matchers::WithinRel(1.0 / fwd.betaV, 1e-14));
    // Gammas are symmetric, not reciprocal.
    CHECK_THAT(rev.gammaT, Catch::Matchers::WithinRel(fwd.gammaT, 1e-14));
    CHECK_THAT(rev.gammaV, Catch::Matchers::WithinRel(fwd.gammaV, 1e-14));
}

TEST_CASE("GERG-2008 changes some reducing parameters relative to GERG-2004", "[GERG]") {
    // Table A8 of the GERG-2008 manuscript revises a subset of pairs.
    // At least one pair must differ, and pairs GERG-2008 did not touch
    // must fall through unchanged.
    bool any_different = false;
    const auto& names = component_names(GERGModel::GERG_2004);
    for (std::size_t i = 0; i < names.size(); ++i) {
        for (std::size_t j = i + 1; j < names.size(); ++j) {
            auto a = get_betasgammas(GERGModel::GERG_2004, names[i], names[j]);
            auto b = get_betasgammas(GERGModel::GERG_2008, names[i], names[j]);
            if (a.betaT != b.betaT || a.gammaT != b.gammaT || a.betaV != b.betaV || a.gammaV != b.gammaV) {
                any_different = true;
            }
        }
    }
    CHECK(any_different);
}

TEST_CASE("GERG reducing parameter lookup rejects unknown fluids", "[GERG]") {
    CHECK_THROWS_AS(get_betasgammas(GERGModel::GERG_2004, "NOT A FLUID", "water"), CoolProp::ValueError);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cmake --build build_catch --target CatchTestRunner -j8
```
Expected: compile failure — `get_betasgammas` undeclared.

- [ ] **Step 3: Transcribe the tables**

Source: teqp GERG.hpp lines 595-766 (`GERG2004::get_betasgammas`) and 1003-1103 (`GERG2008::get_betasgammas`). teqp keys on `std::pair<std::string,std::string>` with `boost::hash`; use `std::map<std::pair<std::string,std::string>, BetasGammas>` instead — `std::map` orders pairs natively and CoolProp should not take a Boost dependency for this.

The GERG-2008 table holds only pairs whose values changed; unlisted pairs fall through to the GERG-2004 table. Preserve teqp's comment structure marking which block is which.

- [ ] **Step 4: Verify the transcription mechanically**

Extract both sides and diff, as in Task 3 Step 4:

```bash
python3 - <<'PY' > /tmp/gergcheck/teqp_bip.txt
import re, pathlib
src = pathlib.Path.home()/'Code/teqp/include/teqp/models/GERG/GERG.hpp'
text = src.read_text()
# rows look like: {{"methane","nitrogen"}, {1.0,1.0,1.0,1.0}},
for m in re.finditer(r'\{\{"([a-z\-]+)","([a-z\-]+)"\},\s*\{([-0-9.e+, ]+)\}\}', text):
    nums = [float(x) for x in m.group(3).split(',') if x.strip()]
    print(m.group(1), m.group(2), ' '.join('%.14g' % v for v in nums))
PY
sort /tmp/gergcheck/teqp_bip.txt | uniq | wc -l
```

The extraction picks up both the GERG-2004 and GERG-2008 blocks in file order; split them at the `namespace GERG2008` line. Emit the same rows from `get_betasgammas` in a temporary `[.gergdump]` test and diff. Delete the dump test before committing.

- [ ] **Step 5: Run tests to verify they pass**

```bash
cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[GERG]"
```
Expected: PASS. The "exists for every binary pair" test covers 153 + 210 lookups in both orders; a single missing row fails it.

- [ ] **Step 6: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/Backends/GERG src/Tests/CoolProp-Tests-GERG.cpp
git commit --no-verify -m "feat(gerg): binary reducing-parameter tables (CoolProp-p8ub)"
```

---

## Task 6: Departure function tables

**Files:**
- Modify: `src/Backends/GERG/GERGData.h`, `src/Backends/GERG/GERGBackend.cpp`
- Test: `src/Tests/CoolProp-Tests-GERG.cpp`

**Interfaces:**
- Produces:
  ```cpp
  struct DepartureCoeffs { std::vector<double> n, t, d, eta, beta, gamma, epsilon; };
  /// Returns false (and leaves F untouched) when the pair has no departure function.
  bool get_Fij(GERGModel model, const std::string& f1, const std::string& f2, double& F);
  DepartureCoeffs get_departurecoeffs(GERGModel model, const std::string& f1, const std::string& f2);
  /// Number of leading terms with eta == 0 and beta == 0 — the polynomial block.
  std::size_t departure_Npower(const DepartureCoeffs& dc);
  ```
  `get_Fij` returns a bool rather than `std::optional<double>` to keep the header free of `<optional>` churn; the semantics match teqp's `ok_missing = true` path.

GERG-2008 reuses GERG-2004's `F_ij` and departure coefficients unchanged (teqp GERG.hpp:999-1001, `using GERG2004::get_Fij; using GERG2004::get_departurecoeffs;`), so both models share one table.

`departure_Npower` exists because CoolProp's `GERG2008DepartureFunction` constructor (`src/Backends/Helmholtz/ExcessHEFunction.h:108`) takes an `Npower` argument splitting the polynomial block from the Gaussian block. teqp stores a single uniform list where polynomial terms simply carry `eta = beta = 0`.

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("GERG has departure functions for exactly the published pairs", "[GERG]") {
    const auto& names = component_names(GERGModel::GERG_2008);
    std::size_t with_departure = 0, with_scaled_F = 0;
    for (std::size_t i = 0; i < names.size(); ++i) {
        for (std::size_t j = i + 1; j < names.size(); ++j) {
            double F = 0;
            if (get_Fij(GERGModel::GERG_2008, names[i], names[j], F)) {
                ++with_departure;
                if (F != 1.0) ++with_scaled_F;
                CAPTURE(names[i], names[j]);
                CHECK_NOTHROW(get_departurecoeffs(GERGModel::GERG_2008, names[i], names[j]));
            } else {
                CHECK_THROWS_AS(get_departurecoeffs(GERGModel::GERG_2008, names[i], names[j]), CoolProp::ValueError);
            }
        }
    }
    // Table A3.6: 15 pairs carry a departure function; 7 of them use the
    // generalised function scaled by F_ij != 1.
    CHECK(with_departure == 15);
    CHECK(with_scaled_F == 7);
}

TEST_CASE("GERG F_ij is symmetric", "[GERG]") {
    double a = 0, b = 0;
    REQUIRE(get_Fij(GERGModel::GERG_2008, "methane", "isobutane", a));
    REQUIRE(get_Fij(GERGModel::GERG_2008, "isobutane", "methane", b));
    CHECK_THAT(a, Catch::Matchers::WithinRel(b, 1e-14));
    CHECK_THAT(a, Catch::Matchers::WithinRel(0.771035405688, 1e-12));
}

TEST_CASE("GERG departure coefficients are internally consistent", "[GERG]") {
    auto dc = get_departurecoeffs(GERGModel::GERG_2008, "methane", "nitrogen");
    REQUIRE(dc.n.size() > 0);
    CHECK(dc.t.size() == dc.n.size());
    CHECK(dc.d.size() == dc.n.size());
    CHECK(dc.eta.size() == dc.n.size());
    CHECK(dc.beta.size() == dc.n.size());
    CHECK(dc.gamma.size() == dc.n.size());
    CHECK(dc.epsilon.size() == dc.n.size());

    // The polynomial block comes first and is contiguous.
    std::size_t np = departure_Npower(dc);
    for (std::size_t k = 0; k < np; ++k) {
        CHECK(dc.eta[k] == 0.0);
        CHECK(dc.beta[k] == 0.0);
    }
    for (std::size_t k = np; k < dc.n.size(); ++k) {
        CHECK((dc.eta[k] != 0.0 || dc.beta[k] != 0.0));
    }
}

TEST_CASE("GERG-2004 and GERG-2008 share departure data", "[GERG]") {
    CHECK(get_departurecoeffs(GERGModel::GERG_2004, "methane", "nitrogen").n
          == get_departurecoeffs(GERGModel::GERG_2008, "methane", "nitrogen").n);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cmake --build build_catch --target CatchTestRunner -j8
```
Expected: compile failure — `get_Fij` undeclared.

- [ ] **Step 3: Transcribe the tables**

Source: teqp GERG.hpp lines 768-802 (`get_Fij`, Table A3.6) and 805-900 (`get_departurecoeffs`).

- [ ] **Step 4: Implement departure_Npower**

```cpp
std::size_t departure_Npower(const DepartureCoeffs& dc) {
    std::size_t np = 0;
    while (np < dc.n.size() && dc.eta[np] == 0.0 && dc.beta[np] == 0.0) {
        ++np;
    }
    // CoolProp's GERG2008DepartureFunction assumes the polynomial block is a
    // contiguous prefix.  If any later term is also polynomial, that
    // assumption is violated and the split would silently drop terms.
    for (std::size_t k = np; k < dc.n.size(); ++k) {
        if (dc.eta[k] == 0.0 && dc.beta[k] == 0.0) {
            throw ValueError("GERG departure coefficients: polynomial terms are not a contiguous prefix");
        }
    }
    return np;
}
```

The throw is the point. A silent miscount here would place Gaussian terms into the polynomial block, producing wrong numbers with no error.

- [ ] **Step 5: Run tests to verify they pass**

```bash
cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[GERG]"
```
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/Backends/GERG src/Tests/CoolProp-Tests-GERG.cpp
git commit --no-verify -m "feat(gerg): departure function tables and F_ij (CoolProp-p8ub)"
```

---

## Task 7: Reference values generated from teqp

Produces the fixtures Tasks 8 and 9 validate against. Nothing here touches CoolProp runtime code.

**Files:**
- Create: `dev/gerg/generate_reference_values.py`
- Create: `dev/gerg/README.md`
- Create: `src/Backends/GERG/GERGReferenceValues.h` (generated, committed)

**Interfaces:**
- Produces, in `namespace CoolProp::GERG::reference`:
  ```cpp
  struct PureRefPoint { const char* name; double T_K, rhomolar, alphar, alphaig, p_Pa, cvmolar, w; };
  struct MixRefPoint { std::vector<const char*> names; std::vector<double> z; double T_K, rhomolar, p_Pa, cvmolar, w; };
  extern const std::vector<PureRefPoint> pure_points_2004, pure_points_2008;
  extern const std::vector<MixRefPoint> mix_points_2004, mix_points_2008;
  ```

- [ ] **Step 1: Confirm the teqp version exposes GERG**

The system Python has teqp 0.15.3, which predates the GERG factory registration (added before tag `v0.18.0`). Set up an isolated environment:

```bash
python3 -m venv /tmp/gergenv
/tmp/gergenv/bin/pip install "teqp>=0.18" numpy
/tmp/gergenv/bin/python -c "
import teqp
m = teqp.make_model({'kind':'GERG2008resid','model':{'names':['methane','nitrogen']}})
print('ok', m.get_Ar00(300.0, 1000.0, __import__('numpy').array([0.9,0.1])))
"
```
Expected: prints `ok <number>`. If `Unknown kind:GERG2008resid`, raise the version bound until it works, and record the working version in `dev/gerg/README.md`.

- [ ] **Step 2: Write the generator**

`dev/gerg/generate_reference_values.py`:

```python
"""Generate GERG reference values from teqp for the CoolProp GERG backend tests.

Run:
    python3 -m venv /tmp/gergenv
    /tmp/gergenv/bin/pip install "teqp>=0.18" numpy
    /tmp/gergenv/bin/python dev/gerg/generate_reference_values.py \
        > src/Backends/GERG/GERGReferenceValues.h

teqp is the reference implementation this backend is validated against; see
docs/superpowers/specs/2026-07-25-gerg-strict-backend-design.md.
"""
import numpy as np
import teqp

R = 8.314472  # J/mol/K, the GERG gas constant

GERG2004_NAMES = ["methane", "nitrogen", "carbondioxide", "ethane", "propane", "n-butane",
                  "isobutane", "n-pentane", "isopentane", "n-hexane", "n-heptane", "n-octane",
                  "hydrogen", "oxygen", "carbonmonoxide", "water", "helium", "argon"]
GERG2008_NAMES = GERG2004_NAMES + ["hydrogensulfide", "n-nonane", "n-decane"]


def models(year, names):
    resid = teqp.make_model({"kind": f"GERG{year}resid", "model": {"names": names}})
    ideal = teqp.make_model({"kind": f"GERG{year}idealgas", "model": {"names": names}})
    return resid, ideal


def point(year, names, z, T, rho):
    """Return (alphar, alphaig, p, cv, w) at (T, rho) for composition z."""
    resid, ideal = models(year, names)
    z = np.asarray(z, dtype=float)

    # teqp's accessor naming is get_Ar<TAU order><DELTA order>.  An earlier
    # draft of this plan had these indices BACKWARDS, which produced values
    # that looked plausible and were wrong.  Verified empirically:
    # get_Ar20 on the ideal-gas model gives -cv_ideal/R (3.303 for methane
    # at 300 K, which is right), while get_Ar02 gives -1.0.
    Ar00 = resid.get_Ar00(T, rho, z)   # alpha^r
    Ar01 = resid.get_Ar01(T, rho, z)   # delta * dalpha^r/ddelta
    Ar02 = resid.get_Ar02(T, rho, z)   # delta^2 * d2alpha^r/ddelta^2
    Ar11 = resid.get_Ar11(T, rho, z)   # delta*tau * d2alpha^r/(ddelta dtau)
    Ar20 = resid.get_Ar20(T, rho, z)   # tau^2 * d2alpha^r/dtau^2
    Aig00 = ideal.get_Ar00(T, rho, z)
    Aig20 = ideal.get_Ar20(T, rho, z)

    p = rho * R * T * (1.0 + Ar01)
    cv = -(Aig20 + Ar20) * R
    Mmix = float(np.dot(z, [molar_mass(year, n) for n in names]))
    w2 = (R * T / Mmix) * (1 + 2 * Ar01 + Ar02 - (1 + Ar01 - Ar11) ** 2 / (Aig20 + Ar20))
    return Ar00, Aig00, p, cv, np.sqrt(w2)
```

These three expressions are verified: against teqp's published validation
point for gas 2 (T = 190.68 K, D = 11.0 mol/L) they reproduce
p = 4.62270367011 MPa to 4e-13 relative.

The `molar_mass` helper and the `Ar20` / `Ar11` accessor names must be checked against the installed teqp — the derivative-accessor naming convention is `get_ArXY` where X is the tau order and Y the delta order, and the speed-of-sound grouping above follows teqp's own test (`~/Code/teqp/src/tests/catch_test_GERG.cxx:675-694`). Read that test and mirror its exact expressions rather than re-deriving them; any disagreement between this script and teqp's test is a bug in this script.

For molar masses, use `get_pure_info` values already transcribed in Task 2 rather than a second source of truth — hardcode them in the script from the same monograph table and add an assertion comparing against `teqp`'s own `MolarMassGERG` if exposed.

- [ ] **Step 3: Choose the sample points**

Pure fluids: for each of the 23 distinct pure EOS, emit points at
`T = {0.7, 1.0, 1.3, 2.0} x Tc` crossed with `rho = {0.01, 0.5, 1.5, 2.5} x rhoc`,
skipping any point where teqp returns a non-finite value. That is up to 16 points per fluid.

Mixtures: transcribe the 21-component gas compositions and state points from teqp's
`mixture_comps` (`~/Code/teqp/src/tests/catch_test_GERG.cxx:119-321`) and
`validation_data` (`:322-...`). Note that `validation_data[i].GasNo - 2` indexes
`mixture_comps`, compositions are in **mole percent** and must be divided by 100,
and `D_molL` is in mol/L so multiply by 1000 for mol/m^3.

**CRITICAL — the columns of `mixture_comps` are NOT in `GERG2008::component_names`
order.** They follow the AGA8 ordering, declared separately as `components` at
`catch_test_GERG.cxx:512`:

```
{"methane","nitrogen","carbondioxide","ethane","propane","isobutane","n-butane",
 "isopentane","n-pentane","n-hexane","n-heptane","n-octane","n-nonane","n-decane",
 "hydrogen","oxygen","carbonmonoxide","water","hydrogensulfide","helium","argon"}
```

Note isobutane BEFORE n-butane, isopentane BEFORE n-pentane, and helium/argon
last — all three differ from `component_names`. Build the teqp model with THIS
vector when evaluating `mixture_comps` rows. Using `component_names` order
instead silently misassigns every composition: it still sums to 1, still
evaluates, and gives answers wrong by a few tenths of a percent. Verified: with
the wrong order gas 2 gives p = 4.6142620561 MPa; with this order it gives
4.62270367011 MPa, matching the published value to 4e-13.

Also note the precision asymmetry in `validation_data`: its `P_MPa` column is
full precision, but `cv_JmolK` and `w_ms` agree with teqp only to about 1e-6
because they were printed from the published AGA8 table. **Generate the
reference values from teqp directly at full double precision**; use
`validation_data` for the state points (GasNo, T, D) and treat its P column as
a cross-check, not as the fixture.

Also emit every binary pair at a single mid-range state
(`T = 250 K`, `rho = 5000 mol/m^3`, `z = [0.4, 0.6]`) so Task 9 exercises all
153 / 210 reducing-parameter rows numerically, not just structurally.

- [ ] **Step 4: Generate and inspect the header**

```bash
/tmp/gergenv/bin/python dev/gerg/generate_reference_values.py > src/Backends/GERG/GERGReferenceValues.h
grep -c "^  {" src/Backends/GERG/GERGReferenceValues.h
```
Expected: several hundred rows. Emit all floating-point values with `%.17g` so the fixtures round-trip exactly.

- [ ] **Step 5: Write dev/gerg/README.md**

Document: which teqp version was used, the exact venv commands, why the file is committed rather than generated at build time (no teqp dependency in CoolProp's build), and that regenerating it is expected to produce a byte-identical file unless teqp's GERG implementation changes.

- [ ] **Step 6: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add dev/gerg src/Backends/GERG/GERGReferenceValues.h
git commit --no-verify -m "test(gerg): teqp-generated reference values (CoolProp-p8ub)"
```

---

## Task 8: Pure-fluid backend construction

**Files:**
- Modify: `src/Backends/GERG/GERGBackend.h`, `src/Backends/GERG/GERGBackend.cpp`
- Test: `src/Tests/CoolProp-Tests-GERG.cpp`

**Interfaces:**
- Consumes: everything from Tasks 2-4, `GERGReferenceValues.h` from Task 7.
- Produces: a working `GERGMixtureBackend` for `N == 1`; helper
  ```cpp
  CoolPropFluid make_gerg_fluid(GERGModel model, const std::string& gerg_name);
  ```

- [ ] **Step 1: Delete the Task 1 skeleton test, then write the failing test**

First delete `TEST_CASE("GERG factory reaches the GERG backend", "[GERG]")`
from `src/Tests/CoolProp-Tests-GERG.cpp`. It asserts the constructor throws
`NotImplementedError`, which this task makes false. Leaving it in place means
this task's own build fails on a stale assertion. The
`"GERG backend families are registered"` test stays — it is still true.

```cpp
#include "src/Backends/GERG/GERGReferenceValues.h"

TEST_CASE("GERG pure fluids reproduce teqp", "[GERG]") {
    for (const auto& pt : CoolProp::GERG::reference::pure_points_2008) {
        CAPTURE(pt.name, pt.T_K, pt.rhomolar);
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{pt.name}));
        AS->update(DmolarT_INPUTS, pt.rhomolar, pt.T_K);
        CHECK_THAT(AS->p(), Catch::Matchers::WithinRel(pt.p_Pa, 1e-10));
        CHECK_THAT(AS->cvmolar(), Catch::Matchers::WithinRel(pt.cvmolar, 1e-10));
        CHECK_THAT(AS->speed_sound(), Catch::Matchers::WithinRel(pt.w, 1e-10));
        CHECK_THAT(AS->alphar(), Catch::Matchers::WithinRel(pt.alphar, 1e-12));
    }
}

TEST_CASE("GERG pure fluid uses the GERG gas constant and reducing state", "[GERG]") {
    std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    CHECK_THAT(AS->gas_constant(), Catch::Matchers::WithinRel(8.314472, 1e-14));
    CHECK_THAT(AS->T_reducing(), Catch::Matchers::WithinRel(190.564, 1e-12));
    CHECK_THAT(AS->rhomolar_reducing(), Catch::Matchers::WithinRel(10.139342719e3, 1e-12));
    CHECK_THAT(AS->molar_mass(), Catch::Matchers::WithinRel(16.042460e-3, 1e-12));
}

TEST_CASE("GERG pure fluid reference state gives h = s = 0 for the ideal gas", "[GERG]") {
    // Consistency check on the whole assembled ideal-gas path, not just the
    // coefficient solve tested in Task 4.
    const double T0 = 298.15, p0 = 101325.0;
    const double rho0 = p0 / (8.314472 * T0);
    for (const auto& name : CoolProp::GERG::component_names(CoolProp::GERGModel::GERG_2008)) {
        CAPTURE(name);
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{name}));
        // Evaluate at a density low enough that the residual part is negligible
        // but use the ideal-gas density directly for an exact statement.
        AS->update(DmolarT_INPUTS, rho0, T0);
        // h_ideal = h - h_residual; compare the ideal-gas contribution only.
        double h_ig_over_RT = 1 + AS->tau() * AS->dalpha0_dTau();
        double s_ig_over_R = AS->tau() * AS->dalpha0_dTau() - AS->alpha0();
        CHECK_THAT(h_ig_over_RT, Catch::Matchers::WithinAbs(0.0, 1e-8));
        CHECK_THAT(s_ig_over_R, Catch::Matchers::WithinAbs(0.0, 1e-8));
    }
}
```

If `alphar()`, `alpha0()`, `dalpha0_dTau()`, or `tau()` are not public on `AbstractState` in this tree, use `keyed_output` equivalents or drop those specific assertions — the `p`/`cvmolar`/`speed_sound` comparisons are the load-bearing ones. Check `include/CoolProp/AbstractState.h` before writing.

- [ ] **Step 2: Run test to verify it fails**

```bash
cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[GERG]"
```
Expected: FAIL — `NotImplementedError: GERG backend is not yet implemented`.

- [ ] **Step 3: Implement make_gerg_fluid**

```cpp
CoolPropFluid make_gerg_fluid(GERGModel model, const std::string& gerg_name) {
    auto info = get_pure_info(model, gerg_name);
    auto pc = get_pure_coeffs(model, gerg_name);
    auto ac = get_alphaig_coeffs(model, gerg_name);

    CoolPropFluid fluid;
    fluid.name = gerg_name;
    fluid.EOSVector.resize(1);
    EquationOfState& EOS = fluid.EOS();

    EOS.R_u = R_GERG;
    EOS.molar_mass = info.M_kgmol;
    EOS.pseudo_pure = false;
    EOS.reduce.T = info.Tc_K;
    EOS.reduce.rhomolar = info.rhoc_molm3;
    fluid.crit = EOS.reduce;

    // GERG-2008 extended range of validity: 60-700 K, p <= 70 MPa.
    EOS.limits.Tmin = 60.0;
    EOS.limits.Tmax = 700.0;
    EOS.limits.pmax = 70e6;
    EOS.limits.rhomax = 1e6;

    // Residual: n_i delta^{d_i} tau^{t_i} exp(-c_i delta^{l_i}), where
    // add_Power sets c = 1 exactly when l > 0 — matching GERG's convention.
    {
        std::vector<CoolPropDbl> n(pc.n.begin(), pc.n.end());
        std::vector<CoolPropDbl> d(pc.d.begin(), pc.d.end());
        std::vector<CoolPropDbl> t(pc.t.begin(), pc.t.end());
        std::vector<CoolPropDbl> l(pc.l.begin(), pc.l.end());
        EOS.alphar.GeneralizedExponential.add_Power(n, d, t, l);
        EOS.alphar.GeneralizedExponential.finish();
    }

    // Ideal gas.  The R*/R ratio (GERG-2008 placement) is folded into the
    // coefficients here so the stock CoolProp terms can be used unmodified.
    {
        const double RR = RSTAR_GERG / R_GERG;
        EOS.alpha0.Lead = IdealHelmholtzLead(RR * ac.n0[1], RR * ac.n0[2]);
        EOS.alpha0.LogTau = IdealHelmholtzLogTau(RR * ac.n0[3]);

        std::vector<CoolPropDbl> n_sinh, th_sinh, n_cosh, th_cosh;
        for (std::size_t k : {std::size_t(4), std::size_t(6)}) {
            if (ac.theta0[k] != 0) { n_sinh.push_back(RR * ac.n0[k]); th_sinh.push_back(ac.theta0[k]); }
        }
        for (std::size_t k : {std::size_t(5), std::size_t(7)}) {
            if (ac.theta0[k] != 0) { n_cosh.push_back(RR * ac.n0[k]); th_cosh.push_back(ac.theta0[k]); }
        }
        if (!n_sinh.empty()) EOS.alpha0.GERG2004Sinh = IdealHelmholtzGERG2004Sinh(n_sinh, th_sinh, info.Tc_K);
        if (!n_cosh.empty()) EOS.alpha0.GERG2004Cosh = IdealHelmholtzGERG2004Cosh(n_cosh, th_cosh, info.Tc_K);
    }

    EOS.validate();
    return fluid;
}
```

Two sign conventions to verify against `src/Helmholtz.cpp:1312-1341` before trusting this:

1. `IdealHelmholtzGERG2004Cosh::all` already applies the minus sign internally
   (`sum00 += -n[i]*log(|cosh(t*tau)|)`), so pass `n` **positive**, matching the
   published table and teqp's `out -= n0[5]*log(|cosh|)`.
2. The term computes `t = theta[i] * Tc/T_red` and evaluates at `t*tau`. For a
   pure fluid `T_red == Tc`, so this reduces to `theta*Tc/T` as GERG requires.
   Confirm `IdealHelmholtzGERG2004Sinh` has the matching **positive** sign.

Check the exact member name for the residual container (`EOS.alphar.GeneralizedExponential` vs. a differently-named member) in `include/CoolProp/fluids/Helmholtz.h` before compiling.

- [ ] **Step 4: Wire the constructor**

```cpp
GERGMixtureBackend::GERGMixtureBackend(GERGModel model, const std::vector<std::string>& names) : m_model(model) {
    if (names.empty()) {
        throw ValueError("GERG backend requires at least one component");
    }
    std::vector<CoolPropFluid> fluids;
    fluids.reserve(names.size());
    for (const auto& user_name : names) {
        fluids.push_back(make_gerg_fluid(model, resolve_component(model, user_name)));
    }
    set_components(fluids);
    if (names.size() == 1) {
        set_mole_fractions(std::vector<CoolPropDbl>(1, 1.0));
    }
}
```

`set_components` calls `set_mixture_parameters()` for `N > 1`; Task 9 overrides that. For `N == 1` it is not called, so pure fluids work at the end of this task.

- [ ] **Step 5: Run tests to verify they pass**

```bash
cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[GERG]"
```
Expected: PASS. If `p` matches but `cvmolar` does not, the ideal-gas second tau-derivative is wrong — check the sinh/cosh split. If `p` is wrong by a constant factor, check `R_u`.

- [ ] **Step 6: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/Backends/GERG src/Tests/CoolProp-Tests-GERG.cpp
git commit --no-verify -m "feat(gerg): pure-fluid backend matching teqp to 1e-10 (CoolProp-p8ub)"
```

---

## Task 9: Mixture wiring

**Files:**
- Modify: `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.h:423`
- Modify: `src/Backends/GERG/GERGBackend.h`, `src/Backends/GERG/GERGBackend.cpp`
- Test: `src/Tests/CoolProp-Tests-GERG.cpp`

**Interfaces:**
- Consumes: Tasks 5, 6, 7.
- Produces: `void GERGMixtureBackend::set_mixture_parameters() override;`

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("GERG binary mixtures reproduce teqp", "[GERG]") {
    for (const auto& pt : CoolProp::GERG::reference::mix_points_2008) {
        CAPTURE(pt.T_K, pt.rhomolar);
        std::vector<std::string> names(pt.names.begin(), pt.names.end());
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", names));
        AS->set_mole_fractions(std::vector<CoolPropDbl>(pt.z.begin(), pt.z.end()));
        AS->update(DmolarT_INPUTS, pt.rhomolar, pt.T_K);
        CHECK_THAT(AS->p(), Catch::Matchers::WithinRel(pt.p_Pa, 1e-10));
        CHECK_THAT(AS->cvmolar(), Catch::Matchers::WithinRel(pt.cvmolar, 1e-10));
        CHECK_THAT(AS->speed_sound(), Catch::Matchers::WithinRel(pt.w, 1e-10));
    }
}

TEST_CASE("GERG-2004 and GERG-2008 disagree where the models disagree", "[GERG]") {
    // Carbon monoxide's pure EOS and reducing state both changed.
    std::shared_ptr<AbstractState> a(AbstractState::factory("GERG2004", std::vector<std::string>{"CarbonMonoxide"}));
    std::shared_ptr<AbstractState> b(AbstractState::factory("GERG2008", std::vector<std::string>{"CarbonMonoxide"}));
    a->update(DmolarT_INPUTS, 5000.0, 200.0);
    b->update(DmolarT_INPUTS, 5000.0, 200.0);
    CHECK(a->p() != b->p());

    // Methane is unchanged between the two models.
    std::shared_ptr<AbstractState> c(AbstractState::factory("GERG2004", std::vector<std::string>{"Methane"}));
    std::shared_ptr<AbstractState> d(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    c->update(DmolarT_INPUTS, 5000.0, 200.0);
    d->update(DmolarT_INPUTS, 5000.0, 200.0);
    CHECK_THAT(c->p(), Catch::Matchers::WithinRel(d->p(), 1e-15));
}

TEST_CASE("GERG mixture state can be copied", "[GERG]") {
    // get_copy() deep-copies the ExcessTerm, which dereferences every
    // off-diagonal departure function pointer.  A null there would crash.
    std::shared_ptr<AbstractState> AS(
      AbstractState::factory("GERG2008", std::vector<std::string>{"Methane", "Water"}));  // pair with F_ij == 0
    AS->set_mole_fractions(std::vector<CoolPropDbl>{0.9, 0.1});
    CHECK_NOTHROW(AS->update(DmolarT_INPUTS, 5000.0, 400.0));
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[GERG]"
```
Expected: FAIL — `MixtureParameters::set_mixture_parameters` throws "Could not match the binary pair", because it is still reading the global BIP library.

- [ ] **Step 3: Make the base method virtual**

In `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.h` line 423:

```cpp
    virtual void set_mixture_parameters();
```

- [ ] **Step 4: Override it**

Declare in `GERGBackend.h`:

```cpp
   protected:
    void set_mixture_parameters() override;
```

and implement, modelled on `MixtureParameters::set_mixture_parameters`
(`src/Backends/Helmholtz/MixtureParameters.cpp:556-644`):

```cpp
void GERGMixtureBackend::set_mixture_parameters() {
    std::vector<CoolPropFluid> comps = get_components();
    const std::size_t N = comps.size();

    STLMatrix beta_v(N, std::vector<CoolPropDbl>(N, 0));
    STLMatrix gamma_v(N, std::vector<CoolPropDbl>(N, 0));
    STLMatrix beta_T(N, std::vector<CoolPropDbl>(N, 0));
    STLMatrix gamma_T(N, std::vector<CoolPropDbl>(N, 0));

    residual_helmholtz->Excess.resize(N);

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = i + 1; j < N; ++j) {
            // comps[k].name holds the GERG component name (set in make_gerg_fluid).
            auto bg = get_betasgammas(m_model, comps[i].name, comps[j].name);

            beta_T[i][j] = bg.betaT;   beta_T[j][i] = 1.0 / bg.betaT;
            gamma_T[i][j] = bg.gammaT; gamma_T[j][i] = bg.gammaT;
            beta_v[i][j] = bg.betaV;   beta_v[j][i] = 1.0 / bg.betaV;
            gamma_v[i][j] = bg.gammaV; gamma_v[j][i] = bg.gammaV;

            double F = 0.0;
            const bool has_departure = get_Fij(m_model, comps[i].name, comps[j].name, F);
            residual_helmholtz->Excess.F[i][j] = has_departure ? F : 0.0;
            residual_helmholtz->Excess.F[j][i] = residual_helmholtz->Excess.F[i][j];

            DepartureFunctionPointer dep;
            if (has_departure) {
                auto dc = get_departurecoeffs(m_model, comps[i].name, comps[j].name);
                dep = std::make_shared<GERG2008DepartureFunction>(dc.n, dc.d, dc.t, dc.eta, dc.epsilon, dc.beta, dc.gamma,
                                                                  departure_Npower(dc));
            } else {
                // ExcessTerm::copy() dereferences every off-diagonal pointer,
                // so a null here would segfault on get_copy().  Install the
                // same zero-valued placeholder MixtureParameters uses.
                std::vector<double> n(1, 0), d(1, 1), t(1, 1), l(1, 0);
                dep = std::make_shared<ExponentialDepartureFunction>(n, d, t, l);
            }
            residual_helmholtz->Excess.DepartureFunctionMatrix[i][j] = dep;
            residual_helmholtz->Excess.DepartureFunctionMatrix[j][i] = dep;
        }
    }

    Reducing = std::make_shared<GERG2008ReducingFunction>(comps, beta_v, gamma_v, beta_T, gamma_T);
}
```

Note the departure function is symmetric (`alphar_ij == alphar_ji`), so sharing one `shared_ptr` across both matrix slots is correct — but confirm nothing mutates a `DepartureFunction` in place per-slot; `ExcessTerm::update` calls `update(tau, delta)` on each entry, and calling it twice on the same object with identical arguments is harmless. If profiling later shows cache thrash, construct two instances.

- [ ] **Step 5: Run tests to verify they pass**

```bash
cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[GERG]"
```
Expected: PASS, including the full AGA8 21-component set.

If pressures are close but not to 1e-10, suspect the beta reciprocal convention: swap `beta_T[i][j]` and `beta_T[j][i]` for one pair and see whether the error grows or shrinks.

- [ ] **Step 6: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/Backends/GERG src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.h src/Tests/CoolProp-Tests-GERG.cpp
git commit --no-verify -m "feat(gerg): mixture reducing and departure wiring (CoolProp-p8ub)"
```

---

## Task 10: Strictness

Every rule here is a guard that must not be able to silently no-op.

**Files:**
- Modify: `src/Backends/GERG/GERGBackend.h`, `src/Backends/GERG/GERGBackend.cpp`
- Test: `src/Tests/CoolProp-Tests-GERG.cpp`

**Interfaces:**
- Produces overrides: `calc_viscosity`, `calc_conductivity`, `calc_surface_tension`, `set_binary_interaction_double`, `set_binary_interaction_string`.

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("GERG rejects components outside the model", "[GERG]") {
    CHECK_THROWS_AS(AbstractState::factory("GERG2008", std::vector<std::string>{"R134a"}), CoolProp::ValueError);
    CHECK_THROWS_AS(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane", "R134a"}), CoolProp::ValueError);
    CHECK_THROWS_AS(AbstractState::factory("GERG2004", std::vector<std::string>{"n-Decane"}), CoolProp::ValueError);
    CHECK_NOTHROW(AbstractState::factory("GERG2008", std::vector<std::string>{"n-Decane"}));
}

TEST_CASE("GERG refuses transport properties", "[GERG]") {
    std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    AS->update(DmolarT_INPUTS, 5000.0, 200.0);
    CHECK_THROWS_AS(AS->viscosity(), CoolProp::NotImplementedError);
    CHECK_THROWS_AS(AS->conductivity(), CoolProp::NotImplementedError);
    CHECK_THROWS_AS(AS->surface_tension(), CoolProp::NotImplementedError);
}

TEST_CASE("GERG refuses binary interaction parameter mutation", "[GERG]") {
    std::shared_ptr<AbstractState> AS(
      AbstractState::factory("GERG2008", std::vector<std::string>{"Methane", "Nitrogen"}));
    CHECK_THROWS_AS(AS->set_binary_interaction_double(0, 1, "betaT", 1.1), CoolProp::ValueError);
    CHECK_THROWS_AS(AS->set_binary_interaction_double(0, 1, "gammaV", 1.1), CoolProp::ValueError);
}

TEST_CASE("GERG fluids carry no superancillary", "[GERG]") {
    // Load-bearing: FlashRoutines::sat_superanc_path_applies returns
    // superancillary values AS THE ANSWER for pure fluids that own one.
    // A borrowed blob would silently return non-GERG saturation values.
    auto* gerg = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(
      AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    std::shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> holder(gerg);
    CHECK(holder->get_superanc() == nullptr);
}

TEST_CASE("GERG respects its range of validity", "[GERG]") {
    std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    CHECK_THROWS(AS->update(PT_INPUTS, 1e5, 40.0));    // below Tmin = 60 K
    CHECK_THROWS(AS->update(PT_INPUTS, 1e5, 900.0));   // above Tmax = 700 K
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[GERG]"
```
Expected: FAIL — transport calls return numbers from CoolProp's models, BIP setters succeed.

- [ ] **Step 3: Implement the overrides**

```cpp
// In GERGBackend.h, protected section:
    CoolPropDbl calc_viscosity() override {
        throw NotImplementedError("Transport properties are not part of the GERG-2004/GERG-2008 models. "
                                  "Use the HEOS backend if you want CoolProp's transport correlations.");
    }
    CoolPropDbl calc_conductivity() override {
        throw NotImplementedError("Transport properties are not part of the GERG-2004/GERG-2008 models. "
                                  "Use the HEOS backend if you want CoolProp's transport correlations.");
    }
    CoolPropDbl calc_surface_tension() override {
        throw NotImplementedError("Surface tension is not part of the GERG-2004/GERG-2008 models.");
    }

   public:
    void set_binary_interaction_double(const std::size_t, const std::size_t, const std::string&, const double) override {
        throw ValueError("GERG binary interaction parameters are fixed by the published model and cannot be modified. "
                         "A mixture with altered beta/gamma is not GERG.");
    }
    void set_binary_interaction_string(const std::size_t, const std::size_t, const std::string&, const std::string&) override {
        throw ValueError("GERG binary interaction parameters are fixed by the published model and cannot be modified.");
    }
```

Check the exact virtual signatures in `include/CoolProp/AbstractState.h` and
`src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.h` before writing — argument
types and constness must match or the override will silently not override.
Compile with `-Woverloaded-virtual` if available and treat any warning here as
a failure, since a non-overriding "override" is exactly the fail-open case this
task exists to prevent.

Also confirm the string-keyed `set_binary_interaction_double(const std::string& CAS1, const std::string& CAS2, ...)`
overload, if one exists on `AbstractState`, is overridden too — otherwise a user
can reach the parameters by CAS and bypass the index-keyed guard.

- [ ] **Step 4: Run tests to verify they pass**

```bash
cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[GERG]"
```
Expected: PASS.

- [ ] **Step 5: Adversarially check each guard**

For every guard above, answer in writing in the commit message or a code
comment: *what input makes this pass when it should fail?* Specifically verify:

- Does `viscosity()` route through `calc_viscosity()`, or does `AbstractState`
  have a caching layer that could return a value without calling it? Read
  `AbstractState::viscosity()`.
- Does the limits check actually fire, or is it gated behind the
  `enable_limits` configuration key? If gated, the range test above passes
  only with default configuration — assert the configuration value in the test
  or drop the claim.
- Can a caller reach `Reducing->set_binary_interaction_double` directly through
  a public member? If `Reducing` is public on the base class, the guard is
  bypassable and should be noted as a known limitation.

- [ ] **Step 6: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/Backends/GERG src/Tests/CoolProp-Tests-GERG.cpp
git commit --no-verify -m "feat(gerg): strictness guards for components, BIPs and transport (CoolProp-p8ub)"
```

---

## Task 11: Ancillaries and VLE

**Files:**
- Create: `dev/gerg/fit_ancillaries.py`
- Modify: `src/Backends/GERG/GERGData.h`, `src/Backends/GERG/GERGBackend.cpp`
- Test: `src/Tests/CoolProp-Tests-GERG.cpp`

**Interfaces:**
- Produces:
  ```cpp
  struct AncillaryCoeffs {
      double T_r, reducing_value, Tmin, Tmax;
      std::vector<double> n, t;
      std::string type;  // "rhoLnoexp" | "rational" | "exponential" etc., matching CoolProp's ANCILLARIES schema
  };
  AncillaryCoeffs get_ancillary(GERGModel model, const std::string& gerg_name, const std::string& which);  // "pV" | "rhoL" | "rhoV"
  /// Evaluate an ancillary at temperature T, returning the same units as
  /// AncillaryCoeffs::reducing_value (Pa for "pV", mol/m^3 for "rhoL"/"rhoV").
  double evaluate_ancillary(const AncillaryCoeffs& anc, double T);
  ```

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("GERG pure fluid saturation converges", "[GERG]") {
    // Ancillaries seed the iteration; the converged answer is pure GERG.
    for (const auto& name : {"Methane", "Ethane", "Propane", "Nitrogen", "CarbonDioxide"}) {
        CAPTURE(name);
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{name}));
        double Tc = AS->T_critical();
        double T = 0.8 * Tc;
        REQUIRE_NOTHROW(AS->update(QT_INPUTS, 0.0, T));
        double p_sat = AS->p();
        double rhoL = AS->rhomolar();

        // Thermodynamic consistency: at (T, rhoL) the pressure must agree.
        std::shared_ptr<AbstractState> chk(AbstractState::factory("GERG2008", std::vector<std::string>{name}));
        chk->update(DmolarT_INPUTS, rhoL, T);
        CHECK_THAT(chk->p(), Catch::Matchers::WithinRel(p_sat, 1e-8));
    }
}

TEST_CASE("GERG ancillaries are close to the converged saturation state", "[GERG]") {
    // The ancillary is only a guess, but a bad guess means slow or failed
    // convergence.  Require 2% or better.
    for (const auto& name : CoolProp::GERG::component_names(CoolProp::GERGModel::GERG_2008)) {
        CAPTURE(name);
        auto anc = CoolProp::GERG::get_ancillary(CoolProp::GERGModel::GERG_2008, name, "rhoL");
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{name}));
        double T = 0.8 * AS->T_critical();
        if (T < anc.Tmin) continue;
        AS->update(QT_INPUTS, 0.0, T);
        double rho_anc = CoolProp::GERG::evaluate_ancillary(anc, T);
        CHECK_THAT(rho_anc, Catch::Matchers::WithinRel(AS->rhomolar(), 0.02));
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[GERG]"
```
Expected: FAIL — no ancillaries, so `QT_INPUTS` has no starting guess.

- [ ] **Step 3: Write the ancillary fitter**

`dev/gerg/fit_ancillaries.py`:

1. For each of the 23 distinct pure EOS, trace the saturation curve with teqp.
   Start just below the GERG reducing temperature and step down to
   `Tmin = max(60 K, 0.3 * Tc)`, using `fastchebpure`'s stepping approach
   (`~/Code/fastchebpure/fastcheb.cpp`, the `build_ancillaries` and VLE-tracing
   sections) as the reference for how to seed each step from the previous one.
2. Regress CoolProp's standard ancillary forms against the traced data.
   Read `dev/fluids/Methane.json`'s `ANCILLARIES` block first and match its
   schema exactly — same `type` strings, same `n`/`t` meaning, same
   `reducing_value` and `T_r` conventions.
3. Emit a C++ table into `GERGData.h` (or a sibling `GERGAncillaries.h` if the
   main header is getting long — split if `GERGData.h` exceeds ~2500 lines).
4. Print the maximum relative deviation per fluid. Any fluid worse than 1%
   needs more terms or a narrower `Tmin`; do not ship a fit that the Step 1
   test would fail.

Record in `dev/gerg/README.md` that these are fitted against the GERG pure EOS
specifically, and must be refitted if a coefficient table changes.

- [ ] **Step 4: Attach ancillaries in make_gerg_fluid**

Populate `fluid.ancillaries` from the table. Follow how `FluidLibrary.h` builds
`Ancillaries` from the JSON block so the same evaluation code path is used.

Also set `EOS.sat_min_liquid`, `EOS.sat_min_vapor`, `fluid.triple_liquid`,
`fluid.triple_vapor`, and `EOS.hs_anchor` — the flash routines read these.
Compute the saturation states at `Tmin` from the traced data rather than
guessing, and set `hs_anchor` at `1.1*Tc`, `0.9*rhoc` as CoolProp does
elsewhere, calling `update_states()` after construction.

- [ ] **Step 5: Run tests to verify they pass**

```bash
cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[GERG]"
```
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add dev/gerg src/Backends/GERG src/Tests/CoolProp-Tests-GERG.cpp
git commit --no-verify -m "feat(gerg): GERG-fitted saturation ancillaries and VLE (CoolProp-p8ub)"
```

---

## Task 12: Documentation and release wiring

**Files:**
- Create: `Web/coolprop/GERG.rst`
- Modify: `Web/coolprop/index.rst` (or the appropriate toctree)
- Modify: `Web/coolprop/changelog.rst`
- Modify: `CoolPropBibTeXLibrary.bib` (verify `Kunz-JCED-2012` is present; add the GERG-2004 monograph if absent)

- [ ] **Step 1: Write the documentation page**

`Web/coolprop/GERG.rst` must state, in this order:

1. What the backends are and how to call them, with a runnable example
   (`PropsSI("Dmolar","T",300,"P",1e6,"GERG2008::METHANE[0.9]&NITROGEN[0.1]")`).
2. The component lists — 18 for GERG-2004, 21 for GERG-2008 — as a table with
   CoolProp names.
3. Range of validity: normal 90-450 K and p <= 35 MPa; extended 60-700 K and
   p <= 70 MPa, which is what the backend enforces.
4. **What these backends deliberately do not provide**, and why: transport
   properties (not part of GERG), superancillaries (would require fitting
   against the GERG pure EOS), and mutable binary interaction parameters.
5. That these are distinct from CoolProp's default HEOS mixture model, which
   uses reference pure-fluid EOS and refitted binary parameters, and will give
   different numbers.
6. A pointer to teqp as the reference implementation and to the design spec.

- [ ] **Step 2: Check the bibliography**

```bash
grep -n "Kunz-JCED-2012" CoolPropBibTeXLibrary.bib
```
If the GERG-2004 monograph has no entry, add one and cite it from the new page.

- [ ] **Step 3: Build the docs**

Follow whatever `Web/` build the repo uses; at minimum confirm the new page is
reachable from a toctree and that Sphinx emits no warning about an orphaned
document.

- [ ] **Step 4: Run the full preflight**

```bash
./dev/ci/preflight.sh
```
Expected: PASS. This is the gate before pushing.

- [ ] **Step 5: Commit and push**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add Web CoolPropBibTeXLibrary.bib
git commit --no-verify -m "docs(gerg): document the strict GERG-2004/GERG-2008 backends (CoolProp-p8ub)"
git push -u origin ihb/gerg-backend
```

- [ ] **Step 6: Adversarial review, then PR**

Per CLAUDE.md this is mandatory and not optional. Review the **full** diff of
`ihb/gerg-backend` against `origin/master`, then:

```bash
gh pr create --title "feat(gerg): strict GERG-2004 / GERG-2008 backends" --body "..."
```

Re-review any commit pushed after the review pass.

- [ ] **Step 7: Close the issue**

```bash
bd close CoolProp-p8ub
bd dolt push
git push
```

---

## Follow-on issues to file

File these with `bd create` when Task 12 completes; they are explicitly out of
scope for this plan.

1. **Superancillaries for the 23 GERG pure EOS.** Requires `fastchebpure` to
   ingest the GERG pure form, `source_eos_hash` stamping, and extension of
   `dev/scripts/check_superanc_freshness.py` and
   `dev/scripts/check_superanc_release_pin.py`. Approximately 580 kB gzip added
   to the shipped data. Lights up the fast saturation path with no backend API
   change.
2. **Wrapper exposure.** Confirm the new backend families appear correctly in
   the Python, Julia, and other wrappers' backend lists, and that
   `wrappers/Python/_nanobind/CoolProp.pyi` needs no update.
3. **Phase envelopes and critical points for GERG mixtures.** Inherited from
   `HelmholtzEOSMixtureBackend` but untested by this plan.

---

## Self-Review

**Spec coverage.**

| Spec section | Task |
|---|---|
| Files / architecture | 1, 2 |
| Upstream `virtual set_mixture_parameters` | 9 |
| Data header layout, 23 distinct pures | 2, 3 |
| Backend registration, two families | 1 |
| Component naming via CoolProp aliases | 2 |
| `R*/R` ratio | 4, 8 |
| Reference state recomputation | 4, 8 |
| Strictness table (all six rows) | 10 |
| Superancillary hazard | 10 (test), follow-on issue |
| Ancillaries fitted not borrowed | 11 |
| Testing: construction, ideal gas, pures, binaries, full mixtures, strictness, cross-model | 2-6, 8, 9, 10 |
| Phasing 1-4 | Tasks 2-8, 9, 11, 10+12 |
| Out of scope items | Follow-on issues section |

No spec requirement is unaddressed.

**Placeholder scan.** No "TBD" or "handle edge cases". Three places instruct
bulk transcription from named teqp source lines rather than inlining thousands
of coefficients — each pairs with a mechanical diff step and a permanent
consistency test, so the engineer is never guessing what to write. Several
steps say "check signature X before compiling"; these are verification
instructions against a named file and line, not deferred decisions.

**Type consistency.** `GERGModel` (Task 1) is used unchanged in Tasks 2-11.
`PureInfo`/`PureCoeffs`/`AlphaigCoeffs`/`BetasGammas`/`DepartureCoeffs` field
names match teqp's structs exactly so rows copy verbatim. `get_Fij` returns
`bool` with an out-parameter in both its definition (Task 6) and its two call
sites (Tasks 6, 9). `departure_Npower` is defined in Task 6 and consumed in
Task 9. `resolve_component` is defined in Task 2 and consumed in Task 8.
`make_gerg_fluid` is defined in Task 8 and extended in Task 11.
