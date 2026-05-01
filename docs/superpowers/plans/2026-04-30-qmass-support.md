# Qmass support across CoolProp backends — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add full Qmass (mass-basis vapor quality) support to CoolProp — 8 new input pairs, a new `iQmass` keyed parameter, an explicit `Qmass()` accessor, with mixture support working from day one in both HEOS and REFPROP backends.

**Architecture:** Free-function conversions in `AbstractState.cpp` (anonymous namespace). Output is explicit. Input is explicit for pure fluids (extends `mass_to_molar_inputs`) and iterative for mixtures via a virtual `update_Qmass_pair` on AbstractState (default secant-on-Qmolar). REFPROP overrides for `QmassT`/`PQmass` to use REFPROP's native `kq=2` flag in `TQFLSHdll`/`PQFLSHdll`; the other 6 mixture pairs fall through to the base iteration. HEOS uses base iteration for all 8 pairs.

**Tech Stack:** C++17, CMake, Catch2 v3 (`catch2/catch_all.hpp`), CoolProp's existing `Brent` solver in `src/Solvers.cpp`.

**Spec:** `docs/superpowers/specs/2026-04-30-qmass-support-design.md`

---

## File Map

**Created:**
- `docs/superpowers/specs/2026-04-30-qmass-support-design.md` (already exists)
- `docs/superpowers/plans/2026-04-30-qmass-support.md` (this file)

**Modified:**
- `include/DataStructures.h` — `iQmass` parameter; 8 new `input_pairs` entries; 8 new `match_pair` branches in `generate_update_pair`
- `src/DataStructures.cpp` — `iQmass` row in `parameter_info_list`; `index_map_insert("Qmass", iQmass)`; 8 new rows in `input_pair_list`
- `include/AbstractState.h` — `_Qmass = cache.next()` cache slot; `Qmass()` public; protected `calc_Qmass`, `calc_phase_molar_masses`, `update_Qmass_pair`, `is_Qmass_pair`, `is_pure_or_pseudopure` helpers
- `src/AbstractState.cpp` — anonymous-namespace free fns `Qmolar_to_Qmass` / `Qmass_to_Qmolar`; `Qmass()`, `calc_Qmass()`, `calc_phase_molar_masses()` (default = pure-fluid trivial); 8 Qmass cases in `mass_to_molar_inputs`; default `update_Qmass_pair` secant; `iQmass` case in `keyed_output`
- `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.h` — declare `calc_phase_molar_masses` override
- `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp` — early-exit hook in `update()`; `calc_phase_molar_masses` override using `components[i]->EOS()->molar_mass`
- `src/Backends/REFPROP/REFPROPMixtureBackend.h` — declare `calc_phase_molar_masses` and `update_Qmass_pair` overrides
- `src/Backends/REFPROP/REFPROPMixtureBackend.cpp` — early-exit hook in `update()`; `calc_phase_molar_masses` override using REFPROP component MMs (`WMMdll` per component); `update_Qmass_pair` override calling `TQFLSHdll`/`PQFLSHdll` with `kq=2` for `QmassT`/`PQmass`, deferring to `AbstractState::update_Qmass_pair` for other 6 pairs
- `src/Tests/CoolProp-Tests.cpp` — add Qmass test cases (pure round-trip, mixture round-trip HEOS + REFPROP, edge cases)

---

## Task 1: Branch setup

**Files:** none (git only)

- [ ] **Step 1: Verify clean tree and branch from latest master**

```bash
git -C /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass status
git -C /Users/ianbell/Code/CoolProp fetch origin master
git -C /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass checkout -b ib/qmass-support origin/master
```

Expected: `On branch ib/qmass-support; Your branch is up to date with 'origin/master'.`

- [ ] **Step 2: Confirm spec and plan files are present**

```bash
ls /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/docs/superpowers/specs/2026-04-30-qmass-support-design.md \
   /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/docs/superpowers/plans/2026-04-30-qmass-support.md
```

If they don't exist on the new branch (because they were committed only to the previous branch), copy them in from the worktree's prior state and stage them — they will be committed in the final docs task.

---

## Task 2: Add `iQmass` parameter and 8 new `input_pairs` enum entries

**Files:**
- Modify: `include/DataStructures.h`

- [ ] **Step 1: Add `iQmass` to the parameters enum**

In `include/DataStructures.h`, around line 91, after `iQ`:

```cpp
    iQ,      ///< Vapor quality (molar; alias for iQmolar in this codebase)
    iQmass,  ///< Mass-basis vapor quality
```

(Update the inline comment on `iQ` for clarity — it is unambiguously molar.)

- [ ] **Step 2: Add 8 new `input_pairs` entries**

In `include/DataStructures.h`, in the `enum input_pairs : int { ... }` block (around line 280), add each Qmass variant immediately after its molar sibling:

```cpp
    QT_INPUTS,        ///< Molar quality, Temperature in K
    QmassT_INPUTS,    ///< Mass-basis quality, Temperature in K
    PQ_INPUTS,        ///< Pressure in Pa, Molar quality
    PQmass_INPUTS,    ///< Pressure in Pa, Mass-basis quality
    QSmolar_INPUTS,   ///< Molar quality, Entropy in J/mol/K
    QmassSmolar_INPUTS, ///< Mass-basis quality, Entropy in J/mol/K
    QSmass_INPUTS,    ///< Molar quality, Entropy in J/kg/K
    QmassSmass_INPUTS, ///< Mass-basis quality, Entropy in J/kg/K
    HmolarQ_INPUTS,   ///< Enthalpy in J/mol, Molar quality
    HmolarQmass_INPUTS, ///< Enthalpy in J/mol, Mass-basis quality
    HmassQ_INPUTS,    ///< Enthalpy in J/kg, Molar quality
    HmassQmass_INPUTS, ///< Enthalpy in J/kg, Mass-basis quality
    DmolarQ_INPUTS,   ///< Density in mol/m^3, Molar quality
    DmolarQmass_INPUTS, ///< Density in mol/m^3, Mass-basis quality
    DmassQ_INPUTS,    ///< Density in kg/m^3, Molar quality
    DmassQmass_INPUTS, ///< Density in kg/m^3, Mass-basis quality
```

- [ ] **Step 3: Verify the enum compiles**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass
mkdir -p build && cd build
cmake .. -DCOOLPROP_STATIC_LIBRARY=ON -DCMAKE_BUILD_TYPE=Release > /dev/null
cmake --build . --target CoolProp -j 4 2>&1 | tail -30
```

Expected: clean build (warnings about unused enums in switch statements are OK at this stage; they will be fixed in subsequent tasks).

- [ ] **Step 4: Commit**

```bash
git add include/DataStructures.h
git commit -m "qmass: add iQmass parameter and 8 new input_pairs enum entries

Adds the enum surface for mass-basis vapor quality. Each new Qmass-input
pair sits immediately next to its molar sibling. Wiring (string maps,
generate_update_pair, conversion logic) follows in subsequent commits."
```

---

## Task 3: Wire enums into `parameter_info_list` and `input_pair_list`

**Files:**
- Modify: `src/DataStructures.cpp`

- [ ] **Step 1: Add `iQmass` row to `parameter_info_list`**

In `src/DataStructures.cpp`, in `parameter_info_list` (around line 34, after the `iQ` entry):

```cpp
  {iQ, "Q", "IO", "mol/mol", "Molar vapor quality", false},
  {iQmass, "Qmass", "IO", "kg/kg", "Mass-basis vapor quality", false},
```

- [ ] **Step 2: Register the alias map entry for `"Qmass"`**

Find the `index_map_insert(...)` calls (around line 132 — `index_map_insert("S", iSmass);` is one of them). Verify `iQ` already has an entry and add a parallel one for `iQmass` if needed. If `iQ` has `index_map_insert("Q", iQ);`, add `index_map_insert("Qmass", iQmass);` next to it. If no aliasing is needed beyond the `parameter_info_list` row, this step is a no-op — skip and move to step 3.

```bash
grep -n 'index_map_insert' /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/src/DataStructures.cpp
```

- [ ] **Step 3: Add 8 rows to `input_pair_list`**

In `src/DataStructures.cpp`, in the `input_pair_list[]` array (around line 508), add each Qmass entry next to its molar sibling. Match the existing convention of using the same short name for mass-and-molar variants where they share one (note `QSmolar_INPUTS` and `QSmass_INPUTS` both use `"QS_INPUTS"` as short name today; preserve that pattern for Qmass variants too):

```cpp
  {QT_INPUTS, "QT_INPUTS", "Molar quality, Temperature in K"},
  {QmassT_INPUTS, "QmassT_INPUTS", "Mass-basis quality, Temperature in K"},
  {QSmolar_INPUTS, "QS_INPUTS", "Molar quality, Entropy in J/mol/K"},
  {QmassSmolar_INPUTS, "QmassS_INPUTS", "Mass-basis quality, Entropy in J/mol/K"},
  {QSmass_INPUTS, "QS_INPUTS", "Molar quality, Entropy in J/kg/K"},
  {QmassSmass_INPUTS, "QmassS_INPUTS", "Mass-basis quality, Entropy in J/kg/K"},
  {HmolarQ_INPUTS, "HQ_INPUTS", "Enthalpy in J/mol, Molar quality"},
  {HmolarQmass_INPUTS, "HQmass_INPUTS", "Enthalpy in J/mol, Mass-basis quality"},
  {HmassQ_INPUTS, "HQ_INPUTS", "Enthalpy in J/kg, Molar quality"},
  {HmassQmass_INPUTS, "HQmass_INPUTS", "Enthalpy in J/kg, Mass-basis quality"},
  {DmassQ_INPUTS, "DmassQ_INPUTS", "Mass density kg/m^3, Molar quality"},
  {DmassQmass_INPUTS, "DmassQmass_INPUTS", "Mass density kg/m^3, Mass-basis quality"},
  {DmolarQ_INPUTS, "DmolarQ_INPUTS", "Molar density in mol/m^3, Molar quality"},
  {DmolarQmass_INPUTS, "DmolarQmass_INPUTS", "Molar density in mol/m^3, Mass-basis quality"},

  {PQ_INPUTS, "PQ_INPUTS", "Pressure in Pa, Molar quality"},
  {PQmass_INPUTS, "PQmass_INPUTS", "Pressure in Pa, Mass-basis quality"},
```

- [ ] **Step 4: Build to verify**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/build
cmake --build . --target CoolProp -j 4 2>&1 | tail -20
```

Expected: clean build.

- [ ] **Step 5: Commit**

```bash
git add src/DataStructures.cpp
git commit -m "qmass: register iQmass parameter and 8 new input_pair_list rows"
```

---

## Task 4: Add 8 new `match_pair` branches to `generate_update_pair`

**Files:**
- Modify: `include/DataStructures.h`

- [ ] **Step 1: Add a Qmass branch next to each existing Q branch**

In `include/DataStructures.h`, in the `generate_update_pair` template (around line 343), add `iQmass` branches:

```cpp
    if (match_pair(key1, key2, iQ, iT, swap)) {
        pair = QT_INPUTS;
    } else if (match_pair(key1, key2, iQmass, iT, swap)) {
        pair = QmassT_INPUTS;
    } else if (match_pair(key1, key2, iP, iQ, swap)) {
        pair = PQ_INPUTS;
    } else if (match_pair(key1, key2, iP, iQmass, swap)) {
        pair = PQmass_INPUTS;
    } else if (match_pair(key1, key2, iP, iT, swap)) {
        pair = PT_INPUTS;
    }
    // ... continue: for each existing iQ-bearing match_pair, add a parallel
    // iQmass branch that returns the new <X>Qmass_INPUTS / Qmass<X>_INPUTS pair.
```

Make sure to add Qmass branches for these 8 combinations (the 6 not yet shown):

| iQmass paired with | Returns |
|---|---|
| `iSmolar` | `QmassSmolar_INPUTS` |
| `iSmass` | `QmassSmass_INPUTS` |
| `iHmolar` | `HmolarQmass_INPUTS` |
| `iHmass` | `HmassQmass_INPUTS` |
| `iDmolar` | `DmolarQmass_INPUTS` |
| `iDmass` | `DmassQmass_INPUTS` |

(`iQmass`+`iT` and `iP`+`iQmass` are already added above.)

- [ ] **Step 2: Build**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/build
cmake --build . --target CoolProp -j 4 2>&1 | tail -10
```

Expected: clean build.

- [ ] **Step 3: Commit**

```bash
git add include/DataStructures.h
git commit -m "qmass: add iQmass branches to generate_update_pair"
```

---

## Task 5: Add `_Qmass` cache slot and `Qmass()` accessor (header)

**Files:**
- Modify: `include/AbstractState.h`

- [ ] **Step 1: Add cache slot next to other two-phase variables**

In `include/AbstractState.h`, around line 149 (after `_rhoLmolar`/`_rhoVmolar`), add:

```cpp
    /// Two-Phase variables
    CAE _rhoLmolar = cache.next(), _rhoVmolar = cache.next();
    CAE _Qmass = cache.next();
```

- [ ] **Step 2: Declare protected helpers**

In the protected section of `AbstractState`, near other `calc_*` virtuals (search for `calc_molar_mass` around line 232):

```cpp
    /// Mass-basis vapor quality. Default implementation uses _Q (molar) and
    /// calc_phase_molar_masses(); override only if a backend has a faster route.
    virtual CoolPropDbl calc_Qmass(void);

    /// Populate phase molar masses (kg/mol) for the current saturated state.
    /// Default base impl handles pure/pseudopure (both = molar_mass()).
    /// Mixture backends must override.
    virtual void calc_phase_molar_masses(double& MM_liquid, double& MM_vapor);

    /// Default iterative Qmass-pair solver (secant on Qmolar). Backends may
    /// override to use a native fast path (e.g. REFPROP TQFLSHdll kq=2).
    virtual void update_Qmass_pair(CoolProp::input_pairs pair, double v1, double v2);
```

- [ ] **Step 3: Declare public `Qmass()` accessor**

In the public section, near `Q()` (search for `return _Q;` around line 1082):

```cpp
    /// Mass-basis vapor quality (kg vapor / kg total). Throws if not two-phase.
    double Qmass(void);
```

- [ ] **Step 4: Compile-check the header (no impl yet → expect link errors only when impl is missing)**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/build
cmake --build . --target CoolProp -j 4 2>&1 | tail -20
```

Expected: header compiles. Linker may complain about missing `calc_Qmass`/`calc_phase_molar_masses`/`update_Qmass_pair`/`Qmass` symbols — that is fine, those are the next task. If the build target is just the static library, missing symbols are deferred until link time of consumers, so the static library build itself should still succeed for now.

- [ ] **Step 5: Commit**

```bash
git add include/AbstractState.h
git commit -m "qmass: declare _Qmass cache slot, Qmass() accessor and helpers"
```

---

## Task 6: Implement free-function conversions and `Qmass()` / `calc_Qmass` / `calc_phase_molar_masses` (default base impl)

**Files:**
- Modify: `src/AbstractState.cpp`

- [ ] **Step 1: Add anonymous-namespace free functions**

At the top of `src/AbstractState.cpp`, after the `#include` lines and `namespace CoolProp {`:

```cpp
namespace {

/// Mass-basis vapor quality from molar quality and phase molar masses.
/// All molar masses in the same units (kg/mol).
inline double Qmolar_to_Qmass(double Qmolar, double MM_liquid, double MM_vapor) {
    const double mv = Qmolar * MM_vapor;
    const double ml = (1.0 - Qmolar) * MM_liquid;
    return mv / (mv + ml);
}

/// Inverse of Qmolar_to_Qmass. Same units convention.
inline double Qmass_to_Qmolar(double Qmass, double MM_liquid, double MM_vapor) {
    const double nv = Qmass / MM_vapor;
    const double nl = (1.0 - Qmass) / MM_liquid;
    return nv / (nv + nl);
}

} // anonymous namespace
```

- [ ] **Step 2: Implement `Qmass()` (public) and `calc_Qmass()` (default)**

Add to `src/AbstractState.cpp`, alongside other definitions like `Q()` (search for `_Q;` to locate the area):

```cpp
double AbstractState::Qmass(void) {
    if (!_Qmass) _Qmass = calc_Qmass();
    return _Qmass;
}

CoolPropDbl AbstractState::calc_Qmass(void) {
    if (!ValidNumber(_Q) || _Q < 0 || _Q > 1) {
        throw ValueError("Qmass requires a two-phase state (0 <= Q <= 1)");
    }
    double MM_l = 0, MM_v = 0;
    calc_phase_molar_masses(MM_l, MM_v);
    return static_cast<CoolPropDbl>(Qmolar_to_Qmass(_Q, MM_l, MM_v));
}
```

- [ ] **Step 3: Implement default `calc_phase_molar_masses` (pure-fluid)**

Append to `src/AbstractState.cpp`:

```cpp
void AbstractState::calc_phase_molar_masses(double& MM_l, double& MM_v) {
    if (_fluid_type == FLUID_TYPE_PURE || _fluid_type == FLUID_TYPE_PSEUDOPURE) {
        const double mm = molar_mass();
        MM_l = mm;
        MM_v = mm;
        return;
    }
    throw NotImplementedError("calc_phase_molar_masses must be overridden by mixture backends");
}
```

- [ ] **Step 4: Wire `iQmass` into `keyed_output`**

Find `keyed_output` in `src/AbstractState.cpp` (search for `case iQ:` — around line 720+). Add a `case iQmass: return Qmass();` next to it.

```cpp
        case iQ:
            return Q();
        case iQmass:
            return Qmass();
```

- [ ] **Step 5: Add a test for pure-fluid `Qmass()`**

In `src/Tests/CoolProp-Tests.cpp`, add at the end of the file (before any closing namespaces):

```cpp
TEST_CASE("Qmass output: pure fluid equals Qmolar", "[Qmass]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    for (double Q : {0.0, 0.1, 0.5, 0.9, 1.0}) {
        AS->update(CoolProp::QT_INPUTS, Q, 350.0);
        CHECK(AS->Qmass() == Approx(Q).epsilon(1e-12));
    }
}
```

- [ ] **Step 6: Build and run the test**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/build
cmake --build . --target CoolProp-Tests -j 4 2>&1 | tail -20
./src/Tests/CoolProp-Tests "[Qmass]"
```

Expected: 1 test case, all assertions PASS.

- [ ] **Step 7: Commit**

```bash
git add src/AbstractState.cpp src/Tests/CoolProp-Tests.cpp
git commit -m "qmass: implement Qmass() output for pure fluids

Adds free-function Qmolar<->Qmass conversions, the Qmass() accessor with
its calc_Qmass default, calc_phase_molar_masses (pure-fluid trivial path),
and iQmass dispatch in keyed_output. Mixture backends must override
calc_phase_molar_masses; that lands in the next task."
```

---

## Task 7: HEOS `calc_phase_molar_masses` override + mixture output test

**Files:**
- Modify: `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.h`
- Modify: `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp`

- [ ] **Step 1: Declare override**

In `HelmholtzEOSMixtureBackend.h`, in the protected section near other `calc_*` overrides:

```cpp
    void calc_phase_molar_masses(double& MM_liquid, double& MM_vapor) override;
```

- [ ] **Step 2: Implement override**

In `HelmholtzEOSMixtureBackend.cpp`, near `calc_molar_mass`:

```cpp
void HelmholtzEOSMixtureBackend::calc_phase_molar_masses(double& MM_l, double& MM_v) {
    if (_fluid_type == FLUID_TYPE_PURE || _fluid_type == FLUID_TYPE_PSEUDOPURE) {
        const double mm = molar_mass();
        MM_l = mm;
        MM_v = mm;
        return;
    }
    const std::vector<CoolPropDbl> x = mole_fractions_liquid();
    const std::vector<CoolPropDbl> y = mole_fractions_vapor();
    if (x.size() != components.size() || y.size() != components.size()) {
        throw ValueError("phase composition vectors do not match component count");
    }
    MM_l = 0;
    MM_v = 0;
    for (std::size_t i = 0; i < components.size(); ++i) {
        const double mm_i = components[i].EOS().molar_mass;
        MM_l += static_cast<double>(x[i]) * mm_i;
        MM_v += static_cast<double>(y[i]) * mm_i;
    }
}
```

(If `components[i].EOS().molar_mass` is the wrong accessor, use `components[i].EOSVector[0].molar_mass` or whatever the surrounding `calc_molar_mass` uses; mirror the existing pattern verbatim.)

- [ ] **Step 3: Add mixture output test**

Append to `src/Tests/CoolProp-Tests.cpp`:

```cpp
TEST_CASE("Qmass output: HEOS mixture round-trips with manual Qmolar->Qmass", "[Qmass][mixture]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "R32&R125"));
    std::vector<CoolPropDbl> z = {0.5, 0.5};
    AS->set_mole_fractions(z);
    for (double Q : {0.1, 0.3, 0.5, 0.7, 0.9}) {
        AS->update(CoolProp::QT_INPUTS, Q, 280.0);
        const auto x = AS->mole_fractions_liquid();
        const auto y = AS->mole_fractions_vapor();
        // Manual MM weighting using the published REFPROP/HEOS component MMs:
        //   R32 ~ 0.05202 kg/mol, R125 ~ 0.12002 kg/mol
        // (We re-read from molar_mass at z to avoid hard-coding):
        const double MM_l =
            static_cast<double>(x[0]) * 0.052024 + static_cast<double>(x[1]) * 0.120021;
        const double MM_v =
            static_cast<double>(y[0]) * 0.052024 + static_cast<double>(y[1]) * 0.120021;
        const double Qmass_expected = (Q * MM_v) / (Q * MM_v + (1.0 - Q) * MM_l);
        CHECK(AS->Qmass() == Approx(Qmass_expected).epsilon(1e-6));
        CHECK(AS->Qmass() != Approx(Q).epsilon(1e-3)); // mixture should differ from molar
    }
}
```

- [ ] **Step 4: Build and run**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/build
cmake --build . --target CoolProp-Tests -j 4 2>&1 | tail -10
./src/Tests/CoolProp-Tests "[Qmass]"
```

Expected: 2 test cases, all PASS.

- [ ] **Step 5: Commit**

```bash
git add src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.h \
        src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp \
        src/Tests/CoolProp-Tests.cpp
git commit -m "qmass: HEOS calc_phase_molar_masses override + mixture output test"
```

---

## Task 8: Pure-fluid path: extend `mass_to_molar_inputs` for 8 Qmass cases

**Files:**
- Modify: `src/AbstractState.cpp`

- [ ] **Step 1: Add a pure-fluid early-exit branch**

`mass_to_molar_inputs` lives at `src/AbstractState.cpp:200`. Add a separate `if` ahead of the existing `switch`, scoped to pure/pseudopure fluids only — for these, `Qmass = Qmolar` exactly, so the rewrite is value-preserving:

```cpp
void AbstractState::mass_to_molar_inputs(CoolProp::input_pairs& input_pair, CoolPropDbl& value1, CoolPropDbl& value2) {
    // Pure / pseudopure: Qmass and Qmolar are identical. Rewrite the pair
    // to its molar sibling without touching values, then continue.
    if (_fluid_type == FLUID_TYPE_PURE || _fluid_type == FLUID_TYPE_PSEUDOPURE) {
        switch (input_pair) {
            case QmassT_INPUTS:        input_pair = QT_INPUTS;        break;
            case PQmass_INPUTS:        input_pair = PQ_INPUTS;        break;
            case QmassSmolar_INPUTS:   input_pair = QSmolar_INPUTS;   break;
            case QmassSmass_INPUTS:    input_pair = QSmass_INPUTS;    break;
            case HmolarQmass_INPUTS:   input_pair = HmolarQ_INPUTS;   break;
            case HmassQmass_INPUTS:    input_pair = HmassQ_INPUTS;    break;
            case DmolarQmass_INPUTS:   input_pair = DmolarQ_INPUTS;   break;
            case DmassQmass_INPUTS:    input_pair = DmassQ_INPUTS;    break;
            default:                   break;
        }
    }
    // ... existing switch on Smass/Hmass/Dmass conversions ...
}
```

- [ ] **Step 2: Add a pure-fluid Qmass-input round-trip test**

Append to `src/Tests/CoolProp-Tests.cpp`:

```cpp
TEST_CASE("Qmass input: pure fluid round-trips for all 8 pairs", "[Qmass][pure]") {
    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    auto sut = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    struct Pair { CoolProp::input_pairs molar; CoolProp::input_pairs mass; CoolProp::parameters partner; };
    const std::vector<Pair> pairs = {
        {CoolProp::QT_INPUTS,        CoolProp::QmassT_INPUTS,      CoolProp::iT},
        {CoolProp::PQ_INPUTS,        CoolProp::PQmass_INPUTS,      CoolProp::iP},
        {CoolProp::QSmolar_INPUTS,   CoolProp::QmassSmolar_INPUTS, CoolProp::iSmolar},
        {CoolProp::QSmass_INPUTS,    CoolProp::QmassSmass_INPUTS,  CoolProp::iSmass},
        {CoolProp::HmolarQ_INPUTS,   CoolProp::HmolarQmass_INPUTS, CoolProp::iHmolar},
        {CoolProp::HmassQ_INPUTS,    CoolProp::HmassQmass_INPUTS,  CoolProp::iHmass},
        {CoolProp::DmolarQ_INPUTS,   CoolProp::DmolarQmass_INPUTS, CoolProp::iDmolar},
        {CoolProp::DmassQ_INPUTS,    CoolProp::DmassQmass_INPUTS,  CoolProp::iDmass},
    };

    for (const auto& p : pairs) {
        // Establish reference state via the existing molar pair.
        // Ordering of (Q, partner) varies per pair; mirror generate_update_pair convention.
        ref->update(CoolProp::QT_INPUTS, 0.4, 350.0);
        const double partner_value = ref->keyed_output(p.partner);
        const double Qmass_value = ref->Qmass();

        // Re-establish via molar pair to extract reference T, p, etc.
        ref->update(CoolProp::QT_INPUTS, 0.4, 350.0);
        const double T_ref = ref->T();
        const double p_ref = ref->p();
        const double Q_ref = ref->Q();

        // Drive the SUT via the new Qmass pair. For pure fluids Qmass==Qmolar.
        // Argument order matches the molar pair's existing argument convention.
        if (p.mass == CoolProp::QmassT_INPUTS) {
            sut->update(p.mass, Qmass_value, 350.0);
        } else if (p.mass == CoolProp::PQmass_INPUTS) {
            sut->update(p.mass, p_ref, Qmass_value);
        } else {
            // For HQ/DQ/QS pairs the partner is value1 or value2 per existing
            // convention in generate_update_pair / update; use partner first
            // for pairs where Q/Qmass is value2 in the molar version, and
            // partner second otherwise. Inspect existing tests / handlers
            // and mirror.
            sut->update(p.mass, partner_value, Qmass_value); // adjust as needed
        }
        CHECK(sut->T() == Approx(T_ref).epsilon(1e-10));
        CHECK(sut->p() == Approx(p_ref).epsilon(1e-10));
        CHECK(sut->Q() == Approx(Q_ref).epsilon(1e-10));
    }
}
```

(If the argument-order convention for some pairs causes a `ValueError`, swap the two values in the `sut->update` call. The molar version of the same pair sets the convention — match it exactly.)

- [ ] **Step 3: Build and run**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/build
cmake --build . --target CoolProp-Tests -j 4 2>&1 | tail -10
./src/Tests/CoolProp-Tests "[Qmass][pure]"
```

Expected: PASS for all 8 pairs.

- [ ] **Step 4: Commit**

```bash
git add src/AbstractState.cpp src/Tests/CoolProp-Tests.cpp
git commit -m "qmass: pure-fluid input pairs route through molar siblings (Qmass == Qmolar)"
```

---

## Task 9: Mixture iteration: default `update_Qmass_pair` + helpers + HEOS hook

**Files:**
- Modify: `include/AbstractState.h` (already declared in Task 5)
- Modify: `src/AbstractState.cpp`
- Modify: `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp`

- [ ] **Step 1: Add `is_Qmass_pair` and `is_pure_or_pseudopure` helpers**

In `include/AbstractState.h` (protected section) and `src/AbstractState.cpp`:

```cpp
// header
inline bool is_Qmass_pair(CoolProp::input_pairs p) {
    switch (p) {
        case CoolProp::QmassT_INPUTS:
        case CoolProp::PQmass_INPUTS:
        case CoolProp::QmassSmolar_INPUTS:
        case CoolProp::QmassSmass_INPUTS:
        case CoolProp::HmolarQmass_INPUTS:
        case CoolProp::HmassQmass_INPUTS:
        case CoolProp::DmolarQmass_INPUTS:
        case CoolProp::DmassQmass_INPUTS:
            return true;
        default:
            return false;
    }
}

// AbstractState method:
bool is_pure_or_pseudopure(void) const {
    return _fluid_type == FLUID_TYPE_PURE || _fluid_type == FLUID_TYPE_PSEUDOPURE;
}
```

`is_Qmass_pair` is a free inline (in `DataStructures.h` is more appropriate — put it there alongside `match_pair`). `is_pure_or_pseudopure` is a member.

- [ ] **Step 2: Implement default `update_Qmass_pair`**

In `src/AbstractState.cpp`:

```cpp
void AbstractState::update_Qmass_pair(CoolProp::input_pairs pair, double v1, double v2) {
    // Identify which slot holds Qmass and what the molar-equivalent pair is.
    // Mapping: each Qmass pair -> molar sibling, with the Q-slot index.
    struct Mapping { CoolProp::input_pairs molar; int qmass_slot; }; // 1 or 2
    Mapping m;
    switch (pair) {
        case CoolProp::QmassT_INPUTS:        m = {CoolProp::QT_INPUTS,        1}; break;
        case CoolProp::PQmass_INPUTS:        m = {CoolProp::PQ_INPUTS,        2}; break;
        case CoolProp::QmassSmolar_INPUTS:   m = {CoolProp::QSmolar_INPUTS,   1}; break;
        case CoolProp::QmassSmass_INPUTS:    m = {CoolProp::QSmass_INPUTS,    1}; break;
        case CoolProp::HmolarQmass_INPUTS:   m = {CoolProp::HmolarQ_INPUTS,   2}; break;
        case CoolProp::HmassQmass_INPUTS:    m = {CoolProp::HmassQ_INPUTS,    2}; break;
        case CoolProp::DmolarQmass_INPUTS:   m = {CoolProp::DmolarQ_INPUTS,   2}; break;
        case CoolProp::DmassQmass_INPUTS:    m = {CoolProp::DmassQ_INPUTS,    2}; break;
        default:
            throw ValueError("update_Qmass_pair called with non-Qmass pair");
    }
    const double Qmass_target = (m.qmass_slot == 1) ? v1 : v2;
    const double partner = (m.qmass_slot == 1) ? v2 : v1;

    if (Qmass_target < 0 || Qmass_target > 1) {
        throw ValueError(format("Qmass out of range [0,1]: %g", Qmass_target));
    }
    // Endpoints bypass iteration.
    if (Qmass_target == 0.0 || Qmass_target == 1.0) {
        const double Qmolar = Qmass_target;
        if (m.qmass_slot == 1) update(m.molar, Qmolar, partner);
        else                   update(m.molar, partner, Qmolar);
        _Qmass = Qmass_target;
        return;
    }

    // Secant on Qmolar. Initial bracket [Qmass_target, Qmass_target+small].
    auto residual = [&](double Qmolar) -> double {
        if (m.qmass_slot == 1) update(m.molar, Qmolar, partner);
        else                   update(m.molar, partner, Qmolar);
        double MM_l = 0, MM_v = 0;
        calc_phase_molar_masses(MM_l, MM_v);
        const double Qmass_calc = Qmolar_to_Qmass(Qmolar, MM_l, MM_v);
        return Qmass_calc - Qmass_target;
    };

    double a = std::max(1e-12, Qmass_target - 1e-3);
    double b = std::min(1.0 - 1e-12, Qmass_target + 1e-3);
    double fa = residual(a);
    double fb = residual(b);
    if (fa * fb > 0) {
        // Widen to full bracket and use bisection.
        a = 1e-12; b = 1.0 - 1e-12;
        fa = residual(a); fb = residual(b);
        if (fa * fb > 0) {
            throw SolutionError(format("update_Qmass_pair: cannot bracket Qmolar for target Qmass=%g", Qmass_target));
        }
    }

    double Qmolar_solution = std::numeric_limits<double>::quiet_NaN();
    for (int iter = 0; iter < 50; ++iter) {
        // Secant step
        const double denom = fb - fa;
        const double c = (std::abs(denom) < 1e-30) ? 0.5 * (a + b) : b - fb * (b - a) / denom;
        const double c_clamped = std::max(1e-12, std::min(1.0 - 1e-12, c));
        const double fc = residual(c_clamped);
        if (std::abs(fc) < 1e-10) {
            Qmolar_solution = c_clamped;
            break;
        }
        // Maintain bracket
        if (fa * fc < 0) { b = c_clamped; fb = fc; }
        else             { a = c_clamped; fa = fc; }
    }
    if (!ValidNumber(Qmolar_solution)) {
        throw SolutionError(format("update_Qmass_pair: did not converge for Qmass=%g", Qmass_target));
    }
    // Ensure final state corresponds to the converged Qmolar.
    if (m.qmass_slot == 1) update(m.molar, Qmolar_solution, partner);
    else                   update(m.molar, partner, Qmolar_solution);
    _Qmass = Qmass_target;
}
```

- [ ] **Step 3: Wire HEOS `update()` early-exit**

In `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp`, at the top of `HelmholtzEOSMixtureBackend::update`, BEFORE `mass_to_molar_inputs` is called:

```cpp
void HelmholtzEOSMixtureBackend::update(CoolProp::input_pairs input_pair, double value1, double value2) {
    if (CoolProp::is_Qmass_pair(input_pair) && !is_pure_or_pseudopure()) {
        update_Qmass_pair(input_pair, value1, value2);
        return;
    }
    // ... existing body ...
}
```

- [ ] **Step 4: Add HEOS mixture round-trip test**

Append to `src/Tests/CoolProp-Tests.cpp`:

```cpp
TEST_CASE("Qmass input: HEOS R32+R125 round-trip via QmassT_INPUTS", "[Qmass][mixture]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "R32&R125"));
    AS->set_mole_fractions({0.5, 0.5});
    AS->update(CoolProp::QT_INPUTS, 0.4, 280.0);
    const double Qmass_observed = AS->Qmass();
    const double p_ref = AS->p();
    const double Q_ref = AS->Q();

    auto AS2 = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "R32&R125"));
    AS2->set_mole_fractions({0.5, 0.5});
    AS2->update(CoolProp::QmassT_INPUTS, Qmass_observed, 280.0);
    CHECK(AS2->p()    == Approx(p_ref).epsilon(1e-8));
    CHECK(AS2->Q()    == Approx(Q_ref).epsilon(1e-8));
    CHECK(AS2->Qmass()== Approx(Qmass_observed).epsilon(1e-10));
}
```

- [ ] **Step 5: Build and run**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/build
cmake --build . --target CoolProp-Tests -j 4 2>&1 | tail -10
./src/Tests/CoolProp-Tests "[Qmass][mixture]"
```

Expected: PASS in 5–8 secant iterations per call.

- [ ] **Step 6: Commit**

```bash
git add include/AbstractState.h include/DataStructures.h src/AbstractState.cpp \
        src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp \
        src/Tests/CoolProp-Tests.cpp
git commit -m "qmass: HEOS mixture support via secant on Qmolar (default base impl)"
```

---

## Task 10: Add HEOS mixture round-trip tests for the other 3 pair families

**Files:**
- Modify: `src/Tests/CoolProp-Tests.cpp`

- [ ] **Step 1: Add round-trip tests for `PQmass`, `HmolarQmass`, `DmolarQmass`**

Append to `src/Tests/CoolProp-Tests.cpp`:

```cpp
TEST_CASE("Qmass input: HEOS mixture round-trip via PQmass / HmolarQmass / DmolarQmass", "[Qmass][mixture]") {
    auto setup = []() {
        auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "R32&R125"));
        AS->set_mole_fractions({0.5, 0.5});
        return AS;
    };

    SECTION("PQmass") {
        auto ref = setup();
        ref->update(CoolProp::PQ_INPUTS, 1.5e6, 0.4);
        const double Qmass = ref->Qmass();
        const double T_ref = ref->T();
        auto sut = setup();
        sut->update(CoolProp::PQmass_INPUTS, 1.5e6, Qmass);
        CHECK(sut->T() == Approx(T_ref).epsilon(1e-8));
        CHECK(sut->Q() == Approx(0.4).epsilon(1e-8));
    }

    SECTION("HmolarQmass") {
        auto ref = setup();
        ref->update(CoolProp::HmolarQ_INPUTS, ref->hmolar() /* placeholder; will set below */, 0.4);
        // To get a sensible Hmolar, do a TQ flash first and read back:
        ref->update(CoolProp::QT_INPUTS, 0.4, 280.0);
        const double H = ref->hmolar();
        const double Qmass = ref->Qmass();
        ref->update(CoolProp::HmolarQ_INPUTS, H, 0.4);
        const double T_ref = ref->T();
        auto sut = setup();
        sut->update(CoolProp::HmolarQmass_INPUTS, H, Qmass);
        CHECK(sut->T() == Approx(T_ref).epsilon(1e-8));
        CHECK(sut->Q() == Approx(0.4).epsilon(1e-6));
    }

    SECTION("DmolarQmass") {
        auto ref = setup();
        ref->update(CoolProp::QT_INPUTS, 0.4, 280.0);
        const double D = ref->rhomolar();
        const double Qmass = ref->Qmass();
        const double T_ref = ref->T();
        auto sut = setup();
        sut->update(CoolProp::DmolarQmass_INPUTS, D, Qmass);
        CHECK(sut->T() == Approx(T_ref).epsilon(1e-8));
        CHECK(sut->Q() == Approx(0.4).epsilon(1e-6));
    }
}
```

- [ ] **Step 2: Build and run**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/build
cmake --build . --target CoolProp-Tests -j 4 2>&1 | tail -10
./src/Tests/CoolProp-Tests "[Qmass][mixture]"
```

Expected: PASS for all sections.

- [ ] **Step 3: Commit**

```bash
git add src/Tests/CoolProp-Tests.cpp
git commit -m "qmass: extend HEOS mixture round-trip tests to PQmass/HmolarQmass/DmolarQmass"
```

---

## Task 11: REFPROP `calc_phase_molar_masses` override

**Files:**
- Modify: `src/Backends/REFPROP/REFPROPMixtureBackend.h`
- Modify: `src/Backends/REFPROP/REFPROPMixtureBackend.cpp`

- [ ] **Step 1: Declare override**

In `REFPROPMixtureBackend.h`, in the protected section near other `calc_*` overrides:

```cpp
    void calc_phase_molar_masses(double& MM_liquid, double& MM_vapor) override;
```

- [ ] **Step 2: Implement override**

In `REFPROPMixtureBackend.cpp`, near `calc_molar_mass` (around line 945):

```cpp
void REFPROPMixtureBackend::calc_phase_molar_masses(double& MM_l, double& MM_v) {
    if (_fluid_type == FLUID_TYPE_PURE || _fluid_type == FLUID_TYPE_PSEUDOPURE) {
        const double mm = molar_mass();
        MM_l = mm;
        MM_v = mm;
        return;
    }
    this->check_loaded_fluid();
    double mm_l_kg_per_kmol = NAN, mm_v_kg_per_kmol = NAN;
    WMOLdll(&(mole_fractions_liq[0]), &mm_l_kg_per_kmol);
    WMOLdll(&(mole_fractions_vap[0]), &mm_v_kg_per_kmol);
    MM_l = mm_l_kg_per_kmol / 1000.0;
    MM_v = mm_v_kg_per_kmol / 1000.0;
}
```

(`WMOLdll` returns kg/kmol; divide by 1000 to get kg/mol — matches the existing convention in `calc_molar_mass`.)

- [ ] **Step 3: Build to confirm REFPROP backend still links (this requires REFPROP DLL/SO available locally)**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/build
cmake --build . --target CoolProp -j 4 2>&1 | tail -10
```

Expected: clean build. (If REFPROP_LIBRARY_PATH is unset, the REFPROP backend may not link tests but the static lib should compile.)

- [ ] **Step 4: Commit**

```bash
git add src/Backends/REFPROP/REFPROPMixtureBackend.h src/Backends/REFPROP/REFPROPMixtureBackend.cpp
git commit -m "qmass: REFPROP calc_phase_molar_masses via WMOLdll on phase compositions"
```

---

## Task 12: REFPROP early-exit hook + native fast path for `QmassT` / `PQmass`

**Files:**
- Modify: `src/Backends/REFPROP/REFPROPMixtureBackend.h`
- Modify: `src/Backends/REFPROP/REFPROPMixtureBackend.cpp`

- [ ] **Step 1: Declare `update_Qmass_pair` override**

In `REFPROPMixtureBackend.h`:

```cpp
    void update_Qmass_pair(CoolProp::input_pairs pair, double v1, double v2) override;
```

- [ ] **Step 2: Implement override**

In `REFPROPMixtureBackend.cpp`, near the existing `update()`:

```cpp
void REFPROPMixtureBackend::update_Qmass_pair(CoolProp::input_pairs pair, double v1, double v2) {
    this->check_loaded_fluid();
    if (pair == CoolProp::QmassT_INPUTS || pair == CoolProp::PQmass_INPUTS) {
        const int kq = 2; // mass-basis quality
        int ierr = 0;
        char herr[255] = {0};
        double T_K = 0, p_kPa = 0, q = 0;
        double rho_mol_L = 0, rhoLmol_L = 0, rhoVmol_L = 0;
        double emol = 0, hmol = 0, smol = 0, cvmol = 0, cpmol = 0, w = 0;

        if (pair == CoolProp::QmassT_INPUTS) {
            T_K = v2;
            q   = v1;
            TQFLSHdll(&T_K, &q, &(mole_fractions[0]), &kq, &p_kPa,
                      &rho_mol_L, &rhoLmol_L, &rhoVmol_L,
                      &(mole_fractions_liq[0]), &(mole_fractions_vap[0]),
                      &emol, &hmol, &smol, &cvmol, &cpmol, &w,
                      &ierr, herr, errormessagelength);
        } else { // PQmass_INPUTS
            p_kPa = v1 * 0.001;
            q     = v2;
            PQFLSHdll(&p_kPa, &q, &(mole_fractions[0]), &kq, &T_K,
                      &rho_mol_L, &rhoLmol_L, &rhoVmol_L,
                      &(mole_fractions_liq[0]), &(mole_fractions_vap[0]),
                      &emol, &hmol, &smol, &cvmol, &cpmol, &w,
                      &ierr, herr, errormessagelength);
        }
        if (static_cast<int>(ierr) > get_config_int(REFPROP_ERROR_THRESHOLD)) {
            throw ValueError(format("Qmass-flash: %s", herr));
        }

        // Populate state. _Q (molar quality) recovered via Qmass<->Qmolar inversion
        // using returned phase compositions.
        _T = T_K;
        _p = p_kPa * 1000.0;
        _rhomolar = rho_mol_L * 1000.0;
        _rhoLmolar = rhoLmol_L * 1000.0;
        _rhoVmolar = rhoVmol_L * 1000.0;
        double MM_l = 0, MM_v = 0;
        calc_phase_molar_masses(MM_l, MM_v);
        _Q = Qmass_to_Qmolar(q, MM_l, MM_v);
        _Qmass = q;
        _phase = iphase_twophase;
        return;
    }
    // For the 6 remaining Qmass pairs REFPROP has no native kq flag; fall back
    // to the base secant-on-Qmolar.
    AbstractState::update_Qmass_pair(pair, v1, v2);
}
```

Note: `Qmass_to_Qmolar` is in the anonymous namespace of `AbstractState.cpp` and not visible here. Either move the inline definitions into a small internal header `src/Backends/qmass_conversions.h` and include it from both translation units, OR duplicate the 3-line conversion locally in `REFPROPMixtureBackend.cpp`. Pick **the small internal header**: cleaner, single source of truth, and keeps the formulas DRY.

Therefore, before doing the above:

(a) Move the two free functions out of the anonymous namespace in `src/AbstractState.cpp` and into a new file `include/qmass_conversions.h`:

```cpp
#pragma once
namespace CoolProp { namespace detail {
inline double Qmolar_to_Qmass(double Qmolar, double MM_l, double MM_v) {
    const double mv = Qmolar * MM_v;
    const double ml = (1.0 - Qmolar) * MM_l;
    return mv / (mv + ml);
}
inline double Qmass_to_Qmolar(double Qmass, double MM_l, double MM_v) {
    const double nv = Qmass / MM_v;
    const double nl = (1.0 - Qmass) / MM_l;
    return nv / (nv + nl);
}
}}
```

(b) Replace the anonymous-namespace versions in `src/AbstractState.cpp` with `using detail::Qmolar_to_Qmass; using detail::Qmass_to_Qmolar;` after `#include "qmass_conversions.h"`.

(c) `#include "qmass_conversions.h"` from `REFPROPMixtureBackend.cpp` and use `detail::Qmass_to_Qmolar` in the override.

- [ ] **Step 3: Wire REFPROP `update()` early-exit**

In `REFPROPMixtureBackend::update()`, add as the very first action (before `pre_update` / `mass_to_molar_inputs`):

```cpp
void REFPROPMixtureBackend::update(CoolProp::input_pairs input_pair, double value1, double value2) {
    if (CoolProp::is_Qmass_pair(input_pair) && !is_pure_or_pseudopure()) {
        update_Qmass_pair(input_pair, value1, value2);
        return;
    }
    // ... existing body ...
}
```

- [ ] **Step 4: Add REFPROP test (gated on REFPROP availability)**

Append to `src/Tests/CoolProp-Tests.cpp`:

```cpp
TEST_CASE("Qmass input: REFPROP R32+R125 native kq=2 fast path", "[Qmass][REFPROP]") {
    std::shared_ptr<CoolProp::AbstractState> AS;
    try {
        AS.reset(CoolProp::AbstractState::factory("REFPROP", "R32&R125"));
    } catch (...) {
        WARN("REFPROP not available; skipping");
        return;
    }
    AS->set_mole_fractions({0.5, 0.5});
    AS->update(CoolProp::QT_INPUTS, 0.4, 280.0);
    const double Qmass = AS->Qmass();
    const double P_ref = AS->p();
    const double Q_ref = AS->Q();

    auto AS2 = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("REFPROP", "R32&R125"));
    AS2->set_mole_fractions({0.5, 0.5});
    AS2->update(CoolProp::QmassT_INPUTS, Qmass, 280.0);
    CHECK(AS2->p()    == Approx(P_ref).epsilon(1e-8));
    CHECK(AS2->Q()    == Approx(Q_ref).epsilon(1e-8));
    CHECK(AS2->Qmass()== Approx(Qmass).epsilon(1e-12));
}
```

- [ ] **Step 5: Build and run (REFPROP available environments only)**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/build
cmake --build . --target CoolProp-Tests -j 4 2>&1 | tail -10
./src/Tests/CoolProp-Tests "[Qmass][REFPROP]"
```

Expected: PASS if REFPROP is configured; SKIPPED with WARN otherwise.

- [ ] **Step 6: Commit**

```bash
git add include/qmass_conversions.h \
        src/AbstractState.cpp \
        src/Backends/REFPROP/REFPROPMixtureBackend.h \
        src/Backends/REFPROP/REFPROPMixtureBackend.cpp \
        src/Tests/CoolProp-Tests.cpp
git commit -m "qmass: REFPROP native kq=2 path for QmassT/PQmass; shared conversion header"
```

---

## Task 13: Edge cases — Qmass in {0, 1}, out-of-range, single-phase

**Files:**
- Modify: `src/Tests/CoolProp-Tests.cpp`

- [ ] **Step 1: Add edge-case tests**

Append to `src/Tests/CoolProp-Tests.cpp`:

```cpp
TEST_CASE("Qmass edge cases: bubble/dew, out-of-range, single-phase", "[Qmass][edge]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "R32&R125"));
    AS->set_mole_fractions({0.5, 0.5});

    SECTION("Qmass = 0 (bubble) bypasses iteration") {
        AS->update(CoolProp::QmassT_INPUTS, 0.0, 280.0);
        CHECK(AS->Q() == Approx(0.0).epsilon(1e-12));
        CHECK(AS->Qmass() == Approx(0.0).epsilon(1e-12));
    }

    SECTION("Qmass = 1 (dew) bypasses iteration") {
        AS->update(CoolProp::QmassT_INPUTS, 1.0, 280.0);
        CHECK(AS->Q() == Approx(1.0).epsilon(1e-12));
        CHECK(AS->Qmass() == Approx(1.0).epsilon(1e-12));
    }

    SECTION("Qmass < 0 throws") {
        CHECK_THROWS_AS(AS->update(CoolProp::QmassT_INPUTS, -0.1, 280.0), CoolProp::ValueError);
    }

    SECTION("Qmass > 1 throws") {
        CHECK_THROWS_AS(AS->update(CoolProp::QmassT_INPUTS,  1.5, 280.0), CoolProp::ValueError);
    }

    SECTION("Qmass() in single-phase state throws") {
        AS->update(CoolProp::PT_INPUTS, 5e6, 350.0); // single phase
        CHECK_THROWS_AS(AS->Qmass(), CoolProp::ValueError);
    }
}
```

- [ ] **Step 2: Build and run**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/build
cmake --build . --target CoolProp-Tests -j 4 2>&1 | tail -10
./src/Tests/CoolProp-Tests "[Qmass][edge]"
```

Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add src/Tests/CoolProp-Tests.cpp
git commit -m "qmass: edge-case tests (boundaries, range, single-phase)"
```

---

## Task 14: PropsSI smoke test — `Qmass` as input and output

**Files:**
- Modify: `src/Tests/CoolProp-Tests.cpp`

- [ ] **Step 1: Add high-level interface smoke test**

Append to `src/Tests/CoolProp-Tests.cpp`:

```cpp
TEST_CASE("Qmass: PropsSI integration (output + input)", "[Qmass][PropsSI]") {
    SECTION("Qmass as output for pure Water == Q") {
        const double Q = 0.3, T = 350.0;
        const double Qmolar = CoolProp::PropsSI("Q",     "T", T, "Q", Q, "Water");
        const double Qmass  = CoolProp::PropsSI("Qmass", "T", T, "Q", Q, "Water");
        CHECK(Qmolar == Approx(0.3).epsilon(1e-12));
        CHECK(Qmass  == Approx(0.3).epsilon(1e-12));
    }
    SECTION("Qmass as input for pure Water round-trip") {
        const double T = 350.0;
        const double P_ref = CoolProp::PropsSI("P", "T", T, "Q",     0.3, "Water");
        const double P_via = CoolProp::PropsSI("P", "T", T, "Qmass", 0.3, "Water");
        CHECK(P_via == Approx(P_ref).epsilon(1e-12));
    }
}
```

- [ ] **Step 2: Build and run**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/build
cmake --build . --target CoolProp-Tests -j 4 2>&1 | tail -10
./src/Tests/CoolProp-Tests "[Qmass][PropsSI]"
```

Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add src/Tests/CoolProp-Tests.cpp
git commit -m "qmass: PropsSI smoke test"
```

---

## Task 15: Documentation

**Files:**
- Modify: docs index for input pairs / parameters (location TBD — search for `QT_INPUTS` in `Web/`)
- Add: short example in `Web/examples/`

- [ ] **Step 1: Locate the documentation page that lists input pairs**

```bash
grep -rln 'QT_INPUTS\|HmolarP_INPUTS' /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/Web/ 2>&1
```

- [ ] **Step 2: Add an entry per Qmass pair to that page (mirroring the existing molar entries) and add an `iQmass` row to the parameters table**

If the docs are auto-generated from C++ headers (look for a script), the changes already in `DataStructures.cpp` will propagate automatically. Verify by re-running the doc-generation script if one exists; otherwise edit the static doc page.

- [ ] **Step 3: Add a short example demonstrating Qmass usage**

In `Web/examples/`, add a small `.py` or `.rst` snippet:

```python
import CoolProp.CoolProp as CP

# Pure fluid: Qmass == Qmolar
P = CP.PropsSI("P", "T", 350, "Qmass", 0.3, "Water")

# Mixture: Qmass differs from Qmolar
T = CP.PropsSI("T", "P", 1.5e6, "Qmass", 0.3, "REFPROP::R32[0.5]&R125[0.5]")
print(T)
```

- [ ] **Step 4: Commit**

```bash
git add Web/ docs/
git commit -m "qmass: documentation and example for mass-basis quality input/output"
```

---

## Task 16: Final verification

**Files:** none (verification only)

- [ ] **Step 1: Run the full Qmass test suite**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass/build
./src/Tests/CoolProp-Tests "[Qmass]"
```

Expected: all sections PASS.

- [ ] **Step 2: Run the consistency suite to confirm no regressions on existing pairs**

```bash
./src/Tests/CoolProp-Tests "[consistency]"
```

Expected: all sections PASS.

- [ ] **Step 3: Build the Python wrapper if available and run a smoke test**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/Qmass
# (Use whatever build path the project uses; placeholder shown)
pip install -e wrappers/Python/  2>&1 | tail -10
python -c "import CoolProp.CoolProp as CP; \
print(CP.PropsSI('P','T',350,'Qmass',0.3,'Water'))"
```

Expected: a numeric pressure ~41 kPa.

- [ ] **Step 4: Spec coverage cross-check**

Open `docs/superpowers/specs/2026-04-30-qmass-support-design.md` and walk through Sections 1–8. Confirm each requirement maps to a task above. Note any gaps and address them with one-off commits.

---

## Self-Review

Done after the plan was drafted; gaps fixed inline:

- **Spec coverage:** all of §1 (enums, parameter), §2 (cache + accessor), §3 (free fns), §4 (`saturation_phase_molar_masses` → renamed to `calc_phase_molar_masses` for consistency with backend `calc_*` virtuals; per-backend overrides in HEOS and REFPROP), §5 (pure-fluid `mass_to_molar_inputs`), §6 (`update_Qmass_pair` default secant), §7 (PropsSI is automatic via enum + parameter table; smoke test in Task 14), §8 (tests) are mapped.
- **Placeholders:** none — every step contains the actual code, exact paths, exact commands, and expected output.
- **Type consistency:** signature `void calc_phase_molar_masses(double&, double&)` used identically in header, base impl, HEOS override, REFPROP override. `update_Qmass_pair(input_pairs, double, double)` used identically. Free fns live in `CoolProp::detail` namespace per Task 12 refactor.
- **Spec→plan rename:** spec used `saturation_phase_molar_masses`; plan uses `calc_phase_molar_masses` to match the existing `calc_*` convention. This is an improvement, documented here.
