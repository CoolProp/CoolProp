# Qmass support across CoolProp backends

**Date:** 2026-04-30
**Status:** Proposed
**Supersedes:** PR #2221, PR #2838 (REFPROP-only attempts)
**Branch strategy:** Fresh branch from `master` (no cherry-pick from PR 2838's `mass-quality-methods`; reference only).

## Problem

CoolProp's vapor quality `Q` is implicitly molar. For mixtures, mass-based vapor quality `Qmass` differs from molar `Qmolar` and is the quantity engineers typically want to specify and read. Today users must convert by hand, which is error-prone and not symmetric with the rest of the API (where `Smass` ↔ `Smolar`, `Hmass` ↔ `Hmolar`, etc. are first-class).

PR 2838 added `QmassT_INPUTS` and `PQmass_INPUTS` for REFPROP only, with conversion logic embedded in the REFPROP backend. Reviewer feedback: cover all 8 Q-input pairs, move conversion to `AbstractState`, and use the existing cache-clearing pattern.

## Goals

1. Full enum coverage: a `Qmass` variant for every existing `Q`-input pair (8 pairs total).
2. A new keyed parameter `iQmass` so `PropsSI("Qmass", …)` and `update(iQmass, …, iT, …)` work.
3. Conversion routines as **free functions** in `AbstractState.cpp` — pure math, no state.
4. Backend-agnostic: works for both HEOS and REFPROP without per-backend overrides for the conversion itself.
5. **Mixture support from day one** — including the iterative Qmass→Qmolar solve required for input pairs.
6. Existing API surface unchanged: `QT_INPUTS`, `PQ_INPUTS`, etc. retain molar-quality semantics. No renames, no aliases.

## Non-goals

- Renaming `QT_INPUTS` → `QmolarT_INPUTS` with backward-compat alias (separate cosmetic PR if ever).
- REFPROP `mass_fractions` member array (mass fractions derived on demand).
- Surface-tension `STNdll` improvements (PR 2838 scope creep — drop).
- `update_with_guesses` Qmass support (follow-up if needed; the molar form already accepts guesses).

## Design

### 1. Enum and parameter additions

**`include/DataStructures.h` — input_pairs enum.** Add 8 new entries, each immediately following its molar sibling:

| Existing | New (added below it) |
|---|---|
| `QT_INPUTS` | `QmassT_INPUTS` |
| `PQ_INPUTS` | `PQmass_INPUTS` |
| `QSmolar_INPUTS` | `QmassSmolar_INPUTS` |
| `QSmass_INPUTS` | `QmassSmass_INPUTS` |
| `HmolarQ_INPUTS` | `HmolarQmass_INPUTS` |
| `HmassQ_INPUTS` | `HmassQmass_INPUTS` |
| `DmolarQ_INPUTS` | `DmolarQmass_INPUTS` |
| `DmassQ_INPUTS` | `DmassQmass_INPUTS` |

Update the `input_pairs` ↔ string map in `src/DataStructures.cpp` accordingly (long name + short description).

**`include/DataStructures.h` — parameters enum.** Add `iQmass` immediately after `iQ` (line 91). Update the parameter table in `src/DataStructures.cpp`: short name `"Qmass"`, units empty (dimensionless), description `"Mass-based vapor quality"`, IO type `OUTPUT_TYPE` (or both, parallel to `iSmass`).

**`include/DataStructures.h` — `generate_update_pair`.** Add 8 new `match_pair` branches keyed on `iQmass`, returning the corresponding `Qmass_INPUTS` pair. Order: place each new branch adjacent to its molar sibling for readability.

**`src/AbstractState.cpp` — `keyed_output` dispatch.** Add `case iQmass: return Qmass();` so `PropsSI("Qmass", …)` and the keyed output table resolve correctly.

### 2. Cached output and accessor

**`include/AbstractState.h`:**

```cpp
// next to other "two-phase" cached values, around line 149
CAE _Qmass = cache.next();
```

Use the `CAE = cache.next()` pattern (auto-cleared by `cache.clear()` — no manual entry in `AbstractState::clear()`).

**`include/AbstractState.h` (public API):**

```cpp
double Qmass(void);              // wraps calc_Qmass with cache
```

**`src/AbstractState.cpp`:**

```cpp
double AbstractState::Qmass(void) {
    if (!_Qmass) _Qmass = calc_Qmass();
    return _Qmass;
}

CoolPropDbl AbstractState::calc_Qmass(void) {
    if (!ValidNumber(_Q) || _Q < 0 || _Q > 1)
        throw ValueError("Qmass requires a two-phase state (0 <= Q <= 1)");
    double MM_l, MM_v;
    saturation_phase_molar_masses(MM_l, MM_v);
    return Qmolar_to_Qmass(_Q, MM_l, MM_v);
}
```

`calc_Qmass` is **non-virtual** — single base implementation, identical across backends.

### 3. Free-function conversions

**`src/AbstractState.cpp` (anonymous namespace at top of file):**

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

} // namespace
```

Pure math, no state, no AbstractState dependency. Same file as the only callers.

### 4. Phase molar-mass helper (private method)

**`include/AbstractState.h` (protected):**

```cpp
/// Populate MM_liquid and MM_vapor (kg/mol) for the current saturated state.
/// Pure/pseudopure: both equal molar_mass(). Mixture: derived from
/// mole_fractions_liquid()/mole_fractions_vapor() and component molar masses.
void saturation_phase_molar_masses(double& MM_liquid, double& MM_vapor);
```

Implementation in `src/AbstractState.cpp`:

```cpp
void AbstractState::saturation_phase_molar_masses(double& MM_l, double& MM_v) {
    if (_fluid_type == FLUID_TYPE_PURE || _fluid_type == FLUID_TYPE_PSEUDOPURE) {
        const double mm = molar_mass();
        MM_l = mm;
        MM_v = mm;
        return;
    }
    const auto x = mole_fractions_liquid();    // virtual, backend-provided
    const auto y = mole_fractions_vapor();
    const auto mm_components = get_mole_fractions_ref();  // not used; need MM per comp
    // Both HEOS and REFPROP expose component molar masses via:
    //   - HEOS: components[i]->EOS()->molar_mass
    //   - REFPROP: WMOLdll on a single-component composition vector
    // Provide a virtual hook calc_component_molar_masses() if needed,
    // OR compute MM_l = molar_mass_of(x), MM_v = molar_mass_of(y) where
    // molar_mass_of(z) is a small helper that does:
    //     save current mole fractions; set to z; call molar_mass(); restore.
    // (See implementation note below.)
    MM_l = molar_mass_of_composition(x);
    MM_v = molar_mass_of_composition(y);
}
```

**Implementation note on per-composition molar mass.** Both backends already expose `calc_molar_mass()` driven by the *current* `mole_fractions` vector. The cleanest route is a tiny helper `molar_mass_of_composition(const std::vector<CoolPropDbl>&)` that:

1. saves the current `mole_fractions` and the cached `_molar_mass`,
2. sets `mole_fractions = z`, clears `_molar_mass`,
3. calls `calc_molar_mass()`,
4. restores both.

This keeps the change confined to `AbstractState` (the helper lives in the base class, reading/writing `mole_fractions` and calling the backend's `calc_molar_mass()` virtual). No new virtuals.

Alternative if save/restore is fragile: add a tiny non-virtual helper that, for each component, computes MM_i = sum-formula via existing component-table lookup. This has backend-specific shape (HEOS has `components[i]->EOS()->molar_mass`; REFPROP has `WMOLdll` per pure call). If the save/restore pattern proves problematic in testing, fall back to a `calc_component_molar_masses()` virtual returning a `vector<double>`. **Default plan: save/restore. Decide during implementation.**

### 5. Input-pair conversion: `mass_to_molar_inputs`

Extend the existing switch in `src/AbstractState.cpp:200`. **Two paths**:

**5a. Pure / pseudopure fluid.** `Qmass = Qmolar` exactly (MM_l = MM_v). Add the 8 cases to the existing switch alongside `SmassT`, `HmassP`, etc. Each case rewrites the input pair to its molar sibling without altering values.

**5b. Mixture.** Algebra is not closed — MM_l and MM_v depend on the equilibrium that the flash itself produces. Iteration required.

For mixture handling, do **not** put iteration inside `mass_to_molar_inputs` (it should remain a pure transform). Instead, add a virtual entry point on the base:

```cpp
// AbstractState.h (protected, virtual with default implementation)
virtual void update_Qmass_pair(input_pairs pair, double v1, double v2);
```

Default implementation in `AbstractState.cpp` does the secant iteration (Section 6). Backends gain one early-exit line at the top of their `update()`:

```cpp
void HelmholtzEOSMixtureBackend::update(input_pairs pair, double v1, double v2) {
    if (is_Qmass_pair(pair) && !is_pure_or_pseudopure()) {
        update_Qmass_pair(pair, v1, v2);
        return;
    }
    // ... existing body ...
}
```

Same one-liner added to `REFPROPMixtureBackend::update`. Pure/pseudopure flows through `mass_to_molar_inputs` as before.

`is_Qmass_pair(input_pairs)` is a free helper returning true for the 8 new pairs. `is_pure_or_pseudopure()` is a small inline (`_fluid_type == FLUID_TYPE_PURE || _fluid_type == FLUID_TYPE_PSEUDOPURE`).

**REFPROP native fast path.** REFPROP's `TQFLSHdll` and `PQFLSHdll` accept a `kq` flag (`kq = 2` → mass-basis quality input). Comment in current code at `src/Backends/REFPROP/REFPROPMixtureBackend.cpp:1788-1791` confirms this. Override `update_Qmass_pair` in `REFPROPMixtureBackend`:

```cpp
void REFPROPMixtureBackend::update_Qmass_pair(input_pairs pair, double v1, double v2) {
    if (pair == QmassT_INPUTS || pair == PQmass_INPUTS) {
        // Direct call to TQFLSHdll/PQFLSHdll with kq=2.
        // After the call, populate _Q via Qmass_to_Qmolar using
        // the returned mole_fractions_liq/vap and component MMs;
        // cache _Qmass = target directly.
        return;
    }
    // 6 remaining Qmass pairs: defer to base-class iteration.
    AbstractState::update_Qmass_pair(pair, v1, v2);
}
```

The two pairs users actually reach for (`QmassT`, `PQmass`) cost the same as their molar siblings in REFPROP — REFPROP swallows the iteration internally. The other 6 fall through to the base secant.

For HEOS, no override: all 8 Qmass pairs use the base-class iteration.

### 6. Mixture iteration: `update_Qmass_pair` (default base implementation)

**Algorithm** (in `src/AbstractState.cpp`):

```
input: pair (a Qmass_INPUTS variant), v1, v2
identify the Qmass slot (one of v1/v2) and the partner property + slot
identify the corresponding molar pair (e.g. QmassT_INPUTS → QT_INPUTS)
target_Qmass = v_qmass_slot

f(Qmolar):
    call this->update(molar_pair, ...) with Qmolar in the Q slot, partner unchanged
    saturation_phase_molar_masses(MM_l, MM_v)
    return Qmolar_to_Qmass(Qmolar, MM_l, MM_v) - target_Qmass

solve f(Qmolar) = 0 on Qmolar ∈ [eps, 1-eps] using secant / Brent
    initial bracket: [target_Qmass, target_Qmass + small] for secant,
                     or [eps, 1-eps] for bracketed Brent
on convergence: state is at the correct equilibrium
set _Qmass = target_Qmass directly in cache
```

**Solver choice:** secant with safety fallback to bisection. Reuse the existing 1D solver in `src/Solvers.cpp` (`Brent`, `Secant`) — already used throughout CoolProp for similar problems.

**Convergence:** monotonic relationship (Qmass strictly increasing in Qmolar for valid two-phase states), so robust convergence in 5–8 iterations to typical 1e-10 tolerance. Exit conditions: `|f| < 1e-10` on Qmass residual, max 30 iterations.

**Edge cases:**
- `target_Qmass == 0`: skip iteration, set Qmolar = 0, recurse to molar pair directly.
- `target_Qmass == 1`: similarly Qmolar = 1.
- `target_Qmass < 0` or `> 1`: throw `ValueError`.
- Secant fails to bracket: fall back to bisection on `[1e-12, 1 - 1e-12]`.

**Cost:** mixture Qmass input is ~5–8× a normal mixture flash. Acceptable; alternative would be co-solving Qmolar and the partner state simultaneously, which is much more invasive.

**Cache hygiene:** each iteration calls `update()` which calls `clear()`, so all caches are fresh on each pass. After the final converged call, `_Q` holds the molar quality; we set `_Qmass` directly so the user can query both without recomputation.

### 7. PropsSI / high-level interface

The high-level layer (`PropsSI`, `PropsSImulti`) routes through `generate_update_pair` → backend `update`. Adding the 8 enum entries plus `iQmass` (registered in the parameter table) is sufficient. No additional plumbing.

`PropsSI("Qmass", "T", 300, "Q", 0.3, "Water")`: `iT` and `iQ` (molar) form `QT_INPUTS`, the flash runs, and the requested output `Qmass` is computed via the new accessor.

`PropsSI("P", "T", 300, "Qmass", 0.3, "Water")`: `iT` and `iQmass` form `QmassT_INPUTS`, dispatched as described above.

### 8. Tests

Add to `src/Tests/Tests.cpp` (or a new `src/Tests/QmassTests.cpp`).

**Pure fluid sanity.** Water, R134a, Propane: for several `(Q, T)` and `(P, Q)` points, confirm:
- `Qmass()` from a `QT_INPUTS` flash equals `Q()` (both 0–1, identical).
- `update(QmassT_INPUTS, 0.3, 300)` reproduces the same state as `update(QT_INPUTS, 0.3, 300)` to `1e-12`.
- All 8 new pairs round-trip on a pure fluid: from a known molar-pair state, read the partner property, re-create via the Qmass-pair, get identical state.

**Mixture: HEOS R32+R125 (50/50 mass).**
- `update(QT_INPUTS, 0.5, 280)` → record `Qmass()`, `MM_l`, `MM_v`, `P`.
- `update(QmassT_INPUTS, recorded_Qmass, 280)` → expect identical `Q()`, `P`, `MM_l`, `MM_v` to `1e-9`.
- Repeat for `PQmass_INPUTS`, `HmolarQmass_INPUTS`, `DmolarQmass_INPUTS`.

**Mixture: REFPROP same fluid pair.** Same round-trip, same tolerances. Skipped if REFPROP unavailable in CI.

**Boundary tests.**
- `Qmass = 0` and `Qmass = 1` (bubble/dew) — bypass iteration, exact match.
- `Qmass = 1e-6` and `Qmass = 1 - 1e-6` — full iteration, tight tolerance.
- `Qmass = -0.1` → `ValueError`. `Qmass = 1.5` → `ValueError`.

**Failure modes.**
- Single-phase state (e.g. supercritical) → `Qmass()` accessor throws `ValueError("requires two-phase state")`.
- `update(QmassT_INPUTS, 0.3, 700)` for Water (T above critical) → propagates the underlying flash error.

### 9. Documentation

- `Web/coolprop/wrappers/Python/index.rst` (or wherever input-pair docs live): add the 8 new pairs and `Qmass` keyed parameter to the table.
- A short example in `Web/examples/`: pure water and a binary mixture round-trip.
- Note in changelog / release notes.

## Architectural summary

```
+------------------------------------------------------------------+
| AbstractState (base)                                             |
|                                                                  |
|  Qmass()                       <- public accessor, cached        |
|    |                                                             |
|    v                                                             |
|  calc_Qmass()                  <- non-virtual                    |
|    |                                                             |
|    +-> saturation_phase_molar_masses(MM_l, MM_v)                 |
|    |     pure: molar_mass()                                      |
|    |     mixture: molar_mass_of_composition(x), (y)              |
|    |                                                             |
|    +-> Qmolar_to_Qmass(_Q, MM_l, MM_v)  <- free fn               |
|                                                                  |
|  mass_to_molar_inputs(pair, v1, v2)                              |
|    + 8 pure-fluid Qmass cases (trivial Qmass = Qmolar)           |
|                                                                  |
|  virtual update_Qmass_pair(pair, v1, v2)                         |
|    default: secant on Qmolar; each iter -> this->update(molar)   |
+------------------------------------------------------------------+
            ^                                       ^
            |                                       |
+-----------+-----------+              +------------+----------------+
| HEOSMixtureBackend     |              | REFPROPMixtureBackend       |
|  update():              |              |  update():                   |
|    if Qmass+mixture     |              |    if Qmass+mixture          |
|      -> update_Qmass    |              |      -> update_Qmass         |
|    else                 |              |    else                      |
|      mass_to_molar      |              |      mass_to_molar           |
|      flash as today     |              |      flash as today          |
|                         |              |                              |
| (no override:           |              | override update_Qmass_pair:  |
|  uses base iteration    |              |   QmassT/PQmass -> TQ/PQFLSH |
|  for all 8 pairs)       |              |     with kq=2 (no iteration) |
|                         |              |   other 6 pairs -> base iter |
+-------------------------+              +------------------------------+
```

## Risk and rollback

- **Risk:** mixture iteration may not converge for some pathological inputs (very narrow two-phase region, near critical point). Mitigation: bracketed Brent fallback; unit tests near critical points.
- **Risk:** `molar_mass_of_composition` save/restore approach interacts badly with cached state in subtle ways. Mitigation: discover during testing; fallback design (component-MM virtual) noted above.
- **Rollback:** All changes are additive — existing enum values, parameter ids, and `update()` signatures unchanged. A revert removes Qmass without breaking existing callers.

## Out of scope / follow-ups

- `update_with_guesses` Qmass variants.
- High-throughput tabular backends (BICUBIC, TTSE) — Qmass mappings can be added when needed; for now, they fall through the standard flash path.
- Renaming `QT_INPUTS` to `QmolarT_INPUTS` with alias.
