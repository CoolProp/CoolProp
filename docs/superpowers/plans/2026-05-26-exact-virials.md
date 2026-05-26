# Phase 1 — Exact Virial Coefficients (dedicated zero-density path) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Beads:** CoolProp-dr3 · **Derivations:** CoolProp-6gf (`docs/superpowers/derivations/virial-axy.md`) · **GitHub:** #2991 (child #1676) · **Phase 2 (separate):** CoolProp-izv

**Goal:** Compute `Bvirial`/`Cvirial`/`dBvirial_dT`/`dCvirial_dT` from the exact analytic zero-density Taylor coefficients (no `delta=1e-12` fudge), for pure fluids **and mixtures**.

**Architecture:** A NEW `virial_coeffs(tau)` path on the residual container that sums per-term contributions to `c1 = d(ar)/d(delta)|_0`, `c2 = (1/2)d^2(ar)/d(delta^2)|_0`, and their `tau`-derivatives — using the closed forms in `virial-axy.md`. It does **not** touch the hot-path `all()` derivative loops (zero regression risk to flash/properties). Each term that hasn't yet got an exact form falls back to evaluating *that term alone* at `delta=1e-12` (the current method, scoped per term), so every fluid stays correct at every commit — virials are never broken mid-conversion.

**Why this is incremental-safe (not an "intermediate state"):** `Bvirial`/`Cvirial` already return correct values for all fluids today via the whole-residual `1e-12` evaluation. Converting one term to exact while others use the per-term `1e-12` fallback keeps every result correct — only accuracy improves. (Contrast Phase 2's `p(rho=0)`, which is binary `0`/`NaN` and must land all-or-nothing.)

**Tech Stack:** C++ (CoolProp core), Catch2, CMake, REFPROP (local), multicomplex (`mcx`).

---

## Mixtures (per user request)

`HelmholtzEOSMixtureBackend` builds the mixture residual `ar(tau,delta,x)` from composition-weighted component departures **plus** binary excess/departure terms (themselves GenExp-type) in the same `residual_helmholtz` container that `calc_alphar_deriv_nocache` reads. Because `virial_coeffs(tau)` operates on that same container at the mixture's reducing state, mixtures are handled by construction — but this MUST be verified against REFPROP (Task 6), not assumed. `B_mix(T,x)` and `C_mix(T,x)` are composition-dependent; the reducing density `rhomolar_reducing()` is already composition-dependent in the backend.

---

## File Structure

- `include/Helmholtz.h` — `struct VirialCoeffs {double c1, c2, dc1_dtau, dc2_dtau;}` (+ `operator+=`); declare `virial_coeffs(tau)` on `BaseHelmholtzContainer`/`ResidualHelmholtzContainer` and a per-term `add_virial_coeffs(tau, VirialCoeffs&)`.
- `src/Helmholtz.cpp` — `ResidualHelmholtzContainer::virial_coeffs` (sum + per-term `1e-12` fallback); exact `add_virial_coeffs` for `GenExp`, `GaoB` (the `1/delta` terms); the finite-at-0 terms (`NonAnalytic`, `SAFT`, cubic, `XiangDeiters`) can use the fallback or a direct `delta=0` eval.
- `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp` — rewrite `calc_Bvirial`/`calc_Cvirial`/`calc_dBvirial_dT`/`calc_dCvirial_dT` (`:1634-1649`) on `virial_coeffs`.
- `src/Tests/CoolProp-tests.cpp` — `[virial]`-tagged tests (pure + mixture + REFPROP + dT).

---

### Task 1: Characterization test — exact path must match the current method

**Files:** Test: `src/Tests/CoolProp-tests.cpp`

Lock in current behavior as a regression oracle before changing anything.

- [ ] **Step 1: Capture current Bvirial/Cvirial for reference fluids**

```cpp
TEST_CASE("Virials match the legacy 1e-12 method (characterization)", "[virial]") {
    struct Row { const char* fluid; double T; };
    for (auto r : {Row{"Nitrogen",300}, Row{"CarbonDioxide",300}, Row{"Methane",250}, Row{"Water",500}}) {
        CAPTURE(r.fluid, r.T);
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", r.fluid));
        AS->update(CoolProp::DmolarT_INPUTS, 1e-10, r.T);
        // Values recorded from current build (fill in from Step-2 run):
        // CHECK(AS->Bvirial() == Approx(<recorded>).epsilon(1e-10));
    }
}
```

- [ ] **Step 2: Run current build, record the values, bake them into the CHECKs**

Run: `cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[virial]" -s`
Action: read the printed `Bvirial()`/`Cvirial()`, substitute as the `Approx(...)` targets, re-run to confirm green.

- [ ] **Step 3: Commit**

```bash
git add src/Tests/CoolProp-tests.cpp
git commit -m "test(virial): characterize current Bvirial/Cvirial as regression oracle (CoolProp-dr3)"
```

---

### Task 2: `VirialCoeffs` type + container path with per-term `1e-12` fallback

**Files:** Modify `include/Helmholtz.h`; `src/Helmholtz.cpp`

- [ ] **Step 1: Add the type and per-term hook**

```cpp
struct VirialCoeffs {
    double c1 = 0, c2 = 0, dc1_dtau = 0, dc2_dtau = 0;
    VirialCoeffs& operator+=(const VirialCoeffs& o) {
        c1 += o.c1; c2 += o.c2; dc1_dtau += o.dc1_dtau; dc2_dtau += o.dc2_dtau; return *this;
    }
};
```

Each residual term declares `virtual bool add_virial_coeffs(CoolPropDbl tau, VirialCoeffs& vc) const;` returning `false` if it has no exact form yet (→ container uses fallback for that term).

- [ ] **Step 2: Container sums exact-or-fallback per term**

```cpp
VirialCoeffs ResidualHelmholtzContainer::virial_coeffs(const CoolPropDbl tau) {
    VirialCoeffs vc;
    auto add = [&](BaseHelmholtzTerm& term) {
        VirialCoeffs t;
        if (!term.add_virial_coeffs(tau, t)) {       // fallback: this term alone at delta=1e-12
            HelmholtzDerivatives d; term.all(tau, 1e-12, d);
            t.c1 = d.dalphar_ddelta;                 // -> d(ar)/d(delta)|_0
            t.c2 = 0.5 * d.d2alphar_ddelta2;
            t.dc1_dtau = d.d2alphar_ddelta_dtau;
            t.dc2_dtau = 0.5 * d.d3alphar_ddelta2_dtau;
        }
        vc += t;
    };
    add(GenExp); add(NonAnalytic); add(SAFT); add(cubic); add(XiangDeiters); add(GaoB);
    return vc;
}
```

> Default base `add_virial_coeffs` returns `false`. So at this point ALL terms use the fallback ⇒ result identical to the legacy whole-residual `1e-12` method (Task 1 oracle stays green).

- [ ] **Step 3: Build** — Run: `cmake --build build_catch --target CatchTestRunner -j8`; Expected: compiles.

- [ ] **Step 4: Commit** — `git commit -m "feat(virial): VirialCoeffs + container path with per-term 1e-12 fallback (CoolProp-dr3)"`

---

### Task 3: Backend uses `virial_coeffs` (still fallback everywhere)

**Files:** Modify `HelmholtzEOSMixtureBackend.cpp:1634-1649`

- [ ] **Step 1: Rewrite the four accessors**

```cpp
CoolPropDbl HelmholtzEOSMixtureBackend::calc_Bvirial() {
    return residual_helmholtz->virial_coeffs(_tau).c1 / rhomolar_reducing();
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_Cvirial() {
    return 2.0 * residual_helmholtz->virial_coeffs(_tau).c2 / pow(rhomolar_reducing(), 2);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_dBvirial_dT() {
    SimpleState r = get_reducing_state(); CoolPropDbl dtau_dT = -r.T / pow(_T, 2);
    return residual_helmholtz->virial_coeffs(_tau).dc1_dtau * dtau_dT / r.rhomolar;
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_dCvirial_dT() {
    SimpleState r = get_reducing_state(); CoolPropDbl dtau_dT = -r.T / pow(_T, 2);
    return 2.0 * residual_helmholtz->virial_coeffs(_tau).dc2_dtau * dtau_dT / pow(r.rhomolar, 2);
}
```

> Confirm `residual_helmholtz->...` mixture plumbing matches `calc_alphar_deriv_nocache` (`:3414`) — it reads the same composition-weighted container, so mixtures route through unchanged.

- [ ] **Step 2: Run characterization oracle** — Run: `./build_catch/CatchTestRunner "[virial]"`; Expected: PASS (fallback reproduces legacy values).

- [ ] **Step 3: Commit** — `git commit -m "feat(virial): route Bvirial/Cvirial through virial_coeffs path (CoolProp-dr3)"`

---

### Task 4: Exact `add_virial_coeffs` for GenExp

**Files:** Modify `src/Helmholtz.cpp` (`ResidualHelmholtzGeneralizedExponential`)

- [ ] **Step 1: Failing test — exact GenExp matches fallback to tight tol**

```cpp
TEST_CASE("GenExp exact virials match 1e-12 path", "[virial]") {
    for (const char* f : {"Methane","Nitrogen","Propane"}) {  // pure GenExp (no nonanalytic)
        CAPTURE(f);
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", f));
        AS->update(CoolProp::DmolarT_INPUTS, 1e-10, 300.0);
        // exact vs recorded legacy value (from Task 1); should now agree to ~1e-12 not ~1e-6
        // CHECK(AS->Bvirial() == Approx(<recorded>).epsilon(1e-12));
    }
}
```

- [ ] **Step 2: Run to verify the *tighter* tolerance fails on fallback** — Expected: the `1e-12`-epsilon CHECK fails (fallback only ~`1e-6`-accurate), motivating the exact form.

- [ ] **Step 3: Implement** per `virial-axy.md` §3. Loop elements; with `u0 = u(tau, delta=0)`:

```cpp
bool ResidualHelmholtzGeneralizedExponential::add_virial_coeffs(CoolPropDbl tau, VirialCoeffs& vc) const {
    for (auto& el : elements) {
        // u0 = -omega*tau^m - beta1*(tau-gam1) - beta2*(tau-gam2)^2 + eta1*eps1 - eta2*eps2^2
        CoolPropDbl u0 = /* ... */, E = el.n * pow(tau, el.t) * exp(u0);
        if (el.d == 1) {
            vc.c1 += E;
            // c2: (-c if l==1) + 2*eta2*eps2 - eta1, all * E
            CoolPropDbl k = (el.l_double == 1 ? -el.c : 0.0) + 2*el.eta2*el.epsilon2 - el.eta1;
            vc.c2 += k * E;
            // dc1_dtau, dc2_dtau from virial-axy.md (du0/dtau factor)
        } else if (el.d == 2) {
            vc.c2 += E;            // dc2_dtau likewise
        } // d>=3: no contribution
    }
    return true;
}
```

> Transcribe `dc1_dtau`/`dc2_dtau` verbatim from `virial-axy.md` §3 (the `du0/dtau = -m*omega*tau^(m-1) - beta1 - 2*beta2*(tau-gam2)` chain). Guard the optional fields (`eta1`/`eta2`/`beta*` may be NaN/unset — mirror the `*_in_u` flags used in `all()`).

- [ ] **Step 4: Run** — Expected: GenExp fluids now match to `1e-12`; characterization oracle still green.

- [ ] **Step 5: Commit** — `git commit -m "feat(virial): exact GenExp virial coefficients (CoolProp-dr3)"`

---

### Task 5: Exact `add_virial_coeffs` for GaoB (the other `1/delta` term)

**Files:** Modify `src/Helmholtz.cpp` (`ResidualHelmholtzGaoB`)

- [ ] **Step 1: Failing test** — a GaoB-using fluid, tighten to `1e-12`.
- [ ] **Step 2: Run to verify fallback-tolerance fails.**
- [ ] **Step 3: Implement** per `virial-axy.md` §5 (`w = eta*eps^2 + 1/(b+beta*(tau-gam)^2)`; `d=1`: `c1=n*tau^t*exp(w)`, `c2=-2*eps*eta*n*tau^t*exp(w)`; `d=2`: `c2=n*tau^t*exp(w)`).
- [ ] **Step 4: Run** — Expected: PASS.
- [ ] **Step 5: Commit** — `git commit -m "feat(virial): exact GaoB virial coefficients (CoolProp-dr3)"`

> NonAnalytic/SAFT/cubic/XiangDeiters have no `1/delta` in source; their per-term `1e-12` fallback is already exact to rounding. Leaving them on the fallback is acceptable — optionally add a direct `delta=0` eval later. (XiangDeiters is GenExp-based; once GenExp `add_virial_coeffs` exists, route XiangDeiters' sub-terms through it for full exactness — a follow-up note, not required for correctness.)

---

### Task 6: REFPROP cross-check — pure AND mixtures

**Files:** Test: `src/Tests/CoolProp-tests.cpp`

- [ ] **Step 1: Pure-fluid REFPROP check**

```cpp
TEST_CASE("Bvirial/Cvirial vs REFPROP (pure)", "[virial][refprop]") {
    for (const char* f : {"Nitrogen","CarbonDioxide","Methane","Propane","Water"}) {
        CAPTURE(f);
        shared_ptr<CoolProp::AbstractState> H(CoolProp::AbstractState::factory("HEOS", f));
        shared_ptr<CoolProp::AbstractState> R(CoolProp::AbstractState::factory("REFPROP", f));
        for (double T : {250.,300.,400.}) {
            H->update(CoolProp::DmolarT_INPUTS, 1e-10, T);
            R->update(CoolProp::DmolarT_INPUTS, 1e-10, T);
            CHECK(H->Bvirial() == Approx(R->Bvirial()).epsilon(1e-6));
            CHECK(H->Cvirial() == Approx(R->Cvirial()).epsilon(1e-5));
        }
    }
}
```

- [ ] **Step 2: MIXTURE REFPROP check (per user request)**

```cpp
TEST_CASE("Bvirial/Cvirial vs REFPROP (mixtures)", "[virial][refprop]") {
    struct M { const char* names; std::vector<double> z; };
    for (auto mix : { M{"Methane&Ethane", {0.7,0.3}},
                      M{"Nitrogen&CarbonDioxide", {0.5,0.5}},
                      M{"Methane&Nitrogen&CarbonDioxide", {0.8,0.1,0.1}} }) {
        CAPTURE(mix.names);
        shared_ptr<CoolProp::AbstractState> H(CoolProp::AbstractState::factory("HEOS", mix.names));
        shared_ptr<CoolProp::AbstractState> R(CoolProp::AbstractState::factory("REFPROP", mix.names));
        H->set_mole_fractions(mix.z); R->set_mole_fractions(mix.z);
        for (double T : {250.,300.,400.}) {
            H->update(CoolProp::DmolarT_INPUTS, 1e-10, T);
            R->update(CoolProp::DmolarT_INPUTS, 1e-10, T);
            CHECK(H->Bvirial() == Approx(R->Bvirial()).epsilon(1e-6));
            CHECK(H->Cvirial() == Approx(R->Cvirial()).epsilon(1e-5));
        }
    }
}
```

> If a `&`-name isn't recognized, build the mixture via the established pair pattern (component names + `set_mole_fractions`). Mixture `C` may need a looser epsilon if REFPROP's `VIRCdll` mixing differs; tighten/loosen empirically and document.

- [ ] **Step 3: Run** — Run: `./build_catch/CatchTestRunner "[virial]"`; Expected: PASS (REFPROP runs locally).

- [ ] **Step 4: Commit** — `git commit -m "test(virial): REFPROP cross-check, pure + mixtures (CoolProp-dr3)"`

---

### Task 7: Temperature-derivative check + multicomplex

**Files:** Test: `src/Tests/CoolProp-tests.cpp`

- [ ] **Step 1: dB/dT, dC/dT vs central finite difference of B(T), C(T)**

```cpp
TEST_CASE("dBvirial_dT/dCvirial_dT match finite difference", "[virial]") {
    for (const char* f : {"Methane","CarbonDioxide","Water"}) {
        CAPTURE(f);
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", f));
        double T = 300, h = 1e-3;
        auto B = [&](double TT){ AS->update(CoolProp::DmolarT_INPUTS, 1e-10, TT); return AS->Bvirial(); };
        double fd = (B(T+h) - B(T-h)) / (2*h);
        AS->update(CoolProp::DmolarT_INPUTS, 1e-10, T);
        CHECK(AS->dBvirial_dT() == Approx(fd).epsilon(1e-6));
    }
}
```

- [ ] **Step 2: (optional) multicomplex absolute check** of `c1` for a GenExp fluid via `one_mcx`.
- [ ] **Step 3: Run** — Expected: PASS.
- [ ] **Step 4: Commit** — `git commit -m "test(virial): dB/dT, dC/dT finite-difference + mcx checks (CoolProp-dr3)"`

---

### Task 8: Pre-PR gate

- [ ] `./dev/ci/preflight.sh` passes (`[virial]`/`[refprop]` run locally).
- [ ] Invoke `superpowers:code-reviewer` on the diff (CLAUDE.md mandatory). Focus: optional-field guards (`eta`/`beta` NaN), `d==1`/`l==1` selection, mixture container plumbing, that `virial_coeffs` is NOT on the flash hot path, fallback removed only where exact is proven.
- [ ] Restore `.beads/issues.jsonl` per CLAUDE.md, `git push`, `gh pr create` referencing #2991.

---

## Self-Review

- **No 1e-12 fudge:** Tasks 4–5 (GenExp/GaoB exact); finite-at-0 terms exact-to-rounding via fallback (Task 5 note).
- **Pure + mixtures:** Task 6 (both, vs REFPROP) — addresses the explicit mixture request.
- **dB/dT, dC/dT:** Task 7.
- **Incremental safety:** Task 2 fallback ⇒ legacy values preserved at every commit; characterization oracle (Task 1) guards it.
- **No hot-path regression:** `virial_coeffs` is a separate method; `all()` loops untouched.
- **Closed forms:** all transcribed from `docs/superpowers/derivations/virial-axy.md` §3/§5; `c1/c2` tables there are the source of truth.
- **Verify-in-tree:** container/mixture plumbing signature vs `calc_alphar_deriv_nocache`; element optional-field names (`l_double`, `eta1`, `epsilon2`, …) vs the GenExp element struct.
