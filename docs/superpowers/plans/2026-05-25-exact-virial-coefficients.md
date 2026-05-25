# Exact Virial Coefficients + Zero-Density-Safe Axy Derivatives Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Beads issue:** CoolProp-izv  **GitHub:** #2991 (child of #1676)

**Goal:** Evaluate virial coefficients *at* zero density (no `delta=1e-12` fudge) and make pressure return exactly `0` at `rho=0` **for every fluid in the database**, by computing the well-behaved `Axy` derivative form analytically in every Helmholtz term and exposing it as `get_Ar(itau,idelta)` / `get_Aig(itau,idelta)`.

**No intermediate state (hard requirement):** the feature must hold for *all* fluids at merge — not "works for Methane, throws for a SAFT/GaoB fluid." Therefore every residual and ideal term type is converted in this single landing; there is **no phasing and no deferral to follow-up issues**. The merge gate is an all-fluids sweep.

**Architecture:** Reduced Helmholtz derivatives are hand-coded per term and stored raw (`dalphar_ddelta`, …). Terms carrying a `1/delta` factor (GenExp's `one_over_delta` at `Helmholtz.cpp:267`/`:164`; GaoB's explicit `/delta` at `:735`; ideal lead `1/delta` at `:1148`) are `Inf`/`NaN` at exactly `delta=0`. We add, in **every** term's `all()`, a parallel analytic accumulation of the dimensionless group `Axy = tau^x * delta^y * d^{x+y}a/(dtau^x ddelta^y)` computed without ever dividing by `delta`, so `Axy` is finite at `delta=0`. Each term's hand-derived `Axy` is cross-checked against its existing multicomplex `one_mcx` autodiff (the correctness oracle). Pressure and virials are then rebuilt on `Ar(0,1)`/`Ar(0,2)` and the analytic zero-density limit.

**Tech Stack:** C++ (CoolProp core), Catch2 (`CatchTestRunner`), CMake, REFPROP (local cross-check), multicomplex (`mcx`, vendored — `one_mcx` exists for GenExp/NonAnalytic/GaoB/XiangDeiters).

---

## Key Math (read before starting)

Residual reduced Helmholtz expanded in `delta` at fixed `tau`:

```
alphar(tau, delta) = c1(tau)*delta + c2(tau)*delta^2 + ...      (no constant; alphar(tau,0)=0)
```

Virial coefficients are the Taylor coefficients:

```
B(T)  = (1/rho_r)   * d(alphar)/d(delta)|_{delta=0}        =  (1/rho_r)   * c1
C(T)  = (1/rho_r^2) * d^2(alphar)/d(delta^2)|_{delta=0}    =  (1/rho_r^2) * 2*c2
dB/dT = (1/rho_r)   * d^2(alphar)/(dtau d delta)|_{0}      * dtau/dT
dC/dT = (1/rho_r^2) * d^3(alphar)/(dtau d delta^2)|_{0}    * dtau/dT
```

The **Axy** form `Axy = tau^x * delta^y * d^{x+y}(a)/(dtau^x ddelta^y)`:

- **Pressure:** `p = rho*R*T*(1 + Ar(0,1))`. At `rho=0`, `Ar(0,1)=0` (finite) so `p=0` exactly — *if* `Ar(0,1)` carries no `1/delta` intermediate. Current code does `rho*R*T*(1 + delta*dalphar_dDelta)` = `0*Inf = NaN`.
- **Ideal-gas lead is the canonical case:** `alpha0` lead is `log(delta)`, so `dalpha0_ddelta=1/delta` (`Inf`) but `Aig(0,1)=delta*(1/delta)=1` exactly. All other ideal-gas terms are functions of `tau` only → their `delta`-derivatives are exactly `0`, so on the ideal side **only `IdealHelmholtzLead` has a `delta`-singularity**; every other ideal term contributes only to `A[itau][0]` (finite, mechanical `tau^itau`-scaling of existing `tau`-derivatives).

**Per-term safe-derivation principle:** `Axy` numerator = `delta^y * tau^x * (raw mixed derivative)`. Re-derive each term so every `1/delta^k` cancels against an explicit `delta^k` *analytically* (in source), never as `delta * (something/delta)` at runtime. The gate that this was done correctly: agreement with `one_mcx` autodiff at small `delta`, and finiteness at `delta=0`.

---

## File Structure

- `include/Helmholtz.h` — `HelmholtzDerivatives` gains the `Axy` store + `getA(itau,idelta)`; `ResidualHelmholtzContainer::all` and the ideal container already sum every term (no structural change, but every summed term now fills `A[][]`).
- `src/Helmholtz.cpp` — safe `Axy` accumulation added to **all** term `all()` methods: `ResidualHelmholtzGeneralizedExponential`, `ResidualHelmholtzNonAnalytic`, `ResidualHelmholtzGaoB`, `ResidualHelmholtzSAFTAssociating`, `ResidualHelmholtzGeneralizedCubic`, `ResidualHelmholtzXiangDeiters`; and ideal: `IdealHelmholtzLead` (delta part) + a mechanical `tau`-scaling sweep across all `IdealHelmholtz*` terms.
- `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.{h,cpp}` — `calc_Ar`/`calc_Aig`; rewrite `calc_pressure`, `calc_Bvirial`/`calc_Cvirial`/`_dT`.
- `include/AbstractState.h` — public `get_Ar`/`get_Aig` + virtual `calc_Ar`/`calc_Aig`.
- `src/Tests/CoolProp-tests.cpp` — `[virial_zero_density]` tests, including the **all-fluids sweep**.

---

### Task 1: Universal failing test — pressure is exactly zero at rho=0 for ALL fluids

**Files:** Test: `src/Tests/CoolProp-tests.cpp`

This is the merge gate that enforces "no intermediate state."

- [ ] **Step 1: Write the failing all-fluids test**

```cpp
TEST_CASE("Pressure is exactly zero at zero density for every fluid", "[virial_zero_density]") {
    std::vector<std::string> fluids;
    {
        std::string list = CoolProp::get_global_param_string("fluids_list");
        std::stringstream ss(list); std::string f;
        while (std::getline(ss, f, ',')) fluids.push_back(f);
    }
    REQUIRE(fluids.size() > 100);  // sanity: the full pure-fluid set
    std::vector<std::string> failures;
    for (const std::string& fluid : fluids) {
        try {
            shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", fluid));
            AS->update(CoolProp::DmolarT_INPUTS, 0.0, 300.0);
            if (AS->p() != 0.0) failures.push_back(fluid + " (p=" + std::to_string(AS->p()) + ")");
        } catch (const std::exception& e) {
            failures.push_back(fluid + " (threw: " + e.what() + ")");
        }
    }
    CAPTURE(failures);
    CHECK(failures.empty());
}
```

> Note: a few fluids may legitimately have `T=300 K` below their triple point; if so, switch that fluid's probe `T` to `1.01 * T_triple` (still zero density). Keep the *fluid set* complete — do not silently drop fluids.

- [ ] **Step 2: Run to verify it fails**

Run: `cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[virial_zero_density]"`
Expected: FAIL — `failures` lists effectively every fluid (`update` throws "p is not a valid number").

- [ ] **Step 3: Commit the failing test**

```bash
git add src/Tests/CoolProp-tests.cpp
git commit -m "test(virial): all fluids must give p==0 at rho=0 (CoolProp-izv)"
```

---

### Task 2: `HelmholtzDerivatives` Axy store + `getA` accessor

**Files:** Modify `include/Helmholtz.h:60-144`

- [ ] **Step 1: Add product-form storage and accessor**

```cpp
// Axy = tau^itau * delta^idelta * d^{itau+idelta} a / (dtau^itau ddelta^idelta),
// accumulated directly by each term WITHOUT dividing by delta -> finite at delta=0.
// Range needed: itau in {0,1}, idelta in {0,1,2} (pressure (0,1)/(0,2); virials (1,1)/(1,2)).
double A[2][3] = {{0,0,0},{0,0,0}};

double getA(std::size_t itau, std::size_t idelta) const {
    if (itau > 1 || idelta > 2) throw ValueError("getA: (itau,idelta) out of supported range");
    return A[itau][idelta];
}
```

- [ ] **Step 2: Build** — Run: `cmake --build build_catch --target CatchTestRunner -j8`; Expected: compiles (zeros default, no behavior change).

- [ ] **Step 3: Commit**

```bash
git add include/Helmholtz.h
git commit -m "feat(virial): add zero-density-safe Axy store to HelmholtzDerivatives (CoolProp-izv)"
```

---

### Task 3: Per-term multicomplex cross-check harness (the correctness oracle)

**Files:** Test: `src/Tests/CoolProp-tests.cpp`

Build the validation tool *first* so every subsequent per-term task (Tasks 4–10) has a concrete gate: hand-derived `A[x][y]` must equal `tau^x*delta^y * (mcx derivative)` at small `delta`, and be finite at `delta=0`.

- [ ] **Step 1: Add a reusable checker (templated over fluid + (x,y) list)**

```cpp
// For a fluid + state, assert get_Ar(x,y) matches tau^x delta^y * (numeric mcx deriv)
// at delta away from 0, and is finite at delta=0.
static void check_Ar_against_mcx(const std::string& fluid, double tau_test, double delta_test) {
    shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", fluid));
    double Tr = AS->T_reducing(), rhor = AS->rhomolar_reducing();
    double T = Tr / tau_test, rho = delta_test * rhor;
    AS->update(CoolProp::DmolarT_INPUTS, rho, T);
    // (x,y) in {(0,1),(0,2),(1,1),(1,2),(1,0)}: compare to delta^y tau^x * raw deriv (existing path)
    CHECK(AS->get_Ar(0,1) == Approx(delta_test            * AS->dalphar_dDelta()      ).epsilon(1e-9));
    CHECK(AS->get_Ar(0,2) == Approx(delta_test*delta_test  * AS->d2alphar_dDelta2()    ).epsilon(1e-9));
    CHECK(AS->get_Ar(1,1) == Approx(delta_test*tau_test    * AS->d2alphar_dDelta_dTau()).epsilon(1e-9));
    // finiteness at zero density:
    AS->update(CoolProp::DmolarT_INPUTS, 0.0, T);
    CHECK(std::isfinite(AS->get_Ar(0,1)));
    CHECK(std::isfinite(AS->get_Ar(0,2)));
    CHECK(std::isfinite(AS->get_Ar(1,1)));
}
```

> The existing fast path (`dalphar_dDelta()` etc.) is the consistency oracle *away from* `delta=0`; `one_mcx` is the independent oracle used in Task 11 for absolute accuracy. The `epsilon(1e-9)` tolerance accommodates ULP-class noise — do **not** bit-compare.

- [ ] **Step 2: Build** — Run: `cmake --build build_catch --target CatchTestRunner -j8`; Expected: compiles (helper unused until per-term tasks call it).

- [ ] **Step 3: Commit**

```bash
git add src/Tests/CoolProp-tests.cpp
git commit -m "test(virial): Ar-vs-fast-path/mcx cross-check harness (CoolProp-izv)"
```

---

### Tasks 4–9: Per-term safe Axy accumulation

Each task follows the identical micro-cycle. **Apply this template to each term**, picking a representative fluid that *activates* that term:

| Task | Term (`src/Helmholtz.cpp`) | Representative fluid | `1/delta`? |
|------|----------------------------|----------------------|------------|
| 4 | `ResidualHelmholtzGeneralizedExponential::all` (`:138`) | Methane | yes (`one_over_delta`) |
| 5 | `ResidualHelmholtzNonAnalytic::all` (`:386`) | Water | no (already finite) |
| 6 | `ResidualHelmholtzGaoB::all` (`:665`) | (fluid using GaoB) | yes (explicit `/delta :735`) |
| 7 | `ResidualHelmholtzSAFTAssociating::all` (`:1077`) | (SAFT fluid) | via `dX_ddelta` |
| 8 | `ResidualHelmholtzGeneralizedCubic::all` (`:615`) | (cubic-backed fluid) | analytic at 0 |
| 9 | `ResidualHelmholtzXiangDeiters::all` (`:800`) | (XD fluid) | inherits GenExp-like |

> Executor: identify each representative fluid by scanning `dev/fluids/*.json` for the term's signature fields (e.g. `gaussian`/`GERG-2008` exponential for GenExp; `nonanalytic` for NonAnalytic; `GaoB` for GaoB; `SAFT`/`association` for SAFT; `alphar` `type":"ResidualHelmholtzGeneralizedCubic"` for cubic; `XiangDeiters` for XD). If a term type has **no** fluid in the database, still implement its safe `Axy` (it can appear via user EOS) and unit-test it with a synthetic term in `[virial_zero_density]`.

**Per-term micro-cycle (template — repeat for each of Tasks 4–9):**

- [ ] **Step 1: Write the failing per-term test** — call the harness for the representative fluid:

```cpp
TEST_CASE("<TermName> Axy is consistent and finite at zero density", "[virial_zero_density]") {
    check_Ar_against_mcx("<RepresentativeFluid>", /*tau*/ 1.2, /*delta*/ 0.7);
}
```

- [ ] **Step 2: Run to verify it fails** — Run: `./build_catch/CatchTestRunner "[virial_zero_density]"`; Expected: FAIL (term does not fill `A[][]` yet, and/or `get_Ar` missing until Task 12 — defer the test body with `[.]` if executing strictly before Task 12, then un-defer).

- [ ] **Step 3: Derive and add the safe `Axy` accumulation** in that term's `all()`, accumulating into `derivs.A[x][y]` using the per-term principle (every `1/delta^k` cancelled analytically against an explicit `delta^k`). Concrete guidance per term:
  - **GenExp (Task 4):** the per-element numerators `ndteu*B_delta`, `ndteu*B_delta2`, `ndteu*B_delta*B_tau`, `ndteu*B_delta2*B_tau`, `ndteu*B_tau` are the `Axy` numerators — but `B_delta` is currently built from `du_ddelta = l*u*one_over_delta`. Rebuild a `Bdelta_safe` where `delta*du_ddelta` is formed directly as a `delta`-power (`-ci*l*pow(delta,l)` for the `delta^l` part; `-2*eta1*delta*(delta-epsilon1)` for the Gaussian part; mirror every active `u`-contribution). Accumulate `derivs.A[0][1]+=ndteu*Bdelta_safe;` etc.
  - **NonAnalytic (Task 5):** raw derivs already finite; accumulate `A[0][1]+=delta*<dalphar_ddelta contribution>`, `A[0][2]+=delta*delta*<d2alphar_ddelta2 contribution>`, `A[1][1]+=delta*tau*<d2alphar_ddelta_dtau>`, `A[1][2]+=delta*delta*tau*<d3alphar_ddelta2_dtau>`, `A[1][0]+=tau*<dalphar_dtau>`.
  - **GaoB (Task 6):** the `:735` form has explicit `/delta` (`n*Ftau*deltadFdeltaddelta/delta`); the `Axy` numerator `delta*dalphar_ddelta = n*Ftau*deltadFdeltaddelta` drops the `/delta` — accumulate that directly. Derive `A[0][2]`/`A[1][*]` analogously from the existing `deltadFdeltaddelta`/`Ftau` building blocks without re-introducing `/delta`.
  - **SAFT (Task 7):** `dalphar_ddelta = m*a*(1/X-0.5)*dX_ddelta` (`:1093`). Form `A[0][1]=delta*` that. Check whether `dX_ddelta`/`X` themselves carry a `1/delta` near `delta=0`; if so, re-derive `delta*dX_ddelta` directly. Gate on the harness.
  - **cubic (Task 8):** delegates to `m_abstractcubic->alphar(tau,delta,z,itau,idelta)`. Cubic EOS are analytic at `rho=0`; accumulate `A[x][y]=tau^x*delta^y*m_abstractcubic->alphar(tau,delta,z,x,y)`. Verify finiteness at `delta=0` (cubics: `alphar ~ a-function * ln((delta+...)/(delta+...))`, derivative finite at 0).
  - **XiangDeiters (Task 9):** sum of GenExp-like sub-terms (`phi0/phi1/phi2`); apply the GenExp treatment to each, scaled by `1`/`acentric`/`theta`.

- [ ] **Step 4: Run the per-term test (+ harness)** — Run: `./build_catch/CatchTestRunner "[virial_zero_density]"`; Expected: PASS (consistent away from 0 within `1e-9`, finite at 0).

- [ ] **Step 5: Commit** — `git commit -m "feat(virial): zero-density-safe Axy for <TermName> (CoolProp-izv)"`

---

### Task 10: Ideal-gas terms — Lead delta-part + tau-scaling sweep

**Files:** Modify `src/Helmholtz.cpp` ideal terms (`IdealHelmholtzLead :1143` and the `IdealHelmholtz*` family `:1154-1430`)

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("Aig(0,1)==1 and Aig is finite at zero density for every fluid", "[virial_zero_density]") {
    std::string list = CoolProp::get_global_param_string("fluids_list");
    std::stringstream ss(list); std::string f; std::vector<std::string> failures;
    while (std::getline(ss, f, ',')) {
        try {
            shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", f));
            AS->update(CoolProp::DmolarT_INPUTS, 0.0, 300.0);
            if (!(AS->get_Aig(0,1) == Approx(1.0))) failures.push_back(f);
        } catch (...) { failures.push_back(f + " (threw)"); }
    }
    CAPTURE(failures); CHECK(failures.empty());
}
```

- [ ] **Step 2: Run to verify it fails** — Expected: FAIL.

- [ ] **Step 3: Implement.** In `IdealHelmholtzLead::all` add `A[0][1]+=1.0; A[0][2]+=-1.0; A[1][0]+=a2*tau; A[1][1]+=0; A[1][2]+=0;` (see Key Math). In **every other** `IdealHelmholtz*::all`, add the `tau`-scaled forms `A[1][0]+=tau*<their dalpha0_dtau contribution>` and (since they have no `delta` dependence) `A[0][1]+=0; A[0][2]+=0; A[1][1]+=0; A[1][2]+=0;` — i.e. only `A[itau][0]` is non-zero for non-Lead ideal terms. `A[0][0]` (contains `log(delta)` in Lead) is intentionally non-finite at 0 and is never read by pressure/virials; comment it.

- [ ] **Step 4: Run** (after Task 12) — Expected: PASS.

- [ ] **Step 5: Commit** — `git commit -m "feat(virial): zero-density-safe Aig for all ideal-gas terms (CoolProp-izv)"`

---

### Task 11: Backend `calc_Ar`/`calc_Aig`

**Files:** Modify `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.{h,cpp}` (near `:598`/`:3414`)

- [ ] **Step 1: Declare** (header):

```cpp
CoolPropDbl calc_Ar(std::size_t itau, std::size_t idelta) override;
CoolPropDbl calc_Aig(std::size_t itau, std::size_t idelta) override;
```

- [ ] **Step 2: Implement** (match the exact residual/ideal container call signatures already used by `calc_alphar_deriv_nocache`/`calc_alpha0_deriv_nocache`):

```cpp
CoolPropDbl HelmholtzEOSMixtureBackend::calc_Ar(std::size_t itau, std::size_t idelta) {
    HelmholtzDerivatives d = residual_helmholtz->all(*this, mole_fractions, _tau, _delta, false);
    return d.getA(itau, idelta);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_Aig(std::size_t itau, std::size_t idelta) {
    HelmholtzDerivatives d = /* ideal container */->all(_tau, _delta);  // mirror calc_alpha0_deriv_nocache plumbing
    return d.getA(itau, idelta);
}
```

- [ ] **Step 3: Build** — Expected: compiles.

- [ ] **Step 4: Commit** — `git commit -m "feat(virial): backend calc_Ar/calc_Aig (CoolProp-izv)"`

---

### Task 12: Public `get_Ar`/`get_Aig` API on AbstractState

**Files:** Modify `include/AbstractState.h` (virtuals near `:301`; public near `:1611`)

- [ ] **Step 1: Base-class virtuals (throw by default)**

```cpp
/// Residual Axy = tau^itau delta^idelta d^{itau+idelta} alphar/(dtau^itau ddelta^idelta); finite at delta=0.
virtual CoolPropDbl calc_Ar(std::size_t itau, std::size_t idelta) { throw NotImplementedError("calc_Ar not implemented"); }
/// Ideal-gas Axy; get_Aig(0,1)==1 (lead term log(delta)).
virtual CoolPropDbl calc_Aig(std::size_t itau, std::size_t idelta) { throw NotImplementedError("calc_Aig not implemented"); }
```

- [ ] **Step 2: Public accessors**

```cpp
CoolPropDbl get_Ar(std::size_t itau, std::size_t idelta) { return calc_Ar(itau, idelta); }
CoolPropDbl get_Aig(std::size_t itau, std::size_t idelta) { return calc_Aig(itau, idelta); }
```

- [ ] **Step 3: Un-defer Tasks 4–10 tests** (remove `[.]`), build + run — Run: `./build_catch/CatchTestRunner "[virial_zero_density]"`; Expected: every per-term test + `Aig(0,1)==1` all-fluids test PASS.

- [ ] **Step 4: Commit** — `git commit -m "feat(virial): public get_Ar/get_Aig API (teqp-style) (CoolProp-izv)"`

---

### Task 13: Rewrite pressure and virials

**Files:** Modify `HelmholtzEOSMixtureBackend.cpp:3031-3049` (pressure), `:1634-1649` (virials)

- [ ] **Step 1: `calc_pressure` on `Ar(0,1)` (+ explicit rho==0 short-circuit as belt-and-suspenders)**

```cpp
CoolPropDbl HelmholtzEOSMixtureBackend::calc_pressure() {
    _delta = _rhomolar / _reducing.rhomolar;
    _tau   = _reducing.T / _T;
    if (_rhomolar == 0) { _p = 0.0; return 0.0; }   // exact: ideal+residual both vanish
    CoolPropDbl R_u = gas_constant();
    _p = _rhomolar * R_u * _T * (1 + calc_Ar(0, 1)); // Ar(0,1) finite everywhere; == delta*dalphar_dDelta for rho>0
    return static_cast<CoolPropDbl>(_p);
}
```

> `calc_Ar(0,1) == delta*dalphar_dDelta` away from 0 (Task 3 gate), so finite-density behavior is unchanged. If `all()` recompute is hot, read `Ar(0,1)` from the cached `HelmholtzDerivatives` rather than recomputing — verify on the flash benchmark.

- [ ] **Step 2: Virials via exact zero-density coefficients (no 1e-12)**

`B*rho_r = c1 = d(alphar)/d(delta)|_0`. Expose `c1,c2,dc1_dtau,dc2_dtau` from the residual container (accumulated in the same Task 4–9 loops as the `d==1`/`d==2` zero-density contributions), then:

```cpp
CoolPropDbl HelmholtzEOSMixtureBackend::calc_Bvirial()  { return c1_zero_density() / rhomolar_reducing(); }
CoolPropDbl HelmholtzEOSMixtureBackend::calc_Cvirial()  { return 2.0*c2_zero_density() / pow(rhomolar_reducing(),2); }
CoolPropDbl HelmholtzEOSMixtureBackend::calc_dBvirial_dT(){ SimpleState r=get_reducing_state(); return dc1_dtau_zero_density()/r.rhomolar * (-r.T/pow(_T,2)); }
CoolPropDbl HelmholtzEOSMixtureBackend::calc_dCvirial_dT(){ SimpleState r=get_reducing_state(); return 2.0*dc2_dtau_zero_density()/pow(r.rhomolar,2) * (-r.T/pow(_T,2)); }
```

- [ ] **Step 3: Run** — Run: `./build_catch/CatchTestRunner "[virial_zero_density]" "[virial]"`; Expected: Task 1 all-fluids pressure test PASSES; existing `[virial]` tests still pass.

- [ ] **Step 4: Commit** — `git commit -m "fix(virial): exact virials + p(rho=0)==0 via safe Axy, all fluids (CoolProp-izv, #2991)"`

---

### Task 14: Absolute validation vs REFPROP + multicomplex

**Files:** Test: `src/Tests/CoolProp-tests.cpp`

- [ ] **Step 1: REFPROP cross-check (runs locally — REFPROP installed)**

```cpp
TEST_CASE("Bvirial/Cvirial match REFPROP", "[virial_zero_density][refprop]") {
    for (const std::string& fluid : {"Nitrogen","CarbonDioxide","Methane","Water","Propane"}) {
        CAPTURE(fluid);
        shared_ptr<CoolProp::AbstractState> H(CoolProp::AbstractState::factory("HEOS", fluid));
        shared_ptr<CoolProp::AbstractState> R(CoolProp::AbstractState::factory("REFPROP", fluid));
        for (double T : {250.0,300.0,400.0}) {
            H->update(CoolProp::DmolarT_INPUTS, 0.0, T);
            R->update(CoolProp::DmolarT_INPUTS, 1e-10, T);
            CHECK(H->Bvirial() == Approx(R->Bvirial()).epsilon(1e-6));
            CHECK(H->Cvirial() == Approx(R->Cvirial()).epsilon(1e-5));
        }
    }
}
```

- [ ] **Step 2: Multicomplex per-term absolute check** — for each term type's representative fluid, compute `delta*d(alphar)/ddelta` via `one_mcx` complex-step at `delta=1e-6` and compare to `get_Ar(0,1)`; confirms each hand-derived form matches autodiff (not just the existing fast path).

- [ ] **Step 3: Run** — Run: `./build_catch/CatchTestRunner "[virial_zero_density]"`; Expected: PASS.

- [ ] **Step 4: Commit** — `git commit -m "test(virial): REFPROP + multicomplex absolute validation (CoolProp-izv)"`

---

### Task 15: Pre-PR gate

- [ ] **Step 1:** `./dev/ci/preflight.sh` passes (auto tag scope; `[virial_zero_density]`/`[refprop]` run locally).
- [ ] **Step 2:** Invoke `superpowers:code-reviewer` against the diff (CLAUDE.md mandatory). Focus: `0*Inf` reintroduction in any term; `NaN` slipping through comparisons; per-term `Axy` vs fast-path consistency away from 0; `getA` index guard; perf of `calc_pressure` (cached vs recompute); the `A[0][0]` non-finite-at-0 comment is present and correct.
- [ ] **Step 3:** Address/justify findings → `git push` → `gh pr create` referencing #2991 and #1676.

---

## Self-Review

**Spec coverage:**
- "Separate path for virial coefficients at zero density" → Task 13 Step 2 (exact Taylor coefficients).
- "General (non-virial) derivative path via Axy" → Tasks 2–12 (`get_Ar`/`get_Aig`).
- "Pressure at rho=0 is zero, not invalid" → Task 1 (all-fluids gate) + Task 13 Step 1.
- "No 1e-12 fudge" → Task 13 Step 2.
- "Similar for ideal-gas terms" → Task 10 + `get_Aig`.
- "No intermediate state / all fluids" → every term converted (Tasks 4–10), all-fluids sweep gates merge (Tasks 1, 10).
- "Better API" → teqp-style `get_Ar`/`get_Aig` (user-selected).
- "Validate" → Task 3 (fast-path), Task 14 (REFPROP + multicomplex).

**Type consistency:** `getA(itau,idelta)`, `A[2][3]`, `calc_Ar`/`calc_Aig`, `get_Ar`/`get_Aig`, `c1_zero_density()`/`c2_zero_density()`/`dc1_dtau_zero_density()`/`dc2_dtau_zero_density()` used consistently across Tasks 2/11/12/13. Index range `itau∈{0,1}, idelta∈{0,1,2}` matches pressure (0,1)/(0,2) and virials (1,1)/(1,2).

**Derivations are pre-computed** — see `docs/superpowers/derivations/virial-axy.md` (and `dev/derivations/virial_axy_derivations.py` to reproduce). Every term's `Axy` closed form (ccode), finiteness-at-δ=0 proof, and virial `c1`/`c2`/`dc1_dtau`/`dc2_dtau` are derived and machine-checked there. Tasks 4–10 transcribe those forms; the §9 implementation-mapping table in that doc says exactly what goes in each `all()`. Remaining verify-in-tree items: backend container call signatures (Task 11) and the `c1/c2` extraction plumbing (Task 13) must match existing `calc_alpha*_deriv_nocache` — shapes given, confirm against the tree.
