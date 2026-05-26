# Phase 2 — Zero-Density-Safe Axy (`get_Ar`/`get_Aig`) + `p(rho=0)==0` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Beads:** CoolProp-izv · **Derivations:** CoolProp-6gf (`docs/superpowers/derivations/virial-axy.md`) · **GitHub:** #2991 (child #1676) · **Phase 1 (separate, virials):** CoolProp-dr3

**Goal:** Expose `get_Ar(itau,idelta)`/`get_Aig(itau,idelta)` (= `tau^x*delta^y*d^{x+y}a/(dtau^x ddelta^y)`, finite at `delta=0`) and make `p(rho=0)` return exactly `0` (not `NaN`) for **every** fluid and mixture.

**Architecture:** Add, in every Helmholtz term's `all()`, a parallel analytic accumulation of the `Axy` product form (no runtime `1/delta`), stored in `HelmholtzDerivatives::A[itau][idelta]`. `get_Ar`/`get_Aig` read it; `p = rho*R*T*(1 + Ar(0,1))`. Each hand-derived `Axy` is gated against `delta^y*tau^x*(raw deriv)` away from 0 and against `one_mcx` autodiff.

**All-or-nothing (hard requirement):** unlike Phase 1's virials, `p(rho=0)` is binary `0`/`NaN` with no valid fallback — so EVERY term is converted in this one landing and the merge gate is an all-fluids + all-mixtures sweep. No phasing, no deferral.

**Independent of Phase 1:** can land before or after CoolProp-dr3 (recommended: after). Phase 1 does not create the `A[][]` store; this phase does.

**Tech Stack:** C++ (CoolProp core), Catch2, CMake, REFPROP (local), multicomplex (`mcx`).

---

## Mixtures (per user request)

`p = rho*R*T*(1 + Ar(0,1))` and `get_Ar` both read the composition-weighted `residual_helmholtz` container, so mixtures route through unchanged — but the `rho=0` gate and a `get_Ar` consistency check MUST be run on mixtures (Tasks 1, 3), not assumed. Binary/ternary excess departure terms are GenExp-type, covered by Task 4.

---

## File Structure

- `include/Helmholtz.h` — `HelmholtzDerivatives` gains `double A[2][3]` + `getA(itau,idelta)`.
- `src/Helmholtz.cpp` — safe `Axy` accumulation in every term's `all()`: `GenExp`, `NonAnalytic`, `GaoB`, `SAFT`, cubic, `XiangDeiters`; ideal `Lead` (+ τ-sweep of other ideal terms).
- `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.{h,cpp}` — `calc_Ar`/`calc_Aig`; rewrite `calc_pressure`.
- `include/AbstractState.h` — public `get_Ar`/`get_Aig` + virtual `calc_Ar`/`calc_Aig`.
- `src/Tests/CoolProp-tests.cpp` — `[axy_zero_density]` tests (all-fluids + mixtures + per-term + mcx).

---

### Task 1: Universal failing test — `p(rho=0)==0` for all fluids AND mixtures

**Files:** Test: `src/Tests/CoolProp-tests.cpp`

The merge gate ("no intermediate state").

- [ ] **Step 1: All-fluids sweep**

```cpp
TEST_CASE("Pressure is exactly zero at zero density for every fluid", "[axy_zero_density]") {
    std::vector<std::string> fluids; {
        std::stringstream ss(CoolProp::get_global_param_string("fluids_list")); std::string f;
        while (std::getline(ss, f, ',')) fluids.push_back(f);
    }
    REQUIRE(fluids.size() > 100);
    std::vector<std::string> fail;
    for (auto& f : fluids) {
        try { shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", f));
              AS->update(CoolProp::DmolarT_INPUTS, 0.0, 300.0);
              if (AS->p() != 0.0) fail.push_back(f + " p=" + std::to_string(AS->p()));
        } catch (std::exception& e) { fail.push_back(f + " threw:" + e.what()); }
    }
    CAPTURE(fail); CHECK(fail.empty());
}
```

> A few fluids may have `T=300 K` below the triple point; switch that fluid's probe `T` to `1.01*T_triple`. Keep the set complete.

- [ ] **Step 2: Mixture sweep (per user request)**

```cpp
TEST_CASE("Pressure is exactly zero at zero density for mixtures", "[axy_zero_density]") {
    struct M { const char* names; std::vector<double> z; };
    std::vector<std::string> fail;
    for (auto mix : { M{"Methane&Ethane",{0.7,0.3}}, M{"Nitrogen&CarbonDioxide",{0.5,0.5}},
                      M{"Methane&Nitrogen&CarbonDioxide",{0.8,0.1,0.1}}, M{"Water&Ammonia",{0.6,0.4}} }) {
        try { shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", mix.names));
              AS->set_mole_fractions(mix.z);
              AS->update(CoolProp::DmolarT_INPUTS, 0.0, 320.0);
              if (AS->p() != 0.0) fail.push_back(std::string(mix.names));
        } catch (std::exception& e) { fail.push_back(std::string(mix.names) + " threw:" + e.what()); }
    }
    CAPTURE(fail); CHECK(fail.empty());
}
```

- [ ] **Step 3: Run** — Expected: FAIL (throws "p is not a valid number" today).
- [ ] **Step 4: Commit** — `git commit -m "test(axy): p==0 at rho=0 for all fluids + mixtures (CoolProp-izv)"`

---

### Task 2: `HelmholtzDerivatives` Axy store + `getA`

**Files:** Modify `include/Helmholtz.h:60-144`

- [ ] **Step 1:**

```cpp
// Axy = tau^itau*delta^idelta*d^{itau+idelta}a/(dtau^itau ddelta^idelta), 1/delta-free -> finite at 0.
// Range: itau in {0,1}, idelta in {0,1,2} (pressure (0,1); cross terms (1,1)/(1,2)).
double A[2][3] = {{0,0,0},{0,0,0}};
double getA(std::size_t itau, std::size_t idelta) const {
    if (itau > 1 || idelta > 2) throw ValueError("getA: (itau,idelta) out of supported range");
    return A[itau][idelta];
}
```

- [ ] **Step 2: Build** — compiles (zeros default). **Step 3: Commit.**

---

### Task 3: Per-term Axy cross-check harness (pure + mixture)

**Files:** Test: `src/Tests/CoolProp-tests.cpp`

Build the validation gate before per-term work.

- [ ] **Step 1:**

```cpp
static void check_Ar(const std::string& fluid, std::vector<double> z, double tau_t, double delta_t) {
    shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", fluid));
    if (z.size() > 1) AS->set_mole_fractions(z);
    double T = AS->T_reducing()/tau_t, rho = delta_t*AS->rhomolar_reducing();
    AS->update(CoolProp::DmolarT_INPUTS, rho, T);
    CHECK(AS->get_Ar(0,1) == Approx(delta_t            * AS->dalphar_dDelta()      ).epsilon(1e-9));
    CHECK(AS->get_Ar(0,2) == Approx(delta_t*delta_t    * AS->d2alphar_dDelta2()    ).epsilon(1e-9));
    CHECK(AS->get_Ar(1,1) == Approx(delta_t*tau_t      * AS->d2alphar_dDelta_dTau()).epsilon(1e-9));
    AS->update(CoolProp::DmolarT_INPUTS, 0.0, T);
    CHECK(std::isfinite(AS->get_Ar(0,1))); CHECK(std::isfinite(AS->get_Ar(0,2)));
}
```

- [ ] **Step 2: Build** (helper unused yet). **Step 3: Commit.**

---

### Tasks 4–9: Per-term safe Axy accumulation

Identical micro-cycle per term; pick a fluid that activates the term. **Transcribe the `Axy` ccode from `docs/superpowers/derivations/virial-axy.md` §2–7** (the source of truth — finiteness already proven there).

| Task | Term (`src/Helmholtz.cpp`) | Rep. fluid | derivations § |
|------|----------------------------|-----------|---------------|
| 4 | `GeneralizedExponential::all` (`:138`) | Methane | §3 |
| 5 | `NonAnalytic::all` (`:386`) | Water | §4 |
| 6 | `GaoB::all` (`:665`) | (GaoB fluid) | §5 |
| 7 | `SAFTAssociating::all` (`:1077`) | (SAFT fluid) | §6 |
| 8 | `GeneralizedCubic::all` (`:615`) | (cubic fluid) | §7 |
| 9 | `XiangDeiters::all` (`:800`) | (XD fluid) | §3 (per sub-term) |

**Micro-cycle (repeat per term):**
- [ ] **Step 1:** failing test calling `check_Ar("<fluid>", {}, 1.2, 0.7);`
- [ ] **Step 2:** run → FAIL (term doesn't fill `A[][]`; `get_Ar` lands Task 11 — defer with `[.]` if strict).
- [ ] **Step 3:** accumulate `derivs.A[x][y] += <form from virial-axy.md §>`. Key points per term:
  - GenExp (§3): compute `Bdelta`/`Bdelta2` via `delta`-powers (no `one_over_delta`); `A01=ndteu*Bdelta`, etc.
  - NonAnalytic (§4): scale existing raw derivs by `delta^y*tau^x` (already finite at 0).
  - GaoB (§5): use the §5 forms (drop the `:735` `/delta`).
  - SAFT (§6): scale existing `X`-derivative-based contributions by `delta^y*tau^x`.
  - cubic (§7): `A[x][y]=tau^x*delta^y*m_abstractcubic->alphar(tau,delta,z,x,y)`.
  - XiangDeiters: combine GenExp sub-terms (`phi0+acentric*phi1+theta*phi2`).
- [ ] **Step 4:** run → PASS (consistency ≤1e-9 away from 0; finite at 0).
- [ ] **Step 5:** commit `feat(axy): zero-density-safe Axy for <Term> (CoolProp-izv)`.

---

### Task 10: Ideal-gas Axy — Lead + τ-sweep

**Files:** Modify ideal terms in `src/Helmholtz.cpp` (`Lead :1143`, family `:1154-1430`)

- [ ] **Step 1: Failing test** — all-fluids `get_Aig(0,1)==1`:

```cpp
TEST_CASE("Aig(0,1)==1 for every fluid", "[axy_zero_density]") {
    std::stringstream ss(CoolProp::get_global_param_string("fluids_list")); std::string f;
    std::vector<std::string> fail;
    while (std::getline(ss, f, ',')) { try {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", f));
        AS->update(CoolProp::DmolarT_INPUTS, 0.0, 300.0);
        if (!(AS->get_Aig(0,1) == Approx(1.0))) fail.push_back(f);
    } catch (...) { fail.push_back(f+" threw"); } }
    CAPTURE(fail); CHECK(fail.empty());
}
```

- [ ] **Step 2: Run** → FAIL.
- [ ] **Step 3: Implement** per `virial-axy.md` §2: `Lead`: `A[0][1]+=1; A[0][2]+=-1; A[1][0]+=a2*tau;` (rest 0). Every other ideal term: `A[1][0]+=tau*<its dalpha0_dtau>` (δ-orders 0). `A00` (with `log(delta)`) left non-finite, never read.
- [ ] **Step 4: Run** (after Task 11) → PASS. **Step 5: Commit.**

---

### Task 11: Backend `calc_Ar`/`calc_Aig` + public `get_Ar`/`get_Aig`

**Files:** `HelmholtzEOSMixtureBackend.{h,cpp}` (`:598`/`:3414`); `include/AbstractState.h` (`:301`/`:1611`)

- [ ] **Step 1: Backend**

```cpp
CoolPropDbl HelmholtzEOSMixtureBackend::calc_Ar(std::size_t it, std::size_t id) {
    return residual_helmholtz->all(*this, mole_fractions, _tau, _delta, false).getA(it, id);
}
CoolPropDbl HelmholtzEOSMixtureBackend::calc_Aig(std::size_t it, std::size_t id) {
    return /* ideal container, mirror calc_alpha0_deriv_nocache */->all(_tau, _delta).getA(it, id);
}
```

- [ ] **Step 2: AbstractState** — virtual `calc_Ar`/`calc_Aig` (throw `NotImplementedError` by default) + public `get_Ar`/`get_Aig` forwarding to them.
- [ ] **Step 3: Un-defer Tasks 4–10 tests**, build + run `[axy_zero_density]` → per-term + `Aig(0,1)==1` PASS.
- [ ] **Step 4: Commit** — `feat(axy): get_Ar/get_Aig public API (teqp-style) (CoolProp-izv)`.

---

### Task 12: Rewrite `calc_pressure`

**Files:** `HelmholtzEOSMixtureBackend.cpp:3031-3049`

- [ ] **Step 1:**

```cpp
CoolPropDbl HelmholtzEOSMixtureBackend::calc_pressure() {
    _delta = _rhomolar / _reducing.rhomolar; _tau = _reducing.T / _T;
    if (_rhomolar == 0) { _p = 0.0; return 0.0; }            // exact short-circuit
    _p = _rhomolar * gas_constant() * _T * (1 + calc_Ar(0,1)); // == old form for rho>0
    return static_cast<CoolPropDbl>(_p);
}
```

> `calc_Ar(0,1)==delta*dalphar_dDelta` away from 0 (Task 3 gate) ⇒ finite-density behavior unchanged. If `all()` recompute is hot, read `Ar(0,1)` from the cached `HelmholtzDerivatives`; verify on the flash benchmark.

- [ ] **Step 2: Run** `[axy_zero_density]` — Task 1 all-fluids + mixture pressure tests PASS.
- [ ] **Step 3: Commit** — `fix(axy): p(rho=0)==0 via safe Ar(0,1), all fluids + mixtures (CoolProp-izv, #2991)`.

---

### Task 13: Multicomplex absolute validation + Task 1 close-out

**Files:** Test: `src/Tests/CoolProp-tests.cpp`

- [ ] **Step 1:** for each term's rep fluid, compare `get_Ar(0,1)` to `one_mcx` complex-step `delta*d(ar)/ddelta` at `delta=1e-6` (absolute oracle, beyond the fast-path consistency).
- [ ] **Step 2:** confirm `[axy_zero_density]` fully green (all fluids + mixtures).
- [ ] **Step 3: Commit** — `test(axy): multicomplex absolute validation of Axy (CoolProp-izv)`.

---

### Task 14: Pre-PR gate

- [ ] `./dev/ci/preflight.sh` passes.
- [ ] `superpowers:code-reviewer` on the diff (CLAUDE.md mandatory). Focus: `0*Inf` reintroduction in any term; per-term `Axy` vs fast-path consistency away from 0; mixture pressure path; `getA` index guard; `calc_pressure` perf (cached vs recompute); `A00` non-finite-at-0 documented.
- [ ] Restore `.beads/issues.jsonl`, `git push`, `gh pr create` referencing #2991, #1676.

---

## Self-Review

- **p(rho=0)==0, all fluids + mixtures:** Task 1 (gate) + Task 12 (fix).
- **get_Ar/get_Aig (teqp-style):** Tasks 2–11.
- **Every term converted (no intermediate state):** Tasks 4–10; all-fluids+mixture sweep gates merge.
- **Ideal-gas:** Task 10 + `get_Aig`.
- **Mixtures:** Tasks 1, 3 (explicit).
- **Closed forms / finiteness:** transcribed from `docs/superpowers/derivations/virial-axy.md` §2–7; per-term gated by Task 3 (fast-path) + Task 13 (mcx).
- **Virials are NOT here** — they're Phase 1 (CoolProp-dr3). This phase may later refactor virials onto `get_Ar`, but that's optional and out of scope.
- **Verify-in-tree:** ideal-container call signature for `calc_Aig` vs `calc_alpha0_deriv_nocache`; GenExp element optional-field names.
