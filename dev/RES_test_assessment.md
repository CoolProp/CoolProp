# RES transport test suite — assessment

Status of `src/Tests/CoolProp-Tests-RES.cpp`, tag `[RES]`, as measured against the current
working tree.

    test cases:  19 |  13 passed |  6 failed
    assertions: 680 | 657 passed | 23 failed

All 13 structural/guard tests pass. All 6 failures are the six sample-data comparisons.
Every failure below has been traced to a cause; **none is an unexplained defect**.

---

## 1. Inventory

### Sample-data comparisons (6) — all currently failing
| # | test | compares | tolerances |
|---|------|----------|-----------|
| 1 | RES viscosity pure fluids (Martinek 2025) | HEOS RES vs `vis_res`, and vs `vis_exp` | 3% / 15% |
| 2 | RES viscosity binary mixtures (Martinek 2025) | HEOS RES vs `vis_res`, and vs `vis_exp` | 5% / 20% |
| 3 | RES conductivity pure fluids (Li 2024) | HEOS RES vs `TC_RES`, and vs `TC_EXP` | 3% / 15% |
| 4 | RES conductivity binary mixtures (Li 2024) | HEOS RES vs `TC_RES`, and vs `TC_EXP` | 5% / 20% |
| 5 | RES viscosity PR backend pure | PR RES vs `vis_exp` only | 15% |
| 6 | RES conductivity PR backend pure | PR RES vs `TC_EXP` only | 15% |

### Guard / API tests (13) — all passing
| # | test | guards |
|---|------|--------|
| 7 | RES throws after alpha change (viscosity) | `n_params_match_alpha` |
| 8 | RES throws after alpha change (conductivity) | same |
| 9 | RES recovers after `set_viscosity_RES_residual_params` | guard reset |
| 10 | RES toggles and setters invalidate memoized values | cache staleness *(new)* |
| 11 | RES conductivity finite with no fitted critical params | NaN from `t_ref==0` *(new)* |
| 12 | RES conductivity doesn't throw a viscosity error | property coupling *(new)* |
| 13 | RES critical enhancement suppressed outside near-critical | `ρ/ρc≥2`, `T/Tc>1.4` *(new)* |
| 14 | RES viscosity residual params round-trip | get/set |
| 15 | RES conductivity residual params round-trip | get/set |
| 16 | HEOS mixture viscosity uses log-mean when RES not enabled | no silent takeover |
| 17 | HEOS pure fluid RES opt-in without params throws | fail-loud *(new)* |
| 18 | HEOS mixture RES opt-in without params throws | fail-loud *(new)* |
| 19 | PR backend without RES params throws NotImplementedError | fail-loud |

---

## 2. Why each sample-data test fails

### (a) Implementation is correct; the model misses experiment — 4 assertions
These fail **only** the `*_exp` check. The RES implementation reproduces the papers' own
published model output to ≤0.34%:

| fluid | vs reference | vs experiment |
|---|---|---|
| R21 (viscosity, 293 K) | +0.0015% | −18% |
| HCL (viscosity, 241 K) | +0.004% | +17% |
| DECANE+METHANE (viscosity, 290 K) | −0.0007% | **+75%** |
| HYDROGEN (viscosity, 150 K) | −0.34% | +25% |

Acceptable model limitation, not a code defect. DECANE+METHANE matching the model to 7
digits while sitting 75% from experiment warrants a look at that data row specifically.

### (b) Reference-data provenance — ARGON+NEON conductivity, −11.6%
Li 2024's SI (`code_SI.py:287`) computes the **mixture** dilute-gas λ₀ from REFPROP directly
(`PropsSI('L','T',T,'P',0.1,'REFPROP::…')`), while the C++ uses the polynomial + Wilke rule.
Measured gap between the two λ₀ routes: ARGON+NEON −11.6%, METHANE+NITROGEN −5.1%,
CO2+NITROGEN −2.3%, BUTANE+METHANE −1.8%. The ARGON+NEON sample is near-ideal-gas
(ρ = 0.87 kg/m³), so λ₀ *is* the whole answer and the 5% tolerance cannot be met.
Pure-fluid conductivity is unaffected — the SI uses the polynomial there.

### (c) Modelling choice in the critical enhancement — PROPANE −4.9%, R143A +1.4%
`s_res` matches REFPROP to ~1e-3 % at both points, so the residual term is fine. The gap is
entirely the viscosity fed into λ_c (λ_c ∝ 1/η): the SI takes η from REFPROP, the C++ uses
`viscosity_RES()`, which carries no critical enhancement of its own and is therefore least
accurate exactly where λ_c matters.

Measured by isolating λ_c (28.3% of λ for PROPANE, 15.5% for R143A) and rescaling it:

| fluid | current (RES η) | with REFPROP η | with HEOS default η |
|---|---|---|---|
| PROPANE | −4.885% | **+0.001%** | −1.004% |
| R143A | +1.432% | **−0.007%** | −1.087% |

REFPROP's η reproduces the reference essentially exactly, confirming the rest of the
implementation. Decision taken: keep RES self-consistency and accept a few-% near-critical
cost rather than depend on a foreign viscosity. R143A is inside the 3% tolerance and passes;
only PROPANE fails.

### (d) PR backend vs experiment — the bulk of the 23 assertions
Two distinct groups:
* **Phase selection** — IHEXANE, PROPYNE, PXYLENE, BENZENE return ≈0.03 W/m/K where
  experiment is 0.08–0.10, i.e. gas-like values for liquid states.
  `transport_pr_pt_with_sat_fallback` is not recovering the right root.
* **Long-chain esters** — MLINOLEA, MLINOLEN, MOLEATE, MPALMITA, MSTEARAT all low by
  21–46%. A genuine cubic-EOS limitation, not a RES problem.

---

## 3. Design defects in the tests themselves

**D1 — Fail-open gate (blocking).** Every sample test wraps its body in
`catch (...) { WARN("Skip " ...) }` and gates on `REQUIRE(ok >= N)` with N of 4–10 against
~100 sample rows. A real regression degrades into a warning as long as N fluids still pass.
This is worse now that excluded fluids throw *by design* — those throws are
indistinguishable from breakage. Split into: expected-throw fluids (assert the throw),
and evaluable fluids (assert the value, no catch-all).

**D2 — No phase handling.** The HEOS tests take whatever phase `PT_INPUTS` returns. This is
what produced BENZENE −41.6% and R41 +266.8%; imposing the reference's phase collapses them
to −0.008% and −1.12%. Impose the phase, or skip straddling points explicitly — don't report
a phase disagreement as a model error.

**D3 — Unreachable tolerance.** Mixture conductivity at 5% cannot pass while the reference
uses a different dilute-gas model (§2b). Either loosen with a comment, restrict to points
where the residual term dominates, or drop the near-ideal-gas rows.

**D4 — Conflated criteria.** One test asserting both "matches the reference implementation"
(tight, a real regression signal) and "matches experiment" (loose, a model-quality
statement) means a failure doesn't say which broke. Separate them.

**D5 — Single sample point per fluid.** Only one published point exists per fluid, and some
are almost pure dilute gas — residual share 4.3% (R161 viscosity), −0.2% (R161
conductivity), 0.3% (VINYLCHLORIDE), 0.9% (R13), 0.2% (R1234YF). At those points the
transport property barely depends on `s_res`, so the test is blind to a systematic error.
Coverage of residual-dominated states is needed to make these tests meaningful.

---

## 4. Suggested changes, highest value first

1. **Fix D1** — replace the catch-all/`ok >= N` pattern with an explicit expected-throw list
   plus hard assertions on the rest. Without this the suite can silently stop testing.
2. **Fix D2** — impose the reference phase in the HEOS sample tests and in
   `transport_pr_pt_with_sat_fallback`. Removes the largest spurious failures in both the
   HEOS and PR groups.
3. **Fix D4** — separate reference-parity from experiment-agreement into distinct test cases.
4. **Fix D3** — restrict or re-tolerance the mixture conductivity test, with the provenance
   reason recorded in a comment so it isn't "fixed" by tightening later.
5. **Address D5** — add residual-dominated state points per fluid, checked against the C++
   implementation for regression rather than against an external reference.
6. Investigate the DECANE+METHANE experimental row (+75% while matching the model exactly).

## 5. Open item

Supporting RES on the REFPROP backend would let deviations be measured at *any* state rather
than only at the published sample points, which is currently the binding constraint on both
D5 and the exclusion decisions. Today the RES API lives on `HelmholtzEOSMixtureBackend`, so
REFPROP cannot evaluate the model at all.
