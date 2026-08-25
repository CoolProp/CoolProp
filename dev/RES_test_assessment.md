# RES transport test suite — assessment

**Status: the overhaul this document called for is done, and it ended somewhere the original
analysis did not anticipate.** D1–D5 are addressed, and the sample-data comparisons that D1–D4
were about have been **removed** rather than repaired. The suite is green:

    test cases:  32 |  32 passed |  0 failed
    assertions: 637 | 637 passed | 0 failed

The original assessment is kept from "Historical assessment" onwards: its per-fluid causes are
still the record of why each deviation exists.

---

## Why the sample comparisons are gone

They were rebuilt first — bucketed, phase-free, parity split from experiment — and that version
passed. It was still the wrong test to keep. Every remaining deviation had a known, explainable
cause, so each needed a paragraph of comment to stop a reader concluding the implementation was
poor; and accommodating those causes required percent-level tolerances on quantities that agree to
1e-7 elsewhere, which is close to no regression signal at all.

What they were standing in for is now measured properly, and better:

| question | where it is answered now | agreement |
|---|---|---|
| is the model implemented correctly? | `dev/RES_reference_check.py` vs Martinek/Li, and the `[RES][REFPROP]` parity tests | median 2.6e-7 / 1.2e-7 |
| do the parameters transfer to HEOS? | `dev/RES_grid_report.py`, 3828 grid points | median 6e-7 |
| does the model agree with experiment? | the published tables themselves; it is a property of RES, not of this code | — |

A unit test asserting "HEOS is within 13% of the reference" adds nothing to a script that measures
6e-7 over 3828 points, and it reads as though HEOS were the problem when the 13% is a documented
dilute-source difference.

## What replaced them

**`RES regression at residual-dominated states`** — 30 pinned values over 10 fluids, on **both**
HEOS and PR, at states where the residual term supplies 55–96% of the property. Golden values from
this implementation; their job is to fail when the model moves. Tolerance 1e-7, against a real
change of 7.5e-5 when the R_D/γ fix landed. This also closes **D5**: the published points give one
state per fluid and several are almost pure dilute gas (residual share 0.2% for R1234YF), where the
dilute polynomial is identical on every backend and the residual term is not exercised at all. It
is also the only RES coverage PR has, so the PR values are pinned separately.

Both backends are evaluated at **(p, T)**, each finding its own density. An earlier draft pinned PR
at HEOS's density, which tests PR outside the frame its coefficients were fitted in — PR's density
is up to 25% from HEOS's at these states, and forcing it onto HEOS's value inflated the apparent
PR-vs-HEOS difference to tens of percent. At (p, T) it is what one would expect of a cubic:

| PR vs HEOS at the same (p, T) | median | p90 | max |
|---|---|---|---|
| viscosity | 3.1% | 10.8% | 13.8% |
| conductivity | 2.7% | 7.9% | 21.8% |
| density | 7.9% | 18.3% | 25.1% |

**`RES withheld fluids report the omission rather than guessing`** — asserts that each of the 14
fluids without a HEOS entry throws, and that a fluid which is *not* withheld still evaluates, so
"everything throws" cannot pass. This is the part of **D1** worth keeping: under the old suite a
designed throw and a broken one were indistinguishable to a `catch (...)`. The list is repeated in
the test on purpose — it is the user-visible contract, so changing it should require editing a test
and justifying it, not just regenerating a data file.

**D2, D3, D4** are moot once the comparisons are gone: no flash, no phase, no conflated criteria,
no unreachable tolerance.

## Two defects found on the way, neither in the original list

**The PR helper was circular.** `transport_pr_pt_with_sat_fallback` evaluated up to three candidate
states and returned whichever came closest to the **experimental value** — choosing the answer that
best matches the expected answer. It is also how the old test could claim "within 15%" per point:
with an honest rule based on density instead, the 90th percentile is 15% and the worst point 57%.
Removed with the tests that used it.

**RES could return a silent NaN.** `s_plus = -s_res/R` is raised to fractional powers throughout,
so a non-positive or non-finite value yields NaN rather than an error, and it left the routine as a
plausible-looking transport value. The PR backend does produce it — a NaN residual entropy for
hydrogen forced onto its liquid root at 150 K / 112 MPa. `check_s_plus()` in
`src/Backends/RES/RESTransport.cpp` now refuses it. The old catch-all had been hiding that fluid
entirely, which is exactly what D1 was about.

## Open item (was §5) — resolved

RES on the REFPROP backend, which §5 identified as the binding constraint on D5 and on the
exclusion decisions, is done. The exclusion lists now come from a 3828-point grid rather than from
single published points; see the Stage 6 entry in `dev/RES_REFPROP_plan.md`.

Still open: the DECANE+METHANE experimental row (+75% while matching the model exactly) remains a
data-provenance question. It no longer fails anything, because nothing asserts against experiment.

---

# Historical assessment

The remainder is the original analysis, retained for the per-fluid causes.

# Historical assessment

The remainder of this document is the original analysis, retained for the per-fluid causes.

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
