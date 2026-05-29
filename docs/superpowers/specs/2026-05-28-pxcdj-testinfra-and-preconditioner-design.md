# P+X test infrastructure + preconditioner (preparation for CoolProp-cdj)

- **Date:** 2026-05-28
- **Issue:** the bd issue filed against this branch (see `bd list --status=open`)
- **Parent effort:** CoolProp-cdj (direct-EOS inner Newton + basin-safety refactor)
- **Status:** Design — autonomous overnight execution

## Scope (what this delivers)

1. A **test bed** for P+X flash refactoring work: broad sweep over ~120 CoolProp pure
   fluids × dense (log-p, T) grid × PH/PS/PU, with watermap instrumentation, named
   regression points, structural tests, and a baseline-fingerprint fixture.
2. A **standalone preconditioner** — per-fluid melt-curve table (p-indexed) plus a
   regime-classified seed function — used to augment the *existing* cached path's
   warm-start ρ. Does NOT change the inner Newton.

The actual direct-EOS inner-Newton refactor (CoolProp-cdj) is the next phase, **out
of scope here.** This branch's role is to make that refactor fast and low-risk.

## Out of scope

- The direct-EOS inner Newton replacing `solver_rho_Tp` (CoolProp-cdj).
- Tolerance / per-fluid perf trade-offs (Ian's calls; surface as open questions).
- Mixture / pseudo-pure paths (the preconditioner is for `is_pure()` only;
  pseudo-pure passes through to the existing cached path unchanged).

## Gates (G1–G5)

| Gate | What it asserts | Tolerance | Applies to |
|---|---|---|---|
| **G1: Correctness (round-trip)** | Every PT-established single-phase state recovers (T, ρ) through P+Y | \|ΔT\|/T ≤ 1e-6 | both paths |
| **G2: Basin equivalence** | Direct path lands on the same root as cached path | \|Δρ\|/ρ ≤ 1e-9, relax to 1e-8 at the TOMS748 floor with rationale | CoolProp-cdj only — *not* required for the preconditioner |
| **G3: Zero-failure vs master** | Every point master succeeds, the path under test must succeed | exact | both paths |
| **G4: Convergence bound** | Inner Newton iters ≤ 5 with preconditioner seed | hard cap | CoolProp-cdj only |
| **G5: Performance regression guard** | Per-fluid mean µs/call ≤ master | bench-asserted | both paths |

This branch enforces **G1, G3, G5** for the preconditioner. G2, G4 are the CoolProp-cdj
phase's gates.

## Test infrastructure

### CI-default tests (`[PXcdj]`)
- The two original Nitrogen/PS supercritical-cold points (86.35 K/4.754 MPa, 90.20 K/7.768 MPa) — must still round-trip after #3020 lands.
- The prototype's known-hard cases: Nitrogen-PS at (121 K, 3.7 MPa), R134a-PH at (250 K, 300 kPa), Water near (277 K, 0.1–200 MPa) for the density anomaly band.
- A representative fluid sample (Water, CO₂, R134a, Propane, Nitrogen, MM) × gas+liquid+supercritical states × PH/PS/PU.
- Structural tests: lazy-init runs once, pseudo-pure passthrough, no-melt-curve fallback path activated, thread-safety of lazy init.

### Opt-in heavy sweeps (`[.]`-tagged)
- **`[PXcdj_sweep][.]`**: all CoolProp pure fluids (`get_global_param_string("fluids_list")`) × dense (log-p, T) single-phase grid (default 40×40, env-tunable). For each (T, p, pair): establish via PT, round-trip via P+Y, record success + (T_solved, ρ_solved, microseconds). Write per-fluid CSV.
- **`[PXcdj_watermap][.]`**: instrumented variant of the above that also records `inner_iter_count`, `basin_rejections`, `fallback_count`, `seed_method` (which preconditioner branch was used), per point.
- **`[PXcdj_bench][.]`**: per-fluid timing distribution comparing preconditioner-on vs preconditioner-off, mean/median/p95 µs/call. Flags any fluid where preconditioner is slower.

### Visualization
- `wrappers/Python/CoolProp/Plots/PXcdjTimingMap.py` (sibling of `TwoPropFlashTimingMap.py`, `HSFlashTimingMap.py`) — renders the watermap CSV as p–T and p–s/p–h heatmaps with overlays (saturation dome, melt curve, critical isobar, supercritical isothermal). Matplotlib output; static PNG/PDF.

### Baseline fixture
- A small CSV committed under `dev/fixtures/pxcdj_master_baseline.csv.gz` (or similar) containing the per-fluid success/failure tally + (T_solved, ρ_solved) statistics at a reference grid, generated against the current master. The CI-default `[PXcdj]` test asserts the current build still matches this fingerprint within tolerance — a regression guard against subtle behavior drift.

## Preconditioner

### Per-fluid melt-curve table (p-indexed, lazy)

```
struct MeltCurveTable {
    std::vector<double> p;          // monotone increasing, log-spaced
    std::vector<double> T_melt;     // T at each p_i
    std::vector<double> rho_melt;   // liquid-side density at each (T_melt, p)
    bool has_data;                  // false if fluid has no melting line
};
```

- Indexed by **pressure** because some fluids (notably water) have melt curves that
  are multi-valued in T but single-valued in p.
- Lazy-initialized once per backend on first `regime_classified_rho_seed` call.
- Generation: sample N ≈ 30 pressures from `p_triple` to `p_max_useful` log-spaced;
  for each, `T_melt_i = calc_melting_line(iT, iP, p_i)`, then
  `ρ_melt_i = update(PT_INPUTS, p_i, T_melt_i).rhomolar()`.
- Lookup: binary search on `p`, linear interpolation between adjacent triples.
- Fluids with `has_melting_line() == false`: `has_data = false`; the seed function
  uses the compressibility-extrapolation fallback.

### Regime-classified seed function

```
double regime_classified_rho_seed(HelmholtzEOSMixtureBackend& H, double T, double p) {
    if (T > Tc) {
        if (T > 0.8*Tmax || p < 0.01*pc) return p/(R*T);          // ideal-gas
        return blend(rho_c, p/(R*T));                             // crit-isochore ↔ ideal
    }
    double Tsat = superanc.get_T_from_p(p);
    if (T < Tsat) {                                                // liquid
        double rhoL = superanc.eval_sat(T, 'D', 0);
        if (H.melt_table.has_data) {
            double rhoM = H.melt_table.rho_at(p);
            double psi = (log(p) - log(p_sat(T))) / (log(p_max_table) - log(p_sat(T)));
            psi = clamp(psi, 0.0, 1.0);
            return rhoL + psi * (rhoM - rhoL);
        }
        return rhoL * (1 + (p - p_sat(T)) / kappa_sat_L);          // compressibility extrapolation
    }
    // vapor
    double rhoV = superanc.eval_sat(T, 'D', 1);
    if (T > 0.5*Tmax || p < 0.1*p_sat(T)) return p/(R*T);          // ideal-gas at high-T / low-p
    return blend(rhoV, p/(R*T));                                   // sat-vapor ↔ ideal on p/p_sat
}
```

Total cost ≈ 50–100 ns per call (superanc reads, log/exp, blend). Always
basin-classified.

### Integration

Wire the preconditioner into `HSU_P_flash_singlephase_Brent`:
- Before the TOMS748 / Halley dispatch, compute `rho_seed = regime_classified_rho_seed(HEOS, T_mid, p)` where `T_mid` is a representative T (e.g. midpoint of the bracket).
- Plumb `rho_seed` into the existing `solver_resid` so the warm path (`update_TP_guessrho`) starts with the preconditioner's seed when iter ≥ 2 (or earlier if confidence is high).
- This is **additive**: if the seed is bad, `update_TP_guessrho` falls back via the existing safety net (which itself falls back to cold `update(PT_INPUTS)`); no new failure modes introduced.

### Toggle

`PXCDJ_PRECONDITIONER` env var (static-const read once); default OFF this branch
so master behavior is preserved. The validation work decides whether to flip the
default ON before merge.

## Validation pass

Final pass before push:
1. Build Release, run `[PXcdj]` — must pass.
2. Run `[PXcdj_sweep]` with `PXCDJ_PRECONDITIONER=0` (master baseline) → generate baseline CSV.
3. Run `[PXcdj_sweep]` with `PXCDJ_PRECONDITIONER=1` (preconditioner on) → comparison CSV.
4. Diff: G1 (round-trip ≤ 1e-6) on both; G3 (no new failures); G5 (per-fluid timing ≤ baseline).
5. Render the watermap plots for spot-check.
6. `./dev/ci/preflight.sh`; clang-format clean; adversarial code review on the diff.
7. Push to a **draft** PR (`--draft`). Update bd issue with the wake-up status block.

## Wake-up status

The draft PR's description will include explicit sections:

- **What works** — gates G1/G3/G5 satisfied, with per-fluid summary stats.
- **What's open** — anything I flagged for Ian's call (e.g., a fluid where preconditioner regresses slightly, a tolerance question, an unexpected master-side bug surfaced by the sweep).
- **What's deferred** — the direct-EOS inner Newton (CoolProp-cdj), any G2/G4 work, any per-fluid opt-outs that need investigation.

## Definition of done

- Branch `ihb/pxcdj-testinfra-precond` pushed as a **draft** PR.
- `[PXcdj]` (CI default) passes; `[PXcdj_sweep]` heavy-but-passes locally.
- Preconditioner is env-gated (`PXCDJ_PRECONDITIONER`), default OFF.
- Baseline fixture committed; watermap plotter committed.
- bd issue updated with PR link + status.
- No master-behavior change when the toggle is OFF.

## Phase 3b findings (integration attempted, not shipped)

Two integration attempts were made; both failed the G3 hard gate. The
preconditioner unit + test bed are shipped; the integration is deferred to
CoolProp-cdj proper (with the inner-Newton replacement).

### Attempt 1 — seed at `T_mid`
- Strategy: compute `regime_classified_rho_seed(H, T_mid, p)` once at function entry; reuse on every Brent probe.
- Result: **10,762 new failures vs baseline.**
- Root cause: Brent probes `Tmin` and `Tmax` first (often very far from `T_mid`); a seed classified for `T_mid` lands the inner Newton in the wrong basin for the bracket-endpoint probes → wrong-root convergence → Brent bisects around a fake solution.

### Attempt 2 — per-call seed re-classification
- Strategy: compute `regime_classified_rho_seed(H, T_probe, p)` on every `solver_resid::call(T)`. Cost: ~50–200 ns × ~9 calls ≈ ~1–2 µs overhead.
- Result: **234 new failures vs baseline** (down from 10,762 — big improvement), but still violates G3 = 0.
- Net effect: **+2,204 newly-successful points** (preconditioner is a clear win overall), but the spec mandates *zero new* failures, not a net gain.
- Failure clustering: 234 failures across just 4 fluids in 2 regimes:
  - **Deep supercritical compressed** (p > ~3·p_c): R12 at p=2×10⁸ Pa, R123 at 2.4–7.6×10⁷, Methane at 4.6–10×10⁸. Classifier picks ideal-gas/critical blend; actual is compressed-liquid density at ~30,000 mol/m³.
  - **Near-T_c handoff** (|T − T_c|/T_c < 0.005): R152A at T ≈ T_c. The sub/supercritical branch switch at exactly T = T_c is discontinuous.

### Attempt 3 — engagement gate (exclude failure regimes)
- Strategy: same as Attempt 2 plus an engagement gate excluding the two failure regimes (`p > 3·p_c`, `|T − T_c|/T_c < 0.005`); fall back to baseline (cold update) otherwise.
- Result: **444 new failures vs baseline** — *worse* than Attempt 2.
- Aggregate ON/OFF timing ratio: **1.257** (~26% slower when engaged).
- Failure clustering: **95% of new failures are pseudo-pure mixtures** (R407C 107, R410A 97, R507A 101, R404A 91, SES36 28) in well-behaved subcritical and moderately-supercritical regimes — the engagement gate's "safe zone" was wrong because `is_pure() == true` returns true for pseudo-pure mixtures whose superancillaries behave differently than true pure-fluid superancillaries.

### Lessons learned

1. **`is_pure()` is not the right gate for the preconditioner.** Pseudo-pure mixtures pass `is_pure()` but their superancillary-derived ρ_sat,L/V values don't seed the cold-path inner Newton correctly. The classifier needs `is_pure() && !is_pseudopure()` or an explicit check that the superancillary behaves as expected.
2. **Brent's bracket-endpoint probes need basin-correct seeds at the actual probe T**, not at any single representative T. Per-call re-classification fixes this for most regimes but exposes (1).
3. **Per-call seed-compute overhead (~200 ns × ~9 outer iters) is not paid back** by faster convergence on most fluids — the cached path's SRK guess + Householder4 + retry is already well-tuned for the bulk regime, and replacing it with the preconditioner's seed adds work without reducing inner-iter count enough to win.
4. **The preconditioner architecture is not enough on its own.** The big speedup hypothesis (preconditioner → fewer inner iters AND cheaper per iter) requires *also* replacing the cached `solver_rho_Tp` + Householder4 path with the direct-EOS inner Newton (CoolProp-cdj proper). The preconditioner alone changes only the seed; the inner-Newton machinery is unchanged and dominates the time.
5. **Net G3-positive but G3-not-zero** is a real outcome of opt-in preconditioning. The spec's hard G3 = 0 is the right gate for a *default-on* change, but if the toggle stays OFF in production, a net-positive ON path is shippable IF the test infrastructure can prove the net gain on representative grids.

### Status

- **Phase 2a** (test bed + named regressions + baseline fixture): ✅ shipped, commit `1b4bc2272`.
- **Phase 3a** (standalone preconditioner unit + unit tests): ✅ shipped, commit `1306d4e24`.
- **Phase 3b** (production integration): ❌ deferred — three attempts failed G3; the architecture needs to be reconsidered in conjunction with the CoolProp-cdj inner-Newton replacement.
- **Plotter** (Python visualization): deferred (lower priority; sweep CSV is directly parseable).
- **Structural tests** (lazy-init, pseudo-pure passthrough, thread-safety): partially covered by the Phase 3a preconditioner tests; full structural suite deferred to integration phase.

### Recommendations for CoolProp-cdj (the next phase)

1. **Re-design the gate**: `is_pure() && !is_pseudopure()` AND the per-fluid melt table actually has data AND the superancillary is consistent (run a startup sanity check).
2. **Co-design the seed with the inner Newton**: the preconditioner's value is in giving a *basin-classified* seed for a *direct-EOS, cache-bypassing* inner Newton. The cached `solver_rho_Tp` doesn't benefit enough from the seed because its iteration cost is already dominated by `CachedElement` and virtual-dispatch overhead.
3. **Consider a stricter ship gate**: instead of G3 = 0, allow "G3 net-positive AND no fluid regresses by more than X%" for an opt-in toggle. This admits net-improvement scenarios that the hard gate rejects.
4. **Address the pseudo-pure mixture handling explicitly** — these are a real and common class (HFCs blends like R407C/R410A/R404A) that need either a separate seed strategy or explicit opt-out.
