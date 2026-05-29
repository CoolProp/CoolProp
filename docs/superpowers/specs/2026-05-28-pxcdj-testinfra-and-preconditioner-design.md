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
