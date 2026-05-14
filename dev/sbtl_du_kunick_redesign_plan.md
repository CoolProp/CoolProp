# SBTL D,U Tables — Kunick Region Redesign Plan

Reference: IAPWS Guideline G13-15 §4 (Spline Functions of (v, u) and Inverse
Functions Based on IAPWS-IF97), figure rendered for any fluid by
`dev/sbtl_du_kunick_diagram.py`.

## Status

D,U code path was **removed** from the SBTL backend in the
PT/PH-focused PR (commit ~`7af349359` on `ihb/SBTL`).  D,U queries now
fall through to the base `TabularBackend` flash, which is correct but
slow (HEOS-iterative).  This doc captures the design for the planned
D,U follow-up PR — picks up where the prior implementation left off.

The previous prototype (commit `94d893941` on the same branch, before
the cleanup) had two rectangular DU tables joined at `D = D_crit` plus
a sibling pair of DT tables.  Water single-phase accuracy was med=3.0e-4
but with a long tail (max 3.5e+2) dominated by dome-edge mis-classification.
The Kunick partition described below should fix that.

## Target (Kunick exact mirror)

Three single-phase regions:

| Region | Bounds                                                       | Coord transform                                                                                              |
| ------ | ------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------ |
| L      | sat-L → p=p_max_L isobar; bounded above by u_c on supercrit | **Normalized v̄(v, u) = (v̄_max − v̄_min) · (v − v_min(u)) / (v_max(u) − v_min(u)) + v̄_min**, v̄_min=1, v̄_max=100 |
| G      | sat-V → p=p_min isobar; bounded below by u_c, above by T_HT  | log v                                                                                                        |
| HT     | T > T_HT, p ≤ p_max_HT                                       | log v                                                                                                        |
| TP     | inside the dome                                              | sat-L/sat-V interpolation (existing path)                                                                    |

Boundary functions of u (for Region L):
  - v_min(u) = v(p = p_max_L, u)              — Cheb1D in u, smooth
  - v_max(u) = v_sat,L(u)                      — direct from rho-superancillary
                                                 (machine-precision, like h_sat
                                                 for PH)

Generalization parameters (water → general):
  - u_c       = u(D_crit, T_crit)              — exact for any fluid
  - T_HT      = T_crit · (1073.15 / 647.096)   — fixed ratio, ~1.658
  - p_max_L   = min(p_max_eos, 4.5 · p_crit)   — water is 100 MPa
  - p_min     = p_triple                       — water is 611.212 Pa

## Implementation steps

### Step 1: `NormalizedDUTable` class (analogous to `NormalizedPHTable`)
Shape:
```cpp
class NormalizedDUTable : public SinglePhaseGriddedTableData {
    enum Region { L, G, HT };
    Region region_;
    Cheb1D v_min_cheb;        // L only: v(p_max_L, u) as a function of u
    // v_max(u) = v_sat,L(u) read from rho-superancillary at lookup
    double xbar_from_v(double v, double u) const;
    double v_from_xbar(double xbar, double u) const;
};
```

Cell-fill iterates over (xbar, u) grid, computes v_grid via the boundary
functions, and queries `AS->update(DmolarUmolar_INPUTS, 1/v_grid, u_grid)`
to populate T, P, S, etc.  Works because the grid is now aligned with the
sat-L curve (xbar=1 ↔ sat-L exactly), so no cells straddle the dome.

### Step 2: Region G as log-v rectangle
Same Hermite-bicubic backbone as current `DUGasTable`, but bounds change:
  - x: log v from log(v_sat,V(u)) to log(v(p_min, u))   — actually rectangular
       in log v over [log v_min, log v_max] of the entire G region; the
       sat-V curve becomes an interior line, but cells crossing the sat-V
       boundary are handled by the new sat-aware cell filter (skip cells
       fully outside G; partial cells use only the G-side corners).

Alternative (cleaner): also use a normalized-v̄ on the G side, with v̄_min = v_sat,V(u),
v̄_max = v(p_min, u).  Same trick as L.  Decide once we know how Kunick's
table sizes compare to current.

### Step 3: Region HT (separate table for T > T_HT)
Rectangular in (log v, u), bounds:
  - x: log v from v at (T_HT, p_max_HT) to v at (T_max_eos, p_min)
  - y: u from u(v_min(T_HT), T_HT) to u(v_min, T_max_eos)

Hermite-bicubic.

### Step 4: Routing in `flash_DmolarUmolar`

```text
state = (D, U);
v = 1 / D;

if (state in dome via 2-phase Newton):
    use TP path (existing, with sat superancillaries)
else if (U > u_c and T_estimated > T_HT):
    use HT
else if (U > u_c):
    use G
else:
    use L  // includes both subcritical liquid and supercritical compressed-liquid
```

T_estimated comes from a quick lookup in the chosen table; if it lands
in the wrong region we re-route once.

### Step 5: Remove the legacy DULiquidTable / DUGasTable
**Already done in the PT/PH PR.**  The old `DULiquidTable`,
`DUGasTable`, `DTLiquidTable`, `DTGasTable` classes plus
`flash_DmolarUmolar`, `flash_DmolarT`, `build_du_tables`,
`build_dt_tables`, the saturation property cache that fed them, and
~830 lines of related code have been removed.  Re-introduce only the
new `NormalizedDUTable` machinery.

## Per-fluid validation matrix

After each step, regenerate the per-fluid (h, P) accuracy figure +
a new (D, U) accuracy figure, and confirm:
  - Water: max < 1e-3 across the full grid (current: 3.5e+2)
  - CO2:   max < 1e-3 (current acceptable but not measured systematically)
  - Helium: max < 1e-2 (the narrow phase diagram makes this hard;
            may need a finer grid in the L region)

## Estimated effort

  - Step 1 (NormalizedDUTable + L): 1 day
  - Steps 2-3 (G, HT separation): 1 day
  - Step 4 (routing): half day
  - Step 5 (cleanup) + validation: half day

Total: ~3 days of focused work for a release-quality D,U path.
