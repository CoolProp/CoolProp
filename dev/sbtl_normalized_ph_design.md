# SBTL: coordinate-aligned PH table — design

Status: **IMPLEMENTED** on `ihb/SBTL`.  This doc remains for the math /
rationale; the code lives in `src/Backends/Tabular/SBTLBackend.{h,cpp}`
under the `NormalizedPHTable` class plus the `_normph_liquid /
_normph_vapor / _normph_super` instances and `_coeffs_normph_*` arrays.
A sibling design captured for the PT path is in
`dev/sbtl_normalized_pt_design.md`.

## Goal

Replace the SBTL backend's uniform `(h, log P)` grid with a partition that puts
the saturation curve on a coordinate axis. Eliminates the cell-straddling
discontinuity errors that otherwise force the existing SBTL backend to perform
no better than BICUBIC on the #1301 window (median 2.9e-2 in the critical region
for both).

## Math

Three subdomains, each with its own grid in `(xnorm, log P)`:

**Region L (subcritical liquid)** — `P ∈ [P_min, P_crit)`, `h ∈ [h_min(P), h_sat,L(P)]`:

    xnorm = (h - h_min(P)) / (h_sat,L(P) - h_min(P))   ∈ [0, 1]

**Region V (subcritical vapor)** — `P ∈ [P_min, P_crit)`, `h ∈ [h_sat,V(P), h_max(P)]`:

    xnorm = (h - h_sat,V(P)) / (h_max(P) - h_sat,V(P))  ∈ [0, 1]

**Region S (supercritical)** — `P ∈ [P_crit, P_max]`, `h ∈ [h_min(P), h_max(P)]`:

    xnorm = (h - h_min(P)) / (h_max(P) - h_min(P))      ∈ [0, 1]

`h_min(P)` and `h_max(P)` are smooth bounds chosen at table-build time so the
sampled (T, ρ) values stay inside the EOS validity envelope. Practical choices:

- `h_min(P)`: enthalpy at the highest density that is in-range at this pressure
  (i.e. compressed liquid limit, or T_min isobar for low P).
- `h_max(P)`: enthalpy at the lowest density that is in-range (vapor /
  ideal-gas limit at T_max).

The inputs to the bicubic evaluation become `(xnorm, log P) ∈ [0,1] × log[P_min, P_max]`
which is rectangular and uniform — no special handling, no straddling.

## Subdomain dispatch

At lookup time, given `(h, P)`:

1. If the user has called `specify_phase()` with a specific single-phase
   value, route to the corresponding table (see "Phase specification" below).
2. Else if `P ≥ P_crit`: region SUPER.
3. Else compute `h_sat,L(P)` and `h_sat,V(P)` via the existing
   `saturation_hmolar_liquid(P)` / `saturation_hmolar_vapor(P)`.
   - `h ≤ h_sat,L`: region L
   - `h ≥ h_sat,V`: region V
   - else: two-phase. Q = (h - h_sat,L) / (h_sat,V - h_sat,L); rho via 1/v
     mixing rule, T = T_sat(P), other props by linear blend.

This dispatch is sub-µs (one comparison + at most two superanc/cache lookups).

## Phase specification (specify_phase)

Identical to the PT path's phase-specification contract — see
`dev/sbtl_normalized_pt_design.md` § "Phase specification" for the full
mapping table and example.  The PH case is slightly less common
(HmassP_INPUTS at h = h_sat,L or h = h_sat,V is rarer than PT at the
sat curve) but the routing follows the same imposed-phase →
LIQUID/VAPOR/SUPER table assignment.

## Data structure

```cpp
class NormalizedPHTable : public SinglePhaseGriddedTableData {
   public:
    enum Region { LIQUID = 0, VAPOR = 1, SUPER = 2 };
    NormalizedPHTable(Region r) : region_(r) { /* logx=false, logy=true */ }

    // xnorm-axis bounds always [0, 1].  P-axis bounds depend on region:
    //   LIQUID, VAPOR: [p_triple, p_crit]
    //   SUPER:         [p_crit, p_max]
    void set_limits() override;

    // Inverse of the normalization: given (xnorm, P), produce h.
    // Used at table-build time to feed HEOS update(HmolarP_INPUTS, ...).
    double h_from_xnorm(double xnorm, double P) const;

    // Forward: given (h, P), produce xnorm.  Used at lookup time.
    // Assumes (h, P) is inside this region.
    double xnorm_from_h(double h, double P) const;

   private:
    Region region_;
    // h_min and h_max as smooth functions of P, sampled on the same log_P
    // grid as the rest of the table.  Cubic-interp between samples.
    std::vector<double> h_min_, h_max_;
};
```

## Lookup path (additions to `SBTLBackend::update`)

For `HmolarP_INPUTS`:

```cpp
const double p = val2_;
const double h = val1_;
double xnorm;
NormalizedPHTable* tbl;

if (p >= p_crit) {
    tbl = &nph_super_;
    xnorm = nph_super_.xnorm_from_h(h, p);
} else {
    const double hL = saturation_hmolar_liquid(p);
    const double hV = saturation_hmolar_vapor(p);
    if (h <= hL) {
        tbl = &nph_liquid_;
        xnorm = nph_liquid_.xnorm_from_h(h, p);
    } else if (h >= hV) {
        tbl = &nph_vapor_;
        xnorm = nph_vapor_.xnorm_from_h(h, p);
    } else {
        // Two-phase.  Compute Q, blend.  No table eval.
        ...
    }
}

std::size_t i, j;
find_native_nearest_good_indices(*tbl, sbtl_coeffs_for(*tbl, ...), xnorm, p, i, j);
_T = evaluate_single_phase(*tbl, ..., iT, xnorm, p, i, j);
_rhomolar = evaluate_single_phase(*tbl, ..., iDmolar, xnorm, p, i, j);
// etc.
```

The existing `find_native_nearest_good_indices` and `evaluate_single_phase`
machinery work unchanged because the table is rectangular in `(xnorm, log P)`.

## Build path

Currently each `SinglePhaseGriddedTableData` derived class implements `set_limits()`
and the base class fills the grid by calling `AS->update(HmolarP_INPUTS, h, p)`
at each grid point (or whatever native input matches). For `NormalizedPHTable`,
override the build to:

1. Compute `h_min(P)` and `h_max(P)` at each P sample.
2. For each `(i, j)` cell, compute `h = h_from_xnorm(xvec[i], yvec[j])`.
3. Call `AS->update(HmolarP_INPUTS, h, P)` and read off T, rho, s, u, etc.

This is the same number of HEOS calls as the existing PH table — just at
different (h, P) points. No build-time regression.

## Storage

Add three new tables to `TabularDataSet` (or hold them as members of
`SBTLBackend`): `single_phase_normph_liquid`, `single_phase_normph_vapor`,
`single_phase_normph_super`. Cache to `~/.CoolProp/Tables/<fluid>/normph_*.bin.z`.

The h_min/h_max bounds need to be persisted alongside the table so that
`xnorm_from_h` produces the same values at lookup as at build. Simplest:
add them to `LIST_OF_SATURATION_VECTORS` (currently has rhoL, hL, etc. — add
hminL_isobar, hmaxV_isobar, etc.). Or store them as additional fields in
`PackablePhaseEnvelopeData`.

## Validation

Re-run the #1301 repro window. Expected results (estimated, to confirm):

- Sub-liquid + sub-vapor regions: median error ≤ 1e-6 (limited by bicubic
  truncation in the smooth coordinate, not by saturation-curve straddling).
- Supercritical region: median error ~1e-3 (still limited by the cell-too-wide
  problem at the critical singularity — same as BICUBIC and unchanged by
  this redesign).
- 99th-percentile dominated by the supercritical region; expect ~5–10 %.

Coordinate alignment fixes the *discontinuity* errors. The supercritical-cusp
errors require *additionally* either:
- Higher table resolution (`TABULAR_NX/NY` knob, already shipped on master)
- Smoother choice of interpolant (e.g. residual against a cubic EOS)
- Adaptive grid refinement near Tc

## Implementation status (post-landing)

All of the below shipped in the PR; cross-referenced for future readers.

1. **`NormalizedPHTable` class** with `set_limits`, `h_from_xnorm`, `xnorm_from_h`,
   `h_lo_isobar`, `h_hi_isobar`, plus `Cheb1D h_lo_cheb / h_hi_cheb` for fast
   xnorm-from-h at lookup.  `SBTLBackend.h` ~80 lines.
2. **Build path**: `populate_normph_bounds()` and `build_normph_table()` in
   `SBTLBackend.cpp` sample HEOS at the normalized grid (using `PQ_INPUTS` at
   sat boundaries, `HmolarP_INPUTS` at interior cells).  Cell coefficients
   built via `build_bspline_coeffs` (2D natural cubic B-spline, Kunick C²).
   ~250 lines.
3. **Serialization**: not implemented.  Tables rebuild from HEOS at every
   construction (~12 s/fluid one-shot).  Future work — see "Outstanding
   build-time" item in `dev/sbtl_pt_outstanding_work.md`.
4. **Lookup path**: subdomain dispatch + `xnorm_from_h(h, p, h_sat)` (with
   the H-superancillary supplying the sat-side bound for machine precision)
   + `evaluate_single_phase` against the active normalized table.
   `SBTLBackend.cpp::update` ~150 lines.
5. **Two-phase queries**: SBTL does not support inside-dome `HmassP_INPUTS`
   queries — those throw `ValueError`, callers should use `PQ_INPUTS` or
   `specify_phase` to disambiguate.
6. **Critical-region HEOS-fallback box** (additional, not in original plan):
   small box around (h_crit, p_crit) routes to HEOS direct because no smooth
   polynomial backbone can reproduce the cusp in ρ(h, P) at the critical
   point.  Box dimensions ±0.1 MPa × ±30% h_crit.

## Risks

- `h_min(P)` / `h_max(P)` boundaries can be subtle near the critical point —
  the (T, ρ) envelope shrinks asymmetrically, and overly aggressive bounds
  produce HEOS validity errors at table-build time. Solution: allow the
  build to detect failure and back off bounds, recording the trimmed range
  per-P-row.
- The supercritical region's xnorm boundaries (`h_min(P)`, `h_max(P)`) are
  not naturally bounded — `h_min` could be the T_min isotherm OR the
  high-density limit, depending on P. For the supercritical table, just
  use `h ∈ [h(T_min, P), h(T_max, P)]`.
- Two-phase blend: properties like cp diverge at the saturation curve (ideal
  cp is bounded but the latent heat term contributes infinity to cp).
  Returning the bounded mixing-rule values is correct for hmolar, smolar,
  rhomolar; cp/cv require special handling (and so does BICUBIC currently).

## Why this is the right design

- **Saturation discontinuity vanishes by construction.** Each cell is
  guaranteed to be on a single side of the curve. No more cell-bumping,
  no more `set_alternate2`-for-vapor-side logic — the table's geometry
  encodes the phase information.
- **No new evaluation primitives.** Bicubic Hermite / 4-point Lagrange
  on `(xnorm, log P)` rectangles uses the existing machinery.
- **Per-fluid tunability.** `h_min`/`h_max` per-region let us push close
  to the EOS limits without globally sacrificing accuracy.
- **Backward compatible.** Existing `(h, log P)` and `(T, log P)` tables
  stay around as fallback. Failure of the normalized lookup (e.g. due to
  table file format mismatch) routes back to the original path.
