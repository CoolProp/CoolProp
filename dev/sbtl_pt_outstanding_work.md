# SBTL — outstanding work after the PT/PH-focused PR

**STATUS**: PT and PH paths are landing-ready for pure fluids.  This
doc lists what's left before the PR can claim broader coverage and
queues the D,U follow-up.  Design rationale for the implemented PT
path lives in `dev/sbtl_normalized_pt_design.md`; same for PH in
`dev/sbtl_normalized_ph_design.md`.

## Scope (this PR)

SBTL is **pure fluids only**, with three input pairs:
`PT_INPUTS`, `HmolarP_INPUTS`, `HmassP_INPUTS`.  Mixtures (real or
pseudo-pure) are rejected at construction; other input pairs throw
cleanly.  Derivative and transport properties route through HEOS
direct after the table-driven update primes the underlying state.

This is the minimum surface area that exercises the coordinate-aligned
normalized PH and PT backbones end-to-end.  Broader coverage (D,U flash;
mixture support; tabular transport) is queued for follow-up PRs.

## Done in this PR

| Path | Coverage | Median | Max | Speedup vs HEOS |
|------|----------|--------|-----|-----------------|
| PT pure fluid | LIQUID + VAPOR + SUPER + critbox | 1e-8 (CO₂), 1.5e-8 (Water) | 1.7e-3 (CO₂), 3.8e-4 (Water) | 8.4× |
| PH pure fluid | LIQUID + VAPOR + SUPER + critbox | (per existing PH validation) | (existing) | (existing) |

The 5 % PT compressed-liquid residue and the triple-corner hang flagged
in the previous version of this doc are both gone — fixed by the
`NormalizedPTTable` (cell boundary aligned with the sat curve) plus the
critical-region HEOS-fallback box.

## Remaining before merge

### 1. Mixtures (follow-up)

Currently rejected at SBTLBackend construction.  Adding mixture support
would replace `saturation_T_LV` (a single `PQ_INPUTS` call) with a
phase-envelope query for bubble/dew-line `T_L(P)` / `T_V(P)`, then route
mixture PT through the existing LIQUID / VAPOR table machinery.  The
SUPER table works as-is.  Estimated 1 day.

### 2. Cryogens / narrow phase diagrams — DONE

Argon sweep ran against the multi-fluid accuracy plots and surfaced two
real bugs that were silently corrupting subcritical PH coverage:

1. **Melting-line ancillary lower-p bound > p_triple** (Argon, Helium).
   `AS->update(PT_INPUTS, p_triple, T)` throws because the melting-line
   polynomial's domain starts just above p_triple.  `safe_h_at_PT`
   swallowed the throw and returned `_HUGE`, which polluted one
   Chebyshev-Lobatto sample point in the h_hi Cheb fit — Clenshaw eval
   at intermediate p then returned NaN through huge-minus-huge
   cancellation.  Every subcritical Argon HmassP query silently fell
   through to a misleading "input pair not supported" terminal throw.
   Fix: `NormalizedPHTable::set_limits` now walks `ymin` upward in 5 %
   steps for LIQUID/VAPOR until HEOS accepts both T_min and T_max
   corner probes.

2. **Probe walk drove SUPER ymin far above p_crit** for the same
   cryogens.  At `(p_crit, T_min)` the fluid is below the melting line,
   so the corner probe always failed and the walk pushed ymin to
   ~4× p_crit before the trial cap stopped it.  This created a missing
   [p_crit, ymin_walked] band in supercritical Argon and CO2 coverage.
   Fix: skip the probe walk for SUPER region entirely (its ymin =
   p_crit is well-defined and doesn't need probing).

The HmassP error path was also tightened: internal failures now re-throw
with `(h, P, xnorm, p_range)` context instead of falling through to the
terminal "input pair not supported" misleading message.

Argon now covers 52 % subcritical points in random sampling (was 0 %);
the multi-fluid PH/PT accuracy plots in `Web/coolprop/SBTL.rst` show
all five fluids (Argon, CO2, Water, R245fa, D6) with consistent
coverage.

### 3. PT accuracy figure (visual confirmation)

`dev/sbtl_acc_fig_pt.py` renders BICUBIC vs SBTL accuracy maps in
(T, P) coords.  Currently 200×200 = 40 000 points × 3 backends = ~25 min
per fluid; too slow for routine regeneration.  Two options:

  - Lower the default to 80×80 (~2 min) and ship the lower-resolution
    figure in the PR description.
  - Add caching keyed on (fluid, grid spec) so subsequent renders are
    instant.

Recommended: just lower the default.  Estimated 5 min.

## Defer-able (low-priority polish)

### 0. Near-critical PH error band — substantially closed in this PR (follow-up for the last decade)

**Status: most of the band is gone.**  Subcritical LIQUID/VAPOR normph
tables now use a two-zone log-uniform yvec (40 % rows in [p_triple,
0.5·p_crit], 60 % rows in [0.5·p_crit, p_max]) instead of one
log-uniform span from p_triple to p_crit.  Cusp-region cell `log p`
spans shrink ~6×; the band's error drops 1000×–10⁵× (random PH at
p=0.93·p_crit, R245fa: 0.59 % → 10⁻⁵ %).

The bullets below describe the residual problem (max ~1–2 % errors
right at the immediate critical-point boundary, swallowed by the
existing critbox) and a follow-up plan to drive it down further.



The multi-fluid PH error plot shows a residual error band on the
subcritical vapor side of the SBTL panel for fluids close to the
critical region (most visible: R245fa, D6).  Not a coordinate or
normalization issue — the H-superancillary delivers `h_sat(p)` at
near-machine precision so xnorm at lookup is exact.

The error is in the *property surface* being interpolated.  Along the
dome boundary (xnorm = 0), `ρ_sat,V(p) ~ ρ_crit - C·(p_crit - p)^β` with
`β = 1/2` (mean-field — CoolProp's analytic Helmholtz EOS in (τ, δ) is
not Ising) has a vertical-tangent cusp as `p → p_crit`.  The cell's
16-coefficient Hermite bicubic in `(xnorm, log p)` cannot reproduce
this cusp — the residual at the cell midpoint is bounded by the cell's
`log p` span times derivatives of `ρ_sat` that diverge near critical.
Probe at p = 0.93 p_crit on R245fa: error peaks on the dome (0.59 % at
η = 0.001) and decays monotonically away (0.25 % at η = 0.95).

#### Plan: adaptive pressure grid (composes with Kunick mirror)

Independent of the property-surface transforms (Kunick), the *grid* in
`log p` can be made adaptive instead of log-uniform.  Concrete proposal:

1. **Seed `yvec` from the H-superancillary breakpoints.**  The
   superancillary is itself piecewise-polynomial in `log p`; its piece
   boundaries are already concentrated near `p_crit` where `h_sat(p)`
   changes character.  We get solved breakpoint placement for free
   instead of log-uniform spacing that ignores the cusp.

2. **Refine to an accuracy spec.**  For each interval
   `(yvec[j], yvec[j+1])` build the cells, then evaluate worst-case
   error at ~5 interior probe points: `p` at log-midpoint;
   `η ∈ {0, 0.25, 0.5, 0.75, 1}`.  If `|ρ_cell − ρ_HEOS| / ρ_HEOS >
   spec` (e.g. `1e-6`), insert `p_mid` into `yvec` and rebuild the
   affected cells.  Recurse until in-spec or depth cap reached
   (the critbox already swallows the immediate cusp neighbourhood
   so depth runaway isn't a concern in practice).

3. **Lookup cost trade-off.**  yvec was log-uniform, so finding the
   cell for a query was `O(1)`.  With adaptive yvec, find `j` via
   `std::lower_bound` on the cached `yvec_log` — `O(log N)`,
   ~5–10 ns for N≈200.  Hot path is currently ~1 µs so this is ~1 %
   overhead.

4. **Persistence.**  No schema change — the on-disk cache already
   serialises `yvec` via msgpack.  Build cost goes from ~25 s to
   ~35 s (maybe 30–50 % more LIQUID/VAPOR rows).

5. **Composes with Kunick mirror below.**  Adaptive p-grid handles
   the `log p` axis; Kunick's transforms handle the *property surface*
   on each cell.  Together: solved breakpoints (adaptive p) +
   no cell contains the cusp (Kunick L/G split at h_crit) +
   smooth surface to interpolate (per-property transforms) +
   HEOS fallback for the immediate cusp neighbourhood (existing critbox).
   Each chips away at a different axis; combined they should hit
   G13-15 spec.

Estimated effort: ~1 week for adaptive p + bisection lookup
infrastructure.  Should be done before or together with the
Kunick mirror so the test harness can attribute residual error.

#### Plan: mirror Kunick's IAPWS G13-15 method (follow-up PR)

Kunick's water tables already solve this problem and ship as the
reference IAPWS G13-15 implementation.  The key insight: don't
interpolate the cusped surface — interpolate a *transformed* surface
that's smooth, then invert the transform at lookup.  Concrete elements
to port to SBTL:

1. **Region split at p_crit and h_crit.**  Kunick's (h, p) plane is
   partitioned into ~6 fixed regions separated by *vertical* lines in p
   (with p_crit as one boundary) and *horizontal* lines in h
   (h_crit ≈ 2087.546 kJ/kg for water).  The critical point sits at a
   corner where 4 regions meet — it isn't *inside* any region.  No
   single cell ever straddles or contains the cusp.

   For SBTL this means promoting our (LIQUID, VAPOR, SUPER) split to a
   finer (LIQUID, VAPOR, SUPER_L, SUPER_G) split with the SUPER_L /
   SUPER_G boundary at the h_crit horizontal.  Compressed-liquid /
   supercritical-liquid extension goes to SUPER_L; vapor /
   supercritical-gas extension goes to SUPER_G.  Each gets its own
   bicubic.

2. **Property-specific coordinate transforms within each region.**
   Within each region, the function being interpolated is transformed
   so the residual of the smooth fit is bounded everywhere:
     - `v̄ = log v` for the gas / supercritical-gas region (kills the
       large dynamic range)
     - `s̄ = √s` near the cusp (absorbs the β=1/2 behaviour)
     - per-property transforms tuned to flatten the third-derivative
       curvature in the transformed coordinate
   The bicubic fits the transformed surface; lookup applies the inverse.
   For our near-critical PH band specifically, fitting
   `1 - ρ/ρ_crit` on the SUPER table with a `(p_crit - p)^0.5`
   transform on the p coordinate of the relevant cells should absorb
   the cusp entirely.

3. **Variable grid density near p_crit.**  Kunick's Tables A8–A10 pack
   hundreds of nodes into the immediate critical-point neighbourhood
   and use sparse spacing elsewhere.  Even with the transform, local
   curvature gets brute-forced down by cell size.  For SBTL this is
   the same as mitigation path 1 from the earlier write-up: piecewise
   yvec with `(p_crit - p)^0.5` spacing in the top 15 % of pressure.

4. **HEOS-fallback box stays.**  The existing critbox handles the
   immediate ±1 bar / ±30 % h_crit neighbourhood that's too cusped for
   any polynomial — works fine, no change.

5. **Property-specific transforms persisted with the table.**  Per
   (region, property) the table stores a transform spec (identity,
   `log`, `√`, `1 - x/x_crit`, etc.) chosen at build time to flatten
   the residual.  Evaluator applies the inverse before returning.
   Storage cost: one tag per (region, property), negligible compared
   to the cell-coefficient payload.

Estimated effort: ~2 weeks for the four region split + transform
infrastructure + per-property build-time transform selection.  This
is also items 1–6 of the existing "IAPWS G13-15 conformance" section
below; landing those would close the cusp band and bring the SBTL
water table to G13-15-grade accuracy simultaneously.

#### Smaller mitigations (still applicable if the full Kunick port slips)

1. **Finer cells in the top 15 % of p.** Replace log-uniform yvec
   spacing with a piecewise refinement that concentrates rows in the
   `[0.85 p_crit, 0.999 p_crit]` band (e.g. spacing in
   `(p_crit - p)^0.5` for that piece).  Likely brings near-critical
   errors down to BICUBIC levels at modest storage cost.  Estimated 1 day.

2. **Non-polynomial correction at xnorm = 0.** Augment the cell's
   Hermite with an explicit `(p_crit - p)^β` term anchored at the
   dome boundary, fit at table-build time.  More invasive (requires
   per-region evaluator changes) but eliminates the cusp residual at
   its source.  Estimated 2-3 days.

PT is unaffected by this mechanism because `T_sat(p)` is smooth where
`ρ_sat(p)` is not — the PT normalization makes T the dome coordinate,
not the density.  Supercritical region is unaffected (no dome, no
cusp).

### A. Supercritical-shoulder corner outliers

Five out of 1 486 random CO₂ points still show err > 1e-3 (max 1.7e-3).
All sit at the supercritical-shoulder edge of the critbox.  Fix is a
Kunick-style L/G split of SUPER along the `T_c` isotherm so each piece
is single-cusp-free.  Estimated 0.5 day.  Defer until someone hits
these in practice.

### B. `evaluate_single_phase_derivative` and `invert_single_phase_x/y`

These throw `ValueError` in the pure-fluid PT/PH scope: derivatives go
through HEOS direct, inversion is unused (no non-PT/PH input pairs).
If broader input-pair support is added, these need real implementations
keyed to the normalized table's `(xnorm, log P)` convention.

### C. Build-time reduction — DONE

On-disk caching for the six normpt/normph regions landed in commit
`5c02e70df`.  msgpack-serialised + zlib-compressed; ~30 MB per file,
six files per fluid (sbtl_norm{ph,pt}_{liquid,vapor,super}.bin.z).
Cold build ~25 s, warm load ~1.8 s, ~14× speedup.  Cache invalidates
on dimension mismatch (TABULAR_NX/NY changed) or msgpack-deserialise
failure.

### D. PT supercritical-only test

The supercritical-far-from-critical region (T > 1.3 × T_c, P > 1.75 × p_c)
shows ~1e-5 errors via the spline backbone — same as BICUBIC.  Worth a
dedicated regression test asserting the spline path doesn't degrade
this in future commits.

## IAPWS G13-15 conformance (follow-up PR)

To propose CoolProp's SBTL backend to IAPWS as a G13-15-conformant
implementation, the present configuration needs substantial work.
A 5 000-point random water sweep against the IAPWS-95 region bounds
from Table A10 of G13-15 shows median / max deviations 4–6 orders of
magnitude past the G13-15 permissible values (10–25 mK in T,
0.001 % in v, 1e-6 kJ/(kg·K) in s, 0.001 % in w / η).  Concrete
items to close the gap:

1. **Hermite bicubic on the normph backbone.**  Currently normph
   tables use cubic B-spline from corner values only.  Switching to
   Hermite bicubic with EOS-supplied first and cross derivatives
   would directly close the SUPER-low-T gap and lift v/T/s accuracy
   by ~3 orders of magnitude.  Estimated 2-3 days.

2. **Multi-piece grids per region.**  Kunick's Tables A8–A10 specify
   piecewise grids with up to ~400 nodes in p and different
   transformations per piece (`p̄ = p`, `p̄ = √p`, `p̄ = ln p`).
   CoolProp's current SBTL uses uniform 200 × 200 in log p × xnorm.
   Build piecewise grids per (region, property) matching Tables A8-A10
   exactly.  Estimated 2-3 days.

3. **Permit EOS extrapolation at table build.**  Kunick's L region
   extends to p = 1100 MPa; CoolProp's water HEOS is valid to 1000 MPa.
   Allow the SBTL builder to call the EOS in its extrapolation regime
   (with a documented validity envelope on the resulting table) so the
   table can cover the full G13-15 bounds.  Estimated 0.5 day.

4. **Property-specific transformations.**  Kunick uses `s̄ = √s` for
   entropy in some regions, `v̄ = log v` for the gas region, etc.  Each
   spline function has its own transformation tuned to reduce
   third-derivative curvature in the transformed coordinate.  Mirror
   these per-property transformations.  Estimated 1 day.

5. **L/G split by h_c at supercritical p.**  Kunick splits L vs G by
   h_c = 2087.546845 kJ/kg (water) for all pressures.  CoolProp
   currently splits by sat-curve at p < p_crit and routes everything
   else through SUPER.  Add an h_c-aware split for p ≥ p_crit so L
   (compressed liquid + supercritical liquid extension) and G (vapor +
   supercritical gas extension) each get their own table.  Estimated
   1 day.

6. **Separate HT region.**  Kunick separates G (T < 1073.15 K) from HT
   (1073.15 < T ≤ 2273.15 K) with different p bounds (HT capped at
   50 MPa) and different grid densities.  Estimated 0.5 day.

7. **Sphinx accuracy page.**  Build a per-region deviation table page
   matching the IAPWS Tables 8–13 layout, computed at doc-build time
   and shown to confirm conformance.  Already drafted in a prior PR
   iteration (commit 2438ca564, since reverted) — drop back in after
   items 1–6 land.  Estimated 0.5 day.

Total: ~2 weeks of focused work for a release-quality IAPWS-conformant
implementation.

## D,U inputs (rescoped to follow-up PR)

D,U code path was removed — `flash_DmolarUmolar`, `flash_DmolarT`, the
DT/DU table classes, and the saturation property cache that fed them
are deleted from SBTL.  D,U queries fall through to the base class
`TabularBackend` flash, which is correct but slow (HEOS-iterative)
for SBTL.

The Kunick-style D,U redesign is documented in
`dev/sbtl_du_kunick_redesign_plan.md` and queued for a follow-up PR
once the PT redesign here is proven.  The (v, u) region diagram lives
in `dev/sbtl_du_kunick_diagram.py`.
