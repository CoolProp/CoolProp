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

### 2. Cryogens / narrow phase diagrams (verify)

`T_c × 1.30` HEOS-fallback box may be meaningless for fluids whose entire
EOS range is small (helium: T_c ≈ 5.2 K, T_max ≈ 2000 K — box would
extend to T = 6.7 K, all fine; but T_min is at least 2.2 K, which is
0.42 × T_c — well below the box).  No actual issue, but worth running
the random PT sweep against a cryogen to confirm.  Estimated 30 min.

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

### C. Build-time reduction

normpt + normph tables together take ~25 s per fluid first-build.  The
existing tabular caching machinery (msgpack of LIST_OF_MATRICES) doesn't
yet cover the normpt / normph data.  Adding it would drop subsequent
loads to <1 s.  Estimated 1 day.

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
