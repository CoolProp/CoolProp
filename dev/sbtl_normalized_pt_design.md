# SBTL: coordinate-aligned PT table — design

Status: **IMPLEMENTED** on `ihb/SBTL`.  Sibling of
`dev/sbtl_normalized_ph_design.md` — same coordinate trick, with `iT` in
place of `iHmolar`.  Code lives in `src/Backends/Tabular/SBTLBackend.{h,cpp}`
under the `NormalizedPTTable` class plus the `_normpt_liquid /
_normpt_vapor / _normpt_super` instances and `_coeffs_normpt_*` arrays.

## Goal

Replace the SBTL backend's uniform `(T, log P)` lookups (which inherit
the legacy `logpT` table's cell-straddling problems near the saturation
curve) with a coordinate-aligned partition that puts the saturation
curve on a coordinate axis.  Eliminates the compressed-liquid 5 % PT
residue that was the centerpiece of `dev/sbtl_pt_outstanding_work.md`.

## Math

Three subdomains, each with its own grid in `(xnorm, log P)` (xnorm
linearly spaced, P log-spaced — same as `NormalizedPHTable`):

**Region L (subcritical liquid)** — `P ∈ [P_min, P_crit)`, `T ∈ [T_lo(P), T_sat,L(P)]`:

    xnorm = (T - T_lo(P)) / (T_sat,L(P) - T_lo(P))   ∈ [0, 1]

**Region V (subcritical vapor)** — `P ∈ [P_min, P_crit)`, `T ∈ [T_sat,V(P), T_hi(P)]`:

    xnorm = (T - T_sat,V(P)) / (T_hi(P) - T_sat,V(P))  ∈ [0, 1]

**Region S (supercritical)** — `P ∈ [P_crit, P_max]`, `T ∈ [T_lo(P), T_hi(P)]`:

    xnorm = (T - T_lo(P)) / (T_hi(P) - T_lo(P))      ∈ [0, 1]

Boundary functions:

  - `T_lo(P) = T_melt(P)` when the fluid has a melting ancillary, else
    `T_min` (= `max(T_triple, T_min_eos)`).  The Helmholtz EOS is valid
    down to the melting curve — including water's negative-slope ice I
    band where `T_melt < T_triple` between ~5 and ~200 MPa.
  - `T_hi(P) = T_max_eos` for VAPOR / SUPER (capped at the published EOS
    upper bound — no extrapolation beyond, which previously produced
    silent ρ errors of 5–50 % at xnorm near 1).  For LIQUID,
    `T_hi(P) = T_sat,L(P)` (the sat curve is the natural upper bound).

The non-saturation boundaries (LIQUID `T_lo`, VAPOR `T_hi`, SUPER both)
get a `Cheb1D` spectral fit at table-build time so `xnorm_from_T` at
lookup is fast and matches the build-time encoding to ≤1e-12.  The
saturation boundary `T_sat,L/V(P)` comes straight from `PQ_INPUTS`
(machine-precision via the IF97 / Helmholtz sat solver).

## Subdomain dispatch

At lookup time, given `(p, T)`:

1. If the user has called `specify_phase()` with a specific single phase,
   route directly to the corresponding table (see "Phase specification"
   below).
2. Else if `P ≥ P_crit`: region SUPER.
3. Else compute `T_sat,L(P) = T_sat,V(P) = T_sat(P)` via a single
   `PQ_INPUTS` HEOS call (cached single-deep in `_Tsat_LV_cache_*` — ~10 µs
   for the cache-miss path, free for the cache-hit).
   - `T ≤ T_sat,L`: region L (LIQUID), `T_sat_active = T_sat,L`
   - `T ≥ T_sat,V`: region V (VAPOR), `T_sat_active = T_sat,V`
   - else: throw `NotImplementedError` — PT_INPUTS inside the two-phase
     dome is ambiguous (this never fires for pure fluids where the dome
     is a 1-D locus, but the throw is the safety net for any case where
     the cache disagrees with the sat solver).

`T_sat_active` is passed into `xnorm_from_T` so the routing decision
and the xnorm computation use the *same* T_sat value (machine-precision
via the IF97/Helmholtz sat solver), avoiding boundary-disagreement
artefacts.

This dispatch is sub-µs after the sat cache warms up.

## Phase specification (specify_phase)

`PT_INPUTS` at exactly `T = T_sat(P)` is formally ambiguous: the query
sits on the saturation curve and could mean either saturated liquid or
saturated vapor.  The user disambiguates with `specify_phase()` *before*
calling `update()`.  The mapping from `imposed_phase_index` to the
SBTL routing decision is:

  | imposed_phase_index             | Action                                                       |
  | ------------------------------- | ------------------------------------------------------------ |
  | `iphase_not_imposed`            | Automatic: route by T vs T_sat (LIQUID at T ≤ T_sat,L; VAPOR at T ≥ T_sat,V; SUPER if p ≥ p_c) |
  | `iphase_liquid`                 | Force LIQUID table at p < p_c; SUPER at p ≥ p_c              |
  | `iphase_gas`                    | Force VAPOR  table at p < p_c; SUPER at p ≥ p_c              |
  | `iphase_supercritical`          | Force SUPER  table regardless of p                            |
  | `iphase_supercritical_liquid`   | Force SUPER  table (caller knows the T < T_c / p ≥ p_c quadrant) |
  | `iphase_supercritical_gas`      | Force SUPER  table (T > T_c, p ≥ p_c)                         |
  | `iphase_critical_point`         | Force SUPER  table; the critbox fallback catches the cusp     |
  | `iphase_twophase`               | Throw — base class TabularBackend will reject PT in dome     |

Implementation: an internal `ForcedRegion` enum derived from
`imposed_phase_index` decides which table to use BEFORE consulting
`T_sat`.  The xnorm computation then uses the sat-T as the xnorm=1
(LIQUID) or xnorm=0 (VAPOR) boundary, so routing through LIQUID with
`T = T_sat,L` evaluates exactly at the saturated-liquid corner of the
table, returning rho_sat,L to machine precision.

The `_phase` member returned to the caller is set to the imposed phase
when one was specified (so `AS->phase()` after the update echoes the
user's intent); when no phase was imposed, it's derived from the active
region (`iphase_liquid` / `iphase_gas` / `iphase_supercritical_*`).

Example:

```cpp
auto AS = AbstractState::factory("SBTL&HEOS", "CO2");
const double P = 15.8e5;                      // CO2: T_sat ≈ 246.2 K, rho_L ≈ 1063, rho_V ≈ 41
// Default: T = T_sat routes to LIQUID by convention.
AS->update(PT_INPUTS, P, T_sat);              // rho() ≈ 1063 (sat-L)
AS->specify_phase(iphase_gas);
AS->update(PT_INPUTS, P, T_sat);              // rho() ≈ 41   (sat-V)
AS->specify_phase(iphase_supercritical_liquid);
AS->update(PT_INPUTS, 100e5, 290);            // SUPER table, T < T_c, p > p_c
AS->unspecify_phase();
```

### Why this matters

Without imposed-phase routing, a user computing properties along the
saturation curve via PT_INPUTS would silently get whichever side the
convention picks (we picked LIQUID).  For workflows that step around
the dome — e.g. integrating from the saturated-vapor curve along an
isobar into superheated vapor — the user must be able to anchor to the
correct side.  Honouring `specify_phase()` is the standard CoolProp
contract for this case, and the SBTL normalized tables make it cheap:
each region has its own table whose xnorm=0 or xnorm=1 boundary IS the
saturation curve, so routing through the right table delivers the
saturated-phase properties exactly.

## Critical-region HEOS-fallback box

Same pattern as the PH path's critical box.  Within a window around
`(T_c, p_c)` no polynomial backbone can faithfully reproduce the cusp
in ρ(T, p) at the critical point and the surrounding shoulder where
ρ drops steeply just above p_crit.  Box dimensions chosen empirically
against a 2 000-pt random CO₂ / Water sweep:

  - `p ∈ [0.75 · p_c, 1.75 · p_c]`
  - `T ∈ [0.90 · T_c, 1.30 · T_c]`

Inside the box: HEOS direct (~10 µs).  Outside: spline lookup (~1.3 µs).

## Out-of-range guards

  - `P < tbl.yvec.front()` or `P > tbl.yvec.back()` → throw clean.
    The base-class fall-through path silently extrapolates the legacy
    `logpT` cell polynomial when the query lies outside the table's
    P range, producing values up to ρ ~ 1e+5 kg/m³ (e.g. the
    notorious `T=223.9 K, P=3.3 bar` for CO₂, where `P < p_triple`
    but `PQ_INPUTS` extrapolates the sat curve along the metastable
    branch instead of throwing).
  - In-dome `PT_INPUTS` (T_sat,L < T < T_sat,V at subcritical P) →
    throw with a "use HmassP_INPUTS or PQ_INPUTS instead" message.

## Property accessors

When a PT query routes through normpt, `_normpt_active_table` and
`_normpt_active_coeffs` are set, and `_normpt_xnorm` caches the
xnorm coordinate.  The property accessors `calc_hmolar / smolar /
umolar / rhomolar` consult these and return values either cached
from the routing step (`_hmolar`, `_smolar`, etc.) or evaluated
from the active normpt cell directly (same `evaluate_single_phase`
the PH path uses).

## Validation

CO₂, 2 000 random `(T, P)` points across `[220, 1000] K × [1, 500] bar`,
single-phase only:

  | metric  | value      |
  |---------|------------|
  | median  | 1.0e-8     |
  | 99th    | 1.7e-4     |
  | max     | 1.7e-3     |
  | # err > 1e-3 | 5 / 1486 |
  | PT timing  | 1.27 µs/call |
  | HEOS PT    | 10.7 µs/call |
  | speedup    | 8.4× |

Water, 2 000 random `(T, P)` points across `[280, 1000] K × [0.01, 500] bar`:

  | metric  | value      |
  |---------|------------|
  | median  | 1.5e-8     |
  | 99th    | 2.0e-5     |
  | max     | 3.8e-4     |
  | # err > 1e-3 | 0 / 2000 |

The CO₂ tail (5 of 1486 points with err 1.1e-3 to 1.7e-3) sits at the
edge of the HEOS-fallback box in the supercritical-shoulder corner.
Tightening the box covers them; alternative is a Kunick-style L/G split
of SUPER along the `T_c` isotherm — see
`dev/sbtl_pt_outstanding_work.md` for the planned half-day follow-up.

## Why this is the right design

- **Mirrors the proven PH path.**  Same `xnorm` × `log p` rectangle,
  same B-spline backbone, same critical-region fallback box.  Lifts
  the existing `evaluate_single_phase` / `find_native_nearest_good_indices`
  machinery unchanged.
- **Saturation curve as a coordinate axis.**  Each cell sits on
  exactly one side of the dome.  No more `set_alternate*`-bumping of
  the LIQUID cell to the GAS side or vice versa, which was the root
  cause of the 5 % compressed-liquid residue.
- **Machine-precision sat boundaries.**  `T_sat,L/V(P)` comes from
  `PQ_INPUTS` (HEOS sat solver), not a polynomial fit to the sat curve.
  Avoids the "sat boundary computed two different ways gives two
  different xnorm at the same query" footgun.
- **Per-fluid correctness via the melting line.**  `T_lo(P) = T_melt(P)`
  preserves the EOS validity range for cryogens with deep melting
  curves and water with its negative-slope ice I band.
