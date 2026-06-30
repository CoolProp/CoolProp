# Research notes: proper mixing models for absorption-cycle fluids

**Status: deferred.** Not part of the current incompressible-backend
zero-division/maintainability cleanup. Captured here so the research doesn't
get lost; revisit as a separate, deliberate project.

## The problem

CoolProp's incompressible solutions (LiBr/H2O, NH3/H2O via the Melinder MAM/
MAM2 fluids, etc.) compute enthalpy and entropy purely by integrating a fitted
`cp(T,x)` surface in temperature. There is no composition-dependent excess
term anywhere in that calculation, so there is no enthalpy of mixing. This is
a deliberate, documented design choice (see the maintainer's own comments on
issue #533 below), not an oversight — but it means absorption-cycle energy
balances across streams of different concentration are wrong by the heat of
mixing.

Related CoolProp issues:
- [#533](https://github.com/CoolProp/CoolProp/issues/533) — MAM/MAM2 enthalpy
  of mixing ~100x too low vs. Melinder/EES reference data; maintainer (jowr)
  confirms the omission is intentional, and later says he regrets the
  decision.
- [#781](https://github.com/CoolProp/CoolProp/issues/781) — jowr's own
  proposal to add a Gibbs-energy-based incompressible class, which would
  naturally produce the missing excess/mixing terms. Never implemented.
- [#1690](https://github.com/CoolProp/CoolProp/issues/1690) — concrete user
  pain: can't write a consistent enthalpy balance across LiBr/H2O streams of
  different concentration in an absorption-chiller generator.
- [#341](https://github.com/CoolProp/CoolProp/issues/341) — long-running
  wishlist thread for a real NH3/H2O *mixture* (not solution) model in
  CoolProp; contains the most useful pointers (below).

## NH3/H2O

Reference formulation: **Tillner-Roth & Friend (1998)**, *A Helmholtz free
energy formulation of the thermodynamic properties of the mixture
{water + ammonia}*, J. Phys. Chem. Ref. Data 27(1):63-96 — adopted as the
IAPWS 2001 guideline. Structure: `a_mix(tau,delta,x) = sum_i x_i * a_i(tau,delta)
+ a^E(tau,delta,x)`, i.e. pure-component Helmholtz energies (water via
IAPWS-95, ammonia via Baehr & Tillner-Roth) plus a composition-dependent
excess/departure term (polynomial + exponential terms in delta, tau, x, x^2).
This is the *same shape* as CoolProp's own HEOS multi-fluid mixture model
(reference + departure function per binary pair) — confirmed CoolProp has
departure-function data for AMMONIA paired with ARGON, BENZENE, BUTANE,
HEXANE, KRYPTON, etc. in `dev/mixtures/mixture_binary_pairs.json`, but **no
AMMONIA/WATER pair**. The architecture to hold this already exists; nobody
has supplied the fit.

Open implementations:
- **teqp** ([usnistgov/teqp](https://github.com/usnistgov/teqp)) — same lead
  author as CoolProp (Ian Bell), NIST-maintained. Has a hardcoded
  `AmmoniaWaterTillnerRoth` class implementing this exact model on a modern
  automatic-differentiation core (C++ with Python bindings). When a
  contributor asked (issue #341, 2023) about implementing Tillner-Roth &
  Friend directly in CoolProp, Bell's own advice was to use teqp instead,
  noting CoolProp's mixture algorithms are "not very refined." **This is the
  strongest single pointer** — both for the NH3/H2O physics and as a worked
  example of how the same author now structures EOS-with-departure-function
  code (AD instead of hand-coded analytic derivatives, which sidesteps the
  whole class of derivative/integral-correctness bugs we hit in
  `Polynomial2DFrac`).
- **iapws** ([jjgomera/iapws](https://github.com/jjgomera/iapws), pure
  Python, **GPL-3.0**) — `iapws/ammonia.py` implements the same
  Tillner-Roth & Friend departure function in full (h, s, u, g, a, cp, cv,
  fugacity all exposed). Readable, working reference math — but GPL-3.0
  means it can be studied, not vendored into MIT-licensed CoolProp.
- A NIST-internal NH3+H2O mixture model (Ian Bell, mentioned in #341,
  ~2017-2023) was put on hold and never published; REFPROP has its own
  (proprietary) implementation.

## LiBr/H2O

Two competing correlations show up in the literature:
- **Pátek & Klomfar (2006)**, *Int. J. Refrigeration* 29:566-578 — simpler,
  polynomial-type fit; this is what EES's `LiBrH2O` library uses, and is
  CoolProp's own current reference (`Patek2006` citation already in
  `Web/fluid_properties/Incompressibles.rst`).
- **Yuan & Herold (2005)** — University of Maryland Sorption Systems
  Consortium — a proper multiproperty **Gibbs-energy** correlation.

Open implementation:
- **LiBr.jl** ([hzgzh/LiBr.jl](https://github.com/hzgzh/LiBr.jl), Julia) —
  implements the Yuan & Herold approach: one Gibbs function `g(T, p, x)`
  fitted to data, with density, enthalpy, entropy, cp, and chemical potential
  (partial molar Gibbs energy) *all* derived by differentiating that single
  function. This is structurally exactly what jowr proposed in #781: the
  enthalpy of mixing isn't a separate fitted term bolted onto a `cp(T)`
  integral, it's `g_mixture(T,x) - sum(x_i * g_i_pure(T))`, automatically
  consistent with every other derived property by construction. License of
  LiBr.jl itself not yet checked — verify before reusing any code; the
  underlying equations are published literature regardless.

## Applied reference

- **openACHP** ([nfette/openACHP](https://github.com/nfette/openACHP),
  Python/Jupyter) — implements and cross-validates *both* NH3/H2O
  (Ibrahim-Klein, plus comparisons against REFPROP/CoolProp) and LiBr/H2O
  (Pátek-Klomfar, EES) specifically for absorption-cycle cycle calculations.
  Useful for seeing how someone actually wires these correlations into
  stream-enthalpy balances once the property layer is right — directly
  addresses the scenario in #1690.

## The common architectural lesson

Every credible implementation here (CoolProp's own HEOS departure functions,
Tillner-Roth & Friend, Yuan & Herold) is the same shape: **one fundamental
thermodynamic potential (Helmholtz or Gibbs energy) as a function of
`(T, p/rho, x)`, with every other property obtained by differentiating it** —
density, cp, enthalpy, entropy, chemical potential all consistent by
construction, and mixing behavior is whatever the potential's `x`-dependence
says it is. CoolProp's current `IncompressibleFluid`/`IncompressibleBackend`
is the opposite shape: five independently-fitted surfaces (`rho(T,x)`,
`cp(T,x)`, `visc(T,x)`, `cond(T,x)`, `psat(T,x)`) with enthalpy/entropy
*reconstructed after the fact* by integrating the cp surface. That's why
there's no mixing term (nothing in that pipeline ever had an `x`-dependent
potential to differentiate).

A future fix likely means giving incompressible *solutions* (not pure
fluids — those have no composition, so no mixing term to begin with) their
own Gibbs-energy-based property class, the way LiBr.jl does, while pure
fluids and simple brines without published excess-Gibbs data could keep the
existing polynomial approach.

## Suggested next steps, when this is picked back up

1. Pick one binary (LiBr/H2O is the best-documented, most-requested
   candidate per #1690/#1331/#2567) and prototype a Gibbs-energy class
   against Yuan & Herold or Pátek & Klomfar, validated against LiBr.jl's
   published values.
2. Decide the JSON schema extension needed to carry a Gibbs-energy fit
   alongside (not replacing) the existing polynomial schema, so pure fluids
   are unaffected.
3. Revisit NH3/H2O only after LiBr/H2O proves the pattern out — it's a
   bigger lift (needs a real departure-function fit, not just a literature
   Gibbs correlation) and teqp's `AmmoniaWaterTillnerRoth` is the right
   reference implementation to validate against.
