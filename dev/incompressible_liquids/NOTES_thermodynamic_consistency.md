# Notes: thermodynamic consistency of the INCOMP backend

Assessment of how consistent the incompressible property model is, what is
already exact, and which options exist for improving it with the math that
is already in the codebase. Companion to `NOTES_mixing_models.md`, which
covers the (bigger) composition/mixing question.

## What the model actually is

Only two fits carry all caloric information:

- `rho(T,x)` -- polynomial, pressure-independent (truly incompressible),
- `c(T,x)` -- polynomial; reported as both cp and cv.

Everything else caloric is *derived* with closed-form identities
(`IncompressibleBackend.cpp`, `raw_calc_hmass`/`raw_calc_smass`):

    h(T,p,x) = Int c dT |_T  +  p * [ (1/rho) * (1 + (T/rho) * drho/dT) ]
    s(T,p,x) = Int (c/T) dT |_T  +  p * [ (1/rho^2) * drho/dT ]
    u = h - p/rho

The `Int ... dT |_T` terms are indefinite integrals evaluated at a point;
reference-state subtraction (`h_ref`, `s_ref` caching) turns them into
definite integrals. The pressure terms are exact because
`(dh/dp)_T = v - T (dv/dT)_p` and `(ds/dp)_T = -(dv/dT)_p` are themselves
pressure-independent when rho does not depend on p.

## What is exact (verified on this branch)

- **T-integrals of the c fit**: `Polynomial2DFrac::integrateCoeffs` (for h)
  and `fracIntCentral`/`fracIntCentralDvector` (for s, the D(j,T) binomial
  expansion documented in `Incompressibles.rst`) are exact analytic
  integrals of the fitted polynomial. No approximation.
- **drho/dT**: exact analytic differentiation via
  `Polynomial2DFrac::deriveCoeffs`. Since the 0/0-at-Tbase fix on this
  branch, the value at `T == Tbase` is the exact limit, not an epsilon-band
  interpolation.
- **Maxwell relation** `(ds/dp)_T = -(dv/dT)_p`: holds identically by
  construction.
- **u = h - p*v**: holds by definition (`calc_umass`).
- **Reference state**: `h`/`s` are exactly `h_ref`/`s_ref` at
  `(T_ref, p_ref, x_ref)`; `set_fractions` re-pins the reference to the
  active composition, so the cancellation is bit-consistent.

There is no remaining numerical approximation in the *forward* h/s path for
the shipped polynomial fits. The only epsilon-band hack left is
`evaluateAwayFromPole` in `IncompressibleFluid.cpp`, guarding the *genuine*
pole of the exponential-type fits -- which affects viscosity/psat/Tfreeze
only and never enters h/s.

## Remaining inconsistencies (in decreasing size)

1. **No enthalpy/entropy of mixing across compositions.** The big one; see
   `NOTES_mixing_models.md`. Not fixable inside the current
   five-independent-fits architecture -- needs a Gibbs-energy-based class.

2. **Reported cp is not (dh/dT)_p.** `calc_dhdTatPx` returns the fitted
   `c(T,x)`, but the model's own enthalpy is `h = H(T,x) + p*g(T,x)` with
   `g = (dh/dp)_T`, so the true derivative is `c + p * dg/dT`. A finite
   difference of `hmass()` over T therefore disagrees with `cpmass()` by
   `p * dg/dT` (small at ambient pressure: g depends on rho and drho/dT
   only, so the term is of order 0.1-1 J/kg/K against c of 2000-4000, but it
   grows linearly with p). Same for `calc_dsdTatPx = c/T`.

3. **(dcp/dp)_T = 0 violates the coupled Maxwell relation.** Consistency
   requires `(dcp/dp)_T = -T (d2v/dT2)_p`, which is nonzero for the fitted
   rho(T). This is the intrinsic price of fitting rho and c independently
   instead of deriving both from one potential.

4. **Non-polynomial density fits cannot produce h/s at all**:
   `drhodTatPx` only implements the polynomial case and throws otherwise.
   No shipped fluid is affected (all density fits are polynomials), but the
   generator does not forbid it either.

## Options with the existing math functions

- **Fix (2) cheaply.** Everything needed for `p * dg/dT` already exists:
  `g' (T) = d/dT [ 1/rho + (T/rho^2) drho/dT ]` requires rho, drho/dT and
  d2rho/dT2 -- the second derivative is one more `deriveCoeffs` pass (or
  `Polynomial2DFrac::derivative` with a second-order call), which the class
  already supports. This would make cpmass consistent with hmass at the cost
  of cpmass no longer exactly reproducing the fitted data at p != p_data.
  Decide deliberately: data-faithful cp (today) vs. model-consistent cp.
  If changed, `calc_dsdTatPx` should get the matching `p * d/dT[(1/rho^2)
  drho/dT]` term so ds/dT stays the exact derivative of smass.

- **(3) is not worth patching.** Any bolt-on correction that forces
  `(dcp/dp)_T = -T v''` re-derives cp from rho, i.e. stops honoring the cp
  data. The clean resolution is the same as for mixing: one fitted
  potential. For a pressure-independent liquid a Gibbs form
  `g(T,p,x) = g0(T,x) + v(T,x) * (p - p0)` fitted to the *same* rho/cp data
  would reproduce today's model almost exactly while making every
  derivative consistent by construction -- and `Polynomial2DFrac` already
  provides all evaluate/derive/integrate machinery a `g0` polynomial needs.
  That is the incremental path to jowr's #781 proposal, and it can be
  prototyped in the Python pipeline before touching C++.

- **Close (4) in the generator, not the backend**: reject non-polynomial
  density/cp fits at JSON-writing time so the constraint is visible where
  fluids are added.

## Suggested order, when picked up

1. Decide the cp question (option above) -- small, self-contained, testable
   against finite differences of hmass.
2. Prototype the `g(T,p,x)` single-potential class for one solution
   (LiBr/H2O per `NOTES_mixing_models.md`), reusing `Polynomial2DFrac`.
3. Only then consider schema/backend changes.
