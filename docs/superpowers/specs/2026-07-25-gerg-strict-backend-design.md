# GERG-2004 / GERG-2008 strict backend

**Date:** 2026-07-25
**Issue:** CoolProp-p8ub
**Status:** design approved, ready for implementation planning

## Goal

Implement the GERG-2004 and GERG-2008 wide-range equations of state for natural
gases as two new CoolProp backend families, derived from
`HelmholtzEOSMixtureBackend` so they inherit the flash algorithms, VLE solvers,
derivative machinery, and property interface.

The backends are *strict*: they admit only the components each model defines,
carry only parameters published with each model, and refuse to answer questions
the model does not cover.  A number obtained from the `GERG2008` backend is a
GERG-2008 number, with no admixture of CoolProp's own correlations.

The reference implementation is teqp's `include/teqp/models/GERG/GERG.hpp`, and
teqp's test values are the acceptance criterion.

## Why the existing CoolProp data cannot be reused

Two independent reasons, either sufficient on its own.

**Binary interaction parameters.**  *(Corrected 2026-07-29 — the original
version of this paragraph was wrong; see below.)*  All 210 GERG-2008 binary
pairs are present in `dev/mixtures/mixture_binary_pairs.json`.  Of those, 194
carry `Kunz-JCED-2012` provenance and are the GERG values themselves; 15 are
`Gernert-Thesis-2013` refits and 1 is `Tkaczuk-JPCRD-2020`.  So 16 of 210
pairs would silently differ from GERG if the default library were used.

The original text claimed 582 `Bell-JCED-2016` refits override GERG pairs.
That was a misreading: 582 is the Bell-2016 count across the *whole* 888-pair
library, and **none** of those rows is a GERG pair — they cover refrigerants
and other non-GERG systems.  The binary-pair argument for a separate backend
is therefore much weaker than first stated: it is 16 pairs, not 582.

**Pure-fluid equations of state.**  CoolProp ships the reference EOS for each
fluid (Setzmann-Wagner for methane, Span-Wagner for CO2, and so on).  GERG uses
its own shortened technical form — 12 to 24 polynomial and exponential terms
with a shared exponent set for most fluids.  These are different equations
producing different numbers.  This is the load-bearing argument and it is
unaffected by the binary-pair correction above.

**The gas constant.**  *(Added 2026-07-29, discovered during implementation.)*
GERG specifies `R = 8.314472 J/mol/K`.  CoolProp's `calc_gas_constant` returns
the CODATA value for any mixture under the default `NORMALIZE_GAS_CONSTANTS`
configuration, which rescales `p`, `alpha^ig`, `c_v` and `w` by about 1.1e-6
while leaving `alpha^r` exact.  A backend that did not override this would be
wrong at the sixth significant figure in a way no structural test would catch.

The backend therefore carries a self-contained parameter set: pure EOS, ideal-gas
coefficients, reducing parameters, departure functions, and ancillaries.

## What CoolProp already provides

The mixture machinery is already GERG-shaped, because CoolProp's multi-fluid
model descends from the same Kunz & Wagner formulation.  This is a
wiring-and-data exercise, not new thermodynamics.

| Need | Existing class | Location |
|---|---|---|
| Reducing function | `GERG2008ReducingFunction(pFluids, beta_v, gamma_v, beta_T, gamma_T)` | `src/Backends/Helmholtz/ReducingFunctions.h:144` |
| Departure function | `GERG2008DepartureFunction(n, d, t, eta, epsilon, beta, gamma)` | `src/Backends/Helmholtz/ExcessHEFunction.h:104` |
| Departure matrix + `F_ij` | `ExcessTerm` | `src/Backends/Helmholtz/ExcessHEFunction.h:200` |
| Ideal-gas sinh/cosh terms | `IdealHelmholtzGERG2004Sinh` / `Cosh` | `include/CoolProp/fluids/Helmholtz.h:1176,1212` |
| Residual polynomial/exponential terms | `ResidualHelmholtzGeneralizedExponential::add_Power` | `include/CoolProp/fluids/Helmholtz.h` |
| Mixture ideal-gas summation | `calc_alpha0_deriv_nocache` mixture branch | `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp:3629` |

The reducing and departure constructors take coefficient matrices directly, so
the GERG backend populates them from its own tables rather than going through
`MixtureParameters::set_mixture_parameters`, which reads the global JSON
library.  The mixture ideal-gas branch already cites *Table B5, Kunz & Wagner,
JCED, 2012* and evaluates each component's `alpha^0` at that component's own
`T_ci` and `rho_ci`, which is exactly GERG's definition.

## Architecture

### Files

```
src/Backends/GERG/GERGData.h          coefficient tables (private to the backend)
src/Backends/GERG/GERGBackend.h       GERGMixtureBackend declaration
src/Backends/GERG/GERGBackend.cpp     construction, registration, strict checks
```

Nothing is installed under `include/CoolProp/`.  The coefficient tables are an
implementation detail of the backend and are not part of the public API.

### Upstream change

`HelmholtzEOSMixtureBackend::set_mixture_parameters()` becomes `virtual`, so
`GERGMixtureBackend` can override it and populate the reducing function and
excess term from GERG tables instead of the global JSON library.
`set_components` is already virtual and needs no change.

*(Corrected 2026-08-07.)*  This originally read "This is the only modification
to existing CoolProp code beyond backend registration."  That is no longer
true.  As shipped, the branch also modifies:

- `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.h` — `set_mixture_parameters()`
  made `virtual` (the change described immediately above; listed here so this
  enumeration is genuinely exhaustive),
- `src/Backends/Helmholtz/ReducingFunctions.cpp` — the removable-singularity
  guard in `GERG2008ReducingFunction::f_Y_ij` and its two first-derivative
  helpers (a shared fix affecting default `HEOS` too; GitHub #1677),
- `src/CoolProp.cpp` — `is_gerg_backend_string()` and the
  `set_reference_stateS` refusal,
- `src/DataStructures.cpp` / `include/CoolProp/DataStructures.h` — the new
  `backend_families` and `backends` enumerators.

Outside `src/`, the branch also touches `CMakeLists.txt`, `dev/ci/preflight.sh`,
`CoolPropBibTeXLibrary.bib` and the `Web/` documentation.

### Data header layout

`GERGData.h` mirrors teqp's namespacing so the two can be diffed against each
other:

```
namespace CoolProp { namespace GERG {
  namespace GERG2004 {
    component_names          // 18 entries
    get_pure_info(name)      // rhoc [mol/m^3], Tc [K], M [kg/mol]  (Table A3.5)
    get_pure_coeffs(name)    // n, t, d, c, l                        (Table A3.2)
    get_alphaig_coeffs(name) // n0, theta0                           (Table A3.1)
    get_betasgammas(f1, f2)  // beta_T, gamma_T, beta_v, gamma_v     (Table A3.8)
    get_Fij(f1, f2)          // optional<double>                     (Table A3.6)
    get_departurecoeffs(f1, f2)
  }
  namespace GERG2008 {
    component_names          // 21 entries
    // overrides only what changed; falls through to GERG2004 otherwise
  }
}}
```

There are **23 distinct pure EOS** across both models.  GERG-2008 changes the
pure EOS for carbon monoxide and isopentane relative to GERG-2004, and adds
hydrogen sulfide, n-nonane, and n-decane.  Everything else falls through.

Reducing parameters (`beta`, `gamma`) exist for **every** binary pair: 153 for
GERG-2004, 210 for GERG-2008.  No pair is unsupported.  **15 pairs** additionally
carry a departure function via a non-null `F_ij`; the remaining pairs take
`F_ij = 0` and contribute no departure term.  Seven of the 15 use the
generalized departure function scaled by `F_ij != 1`.

### Backend registration

Two families in `src/DataStructures.cpp`:

```cpp
{GERG2004_BACKEND_FAMILY, "GERG2004"},
{GERG2008_BACKEND_FAMILY, "GERG2008"},
{GERG2004_BACKEND, "GERG2004Backend", GERG2004_BACKEND_FAMILY},
{GERG2008_BACKEND, "GERG2008Backend", GERG2008_BACKEND_FAMILY},
```

registered through the standard `GeneratorInitializer<...>` path, so
`AbstractState::factory("GERG2008", {...})` and
`PropsSI(..., "GERG2008::METHANE[0.9]&NITROGEN[0.1]")` work without special
cases in `factory()`.

A single `GERGMixtureBackend` class serves both families and both the pure
(N = 1) and mixture cases; the model year is a construction parameter.  GERG is
inherently a mixture model, so there is no pure/mixture class split of the
`HEOS_BACKEND_PURE` / `HEOS_BACKEND_MIX` kind.

### Component naming

Names resolve through CoolProp's normal alias and CAS machinery first, then map
CAS to GERG component via a table in `GERGData.h`.  So `CO2`, `R744`,
`CarbonDioxide`, and `124-38-9` all reach GERG's `carbondioxide`.  A name
CoolProp cannot resolve, and a name CoolProp resolves to a fluid outside the
model's component set, both throw.

## Two quirks that must be replicated exactly

Both are required for bit-level agreement with teqp.

### The `R*/R` ratio

GERG's ideal-gas Helmholtz energy carries the ratio
`R* / R = 8.314510 / 8.314472` multiplying the bracketed sum:

```
alpha^0_oi(rho, T) = ln(rho / rho_ci)
                   + (R*/R) [ n1 + n2 (Tci/T) + n3 ln(Tci/T)
                              + sum_k n_k ln|sinh(theta_k Tci/T)|
                              - sum_k n_k ln|cosh(theta_k Tci/T)| ]
```

teqp's implementation notes that the GERG-2004 monograph places this ratio
incorrectly and the GERG-2008 manuscript corrects it.  We follow the GERG-2008
placement, folding the ratio into the coefficients at table-build time so
CoolProp's existing ideal-gas term classes can be used unmodified.

### Reference state

teqp **discards** the published integration constants `n0[1]` and `n0[2]` and
recomputes them so that `h = 0` and `s = 0` for the ideal gas at
`T0 = 298.15 K`, `p0 = 101325 Pa`, `rho0 = p0 / (R T0)`.  We replicate this
recomputation.  Without it, enthalpy and entropy will not match teqp even
though pressure, `c_v`, and speed of sound will.

*(Corrected 2026-08-07.)*  This originally read "`set_reference_stateS` remains
available to users.  A reference-state change is a pure offset in `alpha^0` and
does not alter the model."  As shipped it **throws** `NotImplementedError` on
the GERG backends: CoolProp applies a reference-state change by writing an
offset into the global fluid-library entry, and the GERG backends do not read
that library, so the call was a silent no-op.  Refusing is better than
pretending.  See `Web/coolprop/GERG.rst`.

## Strictness

| Behaviour | Rule |
|---|---|
| Component outside the model's 18 / 21 | throw |
| Name CoolProp cannot resolve | throw |
| `set_binary_interaction_double` and friends | throw |
| `viscosity()`, `conductivity()`, `surface_tension()` | throw `NotImplementedError` |
| Superancillary | none attached |
| Range | EOS limits set to GERG's extended range (60-700 K, p <= 70 MPa).  A `GERGMixtureBackend::update()` override enforces the **temperature** bound only; `pmax` is deliberately NOT enforced (a valid (T, rho) point inside the two-phase dome legitimately has a pressure far outside the envelope).  *(Clarified 2026-08-07.)* |

*(Corrected 2026-07-29.)*  This row originally said "enforced by CoolProp's
existing limits machinery".  That was false and would have shipped a
fail-open guard: `PT_flash` never compares `T` against `Tmax` at all, so the
upper bound was unenforced, and the lower bound only appeared to work through
an unrelated missing-ancillary accident.  The backend therefore carries its
own `update()` override with a real range check.  *(Updated 2026-08-07:
`update_with_guesses` was originally recorded here as "a known gap"; it is now
overridden too and enforces the same check.)*

Mutating beta, gamma, or `F_ij` and still calling the result GERG is a category
error, so the setters throw rather than silently producing a mutant model.

Transport properties are not part of either GERG model.  CoolProp's transport
correlations would be correct numbers from a different model, and exposing them
through a backend named `GERG2008` invites misattribution.

### The superancillary hazard

This one is load-bearing and worth stating explicitly, because the obvious
shortcut is a silent correctness bug.

`FlashRoutines::sat_superanc_path_applies` (`src/Backends/Helmholtz/FlashRoutines.cpp:558`)
routes saturation calls for pure fluids that own a superancillary straight to
the Chebyshev expansion, and **returns that value as the answer** rather than as
an iteration guess.  Attaching CoolProp's existing superancillary to a GERG
fluid would therefore return Setzmann-Wagner saturation densities labelled
GERG-2008, with no warning and no iteration to correct them.

GERG fluids consequently attach no `SUPERANCILLARY` blob.  `get_superanc()`
returns `nullptr`, `sat_superanc_path_applies` is false, and saturation goes
through the classical ancillary-seeded VLE solver, converging to GERG-consistent
values.

Fitting genuine superancillaries against the 23 GERG pure EOS is a reasonable
follow-on.  It is purely additive: the blob hangs off `EquationOfState`, so
dropping in GERG-fit expansions later lights up the fast path with no backend
API change.  Estimated cost, measured from `dev/fluids/Methane.json`: 65 kB raw
and 25 kB gzip per fluid, so roughly 1.5 MB raw and 580 kB gzip for 23 fluids,
plus a `fastchebpure` run per fluid, `source_eos_hash` stamping, and extension
of the two freshness gates (`dev/scripts/check_superanc_freshness.py`,
`dev/scripts/check_superanc_release_pin.py`).  Filed separately, not part of
this work.

## Ancillaries

Saturation ancillaries are **fitted against the GERG pure EOS**, not borrowed
from CoolProp's fluid library.  Borrowing would work numerically — ancillaries
only seed VLE iteration and the converged answer would still be GERG — but it
would put a non-GERG-derived correlation inside a strict backend, and the cost
of doing it properly is small.

Generation is offline, in `dev/gerg/`:

1. Drive teqp's Python bindings, which already expose both GERG models, to
   trace the saturation curve for each of the 23 distinct pures from near the
   critical point down to `T_min`.  `fastchebpure`'s VLE-tracing helpers
   (`fastcheb.cpp`) provide the stepping logic.
2. Regress CoolProp's standard ancillary forms (`pV`, `rhoL`, `rhoV`) against
   the traced data.
3. Emit a C++ coefficient table into the private header.

At roughly 3.6 kB per fluid in the existing ANCILLARIES form, this is about
80 kB for 23 fluids — negligible against the 23 MB shipped fluid blob.

The generator is a development-time tool.  It is not built as part of CoolProp
and adds no runtime dependency on teqp or fastchebpure.

## Testing

Catch2, tagged `[GERG]`, in `src/Tests/`.

teqp's own tests compare against the AGA8 reference implementation
(`src/tests/GERG2008.cpp`, a translation of the published Fortran) at relative
tolerances of 1e-9 to 1e-10, and only loosely (1e-5, currently commented out)
against the monograph's printed tables.  We therefore port the
reference-implementation-derived values, which are the tighter and more
meaningful check, rather than the printed tables.

Coverage:

- **Construction.** All 18 GERG-2004 and all 21 GERG-2008 components load; all
  153 / 210 binary pairs construct; the 15 departure-carrying pairs return
  coefficients and the rest do not; unknown fluids throw.
- **Ideal gas.** Recomputed integration constants give `h = s = 0` at
  298.15 K / 101325 Pa for every pure, matching teqp's
  `Validate all GERG2008 pures reference states`.
- **Pure thermodynamics.** `p`, `c_v`, and speed of sound for all 23 distinct
  pures against teqp.
- **Binary mixtures.** teqp's `Validate all GERG2008 binaries` set.
- **Full mixtures.** The AGA8 21-component gas compositions from teqp's
  `validation_data` table.
- **Strictness.** Every rule in the strictness table has a test asserting the
  throw. A strictness guard that can silently no-op is a defect, not a nit.
- **Cross-model.** GERG-2004 and GERG-2008 disagree where they should (carbon
  monoxide, isopentane, changed binary pairs) and agree where they should.

## Implementation phasing

Phases 1 and 2 require no ancillaries at all: every teqp validation value is
reachable from `(T, rho)` and `(p, T)` inputs.  This keeps the numerically
delicate work in front and the infrastructure work behind it.

1. **Data tables and pure-fluid thermodynamics.**  `GERGData.h`, pure
   `CoolPropFluid` construction, ideal-gas integration-constant recomputation.
   Validated against teqp for all 23 pures.  No mixtures, no saturation.
2. **Mixture wiring.**  Reducing function and excess term populated from GERG
   tables; `set_mixture_parameters` override.  Validated on binaries and the
   AGA8 gas set.
3. **Ancillaries and VLE.**  Offline generator, fitted tables, saturation tests.
4. **Strictness, registration, documentation.**  Backend families in
   `DataStructures.cpp`, throw paths and their tests, `Web/` documentation
   describing what the backend deliberately does not provide.

## Out of scope

- Superancillaries for GERG pures (separate follow-on; see above).
- Transport properties.
- Mixing GERG components with fluids outside the published sets, or with other
  backends' fluids.
- Any mechanism for user-supplied or modified GERG parameters.
- Extending CoolProp's default HEOS backend to optionally use GERG pure-fluid
  EOS, the way REFPROP's `REFPROP_USE_GERG` flag does.  That is a different
  feature with different semantics.
