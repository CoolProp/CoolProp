# REFPROP-vs-CoolProp functionality gaps

**Date:** 2026-09-06
**Issue:** CoolProp-it6e (epic), CoolProp-it6e.1 .. .8 (individual gaps)
**Status:** assessment complete; each gap needs its own design + plan before implementation

## Goal

Establish, from the source rather than from folklore, which pieces of *core
REFPROP functionality* CoolProp's native `HEOS` backend does not provide, so
that the roadmap argues from evidence and so that non-gaps stop being
re-proposed.

Every claim below is anchored to a file and line in this repository at
`6de0bf8c5`.  Where a claim was checked and turned out to be **wrong**, it is
recorded in "Rejected claims" rather than deleted — the point of the document
is to stop the same wrong ideas recurring.

## Scope

"Core functionality" means the property and equilibrium surface a library user
calls.  Explicitly out of scope: the REFPROP GUI, the Excel add-in, table
generation, and REFPROP's fitting tools.

This document deliberately does **not** become an implementation plan.  It
spans eight independent subsystems whose only common thread is "REFPROP has it";
each needs its own design doc and its own plan, because the work in each is
unrelated to the work in the others.  Bundling them into one plan would produce
a plan nobody can execute.

## Rejected claims

These were considered and are **not** gaps.  Recorded so they are not raised again.

**Multiphase / VLLE / LLE / three-phase flash is NOT a REFPROP gap.**  REFPROP
is two-phase VLE-only, exactly as CoolProp is.  CoolProp's VLE machinery
(`src/Backends/Helmholtz/VLERoutines.cpp`) is structurally two-phase — bubble/dew
plus a single Rachford-Rice split, with no phase-stability-driven phase-count
determination — and REFPROP's flash routines make the same restriction.  Wanting
VLLE is a legitimate goal, but it is a gap against `teqp` or a process
simulator, not against REFPROP, and must not be justified by pointing at
REFPROP.

**Solid-fluid equilibrium is NOT a REFPROP gap.**  REFPROP's `MELTT`/`MELTP`
and `SUBLT`/`SUBLP` return correlated *boundary lines* for pure fluids.  They
are not an SLE flash.  Neither library computes solid-fluid equilibrium, so
"CO2 freeze-out in cryogenic mixtures" is not available from REFPROP either.
What survives from this area is narrow and is recorded as gap 5 below.

**Fugacity, chemical potential, virials, mixture critical points, phase
envelopes, reference-state setting, and BIP get/set are all present natively.**
See `include/CoolProp/AbstractState.h` (`calc_fugacity`,
`calc_fugacity_coefficient`, `calc_chemical_potential`, `calc_Bvirial`,
`calc_Cvirial`), `HelmholtzEOSMixtureBackend.h`
(`_calc_all_critical_points`, `calc_criticality_contour_values`), and
`src/Backends/Helmholtz/PhaseEnvelopeRoutines.cpp`.  Cubics (SRK/PR) and PC-SAFT
already cover the ground REFPROP's `PREOS` toggle addresses.

## The gaps

Ordered by impact, which is not the same as ordered by effort.

### 1. Mixture transport properties are a mole-fraction average, not a model — `CoolProp-it6e.1`, P1

The largest predictive gap.  In
`src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp`:

- `calc_viscosity` (line 832) returns `exp(sum_i x_i * ln eta_i)` for mixtures,
  after `set_warning_string("Mixture model for viscosity is highly approximate")`.
- `calc_conductivity` (line 1086) returns `sum_i x_i * lambda_i`, with the
  matching warning.

Each per-component value is evaluated by constructing a fresh pure-fluid
`HelmholtzEOSBackend` at the *mixture* `(rho, T)`, which is not a state the pure
fluid is necessarily in.  There is no mixture critical enhancement at all.

The underlying routines cannot help: `TransportRoutines.cpp` throws
`"is only for pure and pseudo-pure"` in **15** distinct places, covering every
dilute, initial-density, higher-order, friction-theory, Chung and critical
contribution.

REFPROP's `TRNPRPdll` runs a genuine mixture model — extended corresponding
states with a mixture conformal-state solve, plus friction theory.  CoolProp has
`viscosity_ECS`/`conductivity_ECS` declared in `TransportRoutines.h` (lines
274, 278) with a conformal-state solver, but they are wired for pure fluids
against a reference fluid, not for mixtures.

*Consequence:* for any real mixture, CoolProp's eta and lambda are not
predictive.  The warning string is set, but a caller who does not read
`get_warning_string()` gets a plausible-looking number with no accuracy claim
behind it.

*Effort:* large.  This is the multi-month item.  The existing ECS scaffolding is
the natural starting point, so the first design question is whether the
conformal-state solver generalises to mixtures or needs replacing.

### 2. Ammonia-water is absent entirely — `CoolProp-it6e.2`, P1

`dev/mixtures/mixture_binary_pairs.json` contains 888 binary pairs and
**no Ammonia/Water row**.  Verified programmatically, not by eye.  The
consequence is not "inaccurate" — it is that the mixture cannot be constructed
at all, so every absorption-cycle use case is a hard failure at `factory()`.

REFPROP carries the dedicated Tillner-Roth & Friend Helmholtz model for this
system, which is not expressible in CoolProp's three available departure-function
forms (`Exponential`, `GERG-2008`, `Gaussian+Exponential`; 28 functions total in
`dev/mixtures/mixture_departure_functions.json`).

*Design question to settle first:* fit a BIP row into the existing multi-fluid
framework (cheap, approximate) versus implement TRF as a dedicated model
(faithful, more work).  These have very different costs and the choice is not
obvious.

*Effort:* small if (a), medium if (b).  Highest impact-per-unit-effort on this
list, and the recommended first target.

### 3. Mixture surface tension throws — `CoolProp-it6e.3`, P2

`HelmholtzEOSMixtureBackend.cpp:712` throws
`NotImplementedError("surface tension not implemented for mixtures")`.
REFPROP's `SURFTdll` handles mixtures.  A mixing rule must be chosen before any
code is written; this is a modelling decision, not a plumbing one.

*Effort:* medium, dominated by the model choice and its validation data.

### 4. Four mixture flash input pairs are unimplemented — `CoolProp-it6e.4`, P2

In `src/Backends/Helmholtz/FlashRoutines.cpp`, these throw
`"not ready for mixtures"`:

| Line | Routine |
| ---- | ------- |
| 510  | `DP_flash` |
| 578  | `DQ_flash` |
| 616  | `HQ_flash` |
| 649  | `QS_flash` |

REFPROP ships the complete matrix (`DEFLSH`, `DHFLSH`, `DSFLSH`, `PDFLSH`,
`PEFLSH`, `TDFLSH`, `TEFLSH`, `DQFL2`, ...).

**Sequencing constraint:** `CoolProp-ft05` (P0 — mixture `HSU_P` flash returns
wrong `T`, silently, on master) and `CoolProp-1gth` (P1 — `DmassT`/`DmolarT`
can return a state at a different density than requested) are open defects in
the mixture flash paths that *already exist*.  Adding four more input pairs on
top of a base with known silent-wrong-answer bugs would multiply the surface
area of those bugs.  Fix those first.

*Effort:* medium per pair, and they share machinery.

### 5. No sublimation line, natively or by passthrough — `CoolProp-it6e.5`, P2

Two separate absences:

- **Native.** Zero occurrences of `sublim` in the Helmholtz backend.  The only
  repo-wide hits are the humid-air ice line in `src/HumidAirProp.cpp` and
  `SATT` `kph`-code comments in the REFPROP backend.  Melting *is* covered
  (`src/Backends/Helmholtz/MeltingCaloric.cpp`).
- **Passthrough.** `src/Backends/REFPROP/REFPROPMixtureBackend.cpp` wraps
  `MELTPdll` and `MELTTdll` but **not** `SUBLPdll`/`SUBLTdll`.  So even a user
  who owns REFPROP cannot reach its sublimation curve through CoolProp.

Scope note, to prevent the rejected claim from creeping back: this is boundary
*lines* only.  It is not SLE, and REFPROP does not offer SLE either.

*Effort:* the passthrough is small and mirrors the existing melting-line
wrapping exactly — it is the natural first commit.  The native pure-fluid
ancillary is the larger piece and needs source correlations per fluid.

### 6. No dielectric constant — `CoolProp-it6e.6`, P3

REFPROP exposes `DIELEC`.  CoolProp has no such property; the only `dielectric`
matches under `src/` belong to PC-SAFT's internal electrolyte term
(`src/Backends/PCSAFT/PCSAFTBackend.cpp`), which is not a general-purpose
dielectric constant and is not exposed as a parameter.

Requires a new entry in `DataStructures`, the parameter plumbing, and per-fluid
correlations — the last of which is the real cost.

*Effort:* medium, dominated by data collection.

### 7. No runtime model switching — `CoolProp-it6e.7`, P3

REFPROP's `SETMOD` selects an alternate published EOS or transport model per
fluid at runtime, and can disable the critical enhancement.  CoolProp binds
exactly one EOS per fluid JSON.

This is more relevant than its P3 suggests, because it interacts with work
already open: `CoolProp-9s9u` (REFPROP 10.1 viscosity correlations not
implemented), `CoolProp-3bpg` (how to adopt the 2022 R-134a viscosity given its
ECS reference role) and `CoolProp-f1ez` (nitrogen/argon viscosity supersession
blocked by the conductivity coupling) are all cases where being able to select
between an old and a new correlation would let both ship instead of forcing a
choice.  Worth revisiting the priority if the supersession work stalls.

*Effort:* medium; mostly an architecture change to the fluid-loading path.

### 8. No BIP estimation fallback — `CoolProp-it6e.8`, P3

CoolProp throws for any binary pair outside the 888-pair library unless the
caller explicitly invokes
`apply_simple_mixing_rule(id1, id2, "linear" | "Lorentz-Berthelot")`
(`src/Backends/Helmholtz/MixtureParameters.cpp:286`, rules at :241 and :251).
REFPROP falls back to its own estimation schemes when a pair is missing from
`HMX.BNC`.

**This is a difference in philosophy, not straightforwardly a defect.**
CoolProp's fail-loud default means a user never silently receives an estimated
number believing it to be a fitted one.  That is the better default and should
be preserved.  The design question is therefore whether to offer an *opt-in*
estimation mode, not whether to change the default.

*Effort:* small, but do not start it without agreeing the opt-in framing.

## Non-gap, for completeness

CoolProp ships 137 pure-fluid JSONs (`dev/fluids/*.json`) against REFPROP 10's
~147.  This is the smallest difference on the list and is a data-curation
matter rather than a functionality gap.  No issue filed.

## Recommended sequencing

1. **`CoolProp-it6e.2` (ammonia-water)** first — highest impact per unit effort,
   and it converts a hard `factory()` failure into a working mixture.
2. **`CoolProp-it6e.5` REFPROP sublimation passthrough** — small, self-contained,
   mirrors existing code, good warm-up commit.
3. **`CoolProp-ft05` and `CoolProp-1gth`** (existing flash defects) before
   **`CoolProp-it6e.4`** (new flash pairs).
4. **`CoolProp-it6e.1` (mixture transport)** as its own project, once someone
   can commit to it.  Do not start it as a side quest.

Gaps 3, 6, 7 and 8 are unblocked and can be picked up independently; each needs
its own design doc first, because each turns on a modelling or architecture
decision that this assessment deliberately does not make.
