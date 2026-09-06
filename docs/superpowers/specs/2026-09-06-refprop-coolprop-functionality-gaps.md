# REFPROP-vs-CoolProp functionality gaps

**Date:** 2026-09-06
**Issue:** CoolProp-it6e (epic), CoolProp-it6e.1 .. .10 (individual gaps)
**Status:** assessment complete; each gap needs its own design + plan before implementation

## Goal

Establish, from the source rather than from folklore, which pieces of *core
REFPROP functionality* CoolProp's native `HEOS` backend does not provide, so
that the roadmap argues from evidence and so that non-gaps stop being
re-proposed.

Where a claim was checked and turned out to be **wrong**, it is recorded in
"Rejected claims" rather than deleted — the point of the document is to stop
the same wrong ideas recurring.

## Method, and what it can and cannot support

Two sides, anchored differently.  Read the distinction before quoting anything
here, because it governs how much weight each claim carries.

**CoolProp-side claims** are anchored to a file and line in this repository at
`6de0bf8c5`.  They are checkable by anyone with the tree.

**REFPROP-side claims** are claims about an external product and cannot be
verified from CoolProp's own source.  Where the claim is about REFPROP's
*exported API surface*, it is anchored to `REFPROP_lib.h` from the pinned
`REFPROP-headers` dependency (`cmake/dependencies.cmake:70-73`, tag
`b4faab1b73911c32c4b69c526c7e92f74edb67de`), which lists all 176 exported
routines and is in-tree after a configure at
`build_shared/_deps/refprop_headers-src/`.  That header proves a routine
*exists* or does not.  It says nothing about what a routine *does* internally —
so claims of the form "REFPROP's `TRNPRP` uses mixture ECS plus friction
theory" rest on REFPROP's documentation and are marked as such below.

**How the gap list was selected:** by scanning all 176 exported routines in that
header for CoolProp counterparts, plus the mixture-model and flash coverage
found by reading the HEOS backend.  A reader can therefore distinguish "not a
gap" from "not considered": anything in the header not listed below was checked
and found to have a counterpart.

## Rejected claims

These were considered and are **not** gaps.  Recorded so they are not raised again.

**Multiphase / VLLE / LLE / three-phase flash is NOT a REFPROP gap.**  Neither
library does it.  On the REFPROP side, the header exports 34 `*FLSH`/`*FL1`/
`*FL2` routines and no VLLE or three-phase routine of any name (searched for
`VLLE`, `LLE`, `3PH`, `THREEPH`; zero hits).  On the CoolProp side,
`src/Backends/Helmholtz/VLERoutines.cpp` is two-phase: bubble/dew plus a single
Rachford-Rice split.  Wanting VLLE is a legitimate goal, but it is a gap
against `teqp` or a process simulator, not against REFPROP, and must not be
justified by pointing at REFPROP.

*Correction to an earlier draft of this document:* that draft said CoolProp has
"no phase-stability-driven phase-count determination".  **That is false.**
Michelsen tangent-plane-distance stability analysis is implemented and
load-bearing: `StabilityRoutines::StabilityEvaluationClass`
(`VLERoutines.h:666`, `VLERoutines.cpp:2094` `check_stability_michelsen`,
`:2288` `minimize_tpd`), and `PT_flash_mixtures` gates its split on it —
`FlashRoutines.cpp:179-180` constructs the tester and sets
`do_twophase = !stability_tester.is_stable()`.  It is used at five further
sites.  The accurate statement is that the stability test is *binary* — it
decides one phase versus two and cannot discover a third — not that there is
none.

**Solid-fluid equilibrium is NOT a REFPROP gap.**  REFPROP's `MELTT`/`MELTP`
and `SUBLT`/`SUBLP` return correlated *boundary lines* for pure fluids; they are
not an SLE flash, and the header exports no SLE routine.  So "CO2 freeze-out in
cryogenic mixtures" is not available from REFPROP either.  What survives from
this area is narrow and is gap 5.

**Fugacity, chemical potential, virials, mixture critical points, phase
envelopes, reference-state setting, and BIP get/set are present natively.**
See `include/CoolProp/AbstractState.h` (`calc_fugacity`,
`calc_fugacity_coefficient`, `calc_chemical_potential`, `calc_Bvirial`,
`calc_Cvirial`), `HelmholtzEOSMixtureBackend.h` (`_calc_all_critical_points`,
`calc_criticality_contour_values`), and
`src/Backends/Helmholtz/PhaseEnvelopeRoutines.cpp`.  Cubics (SRK/PR) and PC-SAFT
cover the ground REFPROP's `PREOS` toggle addresses.

## The gaps

Ordered by impact, which is not the same as ordered by effort.

### 1. Mixture transport properties are a mole-fraction average, not a model — `CoolProp-it6e.1`, P1

The largest predictive gap.  In `HelmholtzEOSMixtureBackend.cpp`:

- `calc_viscosity` (:832) returns `exp(sum_i x_i * ln eta_i)` for mixtures,
  after `set_warning_string("Mixture model for viscosity is highly approximate")`.
- `calc_conductivity` (:1086) returns `sum_i x_i * lambda_i`, same warning.

Each per-component value is evaluated by constructing a fresh pure-fluid
`HelmholtzEOSBackend` at the *mixture* `(rho, T)`, which is not a state the pure
fluid is necessarily in.  There is no mixture critical enhancement.

The underlying routines cannot help: `TransportRoutines.cpp` throws
`"is only for pure and pseudo-pure"` in **15** distinct places, covering every
dilute, initial-density, higher-order, friction-theory, Chung and critical
contribution.

CoolProp does have `viscosity_ECS`/`conductivity_ECS` (`TransportRoutines.h:274`,
`:278`) with a conformal-state solver, but both index `HEOS.components[0]` and
are reached only from the pure branch.  REFPROP exports `TRNPRPdll`
(`REFPROP_lib.h:254`), which per REFPROP's documentation runs extended
corresponding states with a mixture conformal-state solve plus friction theory.

*Consequence:* for any real mixture, CoolProp's eta and lambda are not
predictive.  The warning string is set, but a caller who never reads
`get_warning_string()` gets a plausible-looking number with no accuracy claim
behind it.

*Effort:* large — the multi-month item.  First design question: does the
existing conformal-state solver generalise to mixtures, or must it be replaced?

### 2. Ammonia-water is absent entirely — `CoolProp-it6e.2`, P1

`dev/mixtures/mixture_binary_pairs.json` holds 888 binary pairs and **no
Ammonia/Water row** (checked by CAS `7664-41-7`/`7732-18-5` and by name).  The
consequence is not inaccuracy — the mixture cannot be constructed at all:
`set_mixture_parameters()` is called from the constructor
(`HelmholtzEOSMixtureBackend.cpp:131`) and throws at `MixtureParameters.cpp:590`.
So `AbstractState::factory("HEOS", "Ammonia&Water")` is a hard failure, and
every absorption-cycle use case on the HEOS backend fails at construction.
(`factory("REFPROP", ...)` works for a licensed user; this gap is HEOS-only.)

REFPROP carries the dedicated Tillner-Roth & Friend model for this system.

*Design question to settle first:* fit a BIP row into the existing multi-fluid
framework (cheap, approximate) versus implement TRF as a dedicated model
(faithful, more work).  Note the obstruction to (b) is **not** simply the set of
available departure-function forms — it is TRF's composition dependence and its
reducing function.  Establish the real obstruction during design rather than
assuming this one.

*Effort:* small if (a), medium if (b).  Highest impact-per-unit-effort here, and
the recommended first target.

### 3. Mixture surface tension throws — `CoolProp-it6e.3`, P2

`HelmholtzEOSMixtureBackend.cpp:712` throws
`NotImplementedError("surface tension not implemented for mixtures")`.
REFPROP exports `SURFTdll` (`REFPROP_lib.h:238`), documented as handling
mixtures.  A mixing rule must be chosen before any code is written; this is a
modelling decision, not plumbing.

*Effort:* medium, dominated by the model choice and its validation data.

### 4. Four mixture flash input pairs are unimplemented — `CoolProp-it6e.4`, P2

In `src/Backends/Helmholtz/FlashRoutines.cpp`, these throw
`"not ready for mixtures"`:

| Line | Routine | REFPROP counterpart |
| ---- | ------- | ------------------- |
| 510  | `DP_flash` | `PDFLSHdll` (`REFPROP_lib.h:180`) |
| 578  | `DQ_flash` | `DQFL2dll` (`:134`) |
| 616  | `HQ_flash` | `SATHdll` (`:214`), saturation state from `h` with a `kph` code |
| 649  | `QS_flash` | `SATSdll` (`:218`), likewise from `s` |

Note there is **no** `HQFLSH` or `SQFLSH` in REFPROP — HQ and QS are served by
the saturation routines, not by a flash.  An earlier draft claimed REFPROP
"ships the complete matrix" and cited `DEFLSH`/`DHFLSH`/`DSFLSH`/`PEFLSH`/
`TEFLSH`; those are real routines but they correspond to input pairs CoolProp
already supports, so they were not evidence for these four.

There are five `not ready for mixtures` throws in the file; the fifth (:728) is
inside the superancillary helper `resolve_T_via_superancillary`, reachable only
from DQ/HQ/QS and their `_with_guesses` variants, so the gap is four input
pairs, not five.

**Sequencing constraint:** `CoolProp-ft05` (P0 — mixture `HSU_P` flash returns
wrong `T`, silently, on master) and `CoolProp-1gth` (P1 — `DmassT`/`DmolarT`
can return a state at a different density than requested) are open defects in
mixture flash paths that *already exist*.  Adding four more input pairs onto a
base with known silent-wrong-answer bugs multiplies their surface area.  Fix
those first.

*Effort:* medium per pair; they share machinery.

### 5. No sublimation line, natively or by passthrough — `CoolProp-it6e.5`, P2

Two separate absences:

- **Native.** No general per-fluid sublimation ancillary.  There *is* one
  hard-coded ice curve — `psub_Ice` (`src/Ice.cpp:36`, declared
  `include/CoolProp/fluids/Ice.h:4`), reachable via
  `HAProps_Aux("psub_Ice", ...)` (`src/HumidAirProp.cpp:2404`) — so this is a
  missing *general facility*, not a total absence.  (An earlier draft argued
  from "zero occurrences of `sublim`"; that grep misses `psub_Ice` entirely and
  was weak evidence.  Melting, by contrast, is covered generally:
  `src/Backends/Helmholtz/MeltingCaloric.cpp`.)
- **Passthrough.** `src/Backends/REFPROP/REFPROPMixtureBackend.cpp` wraps
  `MELTPdll`/`MELTTdll` but **not** `SUBLPdll`/`SUBLTdll`
  (`REFPROP_lib.h:236-237`).  Even a user who owns REFPROP cannot reach its
  sublimation curve through CoolProp.

Scope note, so the rejected claim does not creep back: this is boundary *lines*
only.  Not SLE; REFPROP has no SLE either.

*Effort:* medium, and larger than "just mirror the melting code" suggests —
there is no sublimation API to mirror onto.  The melting passthrough hangs off
`AbstractState::calc_melting_line` (`AbstractState.h:650`), `melting_line()`
(`:1571`) and `has_melting_line()` (`:1564`); sublimation equivalents must all be
created first.  That new-virtual-plus-parameter-plumbing work is the same shape
as gap 6.

### 6. No dielectric constant — `CoolProp-it6e.6`, P3

REFPROP exports `DIELECdll` (`REFPROP_lib.h:128`).  CoolProp has no such
property; the only matches under `src/` belong to PC-SAFT's internal
electrolyte term (`PCSAFTBackend.cpp`, `PCSAFTBackend.h:37`), which is neither a
general-purpose dielectric constant nor exposed as a parameter.  Also searched
`permittivity` and `epsilon_r`: no hits.

Requires a new `DataStructures` entry, parameter plumbing, and per-fluid
correlations — the last being the real cost.

*Effort:* medium, dominated by data collection.

### 7. No runtime model switching — `CoolProp-it6e.7`, P3

REFPROP's `SETMODdll` (`REFPROP_lib.h:227`) selects an alternate published EOS
or transport model per fluid at runtime.

**The gap is narrower than "CoolProp only has one model per fluid", and the work
is in a different place than you would guess.**  23 of the 137 fluid JSONs
already ship **two** EOS entries (Ammonia carries Gao-2020 *and*
Tillner-Roth-1993; also D4, D5, Helium, HydrogenChloride, MD2M, MD3M, MD4M, …),
and `parse_EOS_listing`
(`src/Backends/Helmholtz/Fluids/FluidLibrary.h:527-531`) loads all of them into
`EOSVector`.  What hard-codes the choice is the accessor:
`include/CoolProp/CoolPropFluid.h:564` and `:567` return `EOSVector[0]`, and
nothing anywhere indexes `EOSVector[i>0]`.  So the alternate models are already
shipped and parsed; the missing piece is **selection plumbing above the loader**,
not loader work.

This interacts with open work: `CoolProp-f1ez` (nitrogen/argon viscosity
supersession blocked by the conductivity coupling) and `CoolProp-3bpg` (how to
adopt the 2022 R-134a viscosity given its ECS reference role) are genuine
forced-choice cases, and 3bpg's option (b) is literally per-fluid model pinning.
`CoolProp-9s9u` is *not* such a case — its items are blocked on missing
manuscript verification points and a NIST HTTP 503, which model switching does
not address.  Worth revisiting this priority if the supersession work stalls.

*Effort:* medium; a selection API plus a policy for what a selected model means
downstream (notably for ECS reference fluids).

### 8. No BIP estimation fallback — `CoolProp-it6e.8`, P3

CoolProp throws for a binary pair outside the 888-pair library unless the caller
supplies parameters by one of three explicit routes:
`apply_simple_mixing_rule(id1, id2, "linear" | "Lorentz-Berthelot")`
(`MixtureParameters.cpp:286`; rules at `:241`, `:251`),
`set_interaction_parameters` with binary-pair JSON (`:396`), or
`set_departure_functions` fed a REFPROP `HMX.BNC` dump
(`parse_HMX_BNC`, `:646`).  What is missing is specifically an *estimation*
scheme that fires automatically, as REFPROP's does when a pair is absent from
`HMX.BNC`.

**This is a difference in philosophy, not straightforwardly a defect.**
CoolProp's fail-loud default means a user never silently receives an estimated
number believing it fitted.  That is the better default and should be preserved.
The design question is whether to offer an *opt-in* estimation mode — not
whether to change the default.

*Effort:* small, but do not start without agreeing the opt-in framing.

### 9. No choked-flow / critical-flow factor — `CoolProp-it6e.9`, P3

REFPROP exports `CSTARdll` (critical flow factor, `REFPROP_lib.h:113`) and
`MASSFLUXdll` (choked mass flux, `:170`).  CoolProp has neither: `grep -rin
'cstar|choked|mass_flux'` over `src/` and `include/` returns nothing outside the
GERG backend.  Relevant to relief-valve and nozzle sizing.

*Effort:* medium; the thermodynamics is standard but needs an isentropic-choking
solver.

### 10. No AGA8 characterization — `CoolProp-it6e.10`, P3

REFPROP exports `SETAGAdll`/`UNSETAGAdll` (`REFPROP_lib.h:222`) to switch a
loaded mixture onto AGA8.  CoolProp carries AGA8 only as GERG *validation
reference data* (`src/Backends/GERG/GERGReferenceValues.h`); there is no AGA8
model a user can select.  Closely related to gap 7 — `SETAGA` is a special case
of runtime model switching — so sequence them together.

*Effort:* medium.

## Non-gap, for completeness

CoolProp ships 137 pure-fluid JSONs (`dev/fluids/*.json`, matching 137 entries in
`dev/all_fluids.json`) against REFPROP 10's roughly 147.  A data-curation matter
rather than a functionality gap.  No issue filed.

## Recommended sequencing

1. **`CoolProp-it6e.2` (ammonia-water)** — highest impact per unit effort, and it
   converts a hard `factory()` failure into a working mixture.
2. **`CoolProp-ft05` and `CoolProp-1gth`** (existing mixture-flash defects)
   before **`CoolProp-it6e.4`** (new flash pairs).
3. **`CoolProp-it6e.7` and `.10` together** — selection plumbing, then AGA8 as
   its first consumer.  Cheaper than their P3 suggests, because the alternate
   EOS data is already parsed and sitting in `EOSVector`.
4. **`CoolProp-it6e.1` (mixture transport)** as its own project, once someone can
   commit to it.  Do not start it as a side quest.

Gaps 3, 5, 6, 8 and 9 are unblocked and can be picked up independently.  Each
needs its own design doc first, because each turns on a modelling or
architecture decision that this assessment deliberately does not make.
