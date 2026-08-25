# Plan: RES on the REFPROP backend, and separating model verification from parameter transfer

## Context

The RES transport work is committed, but every accuracy question we asked came back
ambiguous because two separate questions were entangled in a single measurement:

- **(a) Is the model implemented correctly?** i.e. does the C++ reproduce the published
  REFPROP+RES values.
- **(b) Do the REFPROP-fitted parameters transfer to HEOS?**

Today the only available comparison is "C++ on HEOS vs a published number computed with
REFPROP", which mixes both. Every surprising result so far traced back to that conflation —
PROPANE's −4.9% was the critical-enhancement viscosity source, not `s_res`; BENZENE's −41.6%
was phase selection, not parameters; R161's "clean" sub-1% result was a dilute-dominated
sample point that could not see a 12% `s_res` error.

There is also a second, independent reason to want this: **the parameters were fitted against
REFPROP**, so the REFPROP backend is the configuration in which they are exactly valid rather
than an approximation. Making RES work there is a genuine user-facing feature, not just test
scaffolding.

Once RES runs on REFPROP, the two questions separate cleanly:

- (a) becomes: C++ REFPROP+RES vs the authors' own reference code, on any state points we like.
- (b) becomes: C++ HEOS+RES vs C++ REFPROP+RES — one implementation, two equations of state,
  isolating exactly the parameter-transfer error, over a grid per fluid rather than one point.

Exclusion decisions and Catch2 test scope then follow from (b) data instead of inference.

### Decisions taken

- RES on REFPROP is a **full user-facing feature**, not a test-only hook.
- The dilute-gas term uses the **best available source per backend**: REFPROP-native
  `eta0`/`lambda0` where the reference papers used it, the fitted polynomial elsewhere.

### Key discovery driving the dilute-term decision

The two papers made *opposite* choices, per code path — not per fluid as previously assumed:

| | pure fluid | mixture |
|---|---|---|
| Martinek (viscosity), `main.py:100` | **REFPROP-native eta0** (polynomial fallback) | polynomial + Wilke |
| Li (conductivity), `code_SI.py:287` | polynomial | **REFPROP-native lambda0** (polynomial fallback) |
| CoolProp C++ today | polynomial | polynomial + Wilke |

So the C++ diverges from the reference in exactly two of four paths. This explains both the
ARGON+NEON −11.6% mixture-conductivity gap and the pure-viscosity outliers concentrated in
D2/He/H2/Ne/ETHYLENE (the fluids where the polynomial fits REFPROP's native eta0 worst).
Matching the source per backend is what makes goal (a) *exact* rather than approximate.

---

## Cost assessment (why this is feasible)

`REFPROPMixtureBackend` derives from `AbstractState` directly — a sibling of
`HelmholtzEOSMixtureBackend`, so it inherits none of the RES machinery that cubics get free.

**Already available on REFPROP, no work:** all 13 thermodynamic calls the RES routines make —
`gas_constant`, `T`, `rhomass`, `molar_mass`, `smolar_residual` (properly overridden via
`PHIXdll`), `delta`, `dalphar_dDelta`, `d2alphar_dDelta2`, `T_critical`, `p_critical`,
`rhomass_critical`, `cpmass`, `cvmass`.

**Missing — four couplings, three trivial:**

1. Parameter storage: no `components` vector.
2. `mole_fractions` is zero-padded to 20 on REFPROP; must use `get_mole_fractions()`/`Ncomp`.
3. `is_pure_or_pseudopure` has no `AbstractState` equivalent.
4. `calc_alphar_deriv_nocache` (`TransportRoutines.cpp:1500-1502`) — the only real blocker.
   It evaluates alpha^r derivatives *off-state* at `t_ref`; REFPROP's `PHIXdll` reads the
   cached `_tau`/`_delta` and cannot do that without mutating the object.

---

## Approach

### Stage 1 — Hoist the RES layer to `AbstractState` (enabling refactor, no behaviour change)

The goal is one implementation of the RES API and dispatch, with backend-specific storage.

- `include/CoolProp/AbstractState.h`
  - Move the two enable flags (`viscosity_RES_enabled`, `conductivity_RES_enabled`, currently
    `HelmholtzEOSMixtureBackend.h:81-82`) here as protected members.
  - Add a protected virtual storage accessor, defaulting to "unsupported":
    `virtual std::vector<CoolPropFluid>& get_RES_components()`.
  - Reimplement the eight RES API methods (currently `AbstractState.h:744-786` as throwing
    virtuals, implemented only in `HelmholtzEOSMixtureBackend.cpp:4773-4869`) once here, in
    terms of that accessor. They stay virtual for compatibility but gain a working base body.
  - Add a shared protected helper producing the "RES requested but parameters missing" error,
    so the naming/`ValueError` behaviour is identical on every backend.

- `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.{h,cpp}`
  - Override the accessor to return the existing `components` vector.
  - Delete the eight now-duplicated implementations and the two flag members.
  - `calc_viscosity`/`calc_conductivity` dispatch (`:810-852`, `:1081-1123`) switches to the
    shared helper; behaviour must be unchanged.

- Cubics need no change — they inherit from HEOS and keep working through the same accessor.

**Verification:** `[RES]` results must be bit-identical to the current baseline
(19 cases, 13 passing, 25 failed assertions). This stage is pure refactor.

### Stage 2 — Retarget the routines to `AbstractState&`

- `src/Backends/Helmholtz/TransportRoutines.{h,cpp}`: change the signatures of
  `viscosity_RES` / `conductivity_RES` (`TransportRoutines.h:280,284`) to `AbstractState&`.
  Both must move together — `conductivity_RES` calls `viscosity_RES` internally.
- Replace `HEOS.mole_fractions` with `get_mole_fractions()`.
- Replace `HEOS.is_pure_or_pseudopure` with a composition-size check.
- Replace the `calc_alphar_deriv_nocache` pair with
  `first_partial_deriv(iDmass, iP, iT)` evaluated on a **cached scratch state** at
  `(t_ref, current rho)`. This is available on REFPROP through the inherited generic
  derivative machinery (`AbstractState.cpp:1245-1252`, `get_dT_drho` at `:987`), and is
  exactly what Li's SI does (`PropsSI('d(P)/d(Dmass)|T','T',Tref,'Dmass',rho,...)`).
  Follow the existing cached-reference-state idiom (`viscosity_ecs_reference_state`,
  `HelmholtzEOSMixtureBackend.cpp:858-862`) so the scratch state is built once per instance.

**Risk — call out explicitly:** this changes the numerical route to `dp/drho|t_ref` for HEOS
and cubics. It must be verified as equivalent, not assumed. Pin the current
near-critical conductivity values for a few fluids (PROPANE, R143A, CO2) *before* the change
and compare after; any movement beyond round-off is a defect in the substitution, not an
improvement.

**Verification:** `[RES]` unchanged again; plus the pinned near-critical values.

### Stage 3 — Wire up REFPROP

- `src/Backends/REFPROP/REFPROPMixtureBackend.{h,cpp}`
  - Add a carrier `std::vector<CoolPropFluid>` member and override the Stage-1 accessor.
  - Populate it at setup using the existing free function
    `overlay_RES_transport_by_name(eos_key, fluid, molar_mass_override)`
    (`FluidLibrary.h:1435`, `FluidLibrary.cpp:403-409`), mirroring
    `AbstractCubicBackend::setup()` (`CubicBackend.cpp:38-48`): copy name/aliases onto the
    carrier and pass molar mass explicitly, because the carrier has no EOS
    (`CoolPropFluid::molar_mass()` dereferences `EOSVector[0]` unguarded).
  - Override `calc_viscosity`/`calc_conductivity` (`REFPROPMixtureBackend.cpp:1228-1248`) to
    branch on the RES flags before falling through to the existing `TRNPRPdll` call.
  - Override `get_fluid_constant(i, imolar_mass)`, which `set_*_RES_parameters` needs and
    REFPROP does not currently implement.
  - `phase()` throws for REFPROP mixtures — ensure no RES path depends on it.

**Verification:** `use_viscosity_RES(True)` on a REFPROP state returns a value differing from
the plain `TRNPRPdll` result; the eight API methods round-trip.

### Stage 4 — Dilute-term source per backend

Make the choice explicit and testable rather than buried in a conditional: a protected virtual
on `AbstractState` supplying the dilute term, defaulting to the polynomial (+ Wilke for
mixtures), overridden on REFPROP to use the native call on the paths where the papers did —
native `eta0` for pure viscosity, native `lambda0` for mixture conductivity.

**Verification:** on REFPROP, ARGON+NEON mixture conductivity should collapse from −11.6%
toward the reference, and the pure-viscosity light-fluid outliers (D2, He, H2, Ne, ETHYLENE)
from up to 0.96% toward zero.

### Stage 5 — REFPROP parameter set

`dev/convert_RES_csv_to_json.py` currently applies two filters that must **not** apply to a
REFPROP set: it drops 29 of 151 fluids CoolProp's HEOS cannot construct (REFPROP can), and
drops `HEOS` keys that failed `s_res` vetting (meaningless where the fits are native).

- Emit a `"REFPROP"` key alongside `HEOS`/`PR`/`SRK`, populated for every fluid.
- Change the fluid-level filter to "constructible in *any* backend including REFPROP", and
  omit only the per-backend key that is unusable.
- Regenerate `include/res_transport_parameters_JSON.h`.

### Stage 6 — Restructured verification (the actual objective)

- **(a) Implementation correctness.** Extend `dev/compare_HEOS_vs_REFPROP_RES.py`, or add a
  sibling script, that generates reference values on a grid per fluid using the authors' own
  code — Li `code_SI.py::TC_RES(fluid, MoleFrac, p_Pa, T_K)` and Martinek
  `main.py::viscosity_RES(fluid, MoleFrac, p_Pa, T_K)` — and compares C++ REFPROP+RES against
  them. Target: reproduction to round-off. Any residual gap is an implementation defect.
- **(b) Parameter transfer.** On the same grid, compare C++ HEOS+RES against C++ REFPROP+RES.
  Report per fluid: max/median transport deviation, `s_res` deviation, **residual share** (the
  fraction of the property that actually depends on `s_res`), and phase agreement. Grid points
  must cover residual-dominated states, not just the near-ideal-gas points the published
  samples happen to sit at.
- Re-derive the exclusion lists in `dev/convert_RES_csv_to_json.py` from (b), replacing the
  current single-point rationale. Record the measured basis in the code comment, as now.

### Stage 7 — Test suite overhaul

Driven by `dev/RES_test_assessment.md`, in its stated priority order:

1. Replace the fail-open `catch (...) { WARN }` + `REQUIRE(ok >= N)` gate with an explicit
   expected-throw list plus hard assertions. **Blocking** — a regression currently degrades to
   a warning, and excluded fluids now throw by design, which is indistinguishable from breakage.
2. Impose the reference phase in the HEOS sample tests and in
   `transport_pr_pt_with_sat_fallback`. This is what makes BENZENE fail today at −41.6%
   despite being correct to −0.008% at the right phase.
3. Split reference-parity assertions from experiment-agreement assertions.
4. Re-tolerance or restrict the mixture-conductivity test, with the provenance reason recorded
   so it is not "fixed" by tightening later.
5. Add REFPROP+RES parity tests from Stage 6(a).

---

## Verification summary

| Stage | Check |
|---|---|
| 1 | `[RES]` bit-identical to baseline (19/13/25) |
| 2 | `[RES]` unchanged **plus** pinned near-critical values for PROPANE, R143A, CO2 |
| 3 | REFPROP RES returns a value distinct from `TRNPRPdll`; API round-trips |
| 4 | ARGON+NEON mixture conductivity and light-fluid pure viscosity converge to reference |
| 5 | JSON contains a `REFPROP` key for all 151 fluids; header regenerated |
| 6 | (a) reproduction to round-off; (b) per-fluid grid report with residual share |
| 7 | Suite fails loudly on regression; BENZENE passes |

Build and run throughout with the configuration now documented in CLAUDE.md:

```
cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON \
      -DCOOLPROP_STATIC_LIBRARY=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_catch --target CatchTestRunner --config Release -j8
./build_catch/Release/CatchTestRunner.exe "[RES]"
```

Python-side checks need `pip install .` (the RES API is exposed through
`src/nanobind_interface.cxx`; nanobind is the default interface). Note a full `~[slow]` run
crashes for reasons pre-dating and unrelated to this work — do not treat that as a regression
signal.

## Notes / open risks

- Stage 2's derivative substitution is the main regression risk to existing behaviour; it is
  the one place where "should be equivalent" must be demonstrated numerically.
- REFPROP is a global singleton (`check_loaded_fluid` re-invokes `SETUPdll`); a scratch state
  of the *same* fluid should not thrash, but this needs confirming under Stage 2.
- Stages 1–2 are pure refactor and could be committed independently of any REFPROP work.

---

# REVISIONS (post-approval design review)

A design review completed after approval corrected several points. The goals and the staging
shape stand; the mechanism changes. Where this section conflicts with the text above, **this
section wins**.

## R1. `PHIXdll` is NOT state-coupled — corrects the Cost assessment and Stage 2

The state coupling lives in CoolProp's wrapper, not in REFPROP. `REFPROPMixtureBackend::call_phixdll`
(`.cpp:2283`) reads `_tau`/`_delta` from members and normalises by them, but `PHIXdll` itself takes
`tau` and `del` as **input arguments** and is a pure function. A `call_phixdll(itau, idel, tau, delta)`
overload normalising by the arguments gives REFPROP a full off-state alpha^r evaluator in ~6 lines.

The premise "REFPROP cannot evaluate off-state" was wrong, and with it the entire scratch-state and
SETUPdll-thrashing concern. (Separately, `link_to_loaded_fluids()` at `.cpp:2372` already solves
second-state reuse without `SETUPdll`, with `build_saturation_shim()` at `.cpp:2384` as a working
precedent — not needed here, but it exists.)

## R2. The `first_partial_deriv` substitution is wrong — replaces Stage 2's derivative change

Three defects:

- `first_partial_deriv(iP, iDmass, iT)` returns `(dp/drho_mass)_T` — the reciprocal of what is needed.
- A scratch state at `(T = t_ref, rho)` evaluates alpha^r at `tau = T_reducing/t_ref`, while the
  current code uses `tau_ref = T_critical/t_ref`. These differ for 103 of 119 crit-enhancement
  fluids: mostly ~1e-7 relative, but **CYCLOPRO by 0.1%**. The swap therefore silently changes
  existing HEOS results — precisely what the plan forbids.
- The scratch state needs a real flash on every call because rho changes every call, so the
  `viscosity_ecs_reference_state` idiom does **not** amortise it. That idiom does not fit.

**Replacement** — a narrow, purpose-built virtual on `AbstractState`:

```cpp
virtual CoolPropDbl calc_drhomass_dp_constT_at(double T_eval)
```

- HEOS implementation = the current three lines moved **verbatim** → provably zero numerical change.
- REFPROP implementation = one `THERM2dll` call (already wired at `.cpp:1223`), a pure function of
  its arguments, yielding `dPdrho`/`drhodP` directly. No scratch state, no mutation.

Keep the current-T branch (`TransportRoutines.cpp:1494-1498`) exactly as it is; routing it through
the new virtual would force HEOS off its cached derivatives.

Rejected alternative: a general `calc_alphar_deriv_nocache` virtual on `AbstractState`. Possible via
R1, but it takes a composition argument REFPROP would have to ignore — a silent-wrong-answer trap.

## R3. A pre-existing HEOS bug surfaced — new, isolated Stage 7

`TransportRoutines.cpp:1500` uses `tau_ref = T_critical()/t_ref` while `delta_st = delta()` is
`rho/rho_reducing`. `dp/drho|T` is only correct with `tau = T_reducing/T`. Pre-existing, not
introduced here. Preserve it verbatim through the refactor stages so they remain provable no-ops,
then fix it as a separate measured change with a published diff table.

## R4. Storage: a light store, not a `CoolPropFluid` vector — replaces Stage 1

`AbstractState.h` must not include `CoolPropFluid.h`; that drags in Eigen and the superancillary
machinery. So the approved "protected virtual returning `std::vector<CoolPropFluid>&`" is not viable.

Instead, a new lightweight header `include/CoolProp/RESTransport.h` holding:

- `ViscosityRESData` / `ConductivityRESData`, **moved verbatim** out of `CoolPropFluid.h:314-340`
  (`CoolPropFluid.h` then includes the new header, so `fluid.transport.viscosity_res` keeps working
  and no call site changes);
- `struct RESComponentData { name; viscosity; conductivity; }`;
- `struct RESTransportStore { comps; viscosity_enabled; conductivity_enabled; dilute sources; }`.

`AbstractState` gains one protected `RESTransportStore _RES;`. The eight API methods stop being
virtual and become **concrete** in `src/AbstractState.cpp`: throw `NotImplementedError` when `comps`
is empty (preserving today's behaviour for IF97/Incompressible/Tabular/SBTL/PCSAFT), bound-check,
and clear `_viscosity`/`_conductivity` — which are already `AbstractState` members, so cache
invalidation turns out to be fully generic. Population entry point: protected `set_RES_components(...)`.

`CoolPropFluid::transport.*_res` is demoted to "the shipped record that seeds the store"; the store
is the live per-instance copy. HEOS already holds `components` **by value**, so this names an
existing truth rather than changing semantics.

## R5. Two silent hazards — both must land in the same commit as R4

**(a)** `HelmholtzEOSMixtureBackend::get_copy` (`.cpp:171`) constructs from `components`, so after the
move `SatL`/`SatV`/`TPD_state`/`critical_state` silently lose their RES parameters. Add
`ptr->_RES = this->_RES;`. **Nothing fails to compile** — the highest-risk line in the refactor.

**(b)** `src/Tests/CoolProp-Tests-RES.cpp` lines 447, 448, 473, 539, 540, 559, 560, 575, 576 poke
`get_components()[0].transport.*_res.provided = false` to simulate missing parameters. After the move
these still compile and run but affect nothing, so the guard tests would pass for the wrong reason.
Migrate them to a store accessor in the same commit.

## R6. Move RES out of `TransportRoutines` entirely — refines Stage 2

`TransportRoutines.h:4` includes `HelmholtzEOSMixtureBackend.h`, so leaving RES there forces the
REFPROP backend to include the whole Helmholtz backend. Move `TransportRoutines.cpp:1295-1525` (both
functions, `wilke_mix`, the constants) into a new TU `src/Backends/RES/RESTransport.{h,cpp}`
including only `AbstractState.h` and `RESTransport.h`. Update the six call sites:
`HelmholtzEOSMixtureBackend.cpp:820, 846, 1091, 1117` and `CubicBackend.cpp:858, 873`.

## R7. REFPROP `calc_conductivity` caching trap — Stage 3

`calc_viscosity` (`.cpp:1228-1248`) is one `TRNPRPdll` call populating **both** `_viscosity` and
`_conductivity`; `calc_conductivity` exploits that. A naive `if (RES_enabled) return ...` at the top
of `calc_viscosity` means RES-viscosity on with RES-conductivity off leaves `_conductivity`
unpopulated, and the return throws on an uncached `CachedElement`. Fix: extract a private
`call_TRNPRPdll(double& eta, double& tcx)` that always populates both, and give the two `calc_`
methods independent bodies. Test all four flag combinations.

## R8. Name resolution must not use `NAMEdll` — Stages 3 and 5

`hnam` is `character*12`. With the CoolProp-availability filter removed, seven REFPROP keys exceed
12 characters — `22DIMETHYLBUTANE`, `23DIMETHYLBUTANE`, `3METHYLPENTANE`, `CHLOROBENZENE`,
`ETHYLENEOXIDE`, `PROPYLENEOXIDE`, `VINYLCHLORIDE` — and would silently truncate to non-matching
keys. Instead capture the resolved `.FLD` stems already computed in `set_REFPROP_fluids`
(`.cpp:355-386`, the `resolved_names` local) into a member and use those.

`molar_mass_override` must be `> 0`, taken from `INFOdll` (`wmm/1000`). Predefined `.MIX` mixtures
collapse `fluid_names` to a single joined string — leave RES unsupported there in the first pass,
with a clear throw.

## R9. GERG / Peng-Robinson config guard — Stage 3

`REFPROP_USE_GERG` and `REFPROP_USE_PENGROBINSON` (`.cpp:406-411`, `:520-526`) mean REFPROP is not
evaluating the reference Helmholtz EOS, so the REFPROP-fitted `n_res`/`xita` are invalid. Set
`n_params_match_alpha = false` on both records at setup when either is active — this reuses the
existing fail-loud path with no new machinery.

## R10. Parameter filters apply per key, not per fluid — replaces Stage 5

```
always:                                          dilute, REFPROP, critical_enhancement
if coolprop_has_fluid(fluid):                    PR, SRK
if coolprop_has_fluid(fluid) and not excluded:   HEOS
```

More precise than "drop both filters" — the 29 non-constructible fluids still should not get dead
`PR`/`SRK` entries. Sections go 122 → 151 entries. Add a `REFPROP` entry to `SOURCES` for both
properties. Keep the `EXCLUDE` sets and `--keep-all-heos` as they are.

Note the `REFPROP` and `HEOS` blocks are numerically identical for the 122 overlapping fluids (same
source file). Ship the duplicate; revisit only if size becomes a problem (+~370 KB header).

## R11. The verification gate is a value dump, not "tests pass"

Before Stage 1, dump every computed value from the six sample-data cases to CSV. The gate for
Stages 1, 2 and 4 is that this CSV is **bit-identical**. Stage 5's gate is that the new JSON with the
`REFPROP` key stripped diffs empty against the old one.

## R12. New Stage 0 — purely additive, do first

- `REFPROPMixtureBackend::get_fluid_constant` via `INFOdll`. Needed by `set_*_RES_parameters`, which
  call it at `HelmholtzEOSMixtureBackend.cpp:4803, 4823`; REFPROP does not implement it today.
- The `call_phixdll(itau, idel, tau, delta)` overload per R1.
- Run the off-state purity probe: record `dalphar_dDelta()`, call the new overload and `THERM2dll` at
  a different T, re-read and assert bit-identical; then assert `THERM2dll`'s `dPdrho` equals
  `R*T2*(1 + 2*delta*phi_delta + delta^2*phi_deltadelta)` to ~1e-12. This also settles the
  `T_critical` vs `T_reducing` question of R3 empirically.

## R5 addendum (found during Stage 1 implementation)

R5(a) named only HelmholtzEOSMixtureBackend::get_copy. There is a SECOND copy path and it is
the one that actually crashed: AbstractCubicBackend::copy_internals (CubicBackend.cpp:774).
It constructs the new instance with an EMPTY cubic component list, so setup() skips the RES
seeding block and leaves _RES.comps empty, while set_alpha0_from_components() then pushes
entries into the HEOS components vector. set_cubic_alpha_C bounds-checked _comps.size()
(non-empty) while indexing _RES.comps (empty) -> segmentation fault in SatL/SatV/TPD_state.

Fixes: copy_internals copies _RES to the instance and to each linked state; the guard in
set_cubic_alpha_C bounds-checks _RES.comps.size(), the vector it actually indexes.

Why both copy paths were previously safe: the parameters rode along inside `components`, which
every copy mechanism already duplicated. Centralising the storage is what made the copy paths
load-bearing -- and there are two of them, not one.

## Progress log

Stage 0 DONE (verified): REFPROPMixtureBackend::get_fluid_constant via INFOdll;
  call_phixdll(itau,idel,tau,delta) overload. Probe tests confirm PHIXdll is PURE -- evaluating
  alpha^r at tau/2 left dalphar_dDelta, d2alphar_dDelta2, p, smolar_residual, tau and delta all
  bit-identical. R1 is empirically established.

Stage 1 DONE (gate passed): RESTransportStore on AbstractState; 8 API methods de-virtualised
  and made concrete; all 22 reader sites migrated; get_copy AND copy_internals both carry the
  store; 9 test pokes migrated. Found and fixed a segfault (see R5 addendum) and a cubic
  seeding-order bug that would have disabled RES for every cubic fluid.

Stage 2 DONE (gate passed, bit-identical): RES moved to src/Backends/RES/RESTransport.{h,cpp},
  backend-neutral (includes only AbstractState.h + RESTransport.h). Signatures take
  AbstractState&. mole_fractions -> get_mole_fractions(), is_pure_or_pseudopure -> N==1,
  the two calc_alphar_deriv_nocache calls -> calc_drhomass_dp_constT_at() whose HEOS body is the
  original expression VERBATIM (including the pre-existing T_critical/T_reducing inconsistency,
  preserved on purpose so the refactor stays checkable). RES registered in
  COOLPROP_ENABLED_BACKENDS. Gates: dev/RES_comparison/RES_baseline_stage{0,1,2}.txt

Stage 3 DONE (gate passed) + Stage 5 pulled forward, because Stage 3 is unverifiable without
  it: with no `REFPROP` key in the parameter JSON every fluid's `provided` stays false and every
  RES call on REFPROP throws. R10 was therefore implemented first.

  Parameter set (R10). `dev/convert_RES_csv_to_json.py` now emits the REFPROP-fitted block under
  TWO keys -- `REFPROP` (exact: that is the EOS it was regressed against) and `HEOS` (an
  approximation, gated by HEOS_*_EXCLUDE). Filters moved from per-fluid to per-key: `dilute` and
  `REFPROP` always; `PR`/`SRK` only when a native CoolProp backend can build the fluid; `HEOS`
  additionally not in the exclusion set. 122 -> 151 entries per property, 29 of them REFPROP-only.
  Gate (dev/RES_comparison/gate_json.py logic, run ad hoc): every pre-existing fluid is identical
  once the `REFPROP` key is stripped; every added fluid carries no native-backend key; `REFPROP`
  is present for all 151 and equals `HEOS` wherever `HEOS` survives. PASSED.

  Backend wiring. `resolved_fluid_names` member captures the .FLD stems already computed in
  set_REFPROP_fluids (R8 -- NAMEdll is never used, so the seven >12-char keys survive);
  `setup_RES_transport()` seeds the store from a bare CoolPropFluid carrier with the molar mass
  passed explicitly from INFOdll; `calc_drhomass_dp_constT_at` via one THERM2dll call;
  calc_viscosity/calc_conductivity given independent bodies over a shared `call_TRNPRPdll` (R7).

  Two placement decisions worth keeping:
  - setup_RES_transport() is called from construct(), NOT from set_REFPROP_fluids(). The latter is
    re-invoked by check_loaded_fluid() on every property call, so seeding there would silently
    discard anything the user had set via set_*_RES_parameters().
  - The opportunistic co-caching in calc_viscosity/calc_conductivity is suppressed whenever the
    OTHER property is RES-driven. Without that, a TRNPRPdll value overwrites a RES value and
    leaks out of the next read -- the mirror image of the R7 trap, and invisible unless the two
    properties are read in both orders. The test does read them in both orders.
  - Predefined `.MIX` mixtures clear resolved_fluid_names: SETMIXdll reports Ncomp == 1 for the
    whole mixture, so there is no per-component identity. RES reports unsupported there.

  R9 guard verified through REFPROP_USE_PENGROBINSON rather than REFPROP_USE_GERG, deliberately:
  PREOSdll(0) is called unconditionally on every fresh setup so the flag self-restores, whereas
  GERG04dll is only ever called to ENABLE GERG and would poison every REFPROP test that followed.
  The test asserts the restore rather than assuming it, and matches the exception message
  ("different alpha function") rather than just its type.

### Stage 3 result: implementation correctness is now separated from parameter transfer

The parity sweep (`[RES_refprop_parity]`, a hidden measurement test) runs all 276 published
sample points on REFPROP+RES. **Every one evaluates -- zero fluids drop out**, including the 29
REFPROP-only fluids and all seven >12-character keys. Because the reference columns were computed
with the same EOS, what remains is implementation error alone:

| path | mean abs dev | max abs dev | worst |
|---|---|---|---|
| binary viscosity | 0.0003 % | **0.001 %** | ARGON+NITROGEN |
| pure viscosity | 0.032 % | 0.96 % | D2, ETHYLENE, HELIUM, CO, H2 |
| pure conductivity | 0.061 % | 4.88 % | PROPANE, R143A |
| binary conductivity | 2.04 % | 11.6 % | ARGON+NEON |

Binary viscosity is the one path where both papers and this code use the fitted polynomial plus
Wilke, and it reproduces to 1e-5. That is the clean control: the residual-entropy machinery, the
mixing rule, the units and the component ordering are all correct.

The other three rows are the dilute-source divergence predicted in "Key discovery", now measured
on a single EOS instead of inferred across two:
- pure viscosity outliers are exactly the light fluids where the polynomial fits REFPROP's native
  eta0 worst -- Martinek uses native eta0 on this path;
- binary conductivity is the loosest because Li uses native lambda0 on this path;
- PROPANE / R143A are a separate cause: the Olchowy-Sengers enhancement consumes a viscosity, and
  the reference feeds it REFPROP's native value where this code feeds it the RES viscosity.

**This retires the ambiguity the whole plan was written to resolve.** PROPANE -4.9% and R143A
+1.4% appear at full strength on REFPROP, so they were never parameter transfer. Stage 6(b) can
now attribute whatever HEOS-vs-REFPROP gap remains to the parameters alone.

  Gates run: `[RES]` value dump diffs against RES_baseline_stage2.txt with ZERO removed or
  changed lines (only the two summary counts move) -- masking test-file line numbers as well as
  heap addresses, since inserting tests shifts them. 25 pre-existing failures still exactly 25.
  `[REFPROP]` 1484/1484 (was 1114; +370 new). `[viscosity],[conductivity]` 1131 unchanged.
  cubics 3698 unchanged.

Stage 4 DONE (gate passed, HEOS/cubics bit-identical).

  The plan's premise was right in shape but wrong in two details, both found by reading the
  reference sources rather than inferring from deviations:

  1. Neither paper calls a "native eta0/lambda0 function".  Both flash to a near-zero pressure and
     read the ORDINARY transport property there -- Martinek `PropsSI('V', T, P=1e-9, TheEOS::f)`
     (main.py:100), Li `PropsSI('L', T, P=0.1, 'REFPROP::'+mix)` (code_SI.py:287).  Different
     pressures, and Li hardcodes REFPROP regardless of the backend in use.  Both wrap the call in
     try/except and fall back to the polynomial.
  2. There is a THIRD source choice the plan never mentioned, and it turned out to matter more
     than the dilute term for pure conductivity: **the viscosity inside the Olchowy-Sengers
     enhancement**.  Li feeds it REFPROP's NATIVE viscosity (code_SI.py:116,
     `PropsSI('V', T, Dmass, 'REFPROP::'+f)`), not the RES viscosity being computed.  The existing
     C++ comment already flagged this as a deliberate deviation; it is the ENTIRE remaining
     PROPANE (-4.9%) and R143A (+1.4%) gap.

  Implementation.  One hook rather than two: `AbstractState::calc_transport_native(key,
  rhomolar_eval, value)` -- the backend's own transport model at the current T and composition and
  an ARBITRARY density.  rho = 0 gives the dilute-gas limit; the current density gives the
  enhancement viscosity.  REFPROP implements it in one TRNPRPdll call, which is pure in (T, rho, x),
  so no flash, no scratch state, no mutation -- and it reaches the same numbers the reference gets
  the expensive way.  It returns false (not throws) when REFPROP has no model, and checks
  isfinite && > 0 because REFPROP signals "unavailable" with a large negative sentinel and
  ierr == 0 in some builds.

  Policy lives in one place, `dilute_from_backend()` in RESTransport.cpp, plus the enhancement
  block.  AUTO reproduces the papers per CODE PATH:
      viscosity     pure -> native      mixture -> polynomial + Wilke
      conductivity  pure -> polynomial  mixture -> native (the whole mixture, not per-component)
      enhancement viscosity -> native where available, else the RES viscosity
  AUTO falls back silently when the backend has no native model, mirroring the papers' try/except
  -- and that fallback is exactly what keeps HEOS and the cubics unchanged.  An EXPLICIT
  RES_DILUTE_BACKEND_NATIVE / RES_ENH_VIS_BACKEND_NATIVE does NOT fall back; it throws.  Quietly
  substituting a different model than the caller asked for is how a wrong number ships.

  Three new public setters: set_viscosity_RES_dilute_source, set_conductivity_RES_dilute_source,
  set_conductivity_RES_enhancement_viscosity.  The self-consistent all-RES behaviour is still
  available via RES_DILUTE_POLYNOMIAL + RES_ENH_VIS_RES; it is just no longer the default, because
  the default should reproduce the papers.

### Stage 4 result

| path | Stage 3 max | Stage 4 max | mean |
|---|---|---|---|
| pure viscosity | 0.96 % | **0.011 %** | 0.0005 % |
| pure conductivity | 4.88 % | **0.042 %** | 0.006 % |
| binary viscosity | 0.001 % | 0.001 % (unchanged, correctly) | 0.0003 % |
| binary conductivity | 11.6 % | **0.998 %** | 0.135 % |

Binary viscosity is unchanged BY DESIGN -- both papers and this code use the polynomial there --
which is the control proving the change touched only the intended paths.

  One outlier survives and it is a **feature gap, not an error**: BUTANE+METHANE at -1.0 %.  Li
  applies a critical enhancement to MIXTURES too (`Olchowy_critical_enhancement_mix`,
  code_SI.py:152); this implementation applies it only to pure fluids.  That sample sits at
  T/Tc = 0.974, rho/rhoc = 1.62 -- verified inside Li's near-critical window -- so the missing
  enhancement is worth about that 1 %.  Every other binary is within 0.07 %.  FOLLOW-UP: implement
  the mixture enhancement, then tighten the binary-conductivity bound from 1.2 % to 1e-3.

  Gates: `[RES]~[REFPROP]` value dump vs RES_baseline_stage3.txt -- added lines were the filter
  line and the two summary counts ONLY, i.e. every HEOS and cubic value bit-identical, proving the
  AUTO fallback is a true no-op there.  `[RES]` 33 cases / 1080 assertions, still exactly the same
  6 failing cases and 25 failing assertions.  `[REFPROP]` 1498, `[viscosity],[conductivity]` 1131,
  cubics 3698, `[SBTL]` 5476 -- all pass.  clang-format: every touched file back to its exact HEAD
  violation count (0 new), using dev-only line-range formatting since three of these files carry
  pre-existing violations.

  New baselines: dev/RES_comparison/RES_baseline_stage4.txt, parity_stage4.txt.

Stage 4b DONE -- mixture critical enhancement, and a parameter-sourcing bug it exposed.

  Mixture critical enhancement implemented (Li 2024 `Olchowy_critical_enhancement_mix`,
  code_SI.py:152), gated by a new `RESMixtureEnhancement` policy: AUTO / OFF / ON.

  **AUTO is deliberately OFF on CoolProp's own backends.** The enhancement needs the MIXTURE
  critical point, and HEOS/cubics have to SOLVE for it -- slow, and not reliably convergent --
  while the physical case for a critical enhancement in mixtures is not well established.  The
  policy hangs off a new `AbstractState::calc_has_direct_critical_point()`, false by default and
  true on REFPROP (CRITPdll reports it directly).  So REFPROP reproduces the published mixture
  values out of the box, and nothing on HEOS silently acquires a critical-point solve.  Pure
  fluids are unaffected by the setting; only the mixture path is gated.

  Mixing follows Li exactly: mole-fraction linear on gamma, xi0, Gamma, and on **1/q_D** (the
  table stores qDinv, so the reciprocal is what gets mixed); `t_ref = 1.5*Tc` for mixtures rather
  than the per-fluid table value, which has no mixture entry.  A failed critical-point lookup
  under AUTO returns zero enhancement (Li wraps the whole thing in try/except); under an explicit
  RES_MIX_ENH_ON it propagates, because there the caller asked for it.

### The R_D trap -- worth remembering

  Generalising the pure-fluid block to mixtures meant sourcing R_D and gamma per component instead
  of the hardcoded 1.02 / 1.239.  That made pure conductivity WORSE: 0.042 % -> 0.314 %, all of it
  PROPANE, whose enhancement is ~32 % of lambda at its sample point.

  Our parameter table is byte-identical to Li's and says R_D = 1.030 for PROPANE.  But
  `get_paramters()` **overwrites it with a flat 1.02 for every fluid** (code_SI.py:65) and never
  reads the R_D column at all.  `gamma` on the next line IS read per fluid.  So the shipped R_D
  column is dead data, and reproducing the published values requires ignoring it.

  Correct answer: R_D = 1.02 flat, gamma per fluid.  The original code had R_D right and gamma
  wrong, and the two errors had been cancelling.  Note the trap: the table looked authoritative,
  the deviation looked like a bug in our enhancement formula, and a plausible-looking "fix" (keep
  hardcoding both) would have preserved a real error.  Only running the reference implementation
  and printing what it actually passed settled it.

### Result -- all four paths now reproduce the published values

| path | Stage 3 | Stage 4 | Stage 4b |
|---|---|---|---|
| pure viscosity | 0.96 % | 0.011 % | 0.011 % |
| pure conductivity | 4.88 % | 0.042 % | 0.042 % |
| binary viscosity | 0.001 % | 0.001 % | 0.001 % |
| binary conductivity | 11.6 % | 0.998 % | **0.018 %** |

BUTANE+METHANE went from -0.998 % to +0.001 %, and is now the *best* of the eight binaries.

  HEOS/cubic impact of the gamma fix, measured rather than assumed: exactly three sample values
  move, ETHANE by -0.0075 % and SF6 by +0.0010 %.  No test changes state; [RES] still fails the
  same 6 cases and 25 assertions.  Mixtures on HEOS are untouched, since AUTO leaves the
  enhancement off there.

  Gates: [RES] 36 cases / 1089 assertions, same 6 failures / 25 failing assertions.  [REFPROP]
  1505, [viscosity],[conductivity] 1131, cubics 3698, [SBTL] 5476 -- all pass.  clang-format: no
  new violations anywhere; RESTransport.cpp actually improved 57 -> 43 because the rewritten
  enhancement block replaced non-conforming code.
  Baselines: RES_baseline_stage4b.txt, parity_stage4b.txt.

Stage 6 DONE -- both halves, and the exclusion lists re-derived from measurement.

  Three new dev tools, one shared grid:
    dev/RES_grid_build.py       -> RES_comparison/grid_points.csv  (3828 pts over 147 fluids)
    [RES_grid] Catch2 harness   -> RES_comparison/grid_cpp.csv     (both backends, same points)
    dev/RES_grid_report.py      -> (b) parameter transfer
    dev/RES_reference_check.py  -> (a) vs the authors' own code

  The grid is built in REDUCED coordinates so every fluid is sampled where the residual term
  actually dominates, and phase ambiguity is designed out (supercritical isotherms, compressed
  liquid at >=3x p_sat, superheated vapour at 0.3x p_sat).  Only 2 of 3828 points had the two
  backends disagree on phase, versus the phase artefacts that dominated the old single-point
  comparisons.

### (a) Implementation correctness -- C++ REFPROP+RES vs Martinek/Li reference code

  median **0.000026 %** (viscosity) and **0.000012 %** (conductivity) over 587 sampled points --
  round-off, which was the stated target.

  THE FIRST RUN OF THIS WAS WRONG, and the way it was wrong is the lesson.  The [RES_grid] harness
  pins both backends to the fitted dilute term, which (b) requires; (a) then unknowingly compared
  our POLYNOMIAL against the reference's NATIVE dilute term and reported medians ~100x too large.
  Worse, the gap happened to equal the poly-vs-native gap for the fluids inspected, which read as
  confirmation of a plausible but entirely wrong story ("the reference fell back to the
  polynomial").  It was only disproved by actually calling PropsSI at 1e-9 Pa and watching it
  succeed.  The harness now emits THREE columns -- REFPROP (defaults, for (a)), REFPROP_pinned and
  HEOS (for (b)) -- so neither question can silently read the other's data.

  Two outlier classes remain, both verified rather than inferred:
  - **HELIUM, 4 points, 46-57 %.**  Here the fallback really does fire, and the reason is the
    ENTRY POINT, not the correlation.  We call TRNPRPdll directly at rho = 0, no flash involved.
    The reference asks PropsSI for a property at p = 1e-9 Pa, which must run a PT flash first, and
    at that pressure the flash degenerates -- REFPROP then reports `[TRNPRP error 561] Pure fluid
    correlation produced an erroneous value` and Martinek's try/except falls back to the
    polynomial.  Verified by sweeping the input: at T = 5.2992 K the same query succeeds at every
    DENSITY down to 1e-11 kg/m3 and at every PRESSURE from 1e-3 Pa up, always returning
    1.315073e-6 Pa.s -- which matches our rho = 0 value of 1.315075e-6 to six digits.  Only
    p = 1e-9 Pa fails.  (Li uses 0.1 Pa and never hits this.)  The polynomial the reference falls
    back to gives 3.13e-6 Pa.s, 2.4x the true value, since 5 K is far outside its fit range.  So
    our number is the correct one here.  RES on helium at 4 K is dubious regardless.
  - **The critical singularity, ~30 conductivity points.**  At rho/rhoc = 1, T/Tc = 1.02 the
    enhancement carries 1/(dp/drho) and `arg` is a difference of two large near-equal terms, so
    it is ill-conditioned rather than wrong.  Excluding those 50 points the conductivity max falls
    from 58 % to 6.4 % and the median to 0.000011 %.  Not decomposed further.

### (b) Parameter transfer -- C++ HEOS+RES vs C++ REFPROP+RES

  median 0.00006 % (viscosity) and 0.00005 % (conductivity) over 3243 comparable points: for 133
  of 147 fluids the REFPROP-fitted parameters transfer essentially exactly.

  **New exclusion lists: the SAME 14 fluids for both properties.**
      BENZENE, D5, HEPTANE, MD3M, MD4M, R1123, R1224YDZ, R1233ZDE,
      R1234YF, R1243ZF, R13, R161, R41, VINYLCHLORIDE
  Was 6 and 6, and the new lists are supersets of the old -- as expected, since the old ones were
  derived from single sample points that could not see most of the error.

  Criterion: >1 % transport deviation at residual-dominated states away from the critical
  enhancement, OR >5 % deviation in s_res itself there.  The second clause is load-bearing:
  MD4M's conductivity agrees to 0.27 % while its s_res is off by 27 %, so the sampled states are
  simply insensitive.  Both clauses independently select the same 14 fluids for viscosity and
  conductivity, which is the coherent answer -- s_res is the one input both models share.

  **BENZENE is now excluded, reversing an earlier call.**  It was kept on the strength of one
  sample point agreeing to -0.008 % once the phase was imposed.  On the grid it deviates 1.68 %
  at residual-dominated states with s_res off by 5.6 %; that single point was not diagnostic.

### Three ways this measurement was wrong before it was right -- all silent

  1. **Default source policies.** With the Stage-4 defaults left alone, REFPROP takes the dilute
     term and enhancement viscosity from its own transport model and HEOS cannot, so the
     comparison measured that CHOICE, not the parameters: ~40 fluids past the threshold, with
     s_res agreeing to 0.001 %.  The harness now pins both backends to the fitted model.
  2. **A grid point exactly on the enhancement gate.** rho/rhoc = 2.0 is where Li switches the
     enhancement off.  The backends' rhoc differ slightly, so 172 points had it on for one and
     off for the other -- a pure step artefact that inflated the conductivity list from 12 to 30.
     Grid moved to 1.80 / 2.30, and the report now discards any straddling point regardless.
  3. **Reference units.** viscosity_RES returns Pa.s; only main.py's sample table is in microPa.s.
     Got a uniform 1e6 factor, which at least failed loudly.

  In each case the wrong answer looked entirely plausible.  What caught all three was checking
  that the CAUSE was consistent with the effect -- a 5 % conductivity gap alongside a 0.001 %
  s_res gap cannot be a parameter-transfer failure, whatever the threshold says.

  Test-suite consequence: the R41 case in "RES conductivity does not throw a viscosity error..."
  relied on R41 having conductivity parameters but not viscosity ones.  With the lists now equal
  no such fluid exists, so the case is CONSTRUCTED via RES_data_mutable() instead of found.  Its
  own premise guard is what caught this -- it failed loudly rather than passing vacuously.

  Gates: [RES] 36 cases / 1071 assertions, 6 failing cases and 23 failing assertions (was 25 --
  two fewer because the newly excluded fluids now throw and are skipped by the sample tests
  rather than failing an assertion).  [REFPROP] 1508, [viscosity],[conductivity] 1131, cubics
  3698, [SBTL] 5476 -- all pass.  clang-format unchanged at HEAD's count.

NEXT: Stage 7 -- test-suite overhaul, D1 (the fail-open catch(...) + REQUIRE(ok >= N) gate) first;
  it is now actively misleading, since 14 excluded fluids throw by design and are indistinguishable
  from breakage.  Also outstanding: the converter invariant gate (see the Stage-5 discussion --
  structural checks plus a regenerate-and-diff, NOT cross-key equality, which a future HEOS-specific
  refit would legitimately break).

---

# RESUME HERE (written before context compaction)

## Where things stand

Branch `dev_RES`. Two commits of RES work:

- `d3bf7c7c` feat(transport): the RES model itself (HEOS + cubics)
- `3e10cdcf` refactor(transport): Stages 0-2 below -- RES made backend-neutral

**Stages 0, 1, 2 are DONE and gate-verified. Stage 3 is next.**

Uncommitted: `.claude/settings.json` (unrelated agent config). `dev/RES_comparison/` is
UNTRACKED -- it holds the gate baselines and would be lost by a clean checkout.

## The verification gate (read this before changing anything)

Stages 1, 2 and 4 must be *provable no-ops* for HEOS and cubics. The gate is NOT "tests pass" --
the pass count is unchanged by definition. It is a bit-identical dump of every computed value:

    ./build_tests/Release/CatchTestRunner.exe "[RES]" --success > after.txt 2>&1
    norm(){ grep -v "^Randomness seeded" "$1" | sed -E 's/0x[0-9a-f]+/0xPTR/g; s/[.]cpp\([0-9]+\)/.cpp(LINE)/g'; }
    diff <(norm dev/RES_comparison/RES_baseline_stage6.txt) <(norm after.txt)

Two masks, both required. Heap addresses in `dynamic_cast != nullptr` assertions vary per run.
Source line numbers shift whenever a test is INSERTED above an existing one, which produces
hundreds of spurious diff lines that hide the real ones -- that mask was added in Stage 3 after
exactly that happened. When a stage adds tests, the gate is not "empty diff" but **no removed or
changed lines**: `diff ... | grep "^<"` must return only the two summary counts.

Current gate file: `RES_baseline_stage6.txt`. Current expected [RES] result: **36 cases, 30
passed, 6 failed; 1071 assertions, 1048 passed, 23 failed.** The 6 failing cases are the HEOS
sample-data comparisons, with documented causes -- see `dev/RES_test_assessment.md`. They are NOT
regressions, and the failed-assertion count must stay at exactly 25.

This gate earned its keep: it caught a segfault in Stage 1 that a green build and an unchanged
pass count both hid.

## Environment (non-obvious, cost real time to discover)

- **Build** (MSVC). `-DCOOLPROP_STATIC_LIBRARY=ON` is REQUIRED -- it is the only thing that puts
  `/MD` into the flags; without it `cl` defaults to `/MT` and the link fails against the `/MD`
  Catch2. Documented in CLAUDE.md.

      cmake -B build_tests -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON \
            -DCOOLPROP_STATIC_LIBRARY=ON -DCMAKE_BUILD_TYPE=Release
      cmake --build build_tests --target CatchTestRunner --config Release -j8

  `CMAKE_BUILD_TYPE` is ignored by the VS generator; `--config Release` is what matters. Binary
  lands in `build_tests/Release/CatchTestRunner.exe`. Use `build_tests`, not `build_catch`.

- **Python**: `C:\nospace\miniconda3\envs\devCP\python.exe`. REFPROP is installed and on PATH, so
  `[refprop]` tests RUN rather than skip.

- **Python wrapper**: `pip install .` rebuilds from the working tree. The DEFAULT interface is
  **nanobind** (`src/nanobind_interface.cxx`), NOT the legacy Cython wrapper -- editing only the
  Cython files silently has no effect. Both are bound for the RES API.

## Exit-code and tooling traps (all of these bit at least once)

- Catch2's exit code is NOT the failed-assertion count, and a crashing run can exit non-zero with
  no summary printed. Always redirect to a file and grep for `test cases:` / `assertions:`.
  Exit 139 = segfault. A run that prints "All tests passed" prints NO `test cases:` line.
- NEVER pipe the runner through `tail` -- the pipe's exit status replaces the runner's.
- `grep -c` returning 0 matches exits 1, which fails the whole command. Append `; true`.
- Bash heredocs have failed twice on markdown content in this repo. Use the Write tool for
  markdown/C++ blocks, or write to a scratch file and `cat >>`.
- Bulk string replacement is dangerous here: replacing `transport.conductivity_res` also matched
  inside `transport.conductivity_residual` (an unrelated member) and produced `conductivityidual`.
  Grep for collisions before any bulk edit.

## Stage 3 -- REFPROP wiring (DONE -- see the Progress log above)

Goal: REFPROP carries the RES store and dispatches to `RESTransport::viscosity/conductivity`.
Everything structural is already done; this is wiring. Follow the plan's R7, R8, R9 exactly.

Already in place from Stage 0 (committed, tested):
- `REFPROPMixtureBackend::get_fluid_constant(i, param)` via `INFOdll` -- the RES setters need it.
- `REFPROPMixtureBackend::call_phixdll(itau, idel, tau, delta)` -- public, off-state alpha^r.
  Proven pure by test: evaluating at tau/2 leaves dalphar_dDelta, d2alphar_dDelta2, p,
  smolar_residual, tau, delta bit-identical.

To do:
1. Add a `std::vector<CoolPropFluid>` carrier + seed the store at setup via
   `overlay_RES_transport_by_name(eos_key, carrier, molar_mass_override)`
   (`FluidLibrary.h:1435`, `FluidLibrary.cpp:403`), mirroring `AbstractCubicBackend::setup()`
   (`CubicBackend.cpp:38-58`). `molar_mass_override` MUST be > 0 (from INFOdll, wmm/1000) --
   `CoolPropFluid::molar_mass()` dereferences an empty `EOSVector` otherwise.
   Then call `set_RES_components(...)`.
2. Implement `calc_drhomass_dp_constT_at(T_eval)` for REFPROP via ONE `THERM2dll` call (already
   wired at `REFPROPMixtureBackend.cpp:1223`, in `calc_PIP`); take `dPdrho`/`drhodP` and convert.
   Pure in its arguments -- no scratch state, no mutation.
3. **R7 -- the caching trap.** `calc_viscosity` (`.cpp:1228-1248`) is ONE `TRNPRPdll` call that
   populates BOTH `_viscosity` and `_conductivity`; `calc_conductivity` exploits that by calling
   it. A naive `if (RES_enabled) return ...` at the top means: RES-viscosity ON + RES-conductivity
   OFF leaves `_conductivity` unpopulated and the return THROWS on an uncached CachedElement.
   Fix: extract a private `call_TRNPRPdll(double& eta, double& tcx)` that always populates both,
   then give the two calc_ methods independent bodies. Test all FOUR flag combinations.
4. **R8 -- do NOT resolve names via `NAMEdll`.** `hnam` is `character*12`; seven RES keys exceed
   12 chars (22DIMETHYLBUTANE, 23DIMETHYLBUTANE, 3METHYLPENTANE, CHLOROBENZENE, ETHYLENEOXIDE,
   PROPYLENEOXIDE, VINYLCHLORIDE) and would silently truncate to non-matching keys. Capture the
   resolved .FLD stems already computed in `set_REFPROP_fluids` (`.cpp:355-386`, local
   `resolved_names`) into a member and use those.
   Predefined `.MIX` mixtures collapse `fluid_names` to one joined string -- leave RES unsupported
   there in the first pass, with a clear throw.
5. **R9 -- GERG / Peng-Robinson guard.** `REFPROP_USE_GERG` / `REFPROP_USE_PENGROBINSON`
   (`.cpp:406-411`, `:520-526`) mean REFPROP is not evaluating the reference Helmholtz EOS, so the
   REFPROP-fitted n_res/xita are invalid. Set `n_params_match_alpha = false` on both records at
   setup when either is active -- reuses the existing fail-loud path.
6. `phase()` throws for REFPROP MIXTURES -- ensure no RES path calls it.

Stage 3 has no bit-identical gate (it adds behaviour). Its gate is: REFPROP+RES returns a value
distinct from plain `TRNPRPdll`, the 8 API methods round-trip, and the published `vis_res`/`TC_RES`
columns reproduce to ~1% with the polynomial dilute term (Stage 4 tightens this to <0.05%).

## Remaining stages (detail in this file above)

- **Stage 4** DONE -- see the Progress log.
- **Stage 5** DONE (pulled forward into Stage 3) -- see the Progress log.
- **Stage 6** DONE -- see the Progress log.
- **Stage 7** optional, isolated: fix `tau_ref = T_critical/t_ref` -> `T_reducing/t_ref`. Watch
  CYCLOPRO (0.1%); ~1e-7 for the rest.

## Reference implementations (for Stage 6)

Both are callable at arbitrary (p, T), on the network share:

    ...\Stoffdaten\CoolProp\RES\Li_2024_SI\code_SI\code_SI.py
        TC_RES(fluid, MoleFrac, p_Pa, T_K)
    ...\Stoffdaten\CoolProp\RES\Martinek_2025_SI\supporting_information\code_SI\main.py
        viscosity_RES(fluid, MoleFrac, p_Pa, T_K)

Share root: `\\sccfs.scc.kit.edu\OE\IBPT\Groups\Grohmann\Users\Reichert\Stoffdaten\CoolProp\RES`

The two papers differ PER CODE PATH (not per fluid) on the dilute term -- Martinek uses REFPROP's
native eta0 for PURE fluids, Li uses REFPROP's native lambda0 for MIXTURES; the C++ uses the
polynomial everywhere. This is why ARGON+NEON mixture conductivity is -11.6% and why the
light-fluid pure viscosities (D2, He, H2, Ne, ETHYLENE) are the worst outliers.

## Known pre-existing issues -- NOT regressions, do not chase

- A full `~[slow]` run terminates abnormally (exit 127, no summary). Reproduces with ALL RES tests
  excluded and in the user's own build config, so it predates and is unrelated to this work.
  Localising it is blocked on Catch2's `-f` rejecting test names containing `,` `[` `(`.
- `dev/RES_test_assessment.md` documents 5 design defects in the [RES] suite. D1 (the fail-open
  `catch(...)` + `REQUIRE(ok >= N)` gate) is blocking and should be fixed before relying on the
  suite. D2 (no phase imposition) is why BENZENE fails at -41.6% despite being correct to -0.008%
  at the right phase.
