# Plan: Residual Entropy Scaling (RES) Transport Properties in CoolProp

**Status: implemented.** Phases 1-7 are complete. This document is both the original plan
and the record of where the implementation departed from it - see
[Design decisions](#design-decisions) for the choices that were made deliberately and
[Deviations from the original plan](#deviations-from-the-original-plan) for what changed
during implementation and why.

Companion documents:
- `dev/RES_test_assessment.md` - state of the `[RES]` test suite and every known failure.
- `dev/RES_params/source_documentation.txt` - literature source of each parameter table.

---

## Design decisions

| # | Decision | Status |
|---|----------|--------|
| D1 | **HEOS pure + mixtures**: RES is explicit opt-in via `use_viscosity_RES(true)`. Existing pure-component reference correlations are never silently replaced. | as planned |
| D2 | **Cubic backends (PR/SRK)**: auto-enable RES when parameters exist - there is no existing model to protect. | as planned |
| D3 | **Parameter storage**: one `dev/res_transport_parameters.json`, keyed `fluid x eos_type`, compiled in via `res_transport_parameters_JSON.h`. | as planned |
| D4 | **HEOS parameters** are the REFPROP-fitted values used as an approximation. | **revised** - now vetted per fluid, see D8 |
| D5 | **Alpha-function safety**: `n_res`/`xita` are fitted for a specific alpha function; changing it invalidates them and raises a named-component `ValueError` until refitted values are supplied. | as planned |
| D6 | **Critical-enhancement viscosity source**: the enhancement uses the *RES* viscosity, not a reference correlation, keeping the model self-consistent. Li 2024 instead takes eta from REFPROP. | **new decision** |
| D7 | **Mixture critical enhancement**: not applied (pure fluids only), matching the Li 2024 default (`use_crit_enh_mix=False`). | as planned |
| D8 | **Fluid vetting**: fluids whose HEOS `s_res` moves the transport property too far, and fluids CoolProp cannot construct at all, are omitted from the JSON so the backend raises instead of returning a silently wrong number. | **new decision** |
| D9 | **Missing parameters fail loudly**: requesting RES without parameters raises `ValueError` naming the fluid/component, rather than falling back to the default model. | **new decision** |

### D6 in detail - why the critical enhancement keeps the RES viscosity

The enhancement is proportional to 1/eta, so the choice is measurable. Isolating it and
rescaling:

| fluid | current (RES eta) | with REFPROP eta | with HEOS default eta | share of lambda |
|---|---|---|---|---|
| PROPANE (T/Tc=1.02) | -4.885% | +0.001% | -1.004% | 28.3% |
| R143A (T/Tc=1.37) | +1.432% | -0.007% | -1.087% | 15.5% |

REFPROP's eta reproduces the published reference almost exactly, which confirms the rest of
the implementation (dilute term, residual term, Olchowy-Sengers, the guards) is correct. It is
not reachable from inside the backend without depending on REFPROP. The self-consistent RES
choice was kept; the cost is a few percent confined to the near-critical region, and it is the
reason PROPANE fails its 3% reference tolerance.

### D8 in detail - how a fluid is vetted

Two independent filters in `dev/convert_RES_csv_to_json.py`:

1. **Not constructible in CoolProp** (29 of 151 fluids) - REFPROP's fluid set is much larger
   than CoolProp's. Dropped entirely; the entry could never be looked up.
2. **`s_res` does not transfer** - the HEOS entry is dropped, keeping `dilute`, `PR`, `SRK`
   and `critical_enhancement`.

Filter 2 needs *both* the transport deviation at the published sample point **and** the
`s_res` deviation across a range of states. Neither alone is sufficient:

- There is only **one published sample point per fluid**, and some are almost pure dilute gas,
  where the transport property barely depends on `s_res`. Measured residual share: R161 4.3%
  (viscosity) / -0.2% (conductivity), VINYLCHLORIDE 0.3%, R13 0.9%, R1234YF 0.2%. A sub-1%
  transport deviation at such a point proves nothing - those fluids have systematic `s_res`
  errors elsewhere (R161 12%, VINYLCHLORIDE 32%, R13 17%, R1234YF 6%) and stay excluded.
- HEOS and REFPROP can disagree about the **phase** at a point straddling their near-identical
  saturation curves, which inflates the deviation without implicating the parameters. Imposing
  the reference phase collapses BENZENE from -41.6% to -0.008% (kept) and R41 from +266.8% to
  -1.12% (still over threshold, excluded).

---

## Literature sources

- Viscosity: Martinek 2025 (doi:10.1021/acs.jced.4c00451) - dilute-gas and REFPROP-fitted residual
- Viscosity/conductivity PR + SRK fits: Yang 2025 (doi:10.1021/acsomega.4c10815)
- Conductivity: Li 2024 (doi:10.1021/acs.iecr.4c02946) - dilute-gas, REFPROP-fitted residual, critical enhancement
- Critical enhancement form: Olchowy & Sengers 1989

These are recorded per parameter group in the JSON under a `sources` key.

> **Trap - dilute-gas coefficient order.** The Li 2024 SI ships `Dilute_gas_TC.txt` with header
> `TC0 = n0*T^4 + ... + n4` and columns named `n0;n1;n2;n3;n4`, while the kktoolbox copy in
> `dev/RES_params/` has `TC0 = n4*T^4 + ... + n0` with columns named `n4;n3;n2;n1;n0`. **The
> numbers are identical and the first numeric column is the T^4 coefficient in both**; only the
> naming is inverted. Since the converter indexes by column *name*, feeding it the other file
> would silently reverse the polynomial. The JSON always stores ascending powers,
> `n[0] + n[1]*T + ... + n[4]*T^4`.

---

## Phase 1 - Data structures

`include/CoolProp/CoolPropFluid.h`: `ViscosityRESData` and `ConductivityRESData`, added to
`TransportPropertyData` as `viscosity_res` / `conductivity_res`. Both carry `n_dilute`,
`n_res`, `xita`, `group_num`, `molar_mass`, `n_params_match_alpha` and `provided`; the
conductivity struct additionally carries the Olchowy-Sengers terms (`R_D`, `gamma_uni`,
`Gamma`, `phi0`, `t_ref`, `q_D`) and `crit_provided`.

`HelmholtzEOSMixtureBackend.h`: runtime flags `viscosity_RES_enabled`,
`conductivity_RES_enabled`.

## Phase 2 - Parameter JSON

`dev/res_transport_parameters.json`, generated by `dev/convert_RES_csv_to_json.py` from the
checked-in tables in `dev/RES_params/`. Shape:

```json
{
  "viscosity":    { "CO2": { "dilute": {"n": [...]}, "HEOS": {}, "PR": {}, "SRK": {} } },
  "conductivity": { "CO2": { "dilute": {}, "HEOS": {}, "PR": {}, "SRK": {},
                             "critical_enhancement": {} } },
  "sources":      { "viscosity": {}, "conductivity": {} }
}
```

An `HEOS` key may be absent (see D8); a `critical_enhancement` block is omitted entirely when
the source row is all zeros, which means "not fitted" rather than "zero enhancement".

## Phase 3 - JSON loading

`FluidLibrary.h::load_RES_transport_parameters()` fills the structs, called from `add_one()`
with `eos_type = "HEOS"` and from `AbstractCubicBackend::setup()` with `"PR"`/`"SRK"`.

## Phase 4 - Transport routines

`TransportRoutines::viscosity_RES()` and `conductivity_RES()`. Both guard on `provided` and
`n_params_match_alpha`, compute `s_plus = -smolar_residual()/R`, evaluate the dilute
polynomial per component with Wilke mixing, then the residual term
(`pow = {1.8, 2.4, 2.8}` viscosity, `{1.0, 1.5, 2.0, 2.5}` conductivity).

Conductivity adds the Li 2024 critical enhancement for pure fluids only, suppressed exactly as
the reference implementation does - `rho/rhoc >= 2`, `T/Tc > 1.4`, non-positive `arg`, or
unusable/absent parameters.

## Phase 5 - Backend dispatch

HEOS `calc_viscosity()` / `calc_conductivity()` route to RES when enabled, and raise a
`ValueError` naming the fluid or the offending components when enabled without parameters.
`AbstractCubicBackend` auto-enables when parameters exist and otherwise raises
`NotImplementedError`.

## Phase 6 - AbstractState API

Eight virtuals: `use_viscosity_RES`, `use_conductivity_RES`,
`set_{viscosity,conductivity}_RES_parameters`,
`set_{viscosity,conductivity}_RES_residual_params`,
`get_{viscosity,conductivity}_RES_residual_params`. All eight are exposed to Python.

Toggles and setters clear the memoized `_viscosity`/`_conductivity`, so a change after a
property read takes effect instead of silently returning the previous model's value.

## Phase 6b - Alpha-function invalidation

`set_cubic_alpha_C()` clears `n_params_match_alpha`; the transport routines then raise a
named-component `ValueError` until `set_*_RES_residual_params()` supplies refitted values.

## Phase 7 - Catch2 tests

`src/Tests/CoolProp-Tests-RES.cpp`, tag `[RES]`, registered in `CMakeLists.txt` with
`RES_SAMPLES_DIR` pointing at `dev/RES_samples/`. 19 test cases. See
`dev/RES_test_assessment.md` for status and the known failures.

---

## Deviations from the original plan

| Area | Planned | Actual | Why |
|------|---------|--------|-----|
| HEOS parameters | ship REFPROP-fitted values for every fluid | vetted per fluid; unusable HEOS entries omitted | D8 - silently wrong numbers are worse than a clear exception |
| Missing parameters | unspecified | raises `ValueError` naming fluid/components | D9 - the opt-in previously fell through to the default model with no signal |
| Critical-enhancement viscosity | unspecified | RES viscosity, deliberately not REFPROP's | D6, with the cost measured |
| Critical-enhancement guards | not mentioned | `rho/rhoc >= 2`, `T/Tc > 1.4` ported from the Li 2024 SI | without them the enhancement is applied where the reference returns exactly 0 |
| `crit_provided` | set whenever a block exists | requires `t_ref`, `Gamma`, `phi0`, `q_D` > 0 | all-zero rows mean "not fitted"; treating them as valid divided by zero and returned NaN (D2O, HELIUM, ORTHOHYD) |
| Cache invalidation | not considered | toggles and all four setters clear the memoized values | otherwise a toggle after a read is a silent no-op |
| Test file | `tests/catch/test_RES_transport.cpp` | `src/Tests/CoolProp-Tests-RES.cpp` | matches the existing test layout |
| Python wrapper | deferred | **done**, both interfaces | needed to test the real implementation from Python at all |
| Converter input | read from an external kktoolbox checkout | reads `dev/RES_params/`, now checked in | reproducible without a second repo |
| JSON provenance | not planned | `sources` block with citation + DOI per group | the tables come from three different papers |
| Fluid coverage | all 151 fluids | 122 - 29 dropped as not constructible in CoolProp | those entries could never be looked up |

### Known divergence from the reference implementation

**Mixture dilute-gas conductivity.** Li 2024's SI computes mixture lambda0 from REFPROP
directly (`PropsSI('L','T',T,'P',0.1,'REFPROP::...')`) and only falls back to polynomial +
Wilke when that fails. The C++ always uses polynomial + Wilke, since depending on REFPROP is
not acceptable. Measured gap: ARGON+NEON -11.6%, METHANE+NITROGEN -5.1%, CO2+NITROGEN -2.3%,
BUTANE+METHANE -1.8%. This makes the 5% mixture-conductivity test tolerance unreachable for
near-ideal-gas samples. Pure fluids are unaffected - the SI uses the polynomial there, and the
published pure-fluid reference is reproduced to a median 0.0015% (conductivity) and 0.0001%
(viscosity).

---

## Files

| File | Change |
|------|--------|
| `include/CoolProp/CoolPropFluid.h` | RES structs + fields in `TransportPropertyData` |
| `include/CoolProp/AbstractState.h` | eight virtual API methods |
| `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.{h,cpp}` | flags, dispatch, API impl, cache invalidation |
| `src/Backends/Helmholtz/TransportRoutines.{h,cpp}` | `viscosity_RES`, `conductivity_RES` |
| `src/Backends/Helmholtz/Fluids/FluidLibrary.{h,cpp}` | `load_RES_transport_parameters()` |
| `src/Backends/Cubics/CubicBackend.{h,cpp}` | overrides, setup loading, alpha invalidation |
| `src/nanobind_interface.cxx` | RES API for the default Python interface |
| `wrappers/Python/CoolProp/{cAbstractState.pxd,AbstractState.pxd,AbstractState.pyx,CoolProp.pyx}` | RES API for the legacy Cython interface |
| `src/Tests/CoolProp-Tests-RES.cpp` | Catch2 tests, tag `[RES]` |
| `dev/res_transport_parameters.json` + `include/res_transport_parameters_JSON.h` | parameters |
| `dev/convert_RES_csv_to_json.py` | converter, vetting, `--keep-all-heos` measurement mode |
| `dev/compare_HEOS_vs_REFPROP_RES.py` | per-fluid deviation measurement through the real C++ API |
| `dev/RES_params/`, `dev/RES_samples/` | source tables and published sample data |

---

## Deferred / further work

1. **RES on the REFPROP backend** - currently the API lives on `HelmholtzEOSMixtureBackend`, so
   REFPROP cannot evaluate RES. This is the binding constraint on vetting: deviations can only
   be measured at the one published sample point per fluid. *(new; now the highest-value item)*
2. **Test suite overhaul** - see `dev/RES_test_assessment.md`; the fail-open
   `catch(...)`/`REQUIRE(ok >= N)` gate is the blocking item.
3. **Mixture critical enhancement** - Li 2024 supports it behind `use_crit_enh_mix`; not implemented.
4. **Saturated-state transport** (`viscosity_RES_sat`, `conductivity_RES_sat`).
5. **PR_Twu preset in the JSON** - the runtime setter covers the primary use case.
6. **`xitapower` exponents** are hardcoded, not parameterised.
7. **Existing `rhosr-CS` model** - untouched and independent; coexists with RES.
