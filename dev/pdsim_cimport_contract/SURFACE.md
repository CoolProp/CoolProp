# PDSim ↔ CoolProp Cython contract surface

The exact C-level surface that **ibell/pdsim** consumes from CoolProp via
`cimport` (extracted from `ibell/pdsim @ main`).  Dropping the hand-written
Cython interface in v8 (→ nanobind core + a thin frozen `State` shim)
**must preserve everything below**, or PDSim breaks at compile time — or worse,
silently returns wrong numbers.

`pdsim_surface.pyx` exercises this whole surface; `run_contract.py` checks the
values. If a v8 change drifts the surface, the build fails; if it drifts the
unit convention, the value checks fail.

## cimports PDSim makes

| Statement | Files |
|---|---|
| `from CoolProp.State cimport State` (and `as StateClass`) | recip/_recip, flow/flow{,_models}, core/containers |
| `from CoolProp.CoolProp cimport State as State` | core/state_flooded |
| `from CoolProp.constants_header cimport parameters` | core/containers |
| `cimport CoolProp.constants_header as constants` | core/containers, core/state_flooded |
| `from CoolProp cimport constants_header` | flow/flow_models |

## `State` cdef class — members used

**Construction:** `State(Fluid, {'T':T, 'D':rho})` (dict; keys T/D/P/H/S/Q)

**Methods:** `update(dict)`, `update_Trho`, `update_ph`, `copy() -> State`,
`Props(parameters)`, `set_Fluid(Fluid, backend)`,
`get_T get_p get_h get_rho get_s get_u get_cp get_cp0 get_cv get_MM get_dpdT
get_visc get_cond get_speed_sound get_Q`

**Direct cdef attributes (read at C level):**
`pAS` (the `AbstractState`), `T_`, `p_`, `rho_`, `Fluid` (bytes), `phase` (bytes)

## `AbstractState` (reached through `State.pAS`, 24 call sites)

`keyed_output(parameters)`, `rhomass()`, `cpmass()`, `cvmass()`, `T()`, `p()`,
`update(input_pairs, v1, v2)`, `fluid_names()`,
`first_partial_deriv(parameters, parameters, parameters)`

> PDSim does **not** depend only on `State`; it dips into `AbstractState`
> through `.pAS`. The frozen v8 shim must therefore expose a cimportable `pAS`
> object carrying these methods too — not just `State`.

## `constants_header` — names used

Enum **types:** `parameters`, `input_pairs`
**Values:** `iT iP iHmass iSmass iUmass iDmass iDmolar iCpmass iCp0mass iCvmass
ispeed_sound iconductivity iviscosity imolar_mass iQ iP_critical iT_critical`,
input pair `PT_INPUTS`

## ⚠️ Findings that make this more than "keep the method names"

1. **The legacy `State` class is kPa / kJ, NOT SI.** `get_p() = pAS.p()/1000`
   (kPa), `get_h/get_u = …/1000` (kJ/kg), `get_s/get_cp/get_cv` likewise;
   `get_dpdT = first_partial_deriv(iP,iT,iDmolar)/1000`. But `get_rho`,
   `get_speed_sound`, **`Props()`**, and everything via **`.pAS`** are SI.
   `p_` is kPa; `rho_` is SI. PDSim is written against this split. A shim that
   forwards SI without the `/1000` factors makes PDSim wrong by 1000× **with no
   crash** — the most dangerous possible regression. The value checks lock this.

2. **No `CoolProp/__init__.pxd` ships**, so an external `cimport` of the
   package only resolves when cython's `include_path` points at the cimport
   root (site-packages). v8 should ship `__init__.pxd` (even empty) to make the
   surface cleanly cimportable.

3. **The public surface *leaks* fmt — it is NOT part of PDSim's interface.**
   Nothing in the `State`/`AbstractState` API signatures uses an fmt type;
   it is dragged in purely by transitive `#include`:
     * **fmt** ← `AbstractState.h → CachedElement.h → CoolPropTools.h →
       CPstrings.h`, whose inline `format()` helpers call `fmt::sprintf`.
     * **C++17 `std::filesystem`** ← `CoolPropTools.h → CPfilepaths.h`
       (`write_bytes_atomic`).

   **RapidJSON leak resolved (v8 migration).** The former rapidjson leak via
   `Configuration.h → rapidjson_include.h` was eliminated by the v8
   RapidJSON→nlohmann/json migration: JSON types no longer appear in installed
   headers, and link-time symbol hiding (`coolprop_hide_json_symbols`) keeps
   nlohmann/valijson out of the dynamic export table.

   The remaining fix is header hygiene for fmt: move the `format()` bodies out
   of `CPstrings.h` into a `.cpp` and forward-declare/pimpl the
   `std::filesystem::path` parameter. After that a downstream consumer (and the
   frozen v8 `State` shim) compiles with only `CoolProp.get_include_directory()`
   and **zero** third-party `-I` flags. The fmt `-I` path in `setup_contract.py`
   is a STOPGAP for today's remaining leakage; once the surface is cleaned,
   delete it and this test becomes the regression guard that the leak has not
   returned.

4. **The contract is link-free.** The built module links with no `-lCoolProp`
   (`-undefined dynamic_lookup`); the cdef-class methods resolve at runtime
   through Cython's vtable capsule. This is why the capsule-forwarding design
   works without shipping a linkable library — PDSim never links CoolProp.

## Running

```bash
python dev/pdsim_cimport_contract/run_contract.py          # standalone
RUN_PDSIM_CIMPORT_CONTRACT=1 python -m pytest dev/pdsim_cimport_contract/test_pdsim_contract.py
```

The pytest wrapper is opt-in (skips unless `RUN_PDSIM_CIMPORT_CONTRACT` is set)
so it isn't pulled into default collection — it compiles a shim against the
installed CoolProp. The build auto-locates fmt from any in-repo
`build*/_deps/…`; override with `COOLPROP_FMT_INCLUDE` if needed.
