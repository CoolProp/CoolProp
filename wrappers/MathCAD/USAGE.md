CoolProp Mathcad Wrapper — Function Reference
=============================================

Quick overview of every function the add-in registers with Mathcad Prime. For argument-by-argument detail, worked examples, and the full list of valid parameter/input-pair/phase names, see the canonical reference: [Web/coolprop/wrappers/Mathcad/MathcadWrappers.rst](../../Web/coolprop/wrappers/Mathcad/MathcadWrappers.rst) (rendered at https://coolprop.github.io under Wrappers → Mathcad). The `CoolPropFluidProperties.mcdx` example worksheet (see the main [README](README.md)) has these in live use.


High-Level Functions
====================

One-shot calls — each evaluates a fluid or mixture property directly from its inputs, with no persistent state object to manage. Start here for everyday property lookups; reach for the Low-Level (AbstractState) API below for speed when the same backend fluid/mixture can be used to evaluate properties over many state points.

Fluid & humid-air properties
----------------------------

* ``PropsSI(Output Name, Input Name 1, Input Property 1, Input Name 2, Input Property 2, Fluid Name)`` — returns a fluid-specific parameter that depends on the fluid state (the core property-lookup function).
* ``Props1SI(Fluid, Property Name)`` — returns a fluid-specific parameter that does not depend on the fluid state (e.g. critical temperature).
* ``PropsSImulti(Output Name, Input Name 1, Input Property 1 (Array), Input Name 2, Input Property 2 (Array), Fluid Name)`` — returns a range of fluid-specific parameters for the state ranges defined by the input property arrays.
* ``PhaseSI(Input Name 1, Input Property 1, Input Name 2, Input Property 2, Fluid Name)`` — returns the fluid phase, dependent on the fluid state.
* ``HAPropsSI(Output Name, Input Name 1, Input Property 1, Input Name 2, Input Property 2, Input Name 3, Input Property 3)`` — returns a parameter of humid air, dependent on the fluid state.

Information and Reference State Functions
-----------------------------------------

* ``get_global_param_string(Name)`` — returns the value of a requested CoolProp global parameter (e.g. `"version"`, `"errstring"`).
* ``get_fluid_param_string(Fluid, Name)`` — returns the value of a requested CoolProp fluid-specific string parameter.
* ``set_reference_state(Fluid, Reference State String)`` — sets the reference state to `"IIR"`, `"ASHRAE"`, `"NBP"`, or `"DEF"`.


Mixture / binary-pair settings
------------------------------

* ``get_mixture_binary_pair_data(CAS 1, CAS 2, Name of the parameter to retrieve)`` — returns the value of a requested binary interaction parameter.
* ``set_mixture_binary_pair_data(CAS 1, CAS 2, Parameter Name, Parameter value)`` — sets the value of a binary interaction parameter.
* ``apply_simple_mixing_rule(CAS 1, CAS 2, Mixing Rule)`` — sets a simple mixing rule (`"linear"` or `"Lorentz-Berthelot"`) for a binary pair.
* ``get_predefined_mixture_fluids(MixtureName)`` — returns a semicolon-delimited list of fluid names in a predefined mixture.
* ``get_predefined_mixture_fractions(MixtureName)`` — returns a column vector of mole fractions for a predefined mixture.


Low-Level (AbstractState) API Functions
=======================================

In addition to the high-level, stateless functions above, the wrapper exposes CoolProp's [Low-Level (AbstractState) API](https://coolprop.org/coolprop/LowLevelAPI.html). Note that all of the Low-Level functions implemented begin with a prefix ``AS_``, which stands for ``AbstractState``. Mathcad has no way to hold a C++ object as a worksheet variable, so a persistent low-level fluid/mixture state is instead represented by a plain numeric **handle** (a real scalar), created once with ``AS_factory`` and then passed into the other ``AS_*`` functions:

* ``AS_factory(Backend, Fluids)`` — creates a persistent state and returns a handle.
* ``AS_set_fractions(Handle, Fractions)`` — sets mole/mass/volume fractions for a mixture. Returns ``Handle`` unchanged (reassign it, e.g. ``h := AS_set_fractions(h, x)``, so a downstream low-level call using the new ``h`` is guaranteed to run after this one). **Fraction basis:** mole/mass/volume is auto-detected from the backend — not user-selectable — so just pass fractions in whichever basis the backend expects (mole fractions for HEOS/REFPROP/Cubics/PCSAFT/tabular backends; mass fractions for the Incompressible backend). **Input validation:** checks that ``Fractions`` has one entry per fluid in the handle's mixture, that they sum to 1.0 (within 1e-6), and that the handle actually is a mixture — calling this on a pure-fluid handle is a Custom Error, not a silent no-op.
* ``AS_specify_phase(Handle, Phase)`` / ``AS_unspecify_phase(Handle)`` — impose/remove a fixed phase on a handle for all subsequent updates; call ``AS_specify_phase`` before any ``AS_update``/``AS_props``/``AS_props_multi`` calls on that handle. Both return ``Handle`` unchanged. ``Phase`` is one of ``"phase_liquid"``, ``"phase_gas"``, ``"phase_twophase"``, ``"phase_supercritical"``, ``"phase_supercritical_gas"``, ``"phase_supercritical_liquid"``, ``"phase_critical_point"``, ``"phase_unknown"``, or ``"phase_not_imposed"`` (case sensitive — ``CoolProp::phases`` in ``DataStructures.h``).
* ``AS_free(Handle)`` — releases a handle. Optional — see below.
* ``AS_param_index(Name)`` / ``AS_input_pair_index(Name)`` — resolve an output parameter name (e.g. ``"T"``, ``"Dmolar"``) or input pair name (e.g. ``"PT_INPUTS"``) to the integer index the functions below expect. Each only needs to be called **once anywhere in the worksheet** — the result is an ordinary Mathcad variable, so it can be defined at worksheet scope and reused by any number of program blocks or math regions; it is not limited to setting a local variable within a single Mathcad program structure. A full list of valid input pair names is in the [``CoolProp::input_pairs``](https://coolprop.org/_static/doxygen/html/namespace_cool_prop.html#a85cda1634e1e4c1f76425cfd63edf155) enum in the CoolProp source documentation.
* ``AS_update(Handle, InputPairIdx, Value1, Value2)`` — updates the state to a new input point without reading any output; returns ``Handle`` unchanged. Pair with ``AS_get`` to update once and then read as many outputs as needed with separate calls, as an alternative to ``AS_props``.
* ``AS_get(Handle, ParamIdx)`` — returns one output parameter from the state's *current* point (i.e. whatever ``AS_update`` last set it to).
* ``AS_props(Handle, InputPairIdx, Value1, Value2, ParamIdx)`` — fused update + single output in one call.
* ``AS_props_multi(Handle, InputPairIdx, Value1Array, Value2Array, ParamIdxArray)`` — updates the state for a whole array of input points in one call and returns up to 5 requested outputs as a table (one row per input point, one column per requested output). This is the function to use for evaluating many state points against the same fluid/mixture — it avoids reconstructing the backend per point, which matters most for tabular backends (BICUBIC/TTSE).
* ``AS_build_phase_envelope(Handle, Level)`` — traces the phase envelope (dew/bubble curve) for a Handle. Call once before ``AS_get_phase_envelope_data``. Returns ``Handle`` unchanged. ``Level`` controls how much extra refining is done between traced points: ``"none"`` (**CoolProp's own recommendation** — skips refining), ``"fine"`` (default tolerances — any value other than ``"none"``/``"veryfine"`` behaves the same way), or ``"veryfine"`` (tighter tolerances, more points).
* ``AS_get_phase_envelope_data(Handle, Trigger)`` — returns the traced envelope as a table: one row per point, columns ``T``, ``P``, ``rhomolar_vap``, ``rhomolar_liq``. Raises a Custom Error if ``AS_build_phase_envelope`` hasn't been called yet for this Handle. Per-component compositions (mirroring the Python/C++ interface's ``x``/``y``) aren't included — a meaningfully different, mixture-size-dependent shape from this table; not implemented for now. ``Trigger`` is unused and should be passed as 0 — Mathcad functions **must** pass at least one parameter.
* ``AS_pe_tmax(Handle, Trigger)`` / ``AS_pe_pmax(Handle, Trigger)`` — the cricondentherm (highest-temperature point on the envelope) and cricondenbar (highest-pressure point), each as a 2-element column vector ``[T; P]``. Computed by scanning the same table ``AS_get_phase_envelope_data`` returns for its max — CoolProp tracks these points internally while tracing (and, for most mixtures, solves for their *exact* location rather than settling for the closest traced point — see the canonical reference for detail), but doesn't expose them through the public Low-Level C API this wrapper is built on, so this wrapper reproduces the max-scan on its own side rather than extending that shared surface. Both require ``AS_build_phase_envelope`` first, same as ``AS_get_phase_envelope_data``.
* ``AS_list_handles(Trigger)`` / ``AS_list_states(Trigger)`` — introspection pair, mainly for debugging: each returns every currently-live handle from every ``AS_factory`` call still open anywhere in the worksheet (including ones held only in worksheet-level variables you've since lost track of). ``AS_list_handles`` returns a column vector of the handles; ``AS_list_states`` returns their ``"Backend|Fluids"`` keys as one ``";"``-delimited string, **in the same order** as ``AS_list_handles``'s column vector — Mathcad functions can only return one value each, so the pairing across the two calls is how a single handles-plus-labels listing is reconstructed. Both raise a Custom Error if nothing is currently open. A handle released via ``AS_free`` (or otherwise gone dead) is dropped from the listing automatically — no stale entries. ``Trigger`` is otherwise unused — Mathcad Custom Functions require at least one argument, and there's no natural one to key a whole-registry snapshot on, so this exists to satisfy that requirement (pass a literal ``0`` if nothing better is at hand). Wiring it to a handle already on the sheet instead gives these calls a real dependency edge, so they re-run whenever that handle's defining cell does; otherwise, use **Recalculate Worksheet** to refresh them. Place ``AS_list_handles`` and ``AS_list_states`` adjacent to each other on the sheet — a factory/free call evaluated in between the two, in the same recalculation pass, is the one case where their ordering guarantee could momentarily slip.

**Getting call order right.** Mathcad recalculates by dependency/region order, not top-to-bottom sequential code, so a Handle must always be created before it is used. Two supported patterns:

1. **A Mathcad "program" block** (Programming toolbar): create the handle, make however many ``AS_props``/``AS_props_multi`` calls are needed (or ``AS_update`` followed by as many ``AS_get`` calls as needed), and release it with ``AS_free`` at the end, all as sequential statements in one program region. Recommended when the worksheet just needs one derived result.
2. **One ``AS_factory`` call near the top of a worksheet**, referenced by many downstream low-level function calls/plots. Use **Recalculate Worksheet** (full top-to-bottom recalculation in region order, not a partial/incremental recalc) to guarantee the factory call runs before anything that reads the handle. Suits interactive/exploratory worksheets with many uses of the same fluid state. In this pattern, avoid calling ``AS_free`` from an independent mathcad region — nothing guarantees it runs after every reader of the handle. ``AS_factory`` itself is memoized: recalculating it with the same ``Backend``/``Fluids`` returns the SAME handle rather than rebuilding the backend, so repeatedly recalculating the same ``AS_factory`` call neither leaks state nor pays construction cost again (any phase constraint from a prior ``AS_specify_phase`` call is cleared on reuse, so an edited-away call can't leave it silently in effect; mixture fractions are not reset, since ``AS_set_fractions`` is always re-chained after ``AS_factory`` anyway).

See the ``CoolPropFluidProperties.mcdx`` example worksheet for **Program block** pattern use.
