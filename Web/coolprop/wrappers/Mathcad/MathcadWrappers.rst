.. _mathcadwrappers:

*******************************************
CoolProp Function Implementation in Mathcad
*******************************************

.. contents:: :depth: 3

.. role::green
   :class: green-text

For the most part, the Mathcad wrappers follow the Python implementation and most of the examples on this site can be executed in Mathcad Prime with very little modification.  There are a few key difference:

1. Mathcad Strings always use "double-quotes"
2. Units (*see below*)
3. Functions can be evaluated on their own, assigned to a variable, or both at the same time,

    * PropsSI("D", "T", 295.15, "P", 101325.0, "Water") = 997.773
    * ρ := PropsSI("D", "T", 295.15, "P", 101325.0, "Water")
    * ρ := PropsSI("D", "T", 295.15, "P", 101325.0, "Water") = 997.773

4. Mathcad *can* execute equations in random order as changes are made to the worksheet.  This can sometimes cause unexpected behavior if CoolProp settings are modified non-sequentially.  It is a good idea to press **<Ctrl>-<F9>** periodically to recalculate the entire worksheet from top to bottom.
5. Unfortunately, there is no way to emulate the live Python examples found elsewhere on this web site with Mathcad Prime, so the Mathcad syntax and functionality will be emulated in the equations below.

A majority of these functions and examples of their use are described in the Mathcad file ``CoolPropFluidProperties.mcdx``, found in the  :sfdownloads:`MathcadPrime` folder on SourceForge.

High-Level Functions
====================

PropsSI - State Dependent Fluid Properties
----------------------------------------------

`PropsSI` is the basic, high-level function for returning the scalar value of a specified output property at a fixed state point.::

    PropsSI("Output", "Input1", Val1, "Input2", Val2, "Fluid")

Where,

* "Output" = Requested output property string; see :ref:`Table of Valid Parameters <parameter_table>`.
* "Input1" = First state-point property string.
* Val1 = First state-point property value (scalar variable)
* "Input2" = Second state-point property string.
* Val2 = Second state-point property value (scalar variable)
* "Fluid" = Fluid string (e.g, "Water", "Ammonia", "Air.mix", etc.).

.. note::
    The Fluid string can use a backend prefix (e.g., "HEOS::", "INCOMP::", "REFPROP::", etc.) to specify an alternative EOS; the default being the Helmholtz EOS ("HEOS::") if not provided.

.. note::
    In addition to pure and pseudo-pure fluid strings, the fluid string can be specified as a predefined mixture (e.g., "Air.mix") or an ad-hoc mixture, specifying each pure component and mole fraction (e.g., "O2[0.2096]&N2[0.7812]&AR[0.0092]") where the mole fractions are in square braces [ ] and the components are delimited with "&".

**EXAMPLE:**
    .. math::
       h := PropsSI("H",\ "T",\ 300.0,\ "P",\ 500,\ "Helium") = 1562994.2

|

----


PropsSImulti - Multiple State Dependent Fluid Properties
------------------------------------------------------------

`PropsSImulti` will return a vector/matrix of multiple fluid output properties spanning a range of state points. The return value is an (:math:`m x n`) matrix, where :math:`n` is the number of columns, one for each requested property, and :math:`m` is the number of rows for each state point. For the most part, the parameters of `PropsSImulti` are the same as `PropsSI` with the following exceptions.::

    PropsSImulti("Outputs", "Input1", Vec1, "Input2", Vec2, "Fluid")

Where,

* "Outputs" = Requested output properties string containing a delimited list of one or more valid output property names from the :ref:`Table of Valid Parameters <parameter_table>`.  For this Mathcad wrapper, the delimiter can be any one of (*comma, <space>, colon, semi-colon, or ampersand*), but must be consistent.
* Vec1, Vec2 = State point array pairs corresponding to "Input1" and "Input2".  These are single-column, :math:`m`-element vector arrays and must be the same length or an error will be thrown.

|

**EXAMPLE:**

    Define a fluid:       :math:`fl` := "Water"

    Triple Points:         :math:`T_t := 273.15`     &     :math:`P_t := 611.655`


    Critical Points:       :math:`T_c := 647.096`     &     :math:`P_c := 2.206\cdot10^7`

    Liquid points:       :math:`T_L := mean(T_t,T_c)`     &     :math:`P_L := mean(P_t,P_c)`

    Set Vectors:       :math:`Tvec := \begin{bmatrix} T_t\\ T_L \\ T_L \end{bmatrix}`       :math:`Pvec := \begin{bmatrix} P_L\\ P_L \\ P_c \end{bmatrix}`

    Calc:   :math:`M := PropsSImulti("D\ H",\ "T",\ Tvec,\ "P",\ Pvec,\ fl) = \begin{bmatrix} 1005.334 & 1.115\cdot10^4\\ 886.137 & 7.987\cdot10^5 \\ 893.241 & 8.044\cdot10^5 \end{bmatrix}`

    Extract individual variables from columns:         :math:`\rho = M^{<0>}`     &     :math:`h = M^{<1>}`

|

----

Props1SI - State Independent Fluid Properties
-------------------------------------------------

`Props1SI` returns non-state-dependent properties of a fluid/mixture and does not require state point names or values.  This function only requires the output property and fluid name strings.::

    Props1SI("Output", "Fluid")

Where,

* "Output" = One of the "*trivial*" output-only properties from the :ref:`Table of Valid Parameters <parameter_table>`.
* "Fluid" = Fluid string, as defined above for `PropsSI`


**EXAMPLE:**

    Define a fluid:         :math:`fl` := "Water"

    Triple Point Temperature:         :math:`T_t := Props1SI("Ttriple",\ fl) = 273.15`

    Triple Point Pressure:               :math:`P_t := Props1SI("ptriple",\ fl) = 611.655`


    Critical Temperature:                :math:`T_c := Props1SI("Tcrit",\ fl) = 647.096`

    Critical Pressure:                       :math:`P_c := Props1SI("pcrit",\ fl) = 2.206\cdot10^7`

|

----

PhaseSI - Phase Determination
-----------------------------
The function   is used to find the fluid phase at a specified state point.  The calling structure of the function is the same as `PropsSI` except that no "Output" properties are specified, given that the only output returned is a "string" representing the phase of the fluid.::

    PhaseSI("Input1", Val1, "Input2", Val2, "Fluid")

The input parameters are the same as for `PropsSI`, except there is no "Output" parameter as it is assumed to be "Phase".

.. note::
   The `PropsSI` function can be used directly with the output parameter "Phase", but this returns an enumerated integer value for the phase.  PhaseSI returns a string that represents the phase name for that enumerated value.

**EXAMPLE:**

    Define a fluid:                           :math:`fl` := "Water"

    Triple Point Temperature:         :math:`T_t := Props1SI("Ttriple",\ fl) = 273.15`

    Triple Point Pressure:               :math:`P_t := Props1SI("ptriple",\ fl) = 611.655`


    Critical Temperature:                :math:`T_c := Props1SI("Tcrit",\ fl) = 647.096`

    Critical Pressure:                       :math:`P_c := Props1SI("pcrit",\ fl) = 2.206\cdot10^7`

    Liquid points:             :math:`T_L := mean(T_t,T_c)\ -\ 10.0`     &     :math:`P_L := mean(P_t,P_c)`

    Calc:                           :math:`PhaseSI("T",\ T_L,\ "P",\ P_L,\ fl)` = "Liquid"

|

----

HAPropsSI - Humid Air Fluid Properties
----------------------------------------------

`HA PropsSI`  is used to find fluid properties of humid air.  The physics behind the   function is based on the analysis in ASHRAE RP-1845, which is available online: https://www.tandfonline.com/doi/abs/10.1080/10789669.2009.10390874.  It employs real gas properties for both air and water, as well as the most accurate interaction parameters and enhancement factors.   RP-1845 is based largely on the IAPWS-95 formulation for the properties of water.  The calling structure of the function is as follows: ::

    HAPropsSI("Output", "Input1", Val1, "Input2", Val2, "Input3", Val3)

Where,

* "Output" = Requested output property string; see :ref:`Table of Valid HA Parameters <HAparameter_table>`.
* "Input1" = First state-point property string.
* Val1 = First state-point property value (scalar variable)
* "Input2" = Second state-point property string.
* Val2 = Second state-point property value (scalar variable)
* "Input3" = Third state-point property string.
* Val3 = Third state-point property value (scalar variable)

At least one of the inputs must be "T" (dry bulb temperature), "R" (Relative Humidity between 0.0 and 1.0), "W" (Humidity Ratio), or "Tdp" (dew point).

**EXAMPLE:**

    Enthalpy @ 50% Rel. Humidity:

    :math:`h := HAPropsSI("H",\ "T",\ 298.15,\ "P",\ 101325,\ "R",\ 0.5) = 5.042\cdot10^4`

----

Pseudo-Low-Level Functions
==========================

CoolProp's Low-level functions require the creation of an Abstract State object and then evaluation of properties using that object's member functions.  Mathcad does not have the ability to store objects as variables directly, so these pseudo-Low-Level functions either do not require an abstract state object, or create one temporarily for the purposes of extracting and setting CoolProp data and parameters.  These wrapper functions are listed here.  (A true, persistent Low-Level interface *is* available -- see :ref:`Low-Level (AbstractState) Functions <mathcad_lowlevel_functions>` below, which represents the Abstract State object as a plain numeric handle instead.)

get_global_param_string
-----------------------

This function retrieves global CoolProp parameters that are set and maintained by the CoolProp library.::

    get_global_param_string("GlobalParameter")

Where "GlobalParameter" can be one of the following:

* "version"
* "gitrevision"
* "errstring"
* "warnstring"
* "FluidsList", "fluids_list", "fluidslist" ²
* "incompressible_list_pure" ²
* "incompressible_list_solution" ²
* "mixture_binary_pairs_list" ²
* "parameter_list" - Comma delimited list of valid fluid property strings used by ``PropsSI``
* "predefined_mixtures" ²
* "HOME" - User's $HOME or %HOME% directory
* "REFPROP_version" (if installed)
* "cubic_fluids_schema"
* "cubic_fluids_list"
* "pcsaft_fluids_schema"

.. note::
   The ``"errstring"`` option is *extremely* useful when using CoolProp functions in Mathcad.  While the wrapper functions attempt to trap common errors and display them as meaningful Mathcad error messages, highlighting the offending parameter(s), unknown errors will display as "CoolProp Issue: Use get_global_param_string("errstring") for more info". This is the only way to see the actaul CoolProp error message being thrown, even if the error is already trapped by the Mathcad wrapper.

----

get_fluid_param_string
----------------------

This function retrieves fluid property information for a specific fluid/mixture and its behavior and implementation depend on the backend being used.::

    get_fluid_param_string("Fluid", "FluidParameter")

Where,

* "Fluid" is a fluid name string that follows the rules of of the fluid definition for ``PropsSI``, but is typically called for Pure Fluids.
* "FluidParameter" for the HEOS (default) backend can be any of the following strings:
    * "name" - Primary fluid name
    * "aliases" - Fluid alias names that can be used to reference the fluid
    * "CAS" - Unique, numerical identifier assigned by the Chemical Abstract Service
    * "formula" - Chemical formula of the specified fluid
    * "ASHRAE34" - ASHRAE classification of refrigerant toxicity and flammability
    * "REFPROPname" - Equivalent fluid name in NIST REFPROP
    * "BibTeX-<ref>" -
    * "pure" - returns "true" for pure fluids, "false" for mixtures
    * "INCHI" - International Chemical Identifier representation for chemical structures
    * "INCHI_Key" - 27-character hashed string used for web searching
    * "SMILES" - Simplified Molecular Input Line Entry System; compact, ASCII-based notation for representing 2D/3D chemical structures
    * "CHEMSPIDER_ID" - unique ChemSpider database identifier
    * "JSON" - Returns the full JSON definition of the fluid (*not very useful in Mathcad*)

----

set_reference_state
-------------------

Enthalpy and entropy are relative properties!  Always compare differences rather than absolute values of the enthalpy or entropy to other sources.  That said, if can be useful to set the reference state values for enthalpy and entropy to one of a few standard values. This is done by the use of the low-level ``set_reference_state`` function.::

    set_reference_state("refState")

A number of pre-defined reference states ("refSate") can be used:

* IIR:              h = 200 kJ/kg, s = 1 kJ/kg/K at 0°C
* ASHRAE:     h = 0, s = 0 @ -40°C saturated liquid
* NBP:            h=0, s=0, for saturated liquid at 1 atmosphere
* DEF:              Default reference state from the fluid file

.. warning::
   The changing of the reference state should only be done

   1. at the very beginning of your worksheet, or
   2. at the very beginning of a Mathcad program block, resetting it to "DEF" at the end of the program block

   or unexpected results may occur. It is not recommended to change the reference state during the course of making calculations as done here for demonstration purposes only. Further more, because of Mathcad's top-down calculation order, switching back and forth between reference states can lead to very unexpected results (real or apparent) and is not recommended.

----

get_predefined_mixture_fluids
-----------------------------

This is not a wrapper of a CoolProp functions, but an additional helper function use to assist using predefined mixtures. This function  returns a semicolon delimited string of the predefined mixture component names.::

    get_predefined_mixture_fluids("mixture")

Where "mixture" is a predefined mixture name ending in .mix or .MIX. Mixture names are case sensitive, using the available predefined mixture names retrieved from ``get_global_parameter_string("predefined_mixtures")``.

**EXAMPLES**

    For air:         :math:`get\_predefined\_mixture\_fluids("Air.mix")` = "NITROGEN;ARGON;OXYGEN"

    For R401A:   :math:`get\_predefined\_mixture\_fluids("R401A.mix")` = "R22;R152A;R124"

.. note::
   A few predefined mixtures are missing binary interaction parameters for at least one component pair.  These mixtures are defined, but cannot be used, for now, for property calculations.

----

get_predefined_mixture_fractions
--------------------------------

This is not a wrapper of a CoolProp functions, but an additional helper function use to assist using predefined mixtures. This function returns a Mathcad array (column vector) of the predefined mixture component *mole fractions*.::

    get_predefined_mixture_fractions("mixture")

Where "mixture" is a predefined mixture name ending in .mix or .MIX. Mixture names are case sensitive, using the available predefined mixture names retrieved from ``get_global_parameter_string("predefined_mixtures")``.

**EXAMPLES**

    For air:         :math:`mf_{Air}` := :math:`get\_predefined\_mixture\_fractions("Air.mix")` = :math:`\begin{bmatrix} 0.7812\\ 0.0092 \\ 0.2096 \end{bmatrix}`

    For R401A:   :math:`mf_{R401A}` := :math:`get\_predefined\_mixture\_fractions("R401A.mix")` = :math:`\begin{bmatrix} 0.578854 \\ 0.0.185871 \\ 0.0.235274 \end{bmatrix}`

.. note::
   A few predefined mixtures are missing binary interaction parameters for at least one component pair.  These mixtures are defined, but cannot be used, for now, for property calculations.

|

----

get_mixture_binary_pair_data
----------------------------

Get binary pair interaction parameters and other info for a pair of components::

    get_mixture_binary_pair_data("CAS1", "CAS2", "mix_param")

Where,

* "CAS1", "CAS2" are the CAS identifiers for the two pure fluid components.  These must be CAS numbers of the format "7782-44-7" and cannot be the fluid name (in this case "Oxygen".
* "mix_param" can be any of the following parameters:
    * "name1" - first component name (corresponding to CAS1)
    * "name2" - second component name (corresponding to CAS2)
    * "BibTeX" - Reference for interaction parameters
    * "function" - function for calculating parameters (if not constants)
    * "type" - parameter model used
    * "F" - :math:`F` parameter (if used)
    * "xi" - :math:`\xi` parameter (if used)
    * "betaT", "betaV" - :math:`\beta_{T,ij}` and :math:`\beta_{v,ij}`
    * "gammaT", "gammaV" - :math:`\gamma_{T,ij}` and :math:`\gamma_{v,ij}`
    * "zeta" - :math:`\zeta` parameter

.. note::
   Error message string from CoolProp may indicate that the input CAS numbers need to be reversed to retrieve values.

|

----

apply_simple_mixing_rule
------------------------

The function ``apply_simpl_mixing_rule()`` will take either CAS strings or fluid alias strings as inputs and can be set to either "linear" or the "Lorentz-Berthelot" mixing rules.::

    apply_simple_mixing_rule("Fluid1", "Fluid2", "rule")

Use of this function follows the python example exactly on the Fluid Properties | Mixtures page and will not be repeated here.

|

----

Set_mixture_binary_pair_data
----------------------------

Changes the default parameters with the following calls::

    set_mixture_binary_pair_data("CAS1", "CAS2", "param", value)

Where,

* "param" is any of the mixing parameters as defined above under ``get_mixture_binary_pair_data``
* `value` is the value of the parameter being set.

Use of this function follows the python example exactly on the Fluid Properties | Mixtures page and will not be repeated here.

----

.. _mathcad_lowlevel_functions:

Low-Level (AbstractState) Functions
====================================

CoolProp's `Low-Level (AbstractState) API <https://coolprop.github.io/devdocs/coolprop/LowLevelAPI.html>`_ lets a caller build one persistent fluid/mixture state and reuse it for many flashes/outputs, avoiding the cost of reconstructing the backend for every call -- this matters most for tabular backends (BICUBIC/TTSE), where construction alone can cost 80-140 ms.  Since Mathcad cannot hold a C++ object as a worksheet variable, the state is represented here by a plain numeric **handle** (a real scalar): ``AS_factory`` creates the state and returns the handle; the other ``AS_*`` functions take that handle as their first argument.

.. note::
    All Low-Level functions in the Mathcad wrapper are implemented with the two-letter prefix `AS_` for `AbstractState`.

**Getting call order right.**  Mathcad recalculates by dependency/region order, not top-to-bottom sequential code, so a handle must always be created before it is used.  Two patterns are supported:

1. A Mathcad **program** block (Programming toolbar): create the handle, make however many ``AS_props``/``AS_props_multi`` calls are needed (or ``AS_update`` followed by as many ``AS_get`` calls as needed), and release it with ``AS_free`` at the end, all as sequential statements in one program region.  Recommended when the worksheet just needs one derived result.
2. One ``AS_factory`` call near the top of a worksheet, referenced by many downstream calls/plots.  Use **Recalculate Worksheet** (a full top-to-bottom recalculation in region order, not a partial/incremental recalc) to guarantee the factory call runs before anything that reads the handle.  In this pattern, avoid calling ``AS_free`` from an independent call -- nothing guarantees it runs after every reader of the handle.  ``AS_factory`` itself is memoized: recalculating it with the same ``Backend``/``Fluids`` returns the SAME handle rather than rebuilding the backend, so repeatedly recalculating the same call neither leaks state nor pays construction cost again (any phase constraint from a prior ``AS_specify_phase`` call is cleared on reuse, so an edited-away call can't leave it silently in effect; mixture fractions are not reset, since ``AS_set_fractions`` is always re-chained after ``AS_factory`` anyway).

See the ``CoolPropFluidProperties.mcdx`` example worksheet for both patterns in use.

|

----

AS_factory
----------

Creates a persistent Low-Level fluid/mixture state and returns a handle.::

    AS_factory("Backend", "Fluids")

Where,

* "Backend" is the backend to use, e.g. "HEOS", "REFPROP", "BICUBIC&HEOS".
* "Fluids" is a ``&``-delimited list of fluids, e.g. "Water" or "Methane&Ethane".

.. note::
    Calling this again with the same "Backend"/"Fluids" returns the SAME handle rather than rebuilding the backend -- so recalculating this call repeatedly (every worksheet recalculation re-executes it) neither leaks state nor pays construction cost again.  Any phase constraint set by a prior ``AS_specify_phase`` call is cleared on reuse, so removing/changing that call in the worksheet can't leave a stale constraint in effect.

**EXAMPLE:**

    :math:`h := AS\_factory("HEOS",\ "Water")`

|

----

AS_set_fractions
-----------------

Sets the mole/mass/volume fractions for a mixture handle created by ``AS_factory``.::

    AS_set_fractions(Handle, Fractions)

Where,

* `Handle` is a handle returned by ``AS_factory``.
* `Fractions` is a column vector of mole/mass/volume fractions, one per fluid in the mixture.

.. note::
    **Why it echoes Handle back:** Returns ``Handle`` unchanged.  Reassign it, e.g. ``h := AS_set_fractions(h, x)``, so a downstream call that uses this call's return value as its own ``Handle`` argument is guaranteed to run after this one.

.. note::
    **Fraction basis:** The fraction basis (mole, mass, or volume) is auto-detected from the backend and is **not user-selectable** -- it is a fixed property of the backend, not a runtime setting.  HEOS, REFPROP, Cubics, PCSAFT, and the tabular backends all use mole fractions; the Incompressible backend uses mass fractions.  Just pass fractions in whichever basis the backend you chose in ``AS_factory`` expects.

.. note::
    **Input validation:** This function validates: that ``Fractions`` has exactly one entry per fluid in the handle's mixture; that the entries sum to 1.0 (within 1e-6); and that the handle is actually a mixture in the first place -- calling it on a pure-fluid handle is a Custom Error, not a silent no-op.

|

----

AS_mole_to_mass_fractions / AS_mass_to_mole_fractions
--------------------------------------------------------

Converts an arbitrary composition between mole and mass fractions, using a Low-Level state handle's mixture for component identities and molar masses. Unlike ``AS_set_fractions``, this doesn't read or write the handle's own state at all -- it's a pure unit conversion on the ``MoleFractions``/``MassFractions`` argument, useful as a preprocessing step *before* ``AS_set_fractions`` (e.g. converting a mass-basis composition you have on hand into the mole fractions a HEOS-backed handle actually expects).::

    AS_mole_to_mass_fractions(Handle, MoleFractions)
    AS_mass_to_mole_fractions(Handle, MassFractions)

Where,

* `Handle` is a handle returned by ``AS_factory`` -- only its mixture's component identities and molar masses are used; its own composition/state is untouched.
* `MoleFractions`/`MassFractions` is a column vector of the fractions to convert, one entry per fluid in the mixture, in either basis.

.. note::
    **No new CoolPropLib export:** the C++ API has a direct equivalent of this (``AbstractState::calc_mass_fractions()``, computing ``mass_i = mm_i * mole_i / sum(mm_j * mole_j)`` from whatever mole fractions are already set), but it isn't exposed through the public Low-Level C API this wrapper is built on, and adding it there was deliberately avoided. This function gets the same result a different way: ``AbstractState_fluid_names()`` (already used by ``AS_set_fractions`` above) gives the component names, and ``Props1SI("molar_mass", name)`` -- a plain, handle-independent lookup already used elsewhere in this wrapper -- resolves each one's molar mass. Both are already-public surface; nothing new was added to CoolPropLib.h for this.

.. note::
    **Self-normalizing:** the conversion divides by the actual weighted sum of the input (``sum(mm_j * mole_j)`` or ``sum(mass_j / mm_j)``), not by an assumed 1.0 -- so a composition that doesn't already sum to exactly 1.0 still converts to a correctly-normalized result in the other basis, unlike ``AS_set_fractions``, which requires its input to already sum to 1.0.

|

----

AS_specify_phase
-----------------

Imposes a fixed phase on a Low-Level state handle for all subsequent updates (``AS_update``, ``AS_props``, ``AS_props_multi``).  Call this before any of those, once per handle.  Returns ``Handle`` unchanged, so a downstream Low-Level call that uses this call's return value as its own ``Handle`` argument depends on it.::

    AS_specify_phase(Handle, Phase)

Where,

* `Handle` is a handle returned by ``AS_factory``.
* `Phase` is a phase name (case sensitive): "phase_liquid", "phase_gas", "phase_twophase", "phase_supercritical", "phase_supercritical_gas", "phase_supercritical_liquid", "phase_critical_point", "phase_unknown", or "phase_not_imposed" (``CoolProp::phases`` in ``DataStructures.h``).

|

----

AS_unspecify_phase
--------------------

Removes a phase imposed by ``AS_specify_phase`` from a Low-Level state handle.  Returns ``Handle`` unchanged.::

    AS_unspecify_phase(Handle)

Where,

* `Handle` is a handle returned by ``AS_factory``.

|

----

AS_free
-------

Releases a Low-Level state handle created by ``AS_factory``.  Calling this is optional -- unreleased handles are automatically cleaned up when Mathcad closes -- and is intended for use as the last statement of a Mathcad program block (see above).::

    AS_free(Handle)

|

----

AS_param_index
---------------

Resolves an output parameter name (e.g. "T", "Dmolar", "Hmass") to the integer index ``AS_props``/``AS_props_multi`` expect.  Resolve once and reuse the result, rather than passing the name string on every call.::

    AS_param_index("Name")

.. note::
    This function only needs to be called **once anywhere in the worksheet**, not once per program block.  Its result is an ordinary Mathcad variable, so it can be defined at worksheet scope and referenced from any number of program blocks or independent math regions -- it is not limited to use as a local variable inside a single Mathcad program structure.

|

----

AS_input_pair_index
---------------------

Resolves an input pair name (e.g. "PT_INPUTS", "HmassP_INPUTS") to the integer index ``AS_props``/``AS_props_multi`` expect.::

    AS_input_pair_index("Name")

.. note::
    Like ``AS_param_index``, this only needs to be called **once anywhere in the worksheet** and the result reused throughout -- it is not limited to setting a local variable within a single Mathcad program structure.

    For the full list of valid input pair names, see the `CoolProp::input_pairs <https://coolprop.org/_static/doxygen/html/namespace_cool_prop.html#a85cda1634e1e4c1f76425cfd63edf155>`_ enum in the CoolProp source documentation.

|

----

AS_generate_update_pair
-------------------------

The reverse direction from ``AS_input_pair_index``: given two output-parameter indices, in *either* order, resolves which named input pair they form and returns that name as a string, for further use with ``AS_input_pair_index``/``AS_update``/``AS_props``/``AS_props_multi``.::

    AS_generate_update_pair(ParamIdx1, ParamIdx2)

Where,

* `ParamIdx1`, `ParamIdx2` are output parameter indices from ``AS_param_index``, in either order.

Raises a Custom Error if the two parameters don't form any known input pair.

.. note::
    **No Handle argument:** unlike the other Low-Level functions, this one takes no ``Handle`` -- ``CoolProp::generate_update_pair()`` (the function this wraps) is a pure lookup over the two parameter keys, not tied to any particular fluid/mixture state.

.. note::
    **No value arguments either:** ``generate_update_pair()``'s own signature takes two values alongside the two keys, but its pair-selection logic (a long chain of key-only comparisons) never inspects them -- they exist solely to get copied into its ``out1``/``out2`` parameters in the resolved pair's order, which this function doesn't surface anyway (a Mathcad Custom Function returns one value, and this one returns the resolved name). Passing values through for no purpose would just be dead arguments, so this function only takes the two indices, calling ``generate_update_pair()`` with dummy placeholder values internally. The resolved name itself already answers the ordering question ``out1``/``out2`` exist for: e.g. ``"PT_INPUTS"`` unambiguously means pressure first, temperature second, regardless of which order `ParamIdx1`/`ParamIdx2` were supplied in.

|

----

AS_update
---------

Updates a Low-Level state handle to a new state point without returning any output.  Returns ``Handle`` unchanged, so a downstream Low-Level call that uses this function's return value as its own ``Handle`` argument depends on it.  Pair with ``AS_get`` to update once and then read as many outputs as needed with separate calls, without re-running the flash for each one -- an alternative to ``AS_props``/``AS_props_multi`` when many outputs are wanted from the same point.::

    AS_update(Handle, InputPairIdx, Value1, Value2)

Where,

* `Handle` is a handle returned by ``AS_factory``.
* `InputPairIdx` is an input pair index from ``AS_input_pair_index``.
* `Value1`, `Value2` are the two input property values for that input pair.

|

----

AS_get
------

Returns one output parameter from a Low-Level state handle's *current* point -- i.e. whatever ``AS_update`` (or ``AS_props``) last set it to.::

    AS_get(Handle, ParamIdx)

Where,

* `Handle` is a handle returned by ``AS_factory``.
* `ParamIdx` is an output parameter index from ``AS_param_index``.

**EXAMPLE:**

    Temperature and density at 101325 Pa, 1 kg/kg quality (saturated vapor), updating once and reading two outputs:

    :math:`h := AS\_factory("HEOS",\ "Water")`

    :math:`iPQ := AS\_input\_pair\_index("PQ\_INPUTS")`

    :math:`iT := AS\_param\_index("T")`

    :math:`i\rho := AS\_param\_index("Dmolar")`

    :math:`h := AS\_update(h,\ iPQ,\ 101325,\ 1)`

    :math:`T := AS\_get(h,\ iT) = 373.1`

    :math:`\rho := AS\_get(h,\ i\rho)`

|

----

AS_get_sat_liquid / AS_get_sat_vapor
--------------------------------------

Like ``AS_get`` above, but read the saturated liquid/vapor side of the handle's current point rather than the bulk state -- meaningful when the current point is in the two-phase region, e.g. after a ``Q`` (quality)-based update.::

    AS_get_sat_liquid(Handle, ParamIdx)
    AS_get_sat_vapor(Handle, ParamIdx)

Where,

* `Handle` is a handle returned by ``AS_factory``.
* `ParamIdx` is an output parameter index from ``AS_param_index``.

|

----

AS_mole_fractions_liquid / AS_mole_fractions_vapor
------------------------------------------------------

The saturated liquid/vapor side's mole fractions at the handle's current point, as a column vector.::

    AS_mole_fractions_liquid(Handle, Trigger)
    AS_mole_fractions_vapor(Handle, Trigger)

Where,

* `Handle` is a handle returned by ``AS_factory``.
* `Trigger` is unused -- just pass a dummy integer (``0``), or see the note below for a better choice.

Requires the current point to actually be in the two-phase region (``0 <= quality <= 1``); raises a Custom Error otherwise.

.. note::
    **Why Trigger, when Handle is already an argument:** this is *not* about satisfying Mathcad's one-argument minimum -- ``Handle`` already does that on its own. The real reason is that ``Handle``'s own value never changes when the AbstractState it names is mutated in place: ``AS_update``, ``AS_props``, and ``AS_specify_phase`` all echo ``Handle`` back unchanged, by design (see ``AS_update``'s entry above). So a cell whose only input is ``Handle`` gives Mathcad's dependency graph nothing to key a recalculation on when the underlying point moves. Wire ``Trigger`` to whatever value actually drives the state you want reflected here -- e.g. the quality or mole-fraction value fed into the ``AS_update``/``AS_props`` call that put the state in the two-phase region this function reads -- and this cell re-evaluates whenever that does, instead of needing a full **Recalculate Worksheet**. If this cell already references the freshly-reassigned ``Handle`` from that same update (the normal chaining idiom), that alone may already provide the dependency edge; ``Trigger`` is the explicit fallback for call shapes where it doesn't.

|

----

AS_props
--------

Updates a Low-Level state handle for one input point and returns one output value.::

    AS_props(Handle, InputPairIdx, Value1, Value2, ParamIdx)

Where,

* `Handle` is a handle returned by ``AS_factory``.
* `InputPairIdx` is an input pair index from ``AS_input_pair_index``.
* `Value1`, `Value2` are the two input property values for that input pair.
* `ParamIdx` is an output parameter index from ``AS_param_index``.

**EXAMPLE:**

    Temperature at 101325 Pa, 1 kg/kg quality (saturated vapor):

    :math:`h := AS\_factory("HEOS",\ "Water")`

    :math:`iPQ := AS\_input\_pair\_index("PQ\_INPUTS")`

    :math:`iT := AS\_param\_index("T")`

    :math:`T := AS\_props(h,\ iPQ,\ 101325,\ 1,\ iT) = 373.1`

|

----

AS_props_multi
---------------

Updates a Low-Level state handle for a range of input points and returns up to 5 requested output parameters as a table (one row per input point, one column per requested output) in a single call -- the function to use when evaluating many state points against the same fluid/mixture, since it evaluates the whole array with one native flash loop rather than one Mathcad call per point.::

    AS_props_multi(Handle, InputPairIdx, Value1Array, Value2Array, ParamIdxArray)

Where,

* `Handle` is a handle returned by ``AS_factory``.
* `InputPairIdx` is an input pair index from ``AS_input_pair_index``.
* `Value1Array`, `Value2Array` are column vectors of the two input property values, one row per point (both must be the same length).
* `ParamIdxArray` is a column vector of 1 to 5 output parameter indices from ``AS_param_index``.

|

----

AS_build_phase_envelope
------------------------

Traces the phase envelope (dew/bubble curve) for a Low-Level state handle.  Call once before ``AS_get_phase_envelope_data`` on that handle.  Returns ``Handle`` unchanged, so a downstream Low-Level call that uses this call's return value as its own ``Handle`` argument depends on it.::

    AS_build_phase_envelope(Handle, Level)

Where,

* `Handle` is a handle returned by ``AS_factory``.
* `Level` (string) controls how much extra refining is done between traced points: ``"none"`` (**CoolProp's own recommendation** -- skips refining), ``"fine"`` (default tolerances -- any value other than ``"none"``/``"veryfine"`` behaves the same way), or ``"veryfine"`` (tighter tolerances, more points).

|

----

AS_get_phase_envelope_data
----------------------------

Returns the phase envelope traced by ``AS_build_phase_envelope`` as a table: one row per point, columns ``T``, ``P``, ``rhomolar_vap``, ``rhomolar_liq`` -- matching the fields returned by the C++/Python ``get_phase_envelope_data`` interface, minus the per-component compositions (see the note below).::

    AS_get_phase_envelope_data(Handle, Trigger)

Where,

* `Handle` is a handle returned by ``AS_factory``, after a prior ``AS_build_phase_envelope`` call.
* `Trigger` is unused -- just pass a dummy integer (``0``), or see ``AS_mole_fractions_liquid``'s note above for a better choice.

Raises a Custom Error if ``AS_build_phase_envelope`` hasn't been called yet for this ``Handle``.

.. note::
    **Compositions not included:** the C++/Python interface's ``x``/``y`` per-component compositions (an ``N`` x ``Ncomp`` matrix per phase) are not part of this table -- a meaningfully different, mixture-size-dependent shape. Not implemented for now; a dedicated getter could be added later if needed.

.. note::
    **Cricondentherm / cricondenbar:** see ``AS_pe_tmax``/``AS_pe_pmax`` below.

|

----

AS_pe_tmax
-----------

The cricondentherm -- the point on the phase envelope traced by ``AS_build_phase_envelope`` with the highest temperature -- as a 2-element column vector ``[T; P]``.::

    AS_pe_tmax(Handle, Trigger)

Where,

* `Handle` is a handle returned by ``AS_factory``, after a prior ``AS_build_phase_envelope`` call.
* `Trigger` is unused -- just pass a dummy integer (``0``), or see ``AS_mole_fractions_liquid``'s note above for a better choice.

Raises a Custom Error if ``AS_build_phase_envelope`` hasn't been called yet for this ``Handle``.

.. note::
    **How this is computed:** CoolProp tracks this same point internally while tracing the envelope (``PhaseEnvelopeData::iTsat_max``, set in ``PhaseEnvelopeRoutines::finalize()``), but doesn't expose it through the public Low-Level C API this wrapper is built on -- extending that shared surface (used by every CoolProp wrapper, not just Mathcad's) is out of scope here. Instead, this function fetches the same table ``AS_get_phase_envelope_data`` returns and scans its ``T`` column for the max, entirely on the Mathcad-wrapper side.

.. note::
    **Exactness:** for most mixtures ("Type I", where the traced curve's pressure rises to a single peak then falls), CoolProp doesn't just pick the closest already-traced point for the cricondentherm -- it fits a spline through the nearby points, solves for where :math:`dT_{sat}/dP_{sat} = 0`, and inserts that exact solved point into the envelope. Since that insertion happens before this wrapper ever sees the data, the max-scan above lands on that same exact point. For other mixtures ("Type II"), no such insertion happens, and the result is only as good as how finely the curve was traced -- see ``AS_build_phase_envelope``'s ``Level`` argument to trace more finely if that matters.

|

----

AS_pe_pmax
-----------

The cricondenbar -- the point on the phase envelope traced by ``AS_build_phase_envelope`` with the highest pressure -- as a 2-element column vector ``[T; P]``.::

    AS_pe_pmax(Handle, Trigger)

Where,

* `Handle` is a handle returned by ``AS_factory``, after a prior ``AS_build_phase_envelope`` call.
* `Trigger` is unused -- just pass a dummy integer (``0``), or see ``AS_mole_fractions_liquid``'s note above for a better choice.

Raises a Custom Error if ``AS_build_phase_envelope`` hasn't been called yet for this ``Handle``.

See ``AS_pe_tmax``'s notes above -- both functions work identically, this one scanning ``P`` instead of ``T`` (CoolProp's internal counterpart is ``PhaseEnvelopeData::ipsat_max``).

|

----

AS_list_handles / AS_list_states
---------------------------------

An introspection pair -- mainly useful for debugging -- that lists every Low-Level state currently open anywhere in the worksheet, i.e. every live handle from an ``AS_factory`` call that hasn't been released via ``AS_free`` (or superseded by a later call to ``AS_factory`` with the same Backend/Fluids, per its memoization).::

    AS_list_handles(Trigger)
    AS_list_states(Trigger)

Where,

* `Trigger` is unused by either function -- Mathcad Custom Functions require at least one argument, and there is no argument that naturally belongs to a whole-registry snapshot, so this exists only to satisfy that requirement. Any real scalar works, e.g. a literal ``0``.

``AS_list_handles`` returns a column vector of the currently-live handles; ``AS_list_states`` returns their ``"Backend|Fluids"`` keys (the same string ``AS_factory``'s two arguments were joined into) as one ``";"``-delimited string, **in the same order**. Both raise a Custom Error if no Low-Level states are currently open. A handle released via ``AS_free`` (or otherwise gone dead) is dropped from the listing automatically -- neither function ever reports a stale handle.

.. note::
    **Why two functions:** a Mathcad Custom Function can only return one value -- either a complex array or a string, never both -- so this is the same array-plus-parallel-string pairing already used by ``get_predefined_mixture_fluids``/``get_predefined_mixture_mole_fractions`` above, applied to the Low-Level registry instead of a predefined mixture.

.. note::
    **Ordering guarantee:** the two functions independently snapshot the same underlying registry, ordered by its ``"Backend|Fluids"`` key -- an order that depends only on which keys are *currently* registered, not on when each snapshot was taken. Two calls placed on the same worksheet will therefore agree, unless an ``AS_factory``/``AS_free`` call is evaluated in between them within the same recalculation pass.

.. note::
    **Using Trigger for recalculation:** since ``Trigger``'s value is otherwise ignored, wiring it to a handle already on the sheet -- rather than a bare literal -- gives Mathcad a real dependency edge, so the call re-runs whenever that handle's defining cell does. Without that, use **Recalculate Worksheet** to refresh these two calls, since they otherwise have no dependency edge to anything that changed.

|

----

Applying Mathcad Units to CoolProp Functions
============================================

Mathcad has a built-in units system that allows variables and values to be defined with a specific set of units.  However, Custom Functions provided through DLL add-ins are C++ functions and cannot handle Mathcad's units on inputs or outputs to the functions.  If using values with units, the appropriate unitless values in SI can be provided to the CoolProp function calls and the appropriate units applied ot the results.

.. note::
    To strip units from a Mathcad variable, yet provide the numerical value scaled to a specific unit quantity, the a variable containing units can simply be divided by the desired units expression.  The value will become unitless, but will be scaled to the specified units expression.

    Example:       :math:`P_{psi} / Pa\ =\ P`     (scaled to units of Pascals)

    This is the technique used for plotting input ranges in a specific set of units using Mathcad's 2D Chart Component.

A simple example of a call to ``PropsSI`` using variables with units is,

    Define temperature (with units):       :math:`T\ :=\ 72\ °F`

    Define pressure (with units):             :math:`P\ :=\ 1\ atm`

    Evaluate:       :math:`h\ :=\ PropsSI("H",\ "T",\ \dfrac{T}{K},\ "P",\ \dfrac{P}{Pa},\ "Water")\cdot\dfrac{J}{kg}`

    Show :math:`h` in English Engineering Units:       :math:`h\ =\ 40.133\ \dfrac{BTU}{lb}`

.. note::
    Technically, if input variables are not "stripped" of units, they will be passed as values in Mathcad's Base Units, which are SI.  This is compatible with CoolProp's base units of SI and will work.  However, units still have to be applied to the result and units should be stripped explicitely, as shown above as Mathcad allows the Base Units to be changed.  This will guarantee consistency of units between Mathcad and CoolProp.
