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

* "Outputs" = Requested output properties string containing a delimited list of one or more valid output property names from the :ref:`Table of Valid Parameters <parameter_table>`.  For this Mathcad wrapper, the delimiter can be any one of (*comma, <space>, colon, semi-colon, or ampersand*), but must be consistant.
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

    Define a fluid:         :math:`fl` := "Water"

    Triple Point Temperature:         :math:`T_t := Props1SI("Ttriple",\ fl) = 273.15`

    Triple Point Pressure:               :math:`P_t := Props1SI("ptriple",\ fl) = 611.655`


    Critical Temperature:                :math:`T_c := Props1SI("Tcrit",\ fl) = 647.096`

    Critical Pressure:                       :math:`P_c := Props1SI("pcrit",\ fl) = 2.206\cdot10^7`

    Liquid points:       :math:`T_L := mean(T_t,T_c)\ -\ 10.0`     &     :math:`P_L := mean(P_t,P_c)`

    Calc:                :math:`PhaseSI("T",\ T_L,\ "P",\ P_L,\ fl)` = "Liquid"

|

----

HAPropsSI - Humid Air Fluid Properties
----------------------------------------------

`HA PropsSI`  is used to find fluid properties of humid air.  The physics behind the   function is based on the analysis in ASHRAE RP-1845, which is available online: https://www.tandfonline.com/doi/abs/10.1080/10789669.2009.10390874.  It employs real gas properties for both air and water, as well as the most accurate interaction parameters and enhancement factors.   RP-1845 is based largely on the IAPWS-95 formulation for the properties of water.  The calling structure of the function is as follows: ::

    PropsSI("Output", "Input1", Val1, "Input2", Val2, "Input3", Val3)

Where,

* "Output" = Requested output property string; see :ref:`Table of Valid HA Parameters <HAparameter_table>`.
* "Input1" = First state-point property string.
* Val1 = First state-point property value (scalar variable)
* "Input2" = Second state-point property string.
* Val2 = Second state-point property value (scalar variable)
* "Input3" = Third state-point property string.
* Val2 = Third state-point property value (scalar variable)
At least one of the inputs must be "T" (dry bulb temperature), "R" (Relative Humidity between 0.0 and 1.0), "W" (Humidity Ratio), or "Tdp" (dew point).

**EXAMPLE:**

    .. math::
       h := PropsSI("H",\ "T",\ 298.15,\ "P",\ 101325,\ "R",\ 0.5) = 5.042\cdot10^4

----

Pseudo-Low-Level Functions
==========================

CoolProp's Low-level functions require the creation of an Abstract State object and then evaluation of properties using that object's member functions.  Mathcad does not have the ability to store objects as variables and so the Low-Level interface to CoolProp is not implemented.  However, there are some pseudo-Low-Level functions that to not require an abstract state object, or can at least create one temporarily for the purposes of extracting and setting CoolProp data and parameters.  These wrapper functions are listed here.

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
   The ``"errsting"`` option is *extremely* useful when using CoolProp functions in Mathcad.  While the wrapper functions attempt to trap common errors and display them as meaninful Mathcad error messages, highlighting the offending parameter(s), unknown errors will display as "CoolProp Issue: Use get_global_param_string("errstring") for more info". This is the only way to see the actaul CoolProp error message being thrown, even if the error is already trapped by the Mathcad wrapper.

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

    For air:        :math:`get_predefined_mixture_fluids`("Air.mix") = "NITROGEN;ARGON;OXYGEN"

    For R401A:  :math:`get_predefined_mixture_fluids`("R401A.mix") = "R22;R152A;R124"

.. note::
   A few predefined mixtures are missing binary interaction parameters for at least one component pair.  These mixtures are defined, but cannot be used, for now, for property calculations.

----

get_predefined_mixture_fractions
--------------------------------

This is not a wrapper of a CoolProp functions, but an additional helper function use to assist using predefined mixtures. This function returns a Mathcad array (column vector) of the predefined mixture component *mole fractions*.::

    get_predefined_mixture_fractions("mixture")

Where "mixture" is a predefined mixture name ending in .mix or .MIX. Mixture names are case sensitive, using the available predefined mixture names retrieved from ``get_global_parameter_string("predefined_mixtures")``.

**EXAMPLES**

    For air:        :math:`mf_{Air}` := :math:`get_predefined_mixture_fractions`("Air.mix") = :math:`\begin{bmatrix} 0.7812\\ 0.0092 \\ 0.2096 \end{bmatrix}`

    For R401A:  :math:`mf_{R401A}` := :math:`get_predefined_mixture_fractions'("R401A.mix") = :math:`\begin{bmatrix} 0.578854\\ 0.0.185871 \\ 0.0.235274 \end{bmatrix}`

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
