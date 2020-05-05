.. _high_level_api:

********************
High-Level Interface
********************

.. contents:: :depth: 2

PropsSI function
----------------

For many users, all that is needed is a simple call to the ``PropsSI`` function for pure fluids, pseudo-pure fluids and mixtures.  For humid air properties, see :ref:`Humid air properties <Humid-Air>`.  An example using ``PropsSI``:

.. ipython::

    # Import the PropsSI function
    In [1]: from CoolProp.CoolProp import PropsSI

    # Saturation temperature of Water at 1 atm in K
    In [2]: PropsSI('T','P',101325,'Q',0,'Water')

In this example, the first parameter, :math:`T`, is the *output* property that will be returned from ``PropsSI``.  The second and fourth parameters are the specified *input pair* of properties that determine the state point where the output property will be calculated.  The *output* property and *input pair* properties are text strings and must be quoted.  The third and fifth parameters are the *values* of the *input pair* properties and will determine the state point.  The sixth and last parameter is the fluid for which the *output* property will be calculated; also a quoted string.

More information:

* :ref:`Table of inputs to PropsSI function <parameter_table>`
* :ref:`More examples of the high-level API <Props_Sample>`
* :cpapi:`Documentation for all high-level functions exposed <CoolPropLib.h>`

All :ref:`the wrappers <wrappers>` wrap this function in exactly the same way.

For pure and pseudo-pure fluids, two state variables are required to fix the state.  The equations of state are based on :math:`T` and :math:`\rho` as state variables, so :math:`T, \rho` will always be the fastest inputs.  :math:`P,T` will be a bit slower (3-10 times), followed by input pairs where neither :math:`T` nor :math:`\rho` are specified, like :math:`P,H`; these will be much slower.  If speed is an issue, you can look into table-based interpolation methods using :ref:`TTSE or bicubic interpolation <Tabular_Interpolation>`; or if you are only interested in Water properties, you can look into using the :ref:`IF97 <IF97>` (industrial formulation) backend.

Vapor, Liquid, and Saturation States
------------------------------------

If the input pair (say, :math:`P,T`) defines a state point that lies in the *vapor* region, then the vapor property at that state point will be returned.  Likewise, if the state point lies in the *liquid* region, then the liquid state property at that state point will be returned.  If the state point defined by the input pair lies within 1E-4 % of the saturation pressure, then CoolProp may return an error, because both *liquid* and *vapor* are defined along the saturation curve.

To retrieve either the *vapor* or *liquid* properties along the saturation curve, provide an input pair that includes either the saturation temperature, :math:`T`, or saturation pressure, :math:`p`, along with the *vapor quality*, :math:`Q`.  Use a value of :math:`Q=1` for the saturated *vapor* property or :math:`Q=0` for the saturated *liquid* property.  For example, at a saturation pressure of 1 atm, the *liquid* and *vapor* enthalpies can be returned as follows.

.. ipython::

    # Import the PropsSI function
    In [1]: from CoolProp.CoolProp import PropsSI

    # Saturated vapor enthalpy of Water at 1 atm in J/kg
    In [2]: H_V = PropsSI('H','P',101325,'Q',1,'Water'); print(H_V)

    # Saturated liquid enthalpy of Water at 1 atm in J/kg
    In [3]: H_L = PropsSI('H','P',101325,'Q',0,'Water'); print(H_L)

    # Latent heat of vaporization of Water at 1 atm in J/kg
    In [4]: H_V - H_L

.. note::
   The *latent heat of vaporization* can be calculated using the difference between the vapor and liquid enthalpies at the same point on the saturation curve.

Imposing the Phase (Optional)
-----------------------------

Each call to ``PropsSI()`` requires the phase to be determined based on the provided input pair, and may require a non-trivial flash calculation to determine if the state point is in the single-phase or two-phase region and to generate a sensible initial guess for the solver. For computational efficiency, ``PropsSI()`` allows the phase to be manually imposed through the input key parameters.  If unspecified, PropsSI will attempt to determine the phase automatically.

Depending on the input pair, there may or may not be a speed benefit to imposing a phase.  However, some state points may not be able to find a suitable initial guess for the solver and being able to impose the phase manually may offer a solution if the solver is failing.  Additionally,  with an input pair in the two-phase region, it can be useful to impose a liquid or gas phase to instruct ``PropsSI()`` to return the saturated liquid or saturated gas properties.

To specify the phase to be used, add the "|" delimiter to one (and only one) of the input key strings followed by one of the phase strings in the table below:

+---------------------------------+----------------------------------------------------+
| Phase String                    | Phase Region                                       |
+=================================+====================================================+
| "liquid"                        | p < pcrit & T < Tcrit ; above saturation           |
+---------------------------------+----------------------------------------------------+
| "gas"                           | p < pcrit & T < Tcrit ; below saturation           |
+---------------------------------+----------------------------------------------------+
| "twophase"                      | p < pcrit & T < Tcrit ; mixed liquid/gas           |
+---------------------------------+----------------------------------------------------+
| "supercritical_liquid"          | p > pcrit & T < Tcrit                              |
+---------------------------------+----------------------------------------------------+
| "supercritical_gas"             | p < pcrit & T > Tcrit                              |
+---------------------------------+----------------------------------------------------+
| "supercritical"                 | p > pcrit & T > Tcrit                              |
+---------------------------------+----------------------------------------------------+
| "not_imposed"                   | (Default) CoolProp to determine phase              |
+---------------------------------+----------------------------------------------------+

For example:

.. ipython::

    # Get the density of Water at T = 461.1 K and P = 5.0e6 Pa, imposing the liquid phase
    In [0]: PropsSI('D','T|liquid',461.1,'P',5e6,'Water')

    # Get the density of Water at T = 597.9 K and P = 5.0e6 Pa, imposing the gas phase
    In [0]: PropsSI('D','T',597.9,'P|gas',5e6,'Water')

On each call to ``PropsSI()``, the imposed phase is reset to "not_imposed" as long as no imposed phase strings are used.  A phase string must be appended to an Input key string on each and every call to ``PropsSI()`` to impose the phase.  ``PropsSI()`` will return an error for any of the following syntax conditions:

* If anything other than the pipe, "|", symbol is used as the delimiter
* If the phase string is not one of the valid phase strings in the table above
* If the phase string is applied to more than one of the Input key parameters

In addition, for consistency with the low-level interface, the valid phase strings in the table above may be prefixed with either "phase_" or "iphase_" and still be recognized as a valid phase string.

.. warning::

   When specifying an imposed phase, it is absolutely **critical** that the input pair actually lie within the imposed phase region.  If an incorrect phase is imposed for the given input pair, ``PropsSI()`` may throw unexpected errors or incorrect results may possibly be returned from the property functions.  If the state point phase is not absolutely known, it is best to let CoolProp determine the phase.

Fluid information
-----------------

You can obtain string-encoded information about the fluid from the ``get_fluid_param_string`` function with something like:

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    In [1]: for k in ['formula','CAS','aliases','ASHRAE34','REFPROP_name','pure','INCHI','INCHI_Key','CHEMSPIDER_ID']:
       ...:     item = k + ' --> ' + CP.get_fluid_param_string("R125", k)
       ...:     print(item)
       ...:

.. _trivial_inputs:

Trivial inputs
--------------


In order to obtain trivial inputs that do not depend on the thermodynamic state, in wrappers that support the ``Props1SI`` function, you can obtain the trivial parameter (in this case the critical temperature of water) like:

    Props1SI("Tcrit","Water")

In python, the ``PropsSI`` function is overloaded to also accept two inputs:

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    In [1]: CP.PropsSI("Tcrit","Water")

    In [1]: CP.PropsSI("Tcrit","REFPROP::WATER")

Furthermore, you can in all languages call the ``PropsSI`` function directly using dummy arguments for the other unused parameters:

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    In [1]: CP.PropsSI("Tcrit","",0,"",0,"Water")

PhaseSI function
----------------

It can be useful to know what the phase of a given state point is.  A high-level function called ``PhaseSI`` has been implemented to allow for access to the phase.

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    In [5]: CP.PhaseSI('P',101325,'Q',0,'Water')

The phase index (as floating point number) can also be obtained using the PropsSI function. In python you would do:

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    In [5]: CP.PropsSI('Phase','P',101325,'Q',0,'Water')

where you can obtain the integer indices corresponding to the phase flags using the ``get_phase_index`` function:

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    In [6]: CP.get_phase_index('phase_twophase')

    # Or for liquid
    In [6]: CP.get_phase_index('phase_liquid')

For a given fluid, the phase can be plotted in T-p coordinates:

.. plot::

    import matplotlib
    import numpy as np
    import CoolProp as CP
    import matplotlib.pyplot as plt
    import scipy.interpolate

    Water = CP.AbstractState("HEOS", "Water")
    pc = Water.keyed_output(CP.iP_critical)
    Tc = Water.keyed_output(CP.iT_critical)
    Tmin = 200
    Tmax = 1000
    pmax = Water.keyed_output(CP.iP_max)
    pt = 611.657
    Tt = 273.16
    fillcolor = 'g'

    fig = plt.figure(figsize = (6,6))
    ax = fig.add_subplot(111)
    lw = 3

    # --------------
    # Melting curve
    # --------------
    melt_args = dict(lw = lw, solid_capstyle = 'round')
    TT = []
    PP = list(np.logspace(np.log10(pt), np.log10(pmax),1000))
    for p in PP:
        TT.append(Water.melting_line(CP.iT, CP.iP, p))

    #Zone VI
    for T in np.linspace(max(TT), 355):
        TT.append(T)
        theta = T/273.31
        pi = 1-1.07476*(1-theta**4.6)
        p = pi*632.4e6
        PP.append(p)

    plt.plot(TT,PP,'darkblue',**melt_args)

    # ----------------
    # Saturation curve
    # ----------------
    Ts = np.linspace(273.16, Tc, 1000)
    ps = CP.CoolProp.PropsSI('P','T',Ts,'Q',0,'Water')

    # ------
    # Labels
    # ------

    plt.plot(Ts,ps,'orange',lw = lw, solid_capstyle = 'round')

    # Critical lines
    plt.axvline(Tc, dashes = [2, 2])
    plt.axhline(pc, dashes = [2, 2])

    # Labels
    plt.text(850, 1e8, 'supercritical',ha= 'center')
    plt.text(850, 1e5, 'supercritical_gas', rotation = 90)
    plt.text(450, 1e8, 'supercritical_liquid', rotation = 0, ha = 'center')
    plt.text(350, 3e6, 'liquid', rotation = 45)
    plt.text(450, 5e4, 'gas', rotation = 45)

    plt.ylim(611,1e9)
    plt.gca().set_yscale('log')
    plt.gca().set_xlim(240, 1000)
    plt.ylabel('Pressure [Pa]')
    plt.xlabel('Temperature [K]')
    plt.tight_layout()

.. _partial_derivatives_high_level:

Partial Derivatives
-------------------

First Partial Derivatives for Single-phase States
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For some applications it can be useful to have access to partial derivatives of thermodynamic properties.  A generalized first partial derivative has been implemented into CoolProp, which can be obtained using the ``PropsSI`` function by encoding the desired derivative as a string.  The format of the string is ``d(OF)/d(WRT)|CONSTANT`` which is the same as

.. math::

    \left. \frac{\partial OF}{\partial WRT}\right|_{CONSTANT}

At the low-level, the CoolProp code calls the function :cpapi:`AbstractState::first_partial_deriv`.  Refer to the function documentation to see how the generalized derivative works.

.. warning::

    This derivative formulation is currently only valid for homogeneous (single-phase) states.  Two phase derivatives are not defined, and are for many combinations, invalid.

Here is an example of calculating the constant pressure specific heat, which is defined by the relation

.. math::

    c_p = \left.\frac{\partial h}{\partial T}\right|_{p}

and called through python

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    # c_p using c_p
    In [5]: CP.PropsSI('C','P',101325,'T',300,'Water')

    # c_p using derivative
    In [5]: CP.PropsSI('d(Hmass)/d(T)|P','P',101325,'T',300,'Water')

It is also possible to call the derivatives directly using the :ref:`low-level partial derivatives functionality <partial_derivatives_low_level>`.  The low-level routine is in general faster because it avoids the string parsing.

Second Partial Derivatives for Single-Phase States
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a similar fashion it is possible to evaluate second derivatives.  For instance, the derivative of :math:`c_p` with respect to mass-based specific enthalpy at constant pressure could be obtained by

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    # c_p using derivative
    In [1]: CP.PropsSI('d(d(Hmass)/d(T)|P)/d(Hmass)|P','P',101325,'T',300,'Water')

where the inner part ``d(Hmass)/d(T)|P`` is the definition of :math:`c_p`.

.. warning::

    This derivative formulation is currently only valid for homogeneous (single-phase) states.  Two phase derivatives are not defined, and are for many combinations, invalid.

It is also possible to call the derivatives directly using the :ref:`low-level partial derivatives functionality <partial_derivatives_low_level>`.  The low-level routine is in general faster because it avoids the string parsing.

First Saturation Derivatives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is also possible to retrieve the derivatives along the saturation curves using the high-level interface, encoding the desired derivative as a string just like for the single-phase derivatives.

.. warning::

    This derivative formulation is currently only valid for saturated states where the vapor quality is either 0 or 1.

For instance, to calculate the saturation derivative of enthalpy ALONG the saturated vapor curve, you could do:

.. ipython::

    In [1]: import CoolProp

    In [1]: CoolProp.CoolProp.PropsSI('d(Hmolar)/d(T)|sigma','P',101325,'Q',1,'Water')

It is also possible to call the derivatives directly using the :ref:`low-level partial derivatives functionality <partial_derivatives_low_level>`.  The low-level routine is in general faster because it avoids the string parsing.

.. _predefined_mixtures:

Predefined Mixtures
-------------------

A number of predefined mixtures are included in CoolProp.  You can retrieve the list of predefined mixtures by calling ``get_global_param_string("predefined_mixtures")`` which will return a comma-separated list of predefined mixtures.  In Python, to get the first 6 mixtures, you would do

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    In [1]: CP.get_global_param_string('predefined_mixtures').split(',')[0:6]

and then to calculate the density of air using the mixture model at 1 atmosphere (=101325 Pa) and 300 K, you could do

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    In [1]: CP.PropsSI('D','P',101325,'T',300,'Air.mix')

Exactly the same methodology can be used from other wrappers.

User-Defined Mixtures
---------------------

When using mixtures in CoolProp, you can specify mixture components and composition by encoding the mixture components and mole fractions by doing something like

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    In [1]: CP.PropsSI('D','T',300,'P',101325,'HEOS::R32[0.697615]&R125[0.302385]')

You can handle ternary and multi-component mixtures in the same fashion, just add the other components to the fluid string with a ``&`` separating components and the fraction of the component in ``[`` and ``]`` brackets

.. _high_level_set_reference_state:

Reference States
----------------

Enthalpy and entropy are *relative* properties!  You should always be comparing *differences* in enthalpy rather than absolute values of the enthalpy or entropy.  That said, if can be useful to set the reference state values for enthalpy and entropy to one of a few standard values.  This is done by the use of the ``set_reference_state`` function in python, or the ``set_reference_stateS`` function most everywhere else.  For documentation of the underlying C++ function, see :cpapi:`CoolProp::set_reference_stateS`.

.. warning::

    The changing of the reference state should be part of the initialization of your program, and it is not recommended to change the reference state during the course of making calculations

A number of reference states can be used:

* ``IIR``: h = 200 kJ/kg, s=1 kJ/kg/K at 0C saturated liquid
* ``ASHRAE``: h = 0, s = 0 @ -40C saturated liquid
* ``NBP``: h=0, s=0 for saturated liquid at 1 atmosphere
* ``DEF``: Go back to the default reference state for the fluid

which can be used like

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    In [1]: CP.set_reference_state('n-Propane','ASHRAE')

    # Should be zero (or very close to it)
    In [1]: CP.PropsSI('H', 'T', 233.15, 'Q', 0, 'n-Propane')

    # Back to the original value
    In [1]: CP.set_reference_state('n-Propane','DEF')

    # Should not be zero
    In [1]: CP.PropsSI('H', 'T', 233.15, 'Q', 0, 'n-Propane')

For non-standard reference states, you can specify them directly for the given temperature and molar(!) density.  Here is an example of setting the enthalpy and entropy at 298.15 K and 1 atm to some specified values

.. ipython::
    
    In [1]: import CoolProp.CoolProp as CP

    In [1]: Dmolar = CP.PropsSI("Dmolar", "T", 298.15, "P", 101325, "NH3") 

    In [1]: CP.set_reference_state('NH3', 298.15, Dmolar, -60000.12, 314.159) # fluid, T, D (mol/m^3), h (J/mol), s (J/mol/K)

    In [1]: CP.PropsSI("Hmolar", "T", 298.15, "P", 101325, "NH3")

    In [1]: CP.PropsSI("Smolar", "T", 298.15, "P", 101325, "NH3")

    # Back to the original value
    In [1]: CP.set_reference_state('NH3','DEF')

Calling REFPROP
---------------

If you have the `REFPROP library <http://www.nist.gov/srd/nist23.cfm>`_ installed, you can call REFPROP in the same way that you call CoolProp, but with ``REFPROP::`` preceding the fluid name. For instance, as in python:

.. ipython::

    In [1]: import CoolProp.CoolProp as CP

    # Using properties from CoolProp to get R410A density
    In [2]: CP.PropsSI('D','T',300,'P',101325,'HEOS::R32[0.697615]&R125[0.302385]')

    # Using properties from REFPROP to get R410A density
    In [2]: CP.PropsSI('D','T',300,'P',101325,'REFPROP::R32[0.697615]&R125[0.302385]')

Adding Fluids
-------------

The fluids in CoolProp are all compiled into the library itself, and are given in the `JSON <http://json.org>`_ format.  They are all stored in the ``dev/fluids`` folder relative to the root of the repository.  If you want to obtain the JSON data for a fluid from CoolProp, print out a part of it, and then load it back into CoolProp, you could do:

.. ipython::
    :okexcept:

    In [1]: import CoolProp.CoolProp as CP

    # Get the JSON structure for Water
    In [2]: jj = CP.get_fluid_param_string("Water", "JSON")

    # Now load it back into CoolProp - Oops! This isn't going to work because it is already there.
    In [2]: CP.add_fluids_as_JSON("HEOS", jj)

    # Set the configuration variable allowing for overwriting
    In [2]: CP.set_config_bool(CP.OVERWRITE_FLUIDS, True)

    # Now load it back into CoolProp - Success!
    In [2]: CP.add_fluids_as_JSON("HEOS", jj)

    # Turn overwriting back off
    In [2]: CP.set_config_bool(CP.OVERWRITE_FLUIDS, False)

C++ Sample Code
---------------
.. literalinclude:: snippets/propssi.cxx
   :language: c++

C++ Sample Code Output
----------------------
.. literalinclude:: snippets/propssi.cxx.output


.. _Props_Sample:

Sample Code
-----------

.. ipython::

    In [1]: import CoolProp as CP

    In [1]: print(CP.__version__)

    In [1]: print(CP.__gitrevision__)

    #Import the things you need
    In [1]: from CoolProp.CoolProp import PropsSI

    # Specific heat (J/kg/K) of 20% ethylene glycol as a function of T
    In [2]: PropsSI('C','T',298.15,'P',101325,'INCOMP::MEG-20%')

    # Density of Air at standard atmosphere in kg/m^3
    In [2]: PropsSI('D','T',298.15,'P',101325,'Air')

    # Saturation temperature of Water at 1 atm
    In [2]: PropsSI('T','P',101325,'Q',0,'Water')

    # Saturated vapor enthalpy of R134a at 0C (Q=1)
    In [2]: PropsSI('H','T',273.15,'Q',1,'R134a')

    # Saturated liquid enthalpy of R134a at 0C (Q=0)
    In [2]: PropsSI('H','T',273.15,'Q',0,'R134a')

    # Using properties from CoolProp to get R410A density
    In [2]: PropsSI('D','T',300,'P',101325,'HEOS::R32[0.697615]&R125[0.302385]')

    # Using properties from REFPROP to get R410A density
    In [2]: PropsSI('D','T',300,'P',101325,'REFPROP::R32[0.697615]&R125[0.302385]')

    # Check that the same as using pseudo-pure
    In [2]: PropsSI('D','T',300,'P',101325,'R410A')

    # Using IF97 to get Water saturated vapor density at 100C
    In [2]: PropsSI('D','T',400,'Q',1,'IF97::Water')

    # Check the IF97 result using the default HEOS
    In [2]: PropsSI('D','T',400,'Q',1,'Water')

.. _parameter_table:

Table of string inputs to PropsSI function
------------------------------------------

.. note::

   Please note that any parameter that is indicated as a trivial parameter can be obtained from the ``Props1SI`` function as shown above in :ref:`trivial_inputs`

.. include:: parameter_table.rst.in
