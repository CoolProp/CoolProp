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

More information: 

* :ref:`Table of inputs to PropsSI function <parameter_table>`
* :ref:`More examples of the high-level API <Props_Sample>`
* :cpapi:`Documentation for all high-level functions exposed <CoolPropLib.h>`

All :ref:`the wrappers <wrappers>` wrap this function in exactly the same way.

For pure and pseudo-pure fluids, two state points are required to fix the state.  The equations of state are based on :math:`T` and :math:`\rho` as state variables, so :math:`T, \rho` will always be the fastest inputs.  :math:`P,T` will be a bit slower (3-10 times), and then comes inputs where neither :math:`T` nor :math:`\rho` are given, like :math:`p,h`.  They will be much slower.  If speed is an issue, you can look into table-based interpolation methods using TTSE or bicubic interpolation. 

.. _trivial_inputs:

Trivial inputs
--------------


In order to obtain trivial inputs that do not depend on the thermodynamic state, in wrappers that support the ``Props1SI`` function, you can obtain the trivial parameter (in this case the critical temperature of water) like:

    Props1SI("Tcrit","Water")
    
In python, the ``PropsSI`` function is overloaded to also accept two inputs:

.. ipython::

    In [1]: import CoolProp.CoolProp as CP
    
    In [1]: CP.PropsSI("Tcrit","Water")
    
    In [1]: CP.PropsSI("Tcrit","REFPROP::Water")
    
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
    
Calling REFPROP
---------------

If you have the `REFPROP library <http://www.nist.gov/srd/nist23.cfm>`_ installed, you can call REFPROP in the same way that you call CoolProp, but with ``REFPROP::`` preceding the fluid name. For instance, as in python:

.. ipython::

    In [1]: import CoolProp.CoolProp as CP
    
    # Using properties from CoolProp to get R410A density
    In [2]: CP.PropsSI('D','T',300,'P',101325,'HEOS::R32[0.697615]&R125[0.302385]')
    
    # Using properties from REFPROP to get R410A density
    In [2]: CP.PropsSI('D','T',300,'P',101325,'REFPROP::R32[0.697615]&R125[0.302385]')
    

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
    
    In [1]: print CP.__version__
    
    In [1]: print CP.__gitrevision__
    
    #Import the things you need 
    In [1]: from CoolProp.CoolProp import PropsSI
    
    # Specific heat (J/kg/K) of 20% ethylene glycol as a function of T
    In [2]: PropsSI('C','T',298.15,'P',101325,'INCOMP::MEG-20%')
    
    # Density of Air at standard atmosphere in kg/m^3
    In [2]: PropsSI('D','T',298.15,'P',101325,'Air')
    
    # Saturation temperature of Water at 1 atm
    In [2]: PropsSI('T','P',101325,'Q',0,'Water')
    
    # Saturated vapor density of R134a at 0C
    In [2]: PropsSI('H','T',273.15,'Q',1,'R134a')
    
    # Using properties from CoolProp to get R410A density
    In [2]: PropsSI('D','T',300,'P',101325,'HEOS::R32[0.697615]&R125[0.302385]')
    
    # Using properties from REFPROP to get R410A density
    In [2]: PropsSI('D','T',300,'P',101325,'REFPROP::R32[0.697615]&R125[0.302385]')
    
    # Check that the same as using pseudo-pure
    In [2]: PropsSI('D','T',300,'P',101325,'R410A')

.. _parameter_table:

Table of string inputs to PropsSI function
------------------------------------------

.. note::
   
   Please note that any parameter that is indicated as a trivial parameter can be obtained from the ``Props1SI`` function as shown above in :ref:`trivial_inputs`
   
.. include:: parameter_table.rst.in