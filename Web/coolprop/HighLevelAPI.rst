.. _high_level_api:

********************
High-Level Interface
********************

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

PhaseSI function
----------------

It can be useful to know what the phase of a given state point is.  A high-level function called ``PhaseSI`` has been implemented to allow for access to the phase.

.. ipython::

    In [1]: import CoolProp
    
    In [5]: CoolProp.CoolProp.PhaseSI('P',101325,'Q',0,'Water')

The phase index (as floating point number) can also be obtained using the PropsSI function. In python you would do:

.. ipython::

    In [1]: import CoolProp
    
    In [5]: CoolProp.CoolProp.PropsSI('Phase','P',101325,'Q',0,'Water')
    
where you can obtain the integer indices corresponding to the phase flags using the ``get_phase_index`` function:

.. ipython::

    In [1]: import CoolProp

    In [6]: CoolProp.CoolProp.get_phase_index('phase_twophase')
    
    # Or for liquid
    In [6]: CoolProp.CoolProp.get_phase_index('phase_liquid')
    
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
    ps = CP.CoolProp.PropsSI('P','T',Ts,'Q',[0]*len(Ts),'Water',[1])

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

Code
----
.. literalinclude:: snippets/propssi.cxx
   :language: c++

Output
------
.. literalinclude:: snippets/propssi.cxx.output

.. _parameter_table:

Table of string inputs to PropsSI function
------------------------------------------

.. include:: parameter_table.rst.in

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
