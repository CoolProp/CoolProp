.. _examples:

Examples
========
This page serves as a teaser of the functionality of CoolProp.  These examples are written in the Python programming language.  For more information see:

- :ref:`Pure and Pseudo-Pure fluid properties <Fluid-Properties>`
- :ref:`Mixture properties <mixtures>`
- :ref:`Wrapper-specific code examples <wrappers>`

Sample Props Code
-----------------
To use the ``PropsSI`` function, import it and do some calls, do something like this

.. ipython::

    # Import the things you need
    In [1]: from CoolProp.CoolProp import PropsSI

    # Print some information on the currently used version of coolprop
    In [1]: import CoolProp; print(CoolProp.__version__, CoolProp.__gitrevision__)

    # Density of carbon dioxide at 100 bar and 25C
    In [2]: PropsSI('D', 'T', 298.15, 'P', 100e5, 'CO2')

    # Saturated vapor enthalpy [J/kg] of R134a at 25C
    In [2]: PropsSI('H', 'T', 298.15, 'Q', 1, 'R134a')

All the possible input and output parameters are listed in the
:ref:`High-Level API <high_level_api>` documentation

To understand more about what is going on under the hood, go to :ref:`Fluid-Properties` documentation.


Sample HAPropsSI Code
---------------------
To use the ``HAPropsSI`` function, import it and do some calls, do something like this

.. ipython::

    # import the things you need
    In [1]: from CoolProp.HumidAirProp import HAPropsSI

    # Enthalpy (J per kg dry air) as a function of temperature, pressure,
    #    and relative humidity at STP
    In [2]: h = HAPropsSI('H','T',298.15,'P',101325,'R',0.5); print(h)

    # Temperature of saturated air at the previous enthalpy
    In [2]: T = HAPropsSI('T','P',101325,'H',h,'R',1.0); print(T)

    # Temperature of saturated air - order of inputs doesn't matter
    In [2]: T = HAPropsSI('T','H',h,'R',1.0,'P',101325); print(T)

Or go to the :ref:`Humid-Air` documentation.

Plotting
--------
.. toctree::
    :maxdepth: 1

    python-plotting.rst
    python-cycles.rst
