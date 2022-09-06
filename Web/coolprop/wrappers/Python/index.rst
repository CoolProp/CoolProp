.. _Python:

**************
Python Wrapper
**************

.. contents:: :depth: 2

PyFluids (3-party wrapper)
==========================

It is a simple, full-featured, lightweight CoolProp wrapper for Python.
PyFluids gets published on `PyPI <https://pypi.org/project/pyfluids/>`_,  so you can easily install it using: ::

    pip install pyfluids

All CoolProp features are included: thermophysical properties of pure fluids, mixtures and humid air.
Also you can easily convert the results to a JSON string or Python dict, add new properties or inputs for lookups, and more.

Benefits
--------

* Easy to use: all fluids and properties are at hand, no need to remember CoolProp keys.

* Processes for fluids and humid air are included: there is no need to code it anymore.

* User-friendly interface: writing code is faster.

Examples
--------

To calculate the specific heat of saturated water vapor at *1 atm*: ::

    from pyfluids import Fluid, FluidsList

    water_vapour = Fluid(FluidsList.Water).dew_point_at_pressure(101325)
    print(water_vapour.specific_heat)  # 2079.937085633241

To calculate the dynamic viscosity of propylene glycol aqueous solution with *60 %* mass fraction at *100 kPa* and *-20 °C*: ::

    from pyfluids import Fluid, FluidsList, Input

    propylene_glycol = Fluid(FluidsList.MPG, 60).with_state(
        Input.pressure(100e3), Input.temperature(-20)
    )
    print(propylene_glycol.dynamic_viscosity)  # 0.13907391053938878

To calculate the density of ethanol aqueous solution (with ethanol *40 %* mass fraction) at *200 kPa* and *4 °C*: ::

    from pyfluids import Mixture, FluidsList, Input

    mixture = Mixture([FluidsList.Water, FluidsList.Ethanol], [60, 40]).with_state(
        Input.pressure(200e3), Input.temperature(4)
    )
    print(mixture.density)  # 883.3922771627963

To calculate the wet bulb temperature of humid air at *99 kPa*, *30 °C* and *50 %* relative humidity: ::

    from pyfluids import HumidAir, InputHumidAir

    humid_air = HumidAir().with_state(
        InputHumidAir.pressure(99e3),
        InputHumidAir.temperature(30),
        InputHumidAir.relative_humidity(50),
    )
    print(humid_air.wet_bulb_temperature)  # 21.946578559079228

For any questions or more examples, `see PyFluids on GitHub <https://github.com/portyanikhin/PyFluids>`_.

Automatic installation
======================

Using the ``pip`` installation program, you can install the official release from the pypi server using::

    pip install CoolProp

There are also unofficial `Conda <https://conda.io>`__ packages available from the ``conda-forge`` channel. To
install, use::

    conda install conda-forge::coolprop

If you dare, you can also try the latest nightly release from :sfnightly:`Python`
or get it directly from the development server using::

    pip install -vvv --pre --trusted-host www.coolprop.dreamhosters.com --find-links http://www.coolprop.dreamhosters.com/binaries/Python/ -U --no-cache --force-reinstall CoolProp

Manual installation
===================

Compilation of the python wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

On all platforms, if it is not already there, you need Cython to be installed::

    sudo pip install Cython

Then, follow the commands::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp/wrappers/Python
    # Start the installation
    sudo python setup.py install

If you would like to install CoolProp just for a given version of Python (for
example if ``python`` links to ``python3.4`` and you also have a ``python2.7``
executable), simply use this version of python to execute the ``setup.py``
script::

    sudo python2.7 setup.py install

If you have multiple versions of Visual Studio installed and need to specify the version to use and choice of 32-bit or 64-bit compilation, you can use::

    # 64-bit using VS2008 on Pytnon 2.7
    sudo python setup.py install --cmake-compiler vc9 --cmake-bitness 64

or, equivalently::

    sudo python setup.py install cmake=vc9,64

Omitting the cmake options will use the default (latest) compiler on the machine.


Local installation
------------------

If you prefer not to be sudoer when compiling coolprop on Linux/MacOS, you can
also install it locally using the ``--user`` switch::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp/wrappers/Python
    # Start the installation
    python setup.py install --user

For Pyzo users
--------------

Suppose the directory containing pyzo is on your Desktop in
``~/Desktop/pyzo2014a/``. Then you can install CoolProp to be used within pyzo
by following the same lines as above::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp/wrappers/Python
    # Start the installation (~/Desktop/pyzo2014a/ to be changed according to
    # your effective installation)
    sudo ~/Desktop/pyzo2014a/bin/python setup.py install

Usage
=====

There is example code :ref:`at the end of this page <python_example>`

Once installed, you can use CoolProp for various things:

* Compute special values in SI units::

    import CoolProp.CoolProp as CP
    fluid = 'Water'
    pressure_at_critical_point = CP.PropsSI(fluid,'pcrit')
    # Massic volume (in m^3/kg) is the inverse of density
    # (or volumic mass in kg/m^3). Let's compute the massic volume of liquid
    # at 1bar (1e5 Pa) of pressure
    vL = 1/CP.PropsSI('D','P',1e5,'Q',0,fluid)
    # Same for saturated vapor
    vG = 1/CP.PropsSI('D','P',1e5,'Q',1,fluid)

* Get some nice graphs::

    import CoolProp.Plots as CPP
    ph_plot = CPP.PropertyPlot('Water','Ph')
    ph_plot.savefig('enthalpy_pressure_graph_for_Water.png')

* Solve `thermodynamics exercices`_

* Make your own `more complex graphs`_ if you feel the graphing interface is lacking something

* Make even more complex graphs using `3D stuff`_

.. _thermodynamics exercices: https://github.com/jjfPCSI1/py4phys/blob/master/lib/T6_resolution_cycle_diesel.py
.. _more complex graphs: https://github.com/jjfPCSI1/py4phys/blob/master/lib/T6_diagramme_Ph_coolprop.py
.. _3D stuff: https://github.com/CoolProp/CoolProp/blob/master/dev/TTSE/TTSE_ranges.py

.. _python_example:

Example Code
============

.. literalinclude:: Example.py
   :language: python

Example Code Output
===================

.. literalinclude:: Example.out

Code Warnings
=============

Messages may be issued from the Python CoolProp wrapper via the Python `warnings` module.  This module allows
non-fatal warning messages to be issued to the calling program and stdout to warn of
improper function usage or deprecation of features.  These warnings will, by
default, be issued each and every time a suspect call is made to CoolProp.  While, the best
solution is to correct the calling code according to the message received, sometimes this is
difficult to do in a legacy or third party code and can result in many, many warning messages that obscure
the output and hinder debugging.

Suppressing warning messages
----------------------------

The calling code can suppress or ignore these warning messages by overriding the default
warnings filter and changing the behavior of the warnings module.  As an example, the
following script will result in a `DeprecationWarning` on each call to the deprecated function
Props()::

    from CoolProp.CoolProp import Props
    Rho = Props('D','T',298.15,'P',10000,'R744')
    print("R744 Density at {} K and {} kPa      = {} kg/m³".format(298.15, 10000, Rho))
    H = Props('H','T',298.15,'Q',1,'R134a');
    print("R134a Saturated Liquid Enthalpy at {} K = {} kJ/kg".format(298.15, H))

Example output::

    TestProps.py:14: DeprecationWarning: Props() function is deprecated; Use the PropsSI() function
    Rho = Props('D','T',298.15,'P',10000,'R744')
    R744 Density at 298.15 K and 10000 kPa      = 817.6273812375758 kg/m³
    TestProps.py:16: DeprecationWarning: Props() function is deprecated; Use the PropsSI() function
    H = Props('H','T',298.15,'Q',1,'R134a');
    R134a Saturated Liquid Enthalpy at 298.15 K = 412.33395323186807 kJ/kg

Legacy applications can create a filter override to ignore *all* deprecation warnings by including
the following code just *after* the last import from CoolProp, but *before* any calls to CoolProp::

    import warnings
    warnings.filterwarnings('ignore', category=DeprecationWarning)

To suppress, for example, *only* deprecation warning messages that contain the string "Props()",
the second parameter to filterwarnings() can be a pattern matching regular expression::

    import warnings
    warnings.filterwarnings('ignore', '.*Props()*.', category=DeprecationWarning)

This filter will suppress any `DeprecationWarning` messages that contain the string "Props()" but will
allow all other warning messages to be displayed.  The first parameter, `ignore`, can also be set to
`once`, which will result in a given message to be issued only once and then ignored on further instances.

See `Python >>> Module Warnings <https://docs.python.org/3/library/warnings.html#module-warnings>`_ for more information on using `filterwarnings()`


Module Documentation
====================

.. toctree::

    ../../../apidoc/CoolProp.rst
