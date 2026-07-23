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

Using uv (fast Python package manager)
---------------------------------------

`uv <https://docs.astral.sh/uv/>`__ is an extremely fast Python package and project manager written in Rust.
To use CoolProp with uv::

    # Install CoolProp in the current environment
    uv pip install CoolProp

    # Or create a new project with CoolProp
    uv init my-project
    cd my-project
    uv add CoolProp

    # Run a script with CoolProp (uv will automatically manage the environment)
    uv run python my_script.py

For development from source::

    # Clone the repository
    git clone https://github.com/CoolProp/CoolProp
    cd CoolProp

    # Create a virtual environment (recommended)
    uv venv

    # Install build dependencies
    uv pip install --python .venv/bin/python scikit-build-core cython

    # Install in editable mode with incremental build support
    uv pip install -ve . --python .venv/bin/python --no-build-isolation

**Important**: The ``--no-build-isolation`` flag is required for incremental builds to work properly.
Without it, each build will use a temporary isolated environment with different paths, causing CMake
to reconfigure and rebuild all files even when only one file changed.

Nightly builds
--------------

If you dare, you can also try the latest development build, which is published to
`TestPyPI <https://test.pypi.org/project/CoolProp/>`_::

    pip install --pre --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ CoolProp

The ``--extra-index-url`` lets pip pull CoolProp's dependencies (numpy, etc.) from the
regular PyPI while taking CoolProp itself from TestPyPI.

Manual installation
===================

Compilation of the python wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

CoolProp uses a modern build system based on CMake via scikit-build-core. The build dependencies
(including Cython and CMake) are automatically installed by pip.

To build and install from source::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp
    # Move into the folder you just created
    cd CoolProp
    # Install (pip will automatically handle the build)
    pip install .

Development with Editable Install
++++++++++++++++++++++++++++++++++

An editable install (also known as "development mode") allows you to make changes to the C++
source code and have them take effect after a simple rebuild, without reinstalling the package.

Using pip::

    # Install in editable mode
    pip install -ve .

Using uv (recommended for faster builds and better dependency management)::

    # Create a virtual environment
    uv venv

    # Install build dependencies first (required for --no-build-isolation)
    uv pip install --python .venv/bin/python scikit-build-core cython

    # Install in editable mode with incremental build support
    uv pip install -ve . --python .venv/bin/python --no-build-isolation

When you modify C++ source files, trigger an incremental rebuild::

    # With pip (incremental builds work automatically)
    pip install -ve .

    # With uv (must use --no-build-isolation for incremental builds)
    uv pip install -ve . --python .venv/bin/python --no-build-isolation

The build system (scikit-build-core with CMake/Ninja) will automatically detect which files
have changed and only recompile those files, making the rebuild much faster than a full rebuild.

**Incremental Build Performance**:

* Without ``--no-build-isolation``: Each build uses a fresh isolated environment, causing CMake
  to reconfigure and rebuild all ~60 source files (~50 seconds)
* With ``--no-build-isolation``: CMake reuses the existing build directory and only rebuilds
  changed files (~4 seconds for a single file change)

**Note**: If you use ``--no-build-isolation``, you must manually install the build dependencies
(``scikit-build-core`` and ``cython``) in your virtual environment first.

Local installation
------------------

If you prefer not to install system-wide, you can install locally using the ``--user`` flag::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp
    # Move into the folder you just created
    cd CoolProp
    # Install locally
    pip install --user .

For specific Python versions
-----------------------------

If you have multiple Python versions installed, use the specific Python's pip::

    # For Python 3.11 specifically
    python3.11 -m pip install .

Building with specific compilers
---------------------------------

You can control CMake options via environment variables::

    # Example: Specify a different C++ compiler
    export CXX=/usr/bin/g++-11
    pip install .

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
warnings filter and changing the behavior of the warnings module.  A legacy or third-party
caller can ignore *all* deprecation warnings by including the following code *after* the last
import from CoolProp, but *before* any calls to CoolProp::

    import warnings
    warnings.filterwarnings('ignore', category=DeprecationWarning)

The second parameter to ``filterwarnings()`` can be a regular expression matching only
selected messages, and the first parameter (``ignore``) can be set to ``once`` to issue a
given message only once and then ignore further instances.

.. note::

   The non-SI ``Props()`` and ``HAProps()`` functions — long deprecated — were **removed in
   CoolProp 8.0**.  Use the SI functions ``PropsSI()`` and ``HAPropsSI()`` instead, with all
   inputs and outputs in base SI units (Pa, K, J, kg).  For example the old
   ``Props('D', 'T', 298.15, 'P', 10000, 'R744')`` (pressure in kPa) becomes
   ``PropsSI('D', 'T', 298.15, 'P', 10e6, 'R744')`` (pressure in Pa).

See `Python >>> Module Warnings <https://docs.python.org/3/library/warnings.html#module-warnings>`_ for more information on using `filterwarnings()`

.. _python_wasm_demo:

Running CoolProp in the Browser Using Pyodide
=============================================

The `Pyodide project <https://pyodide.org/en/stable/>`_  
is a Python interpreter compiled to WebAssembly that includes many
of the scientific Python packages, including NumPy, SciPy, and CoolProp.
WebAssembly makes these packages available in the web browser environment, enabling the 
deployment of static websites that can perform Python calculations.

Minimal CoolProp Pyodide example
--------------------------------

The following example shows a minimal web page that calls a CoolProp function using 
Pyodide. See the `Pyodide documentation page <https://pyodide.org/en/stable/usage/index.html>`_ 
for further details on the usage of Pyodide.

.. code-block:: html

    <!doctype html>
    <html>
      <head>
          <!-- Load Pyodide from Pyodide's official CDN-->
          <script src="https://cdn.jsdelivr.net/pyodide/v314.0.2/full/pyodide.js"></script>
      </head>
      <body>
        <script type="text/javascript">
          async function main() {
            // Load Pyodide (the loadPyodide function is provided by the script tag above)
            let pyodide = await loadPyodide();
            
            // Load micropip so that PyPI wheels can be installed
            await pyodide.loadPackage("micropip");
            const micropip = pyodide.pyimport("micropip");

            // Use micropip to install the latest CoolProp
            await micropip.install("CoolProp>=8.0.1");
            
            // Run the Python code, the expression on the last line is returned to javascript
            results = pyodide.runPython(`
                import CoolProp
                import CoolProp.CoolProp as CP

                version = CoolProp.__version__
                temp = CP.PropsSI('T', 'P', 101325, 'Q', 0, 'Water')

                (version, temp)
            `);

            // Popup an alert with the results
            alert(`CoolProp.__version__=${results[0]}\nThe boiling point of water at 1 atm = ${results[1]} K`);
          }
          main();
        </script>
      </body>
    </html>

Getting the latest CoolProp version in Pyodide
----------------------------------------------

Historically, the Pyodide project has bundled CoolProp as part of the Pyodide distribution 
(beginning with Pyodide 0.24.0). This locked the CoolProp version to the version
bundled with Pyodide. Versions of Pyodide 0.28.x and later are able to load compiled Python 
packages from `PyPI <https://pypi.org/>`_ as part of the implementation of 
`PEP 783 <https://peps.python.org/pep-0783/>`_. Starting with 
CoolProp 8.0.1, CoolProp WebAssembly wheels are provided on PyPI. However, for versions 
of Pyodide that bundle CoolProp, the version of CoolProp needs to be specified as 
``>=8.0.1`` in the ``micropip`` call to get the latest CoolProp version since Pyodide 
will use its bundled version first, if available (CoolProp 7.2.0 with Pyodide 314.0.2, 
for example). Future versions of Pyodide will stop including CoolProp as part of the PEP 783 
rollout making it no longer necessary to specify the CoolProp version to get the latest version.

Module Documentation
====================

.. toctree::

    ../../../apidoc/CoolProp.rst
