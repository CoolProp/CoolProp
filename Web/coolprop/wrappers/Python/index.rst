.. _Python:

**************
Python Wrapper
**************

.. contents:: :depth: 2

Automatic installation
======================

Using the ``pip`` installation program, you can install the official release from the pypi server using::

    pip install CoolProp

If you dare, you can also try the latest nightly release from :sfnightly:`Python` 
or get it directly from the development server using::

    pip install -vvv --pre --trusted-host www.coolprop.dreamhosters.com --find-links http://www.coolprop.dreamhosters.com/binaries/Python/ -U --force-reinstall CoolProp
    
.. For those of you who prefer the Anaconda or Miniconda distributions, you can run::
   
       conda install -c https://conda.binstar.org/coolprop coolprop
       
   You can also find our nightly development snapshots on binstar::
    
       conda install -c https://conda.binstar.org/coolprop/label/dev coolprop 


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
    ph_plot = CPP.PropsPlot('Water','Ph')
    ph_plot.savefig('enthalpy_pressure_graph_for_Water.png')

* Solve `thermodynamics exercices`_

* Make you own `more complex graphs`_ if you feel the graphing interface is lacking something

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


Module Documentation
====================

.. toctree::

    ../../../apidoc/CoolProp.rst

