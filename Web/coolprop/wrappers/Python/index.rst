.. _Python:

**************
Python Wrapper
**************

Automatic installation
======================

Using the ``pip`` installation program, you can install CoolProp v5 from the development server using::

    pip install --find-links https://www.coolprop.dreamhosters.com:8010/binaries/Python/ -U --force-reinstall CoolProp
    
Or the official release can be obtained from the pypi server using::

    pip install CoolProp

Manual installation
===================

Special MacOS X requirement
---------------------------

If you never used any command-line installation before, chances are that you 
do not have the compiling utilities needed. Thus you need to first install 
Xcode: see the description on the page http://guide.macports.org/#installing.xcode
After installing, you need to accept the licence by running the following 
command in the Terminal::

   xcodebuild -license
   
and explicitly typing "agree" before closing. Then you can resume to the 
former method.

Windows
-------

To be explained...

Linux and MacOS X
-----------------

If not already there, you need Cython to be installed::

    sudo easy_install Cython

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
also install it locally provided that the `PYTHONPATH` variable is correctly 
set to let Python look into the given directory (here `HOME/.local`). Don't 
worry if it's not the case, the error message should help you to properly set 
your variable::

    python setup.py install --prefix $HOME/.local
    
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

For VS2010:
-----------

http://stackoverflow.com/questions/8044385/64-bit-build-on-microsoft-visual-c-express-2010/8334985#8334985

For VS2008
----------

For MinGW
---------

Create a file called 


Usage
=====

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


