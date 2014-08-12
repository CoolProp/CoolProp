CoolProp Fluid Information
==========================
.. toctree::
    :maxdepth: 1

    Fluids/pseudo_pure_fluids.rst
    Fluids/pure_fluids.rst
    
There are 4 basic classes of fluids that can be used in CoolProp, an example is provided for each one.

Pure Fluids, Pseudo-Pure Fluids
-------------------------------
All the fluids listed in Pure-Fluids and Pseudo-Pure-Fluids sections above can be used.  To use one of these fluids, do something like

.. ipython::

    In [1]: from CoolProp.CoolProp import PropsSI
    
    #Density of dry air at 1 atm. and 25C
    In [1]: PropsSI('D','T',298,'P',101325,'Air')
    
You can also use any of the aliases of the fluids that are the listed on the fluid page.  For instance, R717 is the refrigerant number for ammonia

.. ipython::

    In [1]: from CoolProp.CoolProp import PropsSI
    
    #Density of saturated ammonia vapor at 1 atm.
    In [1]: PropsSI('D','Q',1,'P',101325,'R717')
    
    #Density of saturated ammonia vapor at 1 atm.
    In [1]: PropsSI('D','Q',1,'P',101325,'Ammonia')
    

REFPROP Fluids and Mixtures
---------------------------
If you are on Windows and have REFPROP installed, you can use it with CoolProp.  REFPROP needs to be installed in c:\\Program Files\\REFPROP.  If it is somewhere else, just copy it to this location.

It is also possible to use REFPROP on Linux.  Please follow the instructions from https://github.com/jowr/librefprop.so to install the library from Fortran sources.  Additionally, you also need to copy the fluid and mixture files to /opt/refprop. 

All the pure fluids in REFPROP are used just like the CoolProp fluids except that "REFPROP-" is added at the beginning of the fluid name.  You can use any fluid that is included in REFPROP, but you must use the REFPROP fluid file name.  For CoolProp Fluids, you can use the ``get_REFPROPName()`` function to get the REFPROP name for the fluid.

.. ipython::

    In [1]: from CoolProp.CoolProp import PropsSI
    
    #Saturated isobutane vapor density at 1 atmosphere
    In [1]: PropsSI('D','Q',1,'P',101.325,'REFPROP::ISOBUTAN')
    
You can also use mixtures in REFPROP, there is a special format for the fluid name.  The fluid name is set up like this: ``"REFPROP::R32[0.697615]&R125[0.302385]"`` -  this is R410A.  The numbers within the brackets are the mole fractions of the components.  They must add up to 1.0
    
.. ipython::

    In [1]: from CoolProp.CoolProp import PropsSI
    
    #Saturated R410 vapor density at 1 atmosphere using the mixture properties
    In [1]: PropsSI('D','Q',1,'P',101.325,'REFPROP::R32[0.697615]&R125[0.302385]')
    
In Engineering Equation Solver (EES), you can compose a fluid string using something like ``a$='REFPROP::'||'R32'||'['||yy$||']'||and$||'R125'||'['||xx$||']'`` where ``yy`` and ``xx`` are mole fractions of R32 and R125 respectively.
   
