**************
High-Level API
**************

For many users, all that is needed is a simple call to the ``PropsSI`` function.

Code
----
.. literalinclude:: snippets/propssi.cxx
   :language: c++

Output
------
.. literalinclude:: snippets/propssi.cxx.output


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
    
    In [1]: import timeit
    
    # Specific heat (kJ/kg/K) of 20% ethylene glycol as a function of T
    In [2]: PropsSI('C','T',298.15,'P',101.325,'INCOMP::EG-20%')
    
    # Density of Air at standard atmosphere in kg/m^3
    In [2]: PropsSI('D','T',298.15,'P',101.325,'Air')
    
    # Saturation temperature of Water at 1 atm
    In [2]: PropsSI('T','P',101.325,'Q',0,'Water')
    
    # Saturated vapor density of R134a at 0C
    In [2]: PropsSI('H','T',273.15,'Q',1,'R134a')
    
    # Using properties from REFPROP to get R410A density
    In [2]: PropsSI('D','T',300,'P',100,'REFPROP::R32[0.697615]&R125[0.302385]')
    
    # Check that the same as using pseudo-pure
    In [2]: PropsSI('D','T',300,'P',100,'R410A')