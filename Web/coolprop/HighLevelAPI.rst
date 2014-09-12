.. _high_level_api:

********************
High-Level Interface
********************

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
