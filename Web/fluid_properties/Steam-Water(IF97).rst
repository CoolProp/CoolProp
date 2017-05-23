

.. _Steam-Water(IF97):

IAPWS-IF97 Steam/Water Properties
=================================


General Introduction
--------------------

Steam/Water properties are of critical importance in industry and specifically the power industry.  The standard for water properties is maintained by the International Association for the Properties of Water and Steam (IAPWS_).  The most recent formulation for thermodynamic properties is the IAPWS-95 Helmholtz formulation [R6-95(2016)]_, which is exactly what the CoolProp HEOS backend uses when specifying the fluid as ``Water`` or ``HEOS::Water``.  The equations of state are based on :math:`T` and :math:`\rho` as state variables, so :math:`T, \rho` will always be the fastest inputs.  :math:`p,T` will be a bit slower (3-10 times), and then inputs where neither :math:`T` nor :math:`\rho` are given, like :math:`p,h`.  They will be much slower.

Realizing the need for fast numerical computation of water properties, the IAPWS created a committee to develop a faster formulation for the steam power industry.  The IAPWS release the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam" or IAPWS-IF97.  This formulation replaces the IFC-67 formulation, is pressure based (rather than density based), and uses fundamental equations for the Gibbs free energy and its derivatives.  While IAPWS-IF97 is faster, it is less accurate than the IAPWS-95 formulation, but within 1% uncertainty for all properties except very near the critical point. Details of this formulation can be found in [R7-97(2012)]_, released in 2007 and updated in 2012.

The IF97 backend in CoolProp is based on the latest `IAPWS Releases`_ for thermodynamic and transport properties.  It was originally coded as a fast water property calculation for the HumidAir_ backend.  It has matured to include most of the forward and backward formulations in all phase regions as released by the IAPWS_ and can provide a much faster calculation of water properties than the default HEOS backend (or the REFPROP_ backend, which is also based on IAPWS-95_).  

.. _IAPWS: http://www.iapws.org 
.. _`IAPWS Releases`: http://http://www.iapws.org/release.html  
.. [R6-95(2016)] `Revised Release on the IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use`_
.. [R7-97(2012)] `Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam`: http://www.iapws.org/relguide/IF97-Rev.pdf
.. _`Revised Release on the IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use`: http://http://www.iapws.org/relguide/IAPWS-95.html
.. _`Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam: http://www.iapws.org/relguide/IF97-Rev.html

IF97 Backend
------------

The IF97 backend can be accessed the same way as the other backends, by prefixing the fluid name in calls to ``PropsSI`` and ``Props1SI``.  However, Water is the only fluid name that can be specified, as ``IF97::Water``.  Most output variables and input pairs are accepted as of CoolProp 6.2, with a few exceptions:

1. Generic Partial Derivatives are not yet supported.
2. Mixtures are obviously not supported since this is a Water-only backend.
3. The Reference State cannot be modified for the IF97 formulation.  *The reference state is fixed such that the specific internal energy, :math:`U`, and the specific entropy, :math:`S`, of the saturated liquid at the triple point are set equal to zero.*


IF97 Water Property Examples
----------------------------

A call to the top-level function ``PropsSI`` can provide: temperature, pressure, density, specific heat, internal energy, enthalpy, entropy, speed of sound, viscosity, thermal conductivity, surface tension, and Prandtl Number. Hence, the available output keys are: ``T``, ``P``, ``D``, ``C``, ``Cvmass``, ``U``, ``H``, ``S``, ``A``, ``V``, ``L``, ``I``, and ``Prandtl``.  Molar quantitities can also be returned, but IF97 is mass based.  Trivial outputs, such as ``M``, ``Tmin``, ``Tmax``, ``Pmin``, ``Pmax``, ``Ttriple``, ``Tcrit``, ``ptriple``, and ``pcrit``, are also available.

.. ipython::

    In [1]: from CoolProp.CoolProp import PropsSI
    
    #Specific heat capacity of Water at 500 K and 1 atm
    In [2]: PropsSI('C','T',500,'P',101325,'IF97::Water')

    #Density of Water at 500 K and 1 atm.
    In [3]: PropsSI('D','T',500,'P',101325,'IF97::Water')
    
    #Round trip in thermodynamic properties
    In [4]: T_init = 500.0
    
    In [5]: P_init = 101325
    
    In [6]: D_init = PropsSI('D','T',T_init,'P',P_init,'IF97::Water')
    
    In [7]: S_init = PropsSI('S','D',D_init,'P',P_init,'IF97::Water')
    
    In [8]: H_init = PropsSI('H','S',S_init,'P',P_init,'IF97::Water')
    
    In [9]: T_init = PropsSI('T','H',H_init,'P',P_init,'IF97::Water')
    
    In [10]: T_init

    #Saturation pressure of Water at 500 K
    In [11]: PropsSI('P','T',500,'Q',0,'IF97::Water')

    #Critical temperature for Water
    In [12]: PropsSI('Tcrit','T',0,'P',0,'IF97::Water')

    #Triple Point pressure for Water
    In [13]: PropsSI('ptriple','T',0,'P',0,'IF97::Water')



