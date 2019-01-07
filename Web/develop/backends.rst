.. _backends:

********************
Backends in CoolProp
********************

AbstractState
-------------
The :cpapi:`AbstractState` defines an interface between CoolProp and the rest of the world.  The public methods like :cpapi:`rhomolar<CoolProp::AbstractState::rhomolar>` are meant to be called by other code, while the protected functions like :cpapi:`AbstractState::calc_cpmass(void)` are meant to be implemented by the other backends.

Derived Backends
----------------

The backends in CoolProp provide the implementation of the protocol defined by the abstract base class.  There are a few primary backends:

* :cpapi:`HelmholtzEOSMixtureBackend` : This backend is the backend that provides properties using the CoolProp code.
* :cpapi:`IncompressibleBackend` : This backend provides the thermophysical properties for incompressible pure fluids, incompressible mixtures, and brines
* :cpapi:`REFPROPMixtureBackend` : This backend provides a clean interface between CoolProp and REFPROP
* :cpapi:`IF97Backend` : This backend provides an IAPWS-IF97 implemenation for the properties of steam/water

Example Backend
---------------

Code
----
.. literalinclude:: snippets/ExampleBackend.cxx
   :language: c++

Output
------
.. literalinclude:: snippets/ExampleBackend.cxx.output

