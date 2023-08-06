.. _pcsaft_backend:

************************
PC-SAFT Equations of State
************************

.. contents:: :depth: 2

Introduction
============

CoolProp (as of version 6.4) includes the PC-SAFT equation of state. The PC-SAFT equation of state was originally proposed in 2001 by `Gross and Sadowski <>`. In addition to the hard chain and dispersion terms, the PC-SAFT backend in CoolProp also includes terms for associating, polar, and electrolyte compounds. For the polar term the formulation given by `Gross and Vrabec (2006) <>` was used, and this is sometimes called PCP-SAFT. For electrolyte compounds the equations presented by `Held et al. (2014) <>` were used, and this version of the equation of state is sometimes called electrolyte PC-SAFT (ePC-SAFT).

Caveats
-------

.. warning:: NOT ALL PROPERTIES ARE AVAILABLE AS INPUTS/OUTPUTS

Only a limited subset of properties are available currently. You can do:

* Flash calculations with TP, PQ, DT, QT inputs
* Calculation of some mixture flashes

.. warning:: The flash algorithm for the PC-SAFT backend is not yet as robust as for other backends. For some conditions it may fail to find the correct solution.

Pure Fluids
===========

Usage
-----

Similar to other backends in CoolProp, in the :ref:`high-level interface <high_level_api>`, all that you have to do to evaluate properties using the PC-SAFT equation of state is to change the backend name.

.. ipython::

    In [0]: import CoolProp.CoolProp as CP

    # The multi-parameter Helmholtz backend
    In [0]: CP.PropsSI("T","P",101325,"Q",0,"HEOS::Propane")

    # PC-SAFT
    In [0]: CP.PropsSI("T","P",101325,"Q",0,"PCSAFT::PROPANE")

The same holds for the :ref:`low-level interface <low_level_api>`:

.. ipython::

    In [0]: import CoolProp.CoolProp as CP

    In [0]: AS = CP.AbstractState("PCSAFT", "PROPANE"); AS.update(CP.QT_INPUTS, 0, 300); print(AS.p())

The PC-SAFT equation of state is available for more than 100 fluids for which parameter were available in the literature.

Mixtures
========

Interaction Parameters
----------------------

For mixtures, PC-SAFT generally uses a binary interaction parameter between pairs of fluids. CoolProp does have some of these parameters for the PC-SAFT EOS, and it is possible to add more yourself.

.. ipython::

    In [0]: import CoolProp.CoolProp as CP

    In [0]: CAS_water = CP.get_fluid_param_string("WATER","CAS")

    In [0]: CAS_aacid = "64-19-7"

    In [0]: CP.set_mixture_binary_pair_pcsaft(CAS_water, CAS_aacid, "kij", -0.127)

    In [0]: T = CP.PropsSI("T", "P", 72915.92217342, "Q", 0, "PCSAFT::WATER[0.2691800943]&ACETIC ACID[0.7308199057]")

    In [0]: print(T)
