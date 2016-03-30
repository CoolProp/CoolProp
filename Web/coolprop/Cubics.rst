.. _cubic_backend:

************************
Cubic Equations of State
************************

.. contents:: :depth: 2

Introduction
------------

CoolProp (as of version 6) comes with two standard cubic equations of state: Soave-Redlich-Kwong (SRK) and Peng-Robinson (PR).

Caveats
-------

.. warning:: NOT ALL PROPERTIES ARE AVAILABLE AS INPUTS/OUTPUTS

Only a limited subset of properties are available currently. You can do:

* Flash calculations with TP, PQ, DT, QT inputs
* Calculation of mixture critical point(s)
* Calculation of some mixture flashes

Usage
-----

As a user, in the :ref:`high-level interface <high_level_api>`, all that you have to do to evaluate properties using the SRK or PR backends is to change the backend name.

.. ipython::

    In [0]: import CoolProp.CoolProp as CP

    # The multi-parameter Helmholtz backend
    In [0]: CP.PropsSI("T","P",101325,"Q",0,"HEOS::Propane")

    # SRK
    In [0]: CP.PropsSI("T","P",101325,"Q",0,"SRK::Propane")

    # Peng-Robinson
    In [0]: CP.PropsSI("T","P",101325,"Q",0,"PR::Propane")

The same holds for the :ref:`low-level interface <low_level_api>`:

.. ipython::

    In [0]: import CoolProp.CoolProp as CP

    In [0]: AS = CP.AbstractState("SRK", "Propane"); AS.update(CP.QT_INPUTS, 0, 300); print(AS.p())

    In [0]: AS = CP.AbstractState("PR", "Propane"); AS.update(CP.QT_INPUTS, 0, 300); print(AS.p())

All the fluids available in CoolProp are also available through the cubic equations of state.  The fluids can be extended according to the analysis shown below

Adding your own fluids
----------------------

The cubic fluids in CoolProp are defined based on a JSON format, which could yield something like this for a **FAKE** (illustrative) fluid.  For instance if we had the file ``fake_fluid.json`` (download it here: :download:`fake_fluid.json`):

.. literalinclude:: fake_fluid.json

The JSON-formatted fluid information is validated against the `JSON schema <http://json-schema.org>`_, which can be obtained by calling the ``get_global_param_string`` function:

.. ipython:: 

    In [0]: import CoolProp.CoolProp as CP

    # Just the first bit of the schema (it's a bit large; see below for the complete schema)
    In [0]: CP.get_global_param_string("cubic_fluids_schema")[0:60]

A `JSON schema <http://json-schema.org>`_ defines the structure of the data, and places limits on the values, defines what values must be provided, etc..  There are a number of online tools that can be used to experiment with and validate JSON schemas: http://www.jsonschemavalidator.net , http://jsonschemalint.com/draft4/.

.. ipython:: 

    In [0]: import CoolProp.CoolProp as CP, json

    # Normally you would read it from a file!
    In [0]: fake_fluids = [
       ...:                   { 
       ...:                     "CAS": "000-0-00",
       ...:                     "Tc": 400.0,
       ...:                     "Tc_units": "K",
       ...:                     "acentric": 0.1,
       ...:                     "aliases": [
       ...:                     ],
       ...:                     "molemass": 0.04,
       ...:                     "molemass_units": "kg/mol", 
       ...:                     "name": "FAKEFLUID",
       ...:                     "pc": 4000000.0,
       ...:                     "pc_units": "Pa"
       ...:                   }
       ...:               ]

    # Adds fake fluid for both SRK and PR backends since they need the same parameters from the same internal library
    In [0]: CP.add_fluids_as_JSON("PR", json.dumps(fake_fluids))

    In [0]: CP.PropsSI("T","P",101325,"Q",0,"SRK::FAKEFLUID")

Complete fluid schema
---------------------

For completeness, here is the whole `JSON schema <http://json-schema.org>`_ used for the cubic backends:

.. ipython:: 

    In [0]: import CoolProp.CoolProp as CP

    In [0]: print(CP.get_global_param_string("cubic_fluids_schema"))