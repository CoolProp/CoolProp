.. _REFPROP:

*******************
REFPROP Interface
*******************

.. contents:: :depth: 2

The thermophysical library REFPROP developed by researchers at the National Institute of Standards and Technology in Boulder, Colorado is the gold standard in thermophysical properties.  It is mature and stable, and is currently the library most used in industry and for scientific applications.

Unfortunately, the REFPROP interface is not very user friendly, unless you use FORTRAN.  It is here where CoolProp shines, and CoolProp allows for full interaction with the REFPROP library, while using the nicely abstracted C++ interface for all target languages that CoolProp supports.

The core difference between using CoolProp's internal routines and REFPROP's routines is that you have to change the backend from ``HEOS`` to ``REFPROP``.  In the high-level interface this is demonstrated as:

.. ipython::

    In [0]: from CoolProp.CoolProp import PropsSI
    
    # This one uses CoolProp's internal routines to calculate NBP of water
    In [0]: PropsSI("T","P",101325,"Q",0,"HEOS::Water")

    # This one uses REFPROP's internal routines to calculate NBP of water
    In [0]: PropsSI("T","P",101325,"Q",0,"REFPROP::Water")

or similarly, in the low-level interface, this would be

.. ipython::

    In [0]: import CoolProp

    In [0]: HEOS = CoolProp.AbstractState('HEOS','Water')

    In [0]: REFPROP = CoolProp.AbstractState('REFPROP','Water')    
    
    # This one uses CoolProp's internal routines
    In [0]: HEOS.update(CoolProp.PQ_INPUTS, 101325, 0)

    In [0]: HEOS.T()

    # This one uses REFPROP's internal routines
    In [0]: REFPROP.update(CoolProp.PQ_INPUTS, 101325, 0)

    In [0]: REFPROP.T()

Other flash routines proceed in exactly the same fashion.  Most of the methods available for CoolProp are also available for REFPROP.  For instance, you can generate phase envelopes with both CoolProp and REFPROP, or evaluate the dewpoint of mixtures.  Here are a few examples of the mapping between CoolProp and REFPROP for some more interesting applications:

.. ipython::

    In [0]: import CoolProp

    In [0]: HEOS = CoolProp.AbstractState('HEOS','R32&R2125'); HEOS.set_mole_fractions([0.5,0.5])

    In [0]: REFPROP = CoolProp.AbstractState('REFPROP','R32&R125'); REFPROP.set_mole_fractions([0.5,0.5])
    
    # Here we calculate the critical temperature of this mixture using the routines in both libraries.  The default in REFPROP is to interpolate the phase envelope while CoolProp calculates the thermodynamically correct temperature
    In [0]: HEOS.T_critical()

    In [0]: REFPROP.T_critical()

    # Here we calculate the phase envelope for the mixture using the routines in both libraries.
    In [0]: HEOS.build_phase_envelope(""); PE_HEOS = HEOS.get_phase_envelope_data()

    In [0]: REFPROP.build_phase_envelope(""); PE_REFPROP = REFPROP.get_phase_envelope_data()    

Path Issues
-----------

.. warning::

    In order for REFPROP to be able to be loaded by CoolProp, the default logic for each operating system is used to load the REFPROP shared library.  This means that on windows, the ``PATH`` environmental variable is searched for the ``REFPROP.dll`` (32-bit applications) or ``REFPRP64.dll`` (64-bit applications). On linux/OSX, the default shared library loading protocol is used.  If your REFPROP is installed in a non-standard location (not on the path), make sure that when you run code that uses REFPROP, that you add (temporarily) the location of the REFPROP shared library to your path.

One of the more challenging things to deal with for REFPROP is that REFPROP needs to be able to find the fluid and mixture files at runtime, at a location specified on your computer.  The default logic for where REFPROP looks is quite arcane, and not well documented.  CoolProp allows you to avoid the pains of decoding REFPROP's internal logic for finding these files by explicitly specifying the path that it should tell REFPROP to look for the fluid files.  The configuration key is ``ALTERNATIVE_REFPROP_PATH``, and you can set it doing something like this in python:

.. ipython::

    In [0]: import json, CoolProp.CoolProp as CP

    In [1]: jj = json.loads(CP.get_config_as_json_string())
    
    In [2]: jj['ALTERNATIVE_REFPROP_PATH'] = 'c:\\Program Files\\REFPROP'
    
    In [3]: jj = CP.set_config_as_json_string(json.dumps(jj))

If you do this, internally CoolProp will call the ``SETPATH`` function in REFPROP to tell REFPROP that it should find the ``fluids`` and ``mixtures`` within this directory.  If you don't do this, CoolProp will use whatever default logic REFPROP uses to find the fluid files.

If you are playing around with mixture parameters, you might want to set a different path to the HMX.BNC file which contains the interaction parameters for the mixture.  You can do that by changing the configuration variable ``ALTERNATIVE_REFPROP_HMX_BNC_PATH``

.. ipython::

    In [0]: import json, CoolProp.CoolProp as CP

    In [1]: jj = json.loads(CP.get_config_as_json_string())
    
    In [2]: jj['ALTERNATIVE_REFPROP_HMX_BNC_PATH'] = 'c:\\Program Files\\REFPROP\\fluids\\HMX.BNC'
    
    In [3]: jj = CP.set_config_as_json_string(json.dumps(jj))

If you have set both the ``ALTERNATIVE_REFPROP_PATH`` and ``ALTERNATIVE_REFPROP_HMX_BNC_PATH`` variables, ``ALTERNATIVE_REFPROP_PATH`` "wins", and this path will be used when loading mixture interaction parameters

Advanced Usage
--------------

