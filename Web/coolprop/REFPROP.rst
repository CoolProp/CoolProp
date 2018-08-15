.. _REFPROP:

*******************
REFPROP Interface
*******************

.. contents:: :depth: 2

The thermophysical property library REFPROP developed by researchers at the National Institute of Standards and Technology in Boulder, Colorado is the gold standard in thermophysical properties.  It is mature and stable, and is currently the library most used in industry and for scientific applications.

CoolProp allows for full interaction with the REFPROP library, while using the nicely abstracted C++ interface for all target languages that CoolProp supports. 

The core difference between using CoolProp's internal routines and REFPROP's routines is that you have to change the backend from ``HEOS`` to ``REFPROP``.  In the high-level interface this is demonstrated as:

.. ipython::

    In [0]: from CoolProp.CoolProp import PropsSI
    
    # This one uses CoolProp's internal routines to calculate NBP of water
    In [0]: PropsSI("T","P",101325,"Q",0,"HEOS::Water")

    # This one uses REFPROP's internal routines to calculate NBP of water
    In [0]: PropsSI("T","P",101325,"Q",0,"REFPROP::Water")
    
The list of input keys that can be used in the PropsSI function are given here: :ref:`parameter_table`. 

In the :ref:`low-level interface <low_level_api>` (significantly faster once you have constructed the class!), this would be

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

Other flash routines (see :cpapi:`CoolProp::input_pairs`) proceed in exactly the same fashion.  Most of the methods available for CoolProp are also available for REFPROP.  For instance, you can generate phase envelopes with both CoolProp and REFPROP, or evaluate the dewpoint of mixtures.  Here are a few examples of the mapping between CoolProp and REFPROP for some more interesting applications:

.. ipython::

    In [0]: import CoolProp

    In [0]: HEOS = CoolProp.AbstractState('HEOS','R32&R125'); HEOS.set_mole_fractions([0.5,0.5])

    In [0]: REFPROP = CoolProp.AbstractState('REFPROP','R32&R125'); REFPROP.set_mole_fractions([0.5,0.5])
    
    # CoolProp calculates the thermodynamically correct temperature
    In [0]: HEOS.T_critical()

    # The default in REFPROP is to interpolate the phase envelope
    In [0]: REFPROP.T_critical()

    # Here we calculate and obtain the phase envelope for the mixture using the routines in CoolProp
    In [0]: HEOS.build_phase_envelope(""); PE_HEOS = HEOS.get_phase_envelope_data()

    # Here we calculate and obtain the phase envelope for the mixture using the routines in REFPROP
    In [0]: REFPROP.build_phase_envelope(""); PE_REFPROP = REFPROP.get_phase_envelope_data()    

Path Issues
-----------

.. warning::

    In order for REFPROP to be able to be loaded by CoolProp, the default logic for each operating system is used to load the REFPROP shared library.  This means that on windows, the ``PATH`` environmental variable is searched for the ``REFPROP.dll`` (32-bit applications) or ``REFPRP64.dll`` (64-bit applications). On linux/OSX, the default shared library loading protocol is used.  If your REFPROP is installed in a non-standard location (not on the path), make sure that when you run code that uses REFPROP, that you add (temporarily) the location of the REFPROP shared library to your path or set one of the following configuration variables.

REFPROP needs to be able to find the fluid and mixture files at runtime, at a location specified on your computer.  CoolProp allows you to avoid the pains of decoding REFPROP's internal logic for finding these files by explicitly specifying the path that it should tell REFPROP to look for the fluid files.  

.. warning::

    These configuration variables should be set at the beginning of your script and then not touched again.  Otherwise, you can get some weird behavior!

The configuration key for setting the REFPROP path (see :ref:`configuration`) is ``ALTERNATIVE_REFPROP_PATH``, and you can set it doing something like this in python:

.. ipython::

    In [0]: import json, CoolProp.CoolProp as CP
    
    In [1]: CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'c:\\Program Files\\REFPROP\\')

If you do this, internally CoolProp will construct the full path to the ``fluids`` and ``mixtures`` directories and use them to call REFPROP's `SETUP` function when loading fluids.  If you don't do this, CoolProp will use whatever default logic REFPROP uses to find the fluid files.

If you wish to use a certain shared library, for example to try different REFPROP versions, you explicitly define it via ``ALTERNATIVE_REFPROP_LIBRARY_PATH``. This configuration variable makes CoolProp ignore the ``ALTERNATIVE_REFPROP_PATH`` when loading the shared library and you might have to provide the full path to your shared library here by doing something like this in python:

.. ipython::

    In [0]: import json, CoolProp.CoolProp as CP
    
    In [1]: CP.set_config_string(CP.ALTERNATIVE_REFPROP_LIBRARY_PATH, 'c:\\Program Files\\REFPROP\\REFPRP64.v9.1.dll')

.. warning::
    
    If you use a combination of ``ALTERNATIVE_REFPROP_LIBRARY_PATH`` and ``ALTERNATIVE_REFPROP_PATH``, the shared library gets loaded directly from ``ALTERNATIVE_REFPROP_LIBRARY_PATH`` while the fluid files still will be accessed via ``ALTERNATIVE_REFPROP_PATH``. You can thus have one single folder with fluid files that are used with different shared libraries. Make sure that the fluid files are compatible with all the shared library versions you are using. 

If you are playing around with mixture parameters, you might want to set a different path to the HMX.BNC file which contains the interaction parameters for the mixture.  You can do that by changing the configuration variable  (see :ref:`configuration`) ``ALTERNATIVE_REFPROP_HMX_BNC_PATH``

.. ipython::

    In [0]: import json, CoolProp.CoolProp as CP
    
    In [1]: CP.set_config_string(CP.ALTERNATIVE_REFPROP_HMX_BNC_PATH, 'c:\\Program Files\\REFPROP\\fluids\\HMX.BNC')

If you have set both the ``ALTERNATIVE_REFPROP_PATH`` and ``ALTERNATIVE_REFPROP_HMX_BNC_PATH`` variables, ``ALTERNATIVE_REFPROP_PATH_HMX_BNC_PATH`` "wins", and this path will be used when loading mixture interaction parameters

And now we set them back to their default values

.. ipython::

    In [0]: import json, CoolProp.CoolProp as CP
    
    In [1]: CP.set_config_string(CP.ALTERNATIVE_REFPROP_HMX_BNC_PATH, '')
    
    In [1]: CP.set_config_string(CP.ALTERNATIVE_REFPROP_LIBRARY_PATH, '')

    In [1]: CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, '')

Other Platforms
---------------

On linux and OSX, you can build your own copy of REFPROP shared library using the instructions here: https://github.com/usnistgov/REFPROP-cmake

On linux, here are instructions for adding your shared library to the ``LD_LIBRARY_PATH`` variable: http://stackoverflow.com/a/13428971/1360263

Other Features
--------------

If you want to determine the version of REFPROP that you are actually using, you can do:

.. ipython::

    In [0]: import CoolProp.CoolProp as CP
    
    In [1]: CP.get_global_param_string("REFPROP_version")


If you want to use the GERG-2008 model, you can do this at the beginning of your code:

.. ipython::

    In [0]: import CoolProp.CoolProp as CP
    
    In [1]: CP.set_config_bool(CP.REFPROP_USE_GERG, True)

Subsequently, all calculations will be done with the simplified EOS from the GERG-2008 model

And now we set them back to their default values

.. ipython::

    In [0]: import json, CoolProp.CoolProp as CP
    
    In [1]: CP.set_config_bool(CP.REFPROP_USE_GERG, False)
