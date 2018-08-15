.. _configuration:

***********************
Configuration Variables
***********************

At runtime, there are a several configuration variables that can be used to change the behavior of CoolProp

.. warning:: The adjustment of the internal configuration variables might have side effects that you are not expecting, use with caution!!

From C++ and the SWIG wrappers, the values can be directly set/changed by using the type-specified getter/setter functions :cpapi:`get_config_bool`, :cpapi:`set_config_bool`, :cpapi:`get_config_string`, :cpapi:`set_config_string`, etc., as in something like:

.. ipython::

    In [0]: import CoolProp.CoolProp as CP

    In [1]: current_val = CP.get_config_bool(CP.CRITICAL_WITHIN_1UK)

    In [1]: CP.set_config_bool(CP.CRITICAL_WITHIN_1UK, current_val)

From all languages, the configuration state can obtained by retrieving the configuration state in the form of a  `JSON <http://json.org/>`_ formatted string.  For instance, in python, you can get the default configuration state from 

.. ipython::

    In [0]: import CoolProp.CoolProp as CP, json

    In [1]: json.loads(CP.get_config_as_json_string())
    
Most modern languages have facilities for interfacing with JSON formatted strings and converting them back and forth with language-specific data structures.  For instance, in python, there is the built-in ``json`` package that converts json-formatted strings to python dictionaries, lists, etc.


Configuration Keys
------------------

Here is a list of (alphabetically sorted) configuration keys that can be used, along with a short description of the use of each of the configuration keys:

.. include :: configuration_keys.rst.in