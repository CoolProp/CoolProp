Lua Wrapper for CoolProp
========================

Lua C API and LuaJIT wrappers for CoolProp library (currently implementing the high-level C API).

by Aapo Talvensaari, Helsinki, Finland, January 2015


Hello World with Lua Wrapper
----------------------------

::

  local cp = require "coolprop"
  print(cp.PropsSI("d(Hmass)/d(T)|P", "P", 101325, "T", 300, "Water"))


Installation
------------

This wrapper comes with a support for both LuaJIT FFI and PUC-Lua C-API.

To build a shared library of CoolProp, follow the instructions listed here:
http://www.coolprop.org/coolprop/wrappers/SharedLibrary/index.html


**Installing with Make**

1. Run ``make``
2. Run ``make install``
3. Done.

**Manual Lua (PUC) Installation**

1. First you need to build shared library version of a CoolLib library and install it somewhere in your operating system's library search path.
2. Next you need to compile ``coolprop/cpapi.c`` wrapper. Easiest way to do it is to run ``make`` on Lua wrapper's root folder.
3. Then place the resulting ``coolprop/cpapi.so`` in your Lua's ``package.cpath``, and there under ``coolprop`` directory.
4. Done.

**Manual LuaJIT Installation**

1. First you need to build shared library version of a CoolLib library and install it somewhere in your operating system's library search path.
2. Next you need to place ``coolprop.lua`` and ``coolprop/ffi.lua`` somewhere in your LuaJIT's ``package.path`` so that the LuaJIT can pick them up.
3. Done.


Lua API
-------

The functions will in general return ``nil`` on error.

**number, string coolprop.PropsSI(output, name1, prop1, name2, prop2, ref)**

For many users, all that is needed is a simple call to the ``PropsSI`` function for pure fluids, pseudo-pure fluids
and mixtures. This function will return ``nil`` on error, and additional string that contains the error message.

::

  local cp = require "coolprop"
  print(cp.PropsSI("C", "P", 101325, "T", 300, "Water"))
  print(cp.PropsSI("d(Hmass)/d(T)|P", "P", 101325, "T", 300, "Water"))
  print(cp.PropsSI("D", "P", 101325, "T", 300, "Air.mix"))


**number, string coolprop.Props1SI(fluidname, output)**

::

  local cp = require "coolprop"
  print(cp.Props1SI("Water", "Phase"))


**string coolprop.PhaseSI(name1, prop1, name2, prop2, ref)**

::

  local cp = require "coolprop"
  print(cp.PhaseSI("P", 101325, "Q", 0, "Water"))


**number, string coolprop.HAPropsSI(output, name1, prop1, name2, prop2, name3, prop3)**

::

  local cp = require "coolprop"
  print(cp.HAPropsSI('H','T',298.15,'P',101325,'R',0.5))


**string coolprop.get_global_param_string(param)**

Returns global parameter string. On error, returns ``nil``.

::

  local cp = require "coolprop"
  print(cp.get_global_param_string("predefined_mixtures"))


**number coolprop.get_param_index(param)**

Returns parameter index.

::

  local cp = require "coolprop"
  print(cp.get_param_index("T"))


**number coolprop.F2K(f)**

This function converts fahrenheits to kelvins.

::

  local cp = require "coolprop"
  print(10, "Fahrenheits is", cp.F2K(10) , "Kelvins")


**number coolprop.K2F(f)**

This function converts kelvins to fahrenheits.

::

  local cp = require "coolprop"
  print(cp.F2K(10), "Kelvins is", cp.K2F(cp.F2K(10)), "Fahrenheits")


**string coolprop.error()**

Returns the last error occurred.

::

  local cp = require "coolprop"
  print(cp.error())


**string coolprop.FluidsList()**

Returns the list of available fluids.

::

  local cp = require "coolprop"
  print(cp.FluidsList())


**string coolprop.version()**

Returns the version of the CoolLib library that is installed.

::

  local cp = require "coolprop"
  print(cp.version())


**string coolprop.gitrevision()**

Returns the Git revision of the CoolLib library that is installed.

::

  local cp = require "coolprop"
  print(cp.gitrevision())


**number coolprop.get_debug_level()**

Returns the current debug level.

::

  local cp = require "coolprop"
  print(cp.get_debug_level())


**coolprop.set_debug_level(level)**

Sets the debug level.

::

  local cp = require "coolprop"
  cp.set_debug_level(0)


**boolean coolprop.redirect_stdout(file)**

Sets the output to a file (to given path of the file).

::

  local cp = require "coolprop"
  cp.redirect_stdout("output.log")


Additional APIs (TBD)
---------------------

- **string coolprop.get_parameter_information_string(key)**
- **number coolprop.get_mixture_binary_pair_data(cas1, cas2, key)**
- **string coolprop.get_fluid_param_string(fluid, param)**
- **boolean coolprop.set_reference_stateS(ref, state)**
- **boolean coolprop.set_reference_stateD(ref, t, rho, h0, s0)**
- **number, string coolprop.saturation_ancillary(fluid, output, q, input, value)**
