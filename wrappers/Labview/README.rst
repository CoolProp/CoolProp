Labview wrapper of CoolProp
============================

by Ian Bell, Arnaud Legros, Jan Nolens and Sergei Gusev

University of Liege, Ghent University

October 2013

Available libraries
-------------------
CoolProp.vi: 
Basic Library to get the properties from CoolProp.dll

CoolProp.llb
More advanced library, allowing to compute thermodynamic diagrams, real-time calculation of properties,
measurement processing, etc.

See the Diagrams_coolprop.vi in the LLB to get started.  To use, you must update the path to the CoolProp.dll file in the VI.

To Install
----------
1. Copy the files CoolProp.vi, Coolprop.llb and CoolProp.dll from this folder to somewhere you want
2. Add CoolProp module to your code

To Use
------
Call it using the same sorts of input parameters as Props function : http://coolprop.sourceforge.net/apidoc/CoolProp.html#CoolProp.CoolProp.Props

Notes
-----
Wrapper currently in its infancy, absolutely no error checking is carried out, use at your own risk.

For Developers
--------------

1. To regenerate DLL, run the build script (wrappers/Labview/BuildDLL.bat).  You will need to have Visual Studio 2010 installed, or change the path to vcvarsall.bat in build script
2. Uses __cdecl in combination with extern "C" to make Labview happy
3. For realtime targets, instructions are here: http://www.ni.com/white-paper/5694/en
