Changelog for CoolProp
======================

5.0.3
-----
Bugfix, with some new functionality

The most important fix is for users of Microsoft Excel on windows. It is imperative to download a new CoolProp.dll, there was a serious bug in how Excel and CoolProp interact that has been fixed.

Issues Closed:

* `#293 <http://github.com/CoolProp/CoolProp/issues/293>`_ : Requirement for zipped source code file
* `#292 <http://github.com/CoolProp/CoolProp/issues/292>`_ : Update CycloHexane EOS
* `#289 <http://github.com/CoolProp/CoolProp/issues/289>`_ : Two-phase states don't work for DY flash
* `#288 <http://github.com/CoolProp/CoolProp/issues/288>`_ : Some calls in Excel throw FPU exceptions which throw error messages
* `#287 <http://github.com/CoolProp/CoolProp/issues/287>`_ : Predefined mixtures cannot be used in PropsSI
* `#285 <http://github.com/CoolProp/CoolProp/issues/285>`_ : Cannot solve for conductivity and viscosity
* `#284 <http://github.com/CoolProp/CoolProp/issues/284>`_ : Create build steps on the master that allow us to automate the releasing even more
* `#283 <http://github.com/CoolProp/CoolProp/issues/283>`_ : Change fullclean logic to use git pull to wipe the folder completely
* `#282 <http://github.com/CoolProp/CoolProp/issues/282>`_ : SWIG wrappers not converting errors in PropsSI to exceptions
* `#280 <http://github.com/CoolProp/CoolProp/issues/280>`_ : Describe the predefined mixtures with examples on website

5.0.2
-----
Bugfix.

Issues Closed:

* `#281 <http://github.com/CoolProp/CoolProp/issues/281>`_ : Surface Tension Errors
* `#278 <http://github.com/CoolProp/CoolProp/issues/278>`_ : Add script to generate milestone text automatically
* `#277 <http://github.com/CoolProp/CoolProp/issues/277>`_ : Fix doxygen docs for generalized residual helmholtz term
* `#275 <http://github.com/CoolProp/CoolProp/issues/275>`_ : Logscale densities for consistency plots
* `#274 <http://github.com/CoolProp/CoolProp/issues/274>`_ : P and D as inputs produces some errors
* `#273 <http://github.com/CoolProp/CoolProp/issues/273>`_ : hmolar, smolar etc. are incorrect for HEOS backend with PD inputs
* `#272 <http://github.com/CoolProp/CoolProp/issues/272>`_ : 32bit Pre-compiled Binary for C#
* `#254 <http://github.com/CoolProp/CoolProp/issues/254>`_ : Error : hapropsSI("R";"T";253;"B";252;"P";101325) lead to an error

5.0.1
-----
The first release with the automated release script. No major code changes.

5.0.0
-----
**MAJOR** The new version of CoolProp implementing the new structure based on AbstractState
**MAJOR** Some features have been temporarily (or permanently) deprecated
**MAJOR** CoolProp now supports mixtures
**MAJOR** Buildbot system powered by CMake set up to run builds after every commit