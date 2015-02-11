Changelog for CoolProp
======================


5.0.7
-----

New Features:

* Added a Lua wrapper

Issues Closed:

* `#460 <http://github.com/CoolProp/CoolProp/issues/460>`_ : PropsSI ("Q", "P", valueP, "H", valueH, "REFPROP-R410A") only return 0
* `#459 <http://github.com/CoolProp/CoolProp/issues/459>`_ : PropsSI ("D", "P", valueP, "T", valueT, "R407C") return bad result in L+V Phasis
* `#456 <http://github.com/CoolProp/CoolProp/issues/456>`_ : Slave alert
* `#454 <http://github.com/CoolProp/CoolProp/issues/454>`_ : Add density dependency to entropy and enthalpy of incomprerssible fluids
* `#452 <http://github.com/CoolProp/CoolProp/issues/452>`_ : Allow mixtures to have zero mole fractions
* `#450 <http://github.com/CoolProp/CoolProp/issues/450>`_ : Calling PropsSI to get thermal conductivity throws an exception
* `#448 <http://github.com/CoolProp/CoolProp/issues/448>`_ : Retrieving acentric factor through Props1SI fails
* `#443 <http://github.com/CoolProp/CoolProp/issues/443>`_ : Javascript index.html is missing
* `#437 <http://github.com/CoolProp/CoolProp/issues/437>`_ : REFPROP predefined mixtures no longer work
* `#434 <http://github.com/CoolProp/CoolProp/issues/434>`_ : R404A Refprop value differs from Refprop Value in Excel
* `#432 <http://github.com/CoolProp/CoolProp/issues/432>`_ : All the mixture interaction parameters of Gernert are wrong
* `#431 <http://github.com/CoolProp/CoolProp/issues/431>`_ : REFPROP should not be reloaded after every call to PropsSI
* `#430 <http://github.com/CoolProp/CoolProp/issues/430>`_ : HAPropsSI is missing from the SWIG wrapper
* `#429 <http://github.com/CoolProp/CoolProp/issues/429>`_ : Entropy of Melinder fluids giving wrong results
* `#428 <http://github.com/CoolProp/CoolProp/issues/428>`_ : On windows, do not error out if REFPROP fluid files are not found in c:\Program Files\REFPROP
* `#427 <http://github.com/CoolProp/CoolProp/issues/427>`_ : HapropsSi("W","B", 279.15, "T", 293.15, "P", 101325) lead to a "-1.#IND" value
* `#425 <http://github.com/CoolProp/CoolProp/issues/425>`_ : Incompressible viscosity
* `#419 <http://github.com/CoolProp/CoolProp/issues/419>`_ : HapropSI ("T","B",273.15+37,"D",273.15+36.44,"P",101325) lead to an error ...
* `#416 <http://github.com/CoolProp/CoolProp/issues/416>`_ : Sphinx docs
* `#413 <http://github.com/CoolProp/CoolProp/issues/413>`_ : Incompressible entropy
* `#410 <http://github.com/CoolProp/CoolProp/issues/410>`_ : Phase envelope segfaults for pure fluids
* `#409 <http://github.com/CoolProp/CoolProp/issues/409>`_ : Trivial outputs
* `#408 <http://github.com/CoolProp/CoolProp/issues/408>`_ : HapropsSI function issues
* `#403 <http://github.com/CoolProp/CoolProp/issues/403>`_ : Error in new CoolProp version in the function HAPropsSI (variable combination 'PH' and 'W')
* `#401 <http://github.com/CoolProp/CoolProp/issues/401>`_ : Linux/OSX error with refprop 9.1* and mixtures containing  R1234YF
* `#400 <http://github.com/CoolProp/CoolProp/issues/400>`_ : HAPropsSI(Output, "B",valueB, "R", 1, "P", 101325) lead to an error
* `#398 <http://github.com/CoolProp/CoolProp/issues/398>`_ : HAPropsSI(Output, "B",252.84, "D";250.85, "P", 101325) lead to an infinite value
* `#387 <http://github.com/CoolProp/CoolProp/issues/387>`_ : Vectorised PropSI breaks plotting functions
* `#386 <http://github.com/CoolProp/CoolProp/issues/386>`_ : Bibtex numbering
* `#307 <http://github.com/CoolProp/CoolProp/issues/307>`_ : Transport Properties for Mixtures


5.0.6
-----

New Features:

* Mathematica wrapper finished

Issues Closed:

* `#396 <http://github.com/CoolProp/CoolProp/issues/396>`_ : Initialize fail for HEOS in mixture with Argon and CarbonDioxide (in Matlab)
* `#395 <http://github.com/CoolProp/CoolProp/issues/395>`_ : keyed_output and incompressibles
* `#394 <http://github.com/CoolProp/CoolProp/issues/394>`_ : Python list inputs
* `#391 <http://github.com/CoolProp/CoolProp/issues/391>`_ : release.bsh and source file
* `#390 <http://github.com/CoolProp/CoolProp/issues/390>`_ : Transport properties of water
* `#389 <http://github.com/CoolProp/CoolProp/issues/389>`_ : HAPropsSI("D", "T",273.15+20, "R", 0.8, "P", 101325) lead to an error
* `#384 <http://github.com/CoolProp/CoolProp/issues/384>`_ : Put the example.nb Mathematica file in the main folder
* `#383 <http://github.com/CoolProp/CoolProp/issues/383>`_ : When doing release, force a full build of the docs
* `#382 <http://github.com/CoolProp/CoolProp/issues/382>`_ : Fix up the mathematica docs
* `#379 <http://github.com/CoolProp/CoolProp/issues/379>`_ : After a release is done, delete the release folder
* `#378 <http://github.com/CoolProp/CoolProp/issues/378>`_ : Also integrate the sphinx docs into the binaries/release/unstable folder output
* `#377 <http://github.com/CoolProp/CoolProp/issues/377>`_ : Remove old mathematica files
* `#376 <http://github.com/CoolProp/CoolProp/issues/376>`_ : Add python to list of prerequisites for self-compilation in the docs
* `#329 <http://github.com/CoolProp/CoolProp/issues/329>`_ : Configure buildbot to send emails when we break things

5.0.5
-----

New Features:

* Added Mathematica wrapper
* Added ``Prandtl()`` function to ``AbstractState``
* Added vectorized ``PropsSImulti`` function that can return a matrix of outputs for vectors of state inputs and desired outputs

Removed Features:

* All the ``PropsSI`` overloads.  For all other types of inputs, the ``PropsSImulti`` function is now used

Issues Closed:

* `#375 <http://github.com/CoolProp/CoolProp/issues/375>`_ : If one input and one output to PropsSI, bubble error cleanly
* `#373 <http://github.com/CoolProp/CoolProp/issues/373>`_ : Move predefined mixture parsing to HelmholtzEOS initializer function
* `#372 <http://github.com/CoolProp/CoolProp/issues/372>`_ : Prandtl number is missing from AbstractState
* `#371 <http://github.com/CoolProp/CoolProp/issues/371>`_ : Parse inputs to PropsSI/PropsSI(vectorized) and turn into a vector of inputs
* `#370 <http://github.com/CoolProp/CoolProp/issues/370>`_ : Docs are missing all the fluid files
* `#368 <http://github.com/CoolProp/CoolProp/issues/368>`_ : CoolProp on iOS
* `#367 <http://github.com/CoolProp/CoolProp/issues/367>`_ : Python module architecture
* `#366 <http://github.com/CoolProp/CoolProp/issues/366>`_ : Get value of universal gas constant
* `#365 <http://github.com/CoolProp/CoolProp/issues/365>`_ : REFPROP_lib.h is missed in 5.0.4 source code zip
* `#364 <http://github.com/CoolProp/CoolProp/issues/364>`_ : Liquid and vapor saturation pressures are not the same for some fluids
* `#363 <http://github.com/CoolProp/CoolProp/issues/363>`_ : Revision synchronisation
* `#359 <http://github.com/CoolProp/CoolProp/issues/359>`_ : Add high-level function that allows for multiple outputs
* `#357 <http://github.com/CoolProp/CoolProp/issues/357>`_ : Vector functions and state class
* `#349 <http://github.com/CoolProp/CoolProp/issues/349>`_ : Host v4 docs

5.0.4
-----

BUGFIX: Lots of bugs squashed. 

New features: 

* Julia wrapper added
* Derivatives along the saturation line for pure fluids implemented
* Exposed the configuration getter/setter through SWIG (except for MATLAB)
* Added transport properties for xylenes and Ethylbenzene
* Surface tension for HFC pseudo-pures added

Issues Closed:

* `#355 <http://github.com/CoolProp/CoolProp/issues/355>`_ : In MSVC, too many symbols are exported in SWIG+MATLAB
* `#354 <http://github.com/CoolProp/CoolProp/issues/354>`_ : REFPROP headers
* `#353 <http://github.com/CoolProp/CoolProp/issues/353>`_ : Using HAPropsSI within circular reference on Mac Excel 2011 causes div/0 error!
* `#350 <http://github.com/CoolProp/CoolProp/issues/350>`_ : Python module docs
* `#347 <http://github.com/CoolProp/CoolProp/issues/347>`_ : Implement calc_melting_line for incompressibles
* `#346 <http://github.com/CoolProp/CoolProp/issues/346>`_ : Memory sanitizer is reporting errors with RPVersion function call
* `#344 <http://github.com/CoolProp/CoolProp/issues/344>`_ : skip typeerror in Excel to make 32-bit xlam work in 64-bit excel
* `#342 <http://github.com/CoolProp/CoolProp/issues/342>`_ : Refprop mixture with 4 components error
* `#339 <http://github.com/CoolProp/CoolProp/issues/339>`_ : Some SWIG tests fail due to the inclusion of rapidjson header
* `#337 <http://github.com/CoolProp/CoolProp/issues/337>`_ : ECS not yielding the proper values for eta and lambda
* `#332 <http://github.com/CoolProp/CoolProp/issues/332>`_ : Make the REFPROP wrapper code 1% more sane
* `#331 <http://github.com/CoolProp/CoolProp/issues/331>`_ : Excel wapper shouts errors (in Excel 2013)
* `#330 <http://github.com/CoolProp/CoolProp/issues/330>`_ : Implement ECS model for viscosity of xylenes and ethylbenzene
* `#326 <http://github.com/CoolProp/CoolProp/issues/326>`_ : expose configuration through SWIG
* `#325 <http://github.com/CoolProp/CoolProp/issues/325>`_ : Implement the generalized derivatives for REFPROP as well
* `#324 <http://github.com/CoolProp/CoolProp/issues/324>`_ : SetPath for Refprop
* `#322 <http://github.com/CoolProp/CoolProp/issues/322>`_ : Add method to AbstractState to return mixture component names
* `#321 <http://github.com/CoolProp/CoolProp/issues/321>`_ : Add more R-number aliases
* `#320 <http://github.com/CoolProp/CoolProp/issues/320>`_ : HAPropsSI("T", "V", 0.83, "R", 1, "P", 101325) & lead to infinite value
* `#319 <http://github.com/CoolProp/CoolProp/issues/319>`_ : Error in entropy calculation with TH inputs
* `#314 <http://github.com/CoolProp/CoolProp/issues/314>`_ : Add surface tension reference information to docs
* `#312 <http://github.com/CoolProp/CoolProp/issues/312>`_ : Small examples of the use of derivatives should be in docs
* `#309 <http://github.com/CoolProp/CoolProp/issues/309>`_ : MEG properties
* `#308 <http://github.com/CoolProp/CoolProp/issues/308>`_ : Set maximum states for saturation curves for pseudo-pures properly
* `#306 <http://github.com/CoolProp/CoolProp/issues/306>`_ : Surface Tension for HFC Pseudo-Pure is missing
* `#304 <http://github.com/CoolProp/CoolProp/issues/304>`_ : Develop some docs about hooking up with Julia code
* `#294 <http://github.com/CoolProp/CoolProp/issues/294>`_ : Add the clang sanitize tests to buildbot
* `#247 <http://github.com/CoolProp/CoolProp/issues/247>`_ : Implement thermal conductivity for o-Xylene, m-Xylene, p-Xylene, and Ethylbenzene
* `#238 <http://github.com/CoolProp/CoolProp/issues/238>`_ : add a function to retrieve derivatives along the saturation curve


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