Changelog for CoolProp
======================

6.4.1
-----

Highlights:

* Fixed the nightly builds
* Added documentation for the NuGet packages by MadsKirkFoged

Deprecated:

* Removed the n-dimensional input for Python due to too many bugs

Issues Closed:

* `#1960 <https://github.com/CoolProp/CoolProp/issues/1960>`_ : Docs do not build
* `#1952 <https://github.com/CoolProp/CoolProp/issues/1952>`_ : CoolProp Module not Found
* `#1942 <https://github.com/CoolProp/CoolProp/issues/1942>`_ : Help with DLL hell
* `#1940 <https://github.com/CoolProp/CoolProp/issues/1940>`_ : CoolProp import doesn't work in v6.4.0

Pull Requests merged:

* `#1957 <https://github.com/CoolProp/CoolProp/pull/1957>`_ : Update BaseObjects.py
* `#1949 <https://github.com/CoolProp/CoolProp/pull/1949>`_ : fixed typo with units for sigma
* `#1947 <https://github.com/CoolProp/CoolProp/pull/1947>`_ : Update index.rst
* `#1932 <https://github.com/CoolProp/CoolProp/pull/1932>`_ : Allow nD-array input to PropsSI


6.4.0
-----

Highlights:

* Added a working version of PC-SAFT EOS (big thanks to Zach Baird)
* Updated EOS for neon, helium
* Added Python 3.8 wheels

Deprecated:

* Python 2.7 interface. This will be the last release supporting Python 2.7
* 32-bit builds on OSX

Issues Closed:

* `#1922 <https://github.com/CoolProp/CoolProp/issues/1922>`_ : Bug: The density of air is off by a magnitude of 10 given the P and T inputs specified below
* `#1915 <https://github.com/CoolProp/CoolProp/issues/1915>`_ : Error in Low Level Interface example code on coolprop.org
* `#1857 <https://github.com/CoolProp/CoolProp/issues/1857>`_ : calc_name is not implemented for REFPROP backend
* `#1856 <https://github.com/CoolProp/CoolProp/issues/1856>`_ : List of incompressible fluids/mixtures is missing
* `#1855 <https://github.com/CoolProp/CoolProp/issues/1855>`_ : Python AbstractState has no attribute 'compressibility_factor'
* `#1848 <https://github.com/CoolProp/CoolProp/issues/1848>`_ : PR density calc fail

Pull Requests merged:

* `#1921 <https://github.com/CoolProp/CoolProp/pull/1921>`_ : fix typos in pcsaft json
* `#1916 <https://github.com/CoolProp/CoolProp/pull/1916>`_ : Fix second_partial_deriv example
* `#1913 <https://github.com/CoolProp/CoolProp/pull/1913>`_ : Fix pcsaft flash
* `#1906 <https://github.com/CoolProp/CoolProp/pull/1906>`_ : Update index.rst
* `#1903 <https://github.com/CoolProp/CoolProp/pull/1903>`_ : Update index.rst adding PropiedadesDeFluidos Tool
* `#1901 <https://github.com/CoolProp/CoolProp/pull/1901>`_ : for CDNJS auto-update request
* `#1896 <https://github.com/CoolProp/CoolProp/pull/1896>`_ : Update init.py byte string in split for Py3.8
* `#1891 <https://github.com/CoolProp/CoolProp/pull/1891>`_ : Fix uppercase-only fluid naming
* `#1885 <https://github.com/CoolProp/CoolProp/pull/1885>`_ : fixed typo in getos() and else statement
* `#1881 <https://github.com/CoolProp/CoolProp/pull/1881>`_ : Fixed crash in PTflash_twophase::build_arrays
* `#1878 <https://github.com/CoolProp/CoolProp/pull/1878>`_ : Py38
* `#1877 <https://github.com/CoolProp/CoolProp/pull/1877>`_ : Adding PC-SAFT EOS
* `#1875 <https://github.com/CoolProp/CoolProp/pull/1875>`_ : Document apply_simple_mixing_rule initialisation behaviour
* `#1842 <https://github.com/CoolProp/CoolProp/pull/1842>`_ : fixes in doxygen docu for fundamental derivative
* `#1838 <https://github.com/CoolProp/CoolProp/pull/1838>`_ : Allow list delimiter to be changed
* `#1837 <https://github.com/CoolProp/CoolProp/pull/1837>`_ : Never allow two-phase SOS. 


6.3.0
-----

Highlights:

* The molar mass gets now reset properly - affected mixture calculations with changing concentrations.
* The humid air calculations check the inputs and outputs according to the limits from the publication.
* The isentropic expansion coefficient can now be accessed directly.
* ... and a lot of little bugfixes (see issues)

Issues Closed:

* `#1820 <https://github.com/CoolProp/CoolProp/issues/1820>`_ : Humid air example fails due to new limits
* `#1815 <https://github.com/CoolProp/CoolProp/issues/1815>`_ : molar_mass not getting cleared
* `#1811 <https://github.com/CoolProp/CoolProp/issues/1811>`_ : Humid air properties above 388 K
* `#1786 <https://github.com/CoolProp/CoolProp/issues/1786>`_ : Incompressible docs not building properly
* `#1784 <https://github.com/CoolProp/CoolProp/issues/1784>`_ : Sphinx builder still broken
* `#1782 <https://github.com/CoolProp/CoolProp/issues/1782>`_ : OSX 10.14 builds
* `#1778 <https://github.com/CoolProp/CoolProp/issues/1778>`_ : There are no nightly builds after 2018/11/04
* `#1777 <https://github.com/CoolProp/CoolProp/issues/1777>`_ : Building from the PyPI sdist on Python 3.7 results in compilation errors on macOS
* `#1775 <https://github.com/CoolProp/CoolProp/issues/1775>`_ : Tmin function cannot be evaluated at 0.0 concentration for incomp fluids
* `#1763 <https://github.com/CoolProp/CoolProp/issues/1763>`_ : Mathcad 15 binary builds dropped as of version 6.2
* `#1762 <https://github.com/CoolProp/CoolProp/issues/1762>`_ : IF97 Documentation Page Error
* `#1760 <https://github.com/CoolProp/CoolProp/issues/1760>`_ : Android Wrapper error at 6.2.0 and 6.2.2dev
* `#1759 <https://github.com/CoolProp/CoolProp/issues/1759>`_ : Memory leak in Mathematica interface
* `#1758 <https://github.com/CoolProp/CoolProp/issues/1758>`_ : Build AbstractState object from stored tabular data
* `#1756 <https://github.com/CoolProp/CoolProp/issues/1756>`_ : Issue with incompressible fluid in v6.2.1
* `#1753 <https://github.com/CoolProp/CoolProp/issues/1753>`_ : numpy.core.multiarray failed to import
* `#1752 <https://github.com/CoolProp/CoolProp/issues/1752>`_ : Add fluids to CoolProp if you are using matlab
* `#1748 <https://github.com/CoolProp/CoolProp/issues/1748>`_ : Apostrophe should be escaped in '...' strings or be used in "..." string
* `#1745 <https://github.com/CoolProp/CoolProp/issues/1745>`_ : Surface Tension calculation failing ungracefully
* `#1744 <https://github.com/CoolProp/CoolProp/issues/1744>`_ : Error value Excel on Mac
* `#1742 <https://github.com/CoolProp/CoolProp/issues/1742>`_ : r404a JT valve error
* `#1741 <https://github.com/CoolProp/CoolProp/issues/1741>`_ : Wrapper for OCTAVE-4.x.x 32 bit for Windows
* `#1734 <https://github.com/CoolProp/CoolProp/issues/1734>`_ : Can I access the older 'Props' functions in CoolProp 6
* `#1732 <https://github.com/CoolProp/CoolProp/issues/1732>`_ : Error 53 Excel Wrapper MacBook Issue
* `#1731 <https://github.com/CoolProp/CoolProp/issues/1731>`_ : Will CoolProp contain R513a refrigerant properties in the near future??
* `#1724 <https://github.com/CoolProp/CoolProp/issues/1724>`_ : REFPROP v10.0 enthalpy/pressure look-up bug
* `#1717 <https://github.com/CoolProp/CoolProp/issues/1717>`_ : Coolprop cannot work on local JS
* `#1708 <https://github.com/CoolProp/CoolProp/issues/1708>`_ : calc_phase missing for REFPROP backend.
* `#1707 <https://github.com/CoolProp/CoolProp/issues/1707>`_ : Inconsistency in incompressible backend when getting fluid name()
* `#1700 <https://github.com/CoolProp/CoolProp/issues/1700>`_ : HAPropsSi broken in 6.1.1dev for R(p,T,h)
* `#1618 <https://github.com/CoolProp/CoolProp/issues/1618>`_ : Excel and Python wrapper, HAPropsSI problem with P_w input
* `#1601 <https://github.com/CoolProp/CoolProp/issues/1601>`_ : ValueError: HAProps failed ungracefully with inputs: "W","T",2.7852603934336025e+02,"P",1.0132500000000000e+05,"H",1.4114309647910737e+04
* `#1465 <https://github.com/CoolProp/CoolProp/issues/1465>`_ : Humid air calculations when dry bulb is above freezing but wet bulb is below freezing

Pull Requests merged:

* `#1823 <https://github.com/CoolProp/CoolProp/pull/1823>`_ : Robustify humid air limit checks
* `#1814 <https://github.com/CoolProp/CoolProp/pull/1814>`_ : Feature/docs py3
* `#1813 <https://github.com/CoolProp/CoolProp/pull/1813>`_ : Added limits to the humid air properties, closes #1811
* `#1810 <https://github.com/CoolProp/CoolProp/pull/1810>`_ : Use std::shared_ptr if available in VS
* `#1797 <https://github.com/CoolProp/CoolProp/pull/1797>`_ : Set _phase in REFPROP backend and return it in calc_phase()
* `#1791 <https://github.com/CoolProp/CoolProp/pull/1791>`_ : Add isentropic expansion coefficient
* `#1787 <https://github.com/CoolProp/CoolProp/pull/1787>`_ : Add conda install instructions from conda-forge
* `#1783 <https://github.com/CoolProp/CoolProp/pull/1783>`_ : Feature/recent compilers
* `#1773 <https://github.com/CoolProp/CoolProp/pull/1773>`_ : Make wrapper for LibreOffice an extension
* `#1770 <https://github.com/CoolProp/CoolProp/pull/1770>`_ : Disambiguate IF97 Fluid Names - Water only
* `#1767 <https://github.com/CoolProp/CoolProp/pull/1767>`_ : Add IF97 links to CoolProp homepage and backends page
* `#1765 <https://github.com/CoolProp/CoolProp/pull/1765>`_ : Patch PropsSI for INCOMP fluid with zero mass fraction
* `#1761 <https://github.com/CoolProp/CoolProp/pull/1761>`_ : Free input strings in Mathematica interface
* `#1755 <https://github.com/CoolProp/CoolProp/pull/1755>`_ : Throw sensible error message if single-phase surface tension requested
* `#1751 <https://github.com/CoolProp/CoolProp/pull/1751>`_ : update syntax to julia 1.0

6.2.1
-----

Only a minor fix to the Javascript wrapper


6.2.0
-----

New Features:

* Added a new EoS for heavy water
* Added pre-built wheels for Python 3.6 and Python 3.7
* The MATLAB wrappers have been abandoned in favour of Python-based calls
* Add phase specification to high-level interface
* LabVIEW VIs can now call PropsSI and Props1SI
* Added a wrapper for Android
* ... and a lot of little bugfixes (see issues)

Issues Closed:

* `#1699 <https://github.com/CoolProp/CoolProp/issues/1699>`_ : Pip Install problem with Ubuntu 18.04 and Python 3.6.3
* `#1682 <https://github.com/CoolProp/CoolProp/issues/1682>`_ : Coolprop.Coolprop module does not exist
* `#1672 <https://github.com/CoolProp/CoolProp/issues/1672>`_ : In ODEintegrator, limits are wrong for backwards
* `#1662 <https://github.com/CoolProp/CoolProp/issues/1662>`_ : Wrong results when using INCOMP (incompressible) fluids using high-level interface after importing matplotlib.pyplot
* `#1661 <https://github.com/CoolProp/CoolProp/issues/1661>`_ : install fail with python 3.6 in anaconda for win10
* `#1659 <https://github.com/CoolProp/CoolProp/issues/1659>`_ : More reducing state resetting needed when fractions set
* `#1654 <https://github.com/CoolProp/CoolProp/issues/1654>`_ : Version 6.1 with Python 3.6.1 Anaconda (64-bit)
* `#1652 <https://github.com/CoolProp/CoolProp/issues/1652>`_ : Problem with saturated vapor internal energy calculations with quality/density inputs
* `#1649 <https://github.com/CoolProp/CoolProp/issues/1649>`_ : Cannot cimport CoolProp into cython in python 3.6
* `#1647 <https://github.com/CoolProp/CoolProp/issues/1647>`_ : Parsing of Mixtures depends on LOCALE
* `#1630 <https://github.com/CoolProp/CoolProp/issues/1630>`_ : Predefined mixture cannot have uppercase .MIX
* `#1629 <https://github.com/CoolProp/CoolProp/issues/1629>`_ : Phase envelopes fail for predefined mixtures with REFPROP backend
* `#1607 <https://github.com/CoolProp/CoolProp/issues/1607>`_ : Tabular Backend Fails with HmassP_INPUTS when iphase_twophase Imposed
* `#1604 <https://github.com/CoolProp/CoolProp/issues/1604>`_ : v6.2?
* `#1603 <https://github.com/CoolProp/CoolProp/issues/1603>`_ : Parse out Zero Mass Fraction Components in High-Level Interface
* `#1602 <https://github.com/CoolProp/CoolProp/issues/1602>`_ : hmass() gives strange result after calling update() with PQ_INPUTS in specific case
* `#1582 <https://github.com/CoolProp/CoolProp/issues/1582>`_ : Buildbot update
* `#1563 <https://github.com/CoolProp/CoolProp/issues/1563>`_ : Unify AbstractState's behavior when using HEOS or the tabular interpolations schemes
* `#1551 <https://github.com/CoolProp/CoolProp/issues/1551>`_ : Import of matplotlib.pyplot results in error for mixtures
* `#1530 <https://github.com/CoolProp/CoolProp/issues/1530>`_ : Catch tests failing
* `#1455 <https://github.com/CoolProp/CoolProp/issues/1455>`_ : apply_simple_mixing_rule broken
* `#1439 <https://github.com/CoolProp/CoolProp/issues/1439>`_ : Wrong dam_dtau for Twu
* `#1426 <https://github.com/CoolProp/CoolProp/issues/1426>`_ : UNIFAQ compile errors
* `#1422 <https://github.com/CoolProp/CoolProp/issues/1422>`_ : Ttriple wrong for REFPROP for water
* `#1406 <https://github.com/CoolProp/CoolProp/issues/1406>`_ : StateContainer print fails
* `#1396 <https://github.com/CoolProp/CoolProp/issues/1396>`_ : Formulas are wrong for dichloroethane and ethylene oxide
* `#1393 <https://github.com/CoolProp/CoolProp/issues/1393>`_ : Crash when set_mole_fractions() not called
* `#1381 <https://github.com/CoolProp/CoolProp/issues/1381>`_ : Calling acentric factor with cubic equation
* `#1372 <https://github.com/CoolProp/CoolProp/issues/1372>`_ : inconsistent result with mixture of "Ethylbenzene[0.5]&P-XYLENE[0.5]"<>"ethylbenzene[0.5]&P-XYLENE[0.5]"
* `#1371 <https://github.com/CoolProp/CoolProp/issues/1371>`_ : Get JSON string for fluid at runtime
* `#1369 <https://github.com/CoolProp/CoolProp/issues/1369>`_ : Return 'n/a'  REFPROP version when not loaded or supported.
* `#1368 <https://github.com/CoolProp/CoolProp/issues/1368>`_ : segmentation fault when calling get_global_param_string("REFPROP_version") from Python
* `#1366 <https://github.com/CoolProp/CoolProp/issues/1366>`_ : Allow fluids to be provided (and overwritten) at runtime; closes #1345
* `#1365 <https://github.com/CoolProp/CoolProp/issues/1365>`_ : SMath Wrapper refactoring
* `#1362 <https://github.com/CoolProp/CoolProp/issues/1362>`_ : LabVIEW VIs to call PropsSI and Props1SI
* `#1361 <https://github.com/CoolProp/CoolProp/issues/1361>`_ : Re-enable alpha0 mixture derivative tests for cubics
* `#1359 <https://github.com/CoolProp/CoolProp/issues/1359>`_ : Allow for cubic transformations in HEOS multi-fluid model
* `#1355 <https://github.com/CoolProp/CoolProp/issues/1355>`_ : SMath Wrapper refactoring
* `#1354 <https://github.com/CoolProp/CoolProp/issues/1354>`_ : splined properties, _fluid_type and _phase not cleared in AbstractState.h
* `#1352 <https://github.com/CoolProp/CoolProp/issues/1352>`_ : Faulty state update with DmassT_inputs in HEOS backend with specified phase
* `#1350 <https://github.com/CoolProp/CoolProp/issues/1350>`_ : Simulation error when using ExternalMedia in Dymola
* `#1348 <https://github.com/CoolProp/CoolProp/issues/1348>`_ : Allow alpha0 to be provided for cubic EOS
* `#1347 <https://github.com/CoolProp/CoolProp/issues/1347>`_ : Add ability to ignore setup errors for REFPROP mixtures
* `#1343 <https://github.com/CoolProp/CoolProp/issues/1343>`_ : Call git in the dev folder
* `#1339 <https://github.com/CoolProp/CoolProp/issues/1339>`_ : Set a standard departure function through the AbstractState
* `#1333 <https://github.com/CoolProp/CoolProp/issues/1333>`_ : Make it possible to use x[i]=0 in some alpha0 derivatives
* `#1329 <https://github.com/CoolProp/CoolProp/issues/1329>`_ : DO NOT allow for over-writing of departure functions when loading defaults
* `#1328 <https://github.com/CoolProp/CoolProp/issues/1328>`_ : Dmass wrong for saturated states for REFPROP
* `#1325 <https://github.com/CoolProp/CoolProp/issues/1325>`_ : Also export HAProps_Aux to pybind11 interface
* `#1324 <https://github.com/CoolProp/CoolProp/issues/1324>`_ : Figure out problem with linux wheels
* `#1323 <https://github.com/CoolProp/CoolProp/issues/1323>`_ : Added PQ and QT Input Pairs to provide Saturation Values to IF97 Backend
* `#1322 <https://github.com/CoolProp/CoolProp/issues/1322>`_ : Bigger buffer size for Julia wrapper
* `#1321 <https://github.com/CoolProp/CoolProp/issues/1321>`_ : Finally fix phase envelopes again
* `#1320 <https://github.com/CoolProp/CoolProp/issues/1320>`_ : Figure out why catch runs take forever on "*"nix
* `#1319 <https://github.com/CoolProp/CoolProp/issues/1319>`_ : Fix python windows builds
* `#1318 <https://github.com/CoolProp/CoolProp/issues/1318>`_ : Move entire Emscripten interface into its own file that is included separately
* `#1317 <https://github.com/CoolProp/CoolProp/issues/1317>`_ : Loading HMX.BNC through the DLL yields weird behavior
* `#1316 <https://github.com/CoolProp/CoolProp/issues/1316>`_ : Added configuration options for MSVCRT linking, changed the output di…
* `#1314 <https://github.com/CoolProp/CoolProp/issues/1314>`_ : Android Wrapper
* `#1312 <https://github.com/CoolProp/CoolProp/issues/1312>`_ : First step toward #1310
* `#1309 <https://github.com/CoolProp/CoolProp/issues/1309>`_ : version 6.1.0 not available from pypi
* `#1308 <https://github.com/CoolProp/CoolProp/issues/1308>`_ : Add Trivial Parameter calls to IF97 Backend
* `#1307 <https://github.com/CoolProp/CoolProp/issues/1307>`_ : get_config_string returns nothing in python
* `#1306 <https://github.com/CoolProp/CoolProp/issues/1306>`_ : Typo in CO2+Argon coefficients
* `#1305 <https://github.com/CoolProp/CoolProp/issues/1305>`_ : Fix some warnings in MSVC 2015
* `#1304 <https://github.com/CoolProp/CoolProp/issues/1304>`_ : Parse refprop HMX.BNC file and load coefficients
* `#1303 <https://github.com/CoolProp/CoolProp/issues/1303>`_ : call refprop from coolprop in scilab on linux
* `#1302 <https://github.com/CoolProp/CoolProp/issues/1302>`_ : Export cubic's alpha functions
* `#1300 <https://github.com/CoolProp/CoolProp/issues/1300>`_ : Add criticality_contour_values to pybind11 interface
* `#1299 <https://github.com/CoolProp/CoolProp/issues/1299>`_ : Add config keys to pybind11 interface
* `#1298 <https://github.com/CoolProp/CoolProp/issues/1298>`_ : HAPropsSI H, p, W lookups not working past 5.0.0
* `#1296 <https://github.com/CoolProp/CoolProp/issues/1296>`_ : Add phases enum to pybind11 interface
* `#1295 <https://github.com/CoolProp/CoolProp/issues/1295>`_ : Specify the minimum delta for spinodal tracer as config variable
* `#1294 <https://github.com/CoolProp/CoolProp/issues/1294>`_ : Add parser for HMX.BNC from REFPROP
* `#1292 <https://github.com/CoolProp/CoolProp/issues/1292>`_ : Source zip file on SourceForge is not correct again...
* `#1289 <https://github.com/CoolProp/CoolProp/issues/1289>`_ : Make triple point accessible in HEOS::get_fluid_constant
* `#1285 <https://github.com/CoolProp/CoolProp/issues/1285>`_ : Allow fluids to be overwritten
* `#1279 <https://github.com/CoolProp/CoolProp/issues/1279>`_ : Add phase specification to high-level interface
* `#1253 <https://github.com/CoolProp/CoolProp/issues/1253>`_ : Implement derivatives of alpha0 w.r.t. tau, delta
* `#1249 <https://github.com/CoolProp/CoolProp/issues/1249>`_ : IF97 Error in CoolProp Wrapping for SMath
* `#969 <https://github.com/CoolProp/CoolProp/issues/969>`_ : Support mixtures with component mole fractions of zero

Pull Requests merged:

* `#1675 <https://github.com/CoolProp/CoolProp/pull/1675>`_ : Let DARWIN build with libc++
* `#1666 <https://github.com/CoolProp/CoolProp/pull/1666>`_ : Make string->float conversion aware of the locale
* `#1665 <https://github.com/CoolProp/CoolProp/pull/1665>`_ : Patches PropsSI imposed phase for backends other than HEOS
* `#1660 <https://github.com/CoolProp/CoolProp/pull/1660>`_ : Update PropsSI() to Parse Imposed Phase Strings on Input Keys
* `#1656 <https://github.com/CoolProp/CoolProp/pull/1656>`_ : Mistake in function 'inline_label' in CoolProp/Plots/Common.py
* `#1645 <https://github.com/CoolProp/CoolProp/pull/1645>`_ : Provide return string from PhaseSI() if phase can't be determined.
* `#1609 <https://github.com/CoolProp/CoolProp/pull/1609>`_ : editorconfig
* `#1606 <https://github.com/CoolProp/CoolProp/pull/1606>`_ : Patch PT_flash() to update _phase with imposed phase, in case it changed
* `#1464 <https://github.com/CoolProp/CoolProp/pull/1464>`_ : Fix a few REFPROP functions; closes #1422
* `#1460 <https://github.com/CoolProp/CoolProp/pull/1460>`_ : Greatly improve the stability of REFPROP mixture calls at saturation …
* `#1459 <https://github.com/CoolProp/CoolProp/pull/1459>`_ : Call SATTP properly when QT inputs are given for REFPROP
* `#1458 <https://github.com/CoolProp/CoolProp/pull/1458>`_ : Actually set the Twu parameters if provided
* `#1457 <https://github.com/CoolProp/CoolProp/pull/1457>`_ : Add ierr checks to calls to SETKTV
* `#1450 <https://github.com/CoolProp/CoolProp/pull/1450>`_ : Fix typo in CoolPropLib.h
* `#1449 <https://github.com/CoolProp/CoolProp/pull/1449>`_ : Move F2K into emscripten_interface.cxx
* `#1448 <https://github.com/CoolProp/CoolProp/pull/1448>`_ : Update the ODE integrator to allow it to integrate backwards
* `#1376 <https://github.com/CoolProp/CoolProp/pull/1376>`_ : Update HumidAirProp.cpp

6.1.0
-----

New features:

* Windows installer for Microsoft Excel
* Added VTPR backend
* Twu and Mathias-Copeman attractive parameters can be set for PR and SRK
* Major improvements to Excel wrapper
* Added EOS for MDM of M. Thol
* Implemented first version of PT flash calculations for two-phase states
* Implemented PT flash for mixtures (not finished)
* Added a pybind11 module for CoolProp
* ... and a lot of little bugfixes (see issues)

Contributors to this release:
ibell, JonWel, jowr, babaksamareh, mikekaganski

* `#1290 <https://github.com/CoolProp/CoolProp/issues/1290>`_ : Catch runs should be Release builds
* `#1288 <https://github.com/CoolProp/CoolProp/issues/1288>`_ : Actually check if T < Tmelt for p > pmin
* `#1287 <https://github.com/CoolProp/CoolProp/issues/1287>`_ : Actually commit new pybind11 submodule
* `#1286 <https://github.com/CoolProp/CoolProp/issues/1286>`_ : in phase envelope construction, potential crash
* `#1284 <https://github.com/CoolProp/CoolProp/issues/1284>`_ : Make low-level interface accessible through high-level interface in FORTRAN
* `#1283 <https://github.com/CoolProp/CoolProp/issues/1283>`_ : Add pure fluid check to VTPR
* `#1282 <https://github.com/CoolProp/CoolProp/issues/1282>`_ : Correct typo, see #1270
* `#1281 <https://github.com/CoolProp/CoolProp/issues/1281>`_ : Add ability to add HEOS fluids as JSON at runtime
* `#1272 <https://github.com/CoolProp/CoolProp/issues/1272>`_ : Solves a bug in VTPR
* `#1271 <https://github.com/CoolProp/CoolProp/issues/1271>`_ : Remove possible division by 0, closes #1270
* `#1269 <https://github.com/CoolProp/CoolProp/issues/1269>`_ : SatL and SatV of type VTPR too
* `#1268 <https://github.com/CoolProp/CoolProp/issues/1268>`_ : Implement fluid_names for cubic backend
* `#1267 <https://github.com/CoolProp/CoolProp/issues/1267>`_ : PengRobinson doesn't pass alpha to SatL and SatV
* `#1266 <https://github.com/CoolProp/CoolProp/issues/1266>`_ : Small fixes for VTPR
* `#1264 <https://github.com/CoolProp/CoolProp/issues/1264>`_ : Update initialization for VTPR
* `#1262 <https://github.com/CoolProp/CoolProp/issues/1262>`_ : Set alpha function in JSON
* `#1261 <https://github.com/CoolProp/CoolProp/issues/1261>`_ : Update CMakeLists.txt
* `#1259 <https://github.com/CoolProp/CoolProp/issues/1259>`_ : Methanol-water mixture: strange results
* `#1258 <https://github.com/CoolProp/CoolProp/issues/1258>`_ : Solves a bug with cubic and mixtures
* `#1257 <https://github.com/CoolProp/CoolProp/issues/1257>`_ : Update iPhone compilation docs
* `#1255 <https://github.com/CoolProp/CoolProp/issues/1255>`_ : Allow ability to set Twu parameters for cubic EOS (from JSON)
* `#1252 <https://github.com/CoolProp/CoolProp/issues/1252>`_ : Implement set_double_array2D
* `#1250 <https://github.com/CoolProp/CoolProp/issues/1250>`_ : Implement coefficient derivatives of dYr_dxi in reducing function
* `#1248 <https://github.com/CoolProp/CoolProp/issues/1248>`_ : Problem with OSX compilation
* `#1240 <https://github.com/CoolProp/CoolProp/issues/1240>`_ : Make psi_plus public
* `#1239 <https://github.com/CoolProp/CoolProp/issues/1239>`_ : Shortcut VTPR when pure fluids, solves #1232
* `#1237 <https://github.com/CoolProp/CoolProp/issues/1237>`_ : Create an installer for selected Windows wrappers
* `#1235 <https://github.com/CoolProp/CoolProp/issues/1235>`_ : Excel 2016 Add-In Updates
* `#1234 <https://github.com/CoolProp/CoolProp/issues/1234>`_ : Add the ability to set limits in Consistency plots
* `#1232 <https://github.com/CoolProp/CoolProp/issues/1232>`_ : VTPR components with one group
* `#1230 <https://github.com/CoolProp/CoolProp/issues/1230>`_ : Allow ability to call REFPROP on OSX
* `#1229 <https://github.com/CoolProp/CoolProp/issues/1229>`_ : ConsistencyPlots updates
* `#1227 <https://github.com/CoolProp/CoolProp/issues/1227>`_ : Make all functions in DepartureFunction overrridable
* `#1226 <https://github.com/CoolProp/CoolProp/issues/1226>`_ : More critical point questions
* `#1222 <https://github.com/CoolProp/CoolProp/issues/1222>`_ : Critical point calc failure
* `#1221 <https://github.com/CoolProp/CoolProp/issues/1221>`_ : Take more steps in stability evaluator (at least 100)
* `#1220 <https://github.com/CoolProp/CoolProp/issues/1220>`_ : Add adaptive integrator code
* `#1219 <https://github.com/CoolProp/CoolProp/issues/1219>`_ : Double post_update in update_TP_guessrho
* `#1217 <https://github.com/CoolProp/CoolProp/issues/1217>`_ : Peng-Robinson issue with Hydrogen
* `#1215 <https://github.com/CoolProp/CoolProp/issues/1215>`_ : Vapour QT_INPUT with VTPR
* `#1214 <https://github.com/CoolProp/CoolProp/issues/1214>`_ : Refactor exceptions in CoolPropLib.cpp close #1200
* `#1213 <https://github.com/CoolProp/CoolProp/issues/1213>`_ : Add tests for Poling example with UNIFAC code
* `#1212 <https://github.com/CoolProp/CoolProp/issues/1212>`_ : Add derivatives of a*rho with respect to tau,delta,x
* `#1211 <https://github.com/CoolProp/CoolProp/issues/1211>`_ : Use aii_term and b0_ii from cubic
* `#1209 <https://github.com/CoolProp/CoolProp/issues/1209>`_ : Correct tau derivatives in VTPR
* `#1208 <https://github.com/CoolProp/CoolProp/issues/1208>`_ : Correct derivatives of am and test for VTPR
* `#1206 <https://github.com/CoolProp/CoolProp/issues/1206>`_ : Segmentation fault when calling get_mass_fractions() with SRK and PR
* `#1204 <https://github.com/CoolProp/CoolProp/issues/1204>`_ : Make all functions in reducing function const
* `#1203 <https://github.com/CoolProp/CoolProp/issues/1203>`_ : Allow VTPR to pass only names by setting default R_u value
* `#1202 <https://github.com/CoolProp/CoolProp/issues/1202>`_ : Better error message when UNIFAC component cannot be found
* `#1201 <https://github.com/CoolProp/CoolProp/issues/1201>`_ : Update MixtureDerivatives.cpp
* `#1199 <https://github.com/CoolProp/CoolProp/issues/1199>`_ : dalpha0_dxi is wrong
* `#1198 <https://github.com/CoolProp/CoolProp/issues/1198>`_ : Cubic CP
* `#1197 <https://github.com/CoolProp/CoolProp/issues/1197>`_ : Cubic QT_INPUTS
* `#1196 <https://github.com/CoolProp/CoolProp/issues/1196>`_ : Update CoolPropLib.def
* `#1195 <https://github.com/CoolProp/CoolProp/issues/1195>`_ : Merge VTPR
* `#1193 <https://github.com/CoolProp/CoolProp/issues/1193>`_ : REFPROP backend is missing acentric factor accessor
* `#1192 <https://github.com/CoolProp/CoolProp/issues/1192>`_ : Missing formulas for some HFO
* `#1191 <https://github.com/CoolProp/CoolProp/issues/1191>`_ : Linked states need to be updated in copy_k
* `#1190 <https://github.com/CoolProp/CoolProp/issues/1190>`_ : Problems running the VB.NET and C# wrappers
* `#1189 <https://github.com/CoolProp/CoolProp/issues/1189>`_ : Cubic backend broken for PQ calls
* `#1188 <https://github.com/CoolProp/CoolProp/issues/1188>`_ : Critical state not copying for cubics
* `#1187 <https://github.com/CoolProp/CoolProp/issues/1187>`_ : All critical points destroy density solver
* `#1185 <https://github.com/CoolProp/CoolProp/issues/1185>`_ : Add 4th order solver (Halley+)
* `#1184 <https://github.com/CoolProp/CoolProp/issues/1184>`_ : Add 4th order alphar derivatives to python
* `#1183 <https://github.com/CoolProp/CoolProp/issues/1183>`_ : QT/PQ inputs needs to polish with Newton-Raphson
* `#1182 <https://github.com/CoolProp/CoolProp/issues/1182>`_ : Add function to generate rapidjson instance from JSON string
* `#1181 <https://github.com/CoolProp/CoolProp/issues/1181>`_ : Add warning about T > Tmax for HS inputs
* `#1180 <https://github.com/CoolProp/CoolProp/issues/1180>`_ : CoolProp add-in for Excel not working on re-opened files
* `#1179 <https://github.com/CoolProp/CoolProp/issues/1179>`_ : Add derivatives of vr and Tr with respect to beta and gamma
* `#1178 <https://github.com/CoolProp/CoolProp/issues/1178>`_ : Android Wrapper Undefined Reference error with latest ndk
* `#1176 <https://github.com/CoolProp/CoolProp/issues/1176>`_ : [VTPR] mole fractions must be set before calling set_temperature
* `#1175 <https://github.com/CoolProp/CoolProp/issues/1175>`_ : Impose phase for REFPROP in low-level interface
* `#1174 <https://github.com/CoolProp/CoolProp/issues/1174>`_ : Update PHP module docs
* `#1172 <https://github.com/CoolProp/CoolProp/issues/1172>`_ : Please Help With Java Wrapper
* `#1170 <https://github.com/CoolProp/CoolProp/issues/1170>`_ : Incorrect InChI keys
* `#1169 <https://github.com/CoolProp/CoolProp/issues/1169>`_ : Issue with PropsSI on Methane-Ethane mixtures
* `#1168 <https://github.com/CoolProp/CoolProp/issues/1168>`_ : Volume translation for cubic
* `#1166 <https://github.com/CoolProp/CoolProp/issues/1166>`_ : Thermodynamic Properties of R1233zd(E)
* `#1165 <https://github.com/CoolProp/CoolProp/issues/1165>`_ : Not erroring if T < Tmin and p > ptriple
* `#1164 <https://github.com/CoolProp/CoolProp/issues/1164>`_ : REFPROP doesn't store mole fractions in phase envelope
* `#1161 <https://github.com/CoolProp/CoolProp/issues/1161>`_ : [VTPR] gE/RT needs to be multiplied by RT
* `#1158 <https://github.com/CoolProp/CoolProp/issues/1158>`_ : Retrieve phase envelope through high-level DLL
* `#1150 <https://github.com/CoolProp/CoolProp/issues/1150>`_ : IF97 backend: Wrong results for cvmass
* `#1148 <https://github.com/CoolProp/CoolProp/issues/1148>`_ : Add new EOS for MDM of Thol
* `#1146 <https://github.com/CoolProp/CoolProp/issues/1146>`_ : MEXW32 is actually 64-bit and crashes MATLAB
* `#1145 <https://github.com/CoolProp/CoolProp/issues/1145>`_ : Re-implement fundamental derivative of gas dynamics
* `#1144 <https://github.com/CoolProp/CoolProp/issues/1144>`_ : Repair use of spinodals and cubic backend
* `#1143 <https://github.com/CoolProp/CoolProp/issues/1143>`_ : PT inputs for cubics without phase specification
* `#1142 <https://github.com/CoolProp/CoolProp/issues/1142>`_ : PQ inputs very slow for cubic backends
* `#1141 <https://github.com/CoolProp/CoolProp/issues/1141>`_ : dichloroethane has the wrong CAS #
* `#1137 <https://github.com/CoolProp/CoolProp/issues/1137>`_ : Nonsensical results for mistaken inputs with INCOMP fluids
* `#1122 <https://github.com/CoolProp/CoolProp/issues/1122>`_ : Calculate density with PropsSi in Javascript
* `#1120 <https://github.com/CoolProp/CoolProp/issues/1120>`_ : Allow state generation from backend_name() return values
* `#1118 <https://github.com/CoolProp/CoolProp/issues/1118>`_ : Fix plots for cases with multiple critical points
* `#1114 <https://github.com/CoolProp/CoolProp/issues/1114>`_ : Export set_binary_interaction_double + Julia wrapper improvement
* `#1111 <https://github.com/CoolProp/CoolProp/issues/1111>`_ : Improvements to SMath wrapper error handling and some small tweaks
* `#1109 <https://github.com/CoolProp/CoolProp/issues/1109>`_ : SMath wrapper: update AssemblyInfo.cs.template
* `#1108 <https://github.com/CoolProp/CoolProp/issues/1108>`_ : SMath copyright year outdated
* `#1107 <https://github.com/CoolProp/CoolProp/issues/1107>`_ : Allow conditional build of SMath in source tree (fixes #1110)

Pull Requests merged:

* `#1283 <https://github.com/CoolProp/CoolProp/pull/1283>`_ : Add pure fluid check to VTPR
* `#1282 <https://github.com/CoolProp/CoolProp/pull/1282>`_ : Correct typo, see #1270
* `#1272 <https://github.com/CoolProp/CoolProp/pull/1272>`_ : Solves a bug in VTPR
* `#1271 <https://github.com/CoolProp/CoolProp/pull/1271>`_ : Remove possible division by 0, closes #1270
* `#1269 <https://github.com/CoolProp/CoolProp/pull/1269>`_ : SatL and SatV of type VTPR too
* `#1266 <https://github.com/CoolProp/CoolProp/pull/1266>`_ : Small fixes for VTPR
* `#1262 <https://github.com/CoolProp/CoolProp/pull/1262>`_ : Set alpha function in JSON
* `#1261 <https://github.com/CoolProp/CoolProp/pull/1261>`_ : Update CMakeLists.txt
* `#1258 <https://github.com/CoolProp/CoolProp/pull/1258>`_ : Solves a bug with cubic and mixtures
* `#1257 <https://github.com/CoolProp/CoolProp/pull/1257>`_ : Update iPhone compilation docs
* `#1239 <https://github.com/CoolProp/CoolProp/pull/1239>`_ : Shortcut VTPR when pure fluids, solves #1232
* `#1234 <https://github.com/CoolProp/CoolProp/pull/1234>`_ : Add the ability to set limits in Consistency plots
* `#1214 <https://github.com/CoolProp/CoolProp/pull/1214>`_ : Refactor exceptions in CoolPropLib.cpp close #1200
* `#1211 <https://github.com/CoolProp/CoolProp/pull/1211>`_ : Use aii_term and b0_ii from cubic
* `#1209 <https://github.com/CoolProp/CoolProp/pull/1209>`_ : Correct tau derivatives in VTPR
* `#1208 <https://github.com/CoolProp/CoolProp/pull/1208>`_ : Correct derivatives of am and test for VTPR
* `#1196 <https://github.com/CoolProp/CoolProp/pull/1196>`_ : Update CoolPropLib.def
* `#1195 <https://github.com/CoolProp/CoolProp/pull/1195>`_ : Merge VTPR
* `#1114 <https://github.com/CoolProp/CoolProp/pull/1114>`_ : Export set_binary_interaction_double + Julia wrapper improvement
* `#1111 <https://github.com/CoolProp/CoolProp/pull/1111>`_ : Improvements to SMath wrapper error handling and some small tweaks
* `#1109 <https://github.com/CoolProp/CoolProp/pull/1109>`_ : SMath wrapper: update AssemblyInfo.cs.template
* `#1107 <https://github.com/CoolProp/CoolProp/pull/1107>`_ : Allow conditional build of SMath in source tree (fixes #1110)
* `#1103 <https://github.com/CoolProp/CoolProp/pull/1103>`_ : One small tweak to Props1
* `#1101 <https://github.com/CoolProp/CoolProp/pull/1101>`_ : Add error handling to some functions, see #1096
* `#1100 <https://github.com/CoolProp/CoolProp/pull/1100>`_ : Allow cmake properly build SMath wrapper
* `#1097 <https://github.com/CoolProp/CoolProp/pull/1097>`_ : Set error string in get_parameter_information_string() and fix SMath wrapper : fixes #1096
* `#1093 <https://github.com/CoolProp/CoolProp/pull/1093>`_ : Revert part of 763d4ce to solve #1077

6.0.0
-----

New features:

* MathCAD wrapper working again (thanks to Jeff Henning)
* Added binary interaction parameters for more than 400 mixtures 
* Added a cubic backend supporting PR and SRK for some calculations
* Added new non-iterative viscosity model for a few refrigerants (especially R32 and R245fa)
* Implemented EOS for HCl, D4, ethylene oxide, and dichloroethane from M. Thol
* ... and a lot of little bugfixes (see issues)

Contributors to this release:
ibell, jowr, henningjp, bilderbuchi, dinojr, mapipolo, Mol3culo, stefann82, arashsk, pypamart, milesabarr, wahlenkus, saha84, EmiCas, Heathckliff, Tom0310, dizzux, davideziviani, paarfi

Issues Closed:

* `#1056 <http://github.com/CoolProp/CoolProp/issues/1056>`_ : Added "set_reference_state" wrapper for Mathcad and Updated Example Worksheets
* `#1053 <http://github.com/CoolProp/CoolProp/issues/1053>`_ : Align Tmax with REFPROP values
* `#1049 <http://github.com/CoolProp/CoolProp/issues/1049>`_ : apply_simple_mixing_rule should be implemented for HEOS instances
* `#1048 <http://github.com/CoolProp/CoolProp/issues/1048>`_ : Calling set_binary_interaction_double on AbstractState instance has no effect
* `#1047 <http://github.com/CoolProp/CoolProp/issues/1047>`_ : Mathcad Wrapper Updates for CoolProp 5.x and 6
* `#1044 <http://github.com/CoolProp/CoolProp/issues/1044>`_ : Manylinux build integration
* `#1041 <http://github.com/CoolProp/CoolProp/issues/1041>`_ : Fixed Minor MSVC Compiler Warnings
* `#1034 <http://github.com/CoolProp/CoolProp/issues/1034>`_ : Strange behaviour of densities at critical point
* `#1033 <http://github.com/CoolProp/CoolProp/issues/1033>`_ : Python builder issues
* `#1032 <http://github.com/CoolProp/CoolProp/issues/1032>`_ : Rewrite mixture derivatives tests to use new format
* `#1031 <http://github.com/CoolProp/CoolProp/issues/1031>`_ : Fixes STRING conflict between Mathcad library and cppformat
* `#1030 <http://github.com/CoolProp/CoolProp/issues/1030>`_ : Add pass-throughs for testing derivatives
* `#1029 <http://github.com/CoolProp/CoolProp/issues/1029>`_ : Sphinx builder
* `#1028 <http://github.com/CoolProp/CoolProp/issues/1028>`_ : ALTERNATIVE_REFPROP_PATH ignored for predefined mixtures
* `#1026 <http://github.com/CoolProp/CoolProp/issues/1026>`_ : Add REFPROP version to REFPROP comparison script
* `#1025 <http://github.com/CoolProp/CoolProp/issues/1025>`_ : Phase envelopes construction failing for example in docs 
* `#1024 <http://github.com/CoolProp/CoolProp/issues/1024>`_ : VLE calcs failing for SRK & PR backends
* `#1023 <http://github.com/CoolProp/CoolProp/issues/1023>`_ : AbstractState.update fails for mixtures containing specific refrigerants using REFPROP backend
* `#1020 <http://github.com/CoolProp/CoolProp/issues/1020>`_ : Add target_link_libraries to CMakeLists.txt
* `#1014 <http://github.com/CoolProp/CoolProp/issues/1014>`_ : Figure out how to make coolprop static library a clean cmake dependency
* `#1012 <http://github.com/CoolProp/CoolProp/issues/1012>`_ : Residual Helmholtz energy not work
* `#1011 <http://github.com/CoolProp/CoolProp/issues/1011>`_ : Update references
* `#1010 <http://github.com/CoolProp/CoolProp/issues/1010>`_ : Derivative of residual Helmholtz energy with delta
* `#1009 <http://github.com/CoolProp/CoolProp/issues/1009>`_ : Can't compute densities at the triple point
* `#1007 <http://github.com/CoolProp/CoolProp/issues/1007>`_ : 'error: key [Ar] was not found in string_to_index'
* `#1006 <http://github.com/CoolProp/CoolProp/issues/1006>`_ : Use c++14 when building on MINGW
* `#1005 <http://github.com/CoolProp/CoolProp/issues/1005>`_ : Derivative of the saturation enthalpy cair_sat = d(hsat)/dT
* `#1003 <http://github.com/CoolProp/CoolProp/issues/1003>`_ : Fix bug in Chung estimation model
* `#1002 <http://github.com/CoolProp/CoolProp/issues/1002>`_ : Add python 3.5 wheel
* `#1001 <http://github.com/CoolProp/CoolProp/issues/1001>`_ : DmolarP broken for Air
* `#1000 <http://github.com/CoolProp/CoolProp/issues/1000>`_ : Fix setting of BIP function
* `#999 <http://github.com/CoolProp/CoolProp/issues/999>`_ : Abbreviate all journal names
* `#998 <http://github.com/CoolProp/CoolProp/issues/998>`_ : Refine phase envelope better on liquid side
* `#997 <http://github.com/CoolProp/CoolProp/issues/997>`_ : Abbreviate IECR in CoolProp reference
* `#996 <http://github.com/CoolProp/CoolProp/issues/996>`_ : Update references for R245fa and R1234ze(E)
* `#995 <http://github.com/CoolProp/CoolProp/issues/995>`_ : Check double_equal in CPnumerics.h
* `#994 <http://github.com/CoolProp/CoolProp/issues/994>`_ : Find a way to simplify includes
* `#993 <http://github.com/CoolProp/CoolProp/issues/993>`_ : Test/Add example for DLL calling from C
* `#992 <http://github.com/CoolProp/CoolProp/issues/992>`_ : Fix reference for R1234ze(E) again
* `#987 <http://github.com/CoolProp/CoolProp/issues/987>`_ : Multiple EOS paper refs run together
* `#986 <http://github.com/CoolProp/CoolProp/issues/986>`_ : Air lookup in Excel v5.1.2
* `#982 <http://github.com/CoolProp/CoolProp/issues/982>`_ : Reorganize CoolPropTools.h into smaller modules
* `#981 <http://github.com/CoolProp/CoolProp/issues/981>`_ : Saturation states
* `#976 <http://github.com/CoolProp/CoolProp/issues/976>`_ : Add high-level functions to Julia wrapper
* `#975 <http://github.com/CoolProp/CoolProp/issues/975>`_ : Correct get_parameter_information_string, fixes #974
* `#973 <http://github.com/CoolProp/CoolProp/issues/973>`_ : Remove warnings when using Julia 0.4 realease
* `#971 <http://github.com/CoolProp/CoolProp/issues/971>`_ : Fix bug in PhaseEnvelopeRoutines::evaluate
* `#970 <http://github.com/CoolProp/CoolProp/issues/970>`_ : Props1SI function missing in Mathematica wrapper on OSX
* `#968 <http://github.com/CoolProp/CoolProp/issues/968>`_ : Update index.rst
* `#967 <http://github.com/CoolProp/CoolProp/issues/967>`_ : SO2 ancillaries broken
* `#964 <http://github.com/CoolProp/CoolProp/issues/964>`_ : Update index.rst
* `#963 <http://github.com/CoolProp/CoolProp/issues/963>`_ : Update index.rst
* `#962 <http://github.com/CoolProp/CoolProp/issues/962>`_ : Update sample.sce
* `#960 <http://github.com/CoolProp/CoolProp/issues/960>`_ : Update index.rst
* `#953 <http://github.com/CoolProp/CoolProp/issues/953>`_ : Remap CoolPropDbl to double
* `#952 <http://github.com/CoolProp/CoolProp/issues/952>`_ : Switch string formatting to use the cppformat library; see #907
* `#951 <http://github.com/CoolProp/CoolProp/issues/951>`_ : Allow gibbs as input to first_partial_deriv()
* `#950 <http://github.com/CoolProp/CoolProp/issues/950>`_ : Wrong units for residual entropy
* `#949 <http://github.com/CoolProp/CoolProp/issues/949>`_ : Fix {} in bibtex to protect title capitalization
* `#948 <http://github.com/CoolProp/CoolProp/issues/948>`_ : Update reference for  EOS-CG
* `#947 <http://github.com/CoolProp/CoolProp/issues/947>`_ : Add Fij to REFPROPMixtureBackend::get_binary_interaction_double
* `#945 <http://github.com/CoolProp/CoolProp/issues/945>`_ : Add EOS for R245ca
* `#944 <http://github.com/CoolProp/CoolProp/issues/944>`_ : Update reference for R1233ze(E)
* `#941 <http://github.com/CoolProp/CoolProp/issues/941>`_ : CoolProp returns same value for p_critical and p_triple for R503
* `#937 <http://github.com/CoolProp/CoolProp/issues/937>`_ : Allow ability to get refprop version
* `#934 <http://github.com/CoolProp/CoolProp/issues/934>`_ : Memory access violation on mixture update at very low pressures using tabular backend
* `#933 <http://github.com/CoolProp/CoolProp/issues/933>`_ : ValueError: Bad phase to solver_rho_Tp_SRK (CoolProp 5.1.2)
* `#932 <http://github.com/CoolProp/CoolProp/issues/932>`_ : Fix EOS reference for oxygen
* `#931 <http://github.com/CoolProp/CoolProp/issues/931>`_ : Remap CoolPropDbl to double permanently
* `#930 <http://github.com/CoolProp/CoolProp/issues/930>`_ : Phase envelopes should be begin at much lower pressure
* `#929 <http://github.com/CoolProp/CoolProp/issues/929>`_ : PT should start with Halley's method everywhere
* `#928 <http://github.com/CoolProp/CoolProp/issues/928>`_ : Add EOS for HCl, D4, ethylene oxide, and dichloroethane
* `#927 <http://github.com/CoolProp/CoolProp/issues/927>`_ : Add ability to use Henry's Law to get guesses for liquid phase composition
* `#926 <http://github.com/CoolProp/CoolProp/issues/926>`_ : hydrogen formula is wrong
* `#925 <http://github.com/CoolProp/CoolProp/issues/925>`_ : Fix HS inputs 
* `#921 <http://github.com/CoolProp/CoolProp/issues/921>`_ : Tabular calcs with mixtures often return Dew T< Bubble T using PQ input pair
* `#920 <http://github.com/CoolProp/CoolProp/issues/920>`_ : Can't find temperature at pressure and entropy
* `#917 <http://github.com/CoolProp/CoolProp/issues/917>`_ : Fix errors in docs
* `#907 <http://github.com/CoolProp/CoolProp/issues/907>`_ : Replace string formatting with C++ format library
* `#905 <http://github.com/CoolProp/CoolProp/issues/905>`_ : Using conda recipes
* `#885 <http://github.com/CoolProp/CoolProp/issues/885>`_ : Duplicate critical points found 
* `#854 <http://github.com/CoolProp/CoolProp/issues/854>`_ : Coolprop R448A, R449A or R450A
* `#816 <http://github.com/CoolProp/CoolProp/issues/816>`_ : Issue with viscosity of R245FA
* `#808 <http://github.com/CoolProp/CoolProp/issues/808>`_ : Implement tangent plane distance
* `#665 <http://github.com/CoolProp/CoolProp/issues/665>`_ : Viscosity convergence issue
* `#279 <http://github.com/CoolProp/CoolProp/issues/279>`_ : Rebuild MathCAD wrapper with v5 support
* `#186 <http://github.com/CoolProp/CoolProp/issues/186>`_ : Convert cubics to HE

Pull Requests merged:

* `#1062 <http://github.com/CoolProp/CoolProp/pull/1062>`_ : Export first_partial_deriv, see #946 #1062
* `#1056 <http://github.com/CoolProp/CoolProp/pull/1056>`_ : Added "set_reference_state" wrapper for Mathcad and Updated Example Worksheets
* `#1053 <http://github.com/CoolProp/CoolProp/pull/1053>`_ : Align Tmax with REFPROP values
* `#1047 <http://github.com/CoolProp/CoolProp/pull/1047>`_ : Mathcad Wrapper Updates for CoolProp 5.x and 6
* `#1041 <http://github.com/CoolProp/CoolProp/pull/1041>`_ : Fixed Minor MSVC Compiler Warnings
* `#1031 <http://github.com/CoolProp/CoolProp/pull/1031>`_ : Fixes STRING conflict between Mathcad library and cppformat
* `#1020 <http://github.com/CoolProp/CoolProp/pull/1020>`_ : Add target_link_libraries to CMakeLists.txt
* `#982 <http://github.com/CoolProp/CoolProp/pull/982>`_ : Reorganize CoolPropTools.h into smaller modules
* `#981 <http://github.com/CoolProp/CoolProp/pull/981>`_ : Saturation states
* `#976 <http://github.com/CoolProp/CoolProp/pull/976>`_ : Add high-level functions to Julia wrapper
* `#975 <http://github.com/CoolProp/CoolProp/pull/975>`_ : Correct get_parameter_information_string, fixes #974
* `#973 <http://github.com/CoolProp/CoolProp/pull/973>`_ : Remove warnings when using Julia 0.4 realease
* `#968 <http://github.com/CoolProp/CoolProp/pull/968>`_ : Update index.rst
* `#964 <http://github.com/CoolProp/CoolProp/pull/964>`_ : Update index.rst
* `#963 <http://github.com/CoolProp/CoolProp/pull/963>`_ : Update index.rst
* `#962 <http://github.com/CoolProp/CoolProp/pull/962>`_ : Update sample.sce
* `#960 <http://github.com/CoolProp/CoolProp/pull/960>`_ : Update index.rst
* `#953 <http://github.com/CoolProp/CoolProp/pull/953>`_ : Remap CoolPropDbl to double
* `#952 <http://github.com/CoolProp/CoolProp/pull/952>`_ : Switch string formatting to use the cppformat library; see #907

5.1.2
-----

New features:

* Android wrapper available
* Javascript interface extended to export AbstractState and some functions
* Fixed a wide range of issues with tables
* ... and a lot of little bugfixes (see issues)

Issues Closed:

* `#914 <http://github.com/CoolProp/CoolProp/issues/914>`_ : Tabular ammonia calc yields very different results using TTSE vs. bicubic, including non-physical and NaN quantities
* `#909 <http://github.com/CoolProp/CoolProp/issues/909>`_ : Fortran wrapper on Win...still unable to run it!
* `#906 <http://github.com/CoolProp/CoolProp/issues/906>`_ : Add DOI for Novec649
* `#904 <http://github.com/CoolProp/CoolProp/issues/904>`_ : Deuterium reference has wrong year
* `#903 <http://github.com/CoolProp/CoolProp/issues/903>`_ : Some BibTeX keys need updating
* `#902 <http://github.com/CoolProp/CoolProp/issues/902>`_ : Trap errors in get_BibTeXKey and throw
* `#901 <http://github.com/CoolProp/CoolProp/issues/901>`_ : Viscosity of some incompressibles off by a factor of 100 and 1000
* `#899 <http://github.com/CoolProp/CoolProp/issues/899>`_ : Cp, Cv, speed_sound cannot be calculated with QT inputs (Q=0 or 1) and tabular backends
* `#897 <http://github.com/CoolProp/CoolProp/issues/897>`_ : Update DEF for new AbstractState functions
* `#896 <http://github.com/CoolProp/CoolProp/issues/896>`_ : Tabular refactor
* `#894 <http://github.com/CoolProp/CoolProp/issues/894>`_ : License on homepage
* `#889 <http://github.com/CoolProp/CoolProp/issues/889>`_ :  MSVCP100.dll and MSVCR100.dll dependency issue...
* `#888 <http://github.com/CoolProp/CoolProp/issues/888>`_ : Multi-output library function
* `#886 <http://github.com/CoolProp/CoolProp/issues/886>`_ : ALTERNATE_REFPROP_PATH ignored in low-level interface
* `#882 <http://github.com/CoolProp/CoolProp/issues/882>`_ : Tabular backends and phase specification
* `#880 <http://github.com/CoolProp/CoolProp/issues/880>`_ : low-level interface MATLAB using shared library
* `#871 <http://github.com/CoolProp/CoolProp/issues/871>`_ : Issues with Cp, Cv, u, and viscosity with QT_INPUTS where Q=0 or 1 (xxx&REFPROP backend)
* `#869 <http://github.com/CoolProp/CoolProp/issues/869>`_ : Fix javascript builder on buildbot
* `#868 <http://github.com/CoolProp/CoolProp/issues/868>`_ : Fix fortran builds on buildbot
* `#865 <http://github.com/CoolProp/CoolProp/issues/865>`_ : Hide tabular generation outputs when debug_level=0
* `#859 <http://github.com/CoolProp/CoolProp/issues/859>`_ : Windows wrapper for Octave not working for v 4.0
* `#853 <http://github.com/CoolProp/CoolProp/issues/853>`_ : Problem with linking shared libraries using Code::Blocks and CoolProp
* `#849 <http://github.com/CoolProp/CoolProp/issues/849>`_ : Tidy up references in online docs
* `#848 <http://github.com/CoolProp/CoolProp/issues/848>`_ : PropsSImulti in Python
* `#845 <http://github.com/CoolProp/CoolProp/issues/845>`_ : Tabular calculations fail with message "Unable to bisect segmented vector slice..."
* `#844 <http://github.com/CoolProp/CoolProp/issues/844>`_ : failure in calculation enthalpy for water
* `#843 <http://github.com/CoolProp/CoolProp/issues/843>`_ : Calling AbstractState.update() using Dmass_P input pair causes stack overflow in tabular backends
* `#842 <http://github.com/CoolProp/CoolProp/issues/842>`_ : Wrong enthalpy calculation for SES36
* `#841 <http://github.com/CoolProp/CoolProp/issues/841>`_ : R1233zd(E) reference
* `#840 <http://github.com/CoolProp/CoolProp/issues/840>`_ : Failure to calculate any state using input pair QT_INPUTS with backend TTSE&REFPROP
* `#838 <http://github.com/CoolProp/CoolProp/issues/838>`_ : Request: implement a configuration variable to specify directory for tabular interpolation data
* `#837 <http://github.com/CoolProp/CoolProp/issues/837>`_ : Exceptions thrown when getting/setting MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB configuration setting
* `#835 <http://github.com/CoolProp/CoolProp/issues/835>`_ : Request: CoolProp.AbstractState.first_saturation_deriv wrapped in CoolPropLib.h
* `#831 <http://github.com/CoolProp/CoolProp/issues/831>`_ : Predefined mixtures fail for BICUBIC&REFPROP backend
* `#826 <http://github.com/CoolProp/CoolProp/issues/826>`_ : Unit conversion problem somewhere in Bicubic backend for enthalpy
* `#825 <http://github.com/CoolProp/CoolProp/issues/825>`_ : PQ_with_guesses assumes bubble point
* `#824 <http://github.com/CoolProp/CoolProp/issues/824>`_ : C-Sharp Wrapper AbstractState mole_fractions_liquid
* `#823 <http://github.com/CoolProp/CoolProp/issues/823>`_ : Documentation for use of static libraries is unclear
* `#822 <http://github.com/CoolProp/CoolProp/issues/822>`_ : Request: PropsSI Inputs of D and Q
* `#821 <http://github.com/CoolProp/CoolProp/issues/821>`_ : Fix pip command for nightly
* `#820 <http://github.com/CoolProp/CoolProp/issues/820>`_ : Add cmake option to generate Android .so library
* `#819 <http://github.com/CoolProp/CoolProp/issues/819>`_ : Expose phase envelope calculations in javascript
* `#814 <http://github.com/CoolProp/CoolProp/issues/814>`_ : saturated_liquid/vapor_keyed_output for tabular backend
* `#812 <http://github.com/CoolProp/CoolProp/issues/812>`_ : Add ability to retrieve mass fractions
* `#810 <http://github.com/CoolProp/CoolProp/issues/810>`_ : Python builds crash on Windows
* `#809 <http://github.com/CoolProp/CoolProp/issues/809>`_ : Implement fluid_param_string in python
* `#807 <http://github.com/CoolProp/CoolProp/issues/807>`_ : Return all critical points
* `#805 <http://github.com/CoolProp/CoolProp/issues/805>`_ : Coolprop function like Refprop Excel Fluidstring Function for mixtures
* `#804 <http://github.com/CoolProp/CoolProp/issues/804>`_ : Allow disabling parameter estimation in REFPROP
* `#802 <http://github.com/CoolProp/CoolProp/issues/802>`_ : Error with two-phase DT inputs for R134a
* `#800 <http://github.com/CoolProp/CoolProp/issues/800>`_ : Add access to contributions to viscosity and conductivity
* `#799 <http://github.com/CoolProp/CoolProp/issues/799>`_ : Add access to conformal state solver in AbstractState
* `#798 <http://github.com/CoolProp/CoolProp/issues/798>`_ : Add linear and Lorentz-Berthelot mixing rules
* `#796 <http://github.com/CoolProp/CoolProp/issues/796>`_ : Add SATTP guess implementation
* `#795 <http://github.com/CoolProp/CoolProp/issues/795>`_ : Provide swigged MATLAB wrapper code
* `#793 <http://github.com/CoolProp/CoolProp/issues/793>`_ : Set interaction parameters in REFPROP through CoolProp
* `#792 <http://github.com/CoolProp/CoolProp/issues/792>`_ : Allow possibility to set interaction parameters even if the mixture isn't already included
* `#789 <http://github.com/CoolProp/CoolProp/issues/789>`_ : Make sure all phases are calculated correctly for BICUBIC&HEOS backend
* `#788 <http://github.com/CoolProp/CoolProp/issues/788>`_ : Make sure all phases are calculated correctly for HEOS backend
* `#786 <http://github.com/CoolProp/CoolProp/issues/786>`_ : Implement conductivity for pentanes
* `#785 <http://github.com/CoolProp/CoolProp/issues/785>`_ : Implement viscosity for Toluene
* `#784 <http://github.com/CoolProp/CoolProp/issues/784>`_ : Add docs for get/set config functions
* `#783 <http://github.com/CoolProp/CoolProp/issues/783>`_ : Failure in PsychScript
* `#777 <http://github.com/CoolProp/CoolProp/issues/777>`_ : No input passed with PT_INPUTS and tabular backed
* `#776 <http://github.com/CoolProp/CoolProp/issues/776>`_ : Fix docs for IF97 backend
* `#773 <http://github.com/CoolProp/CoolProp/issues/773>`_ : Missing files in LabVIEW wrapper folder or documentation needed
* `#772 <http://github.com/CoolProp/CoolProp/issues/772>`_ : Acentric factor of air
* `#770 <http://github.com/CoolProp/CoolProp/issues/770>`_ : Make clear() overridable / clear Helmholtz cache
* `#769 <http://github.com/CoolProp/CoolProp/issues/769>`_ : Improve docs for second partial derivatives
* `#768 <http://github.com/CoolProp/CoolProp/issues/768>`_ : Fix solver for first criticality contour crossing
* `#767 <http://github.com/CoolProp/CoolProp/issues/767>`_ : When tracing criticality contour, make sure that delta is always increasing
* `#764 <http://github.com/CoolProp/CoolProp/issues/764>`_ : Add `calc_speed_sound` to tabular backend
* `#763 <http://github.com/CoolProp/CoolProp/issues/763>`_ : Add and implement all phase functions to tabular backends
* `#762 <http://github.com/CoolProp/CoolProp/issues/762>`_ : Temperature with `HmassP_INPUTS` with twophase fluid and tabular
* `#761 <http://github.com/CoolProp/CoolProp/issues/761>`_ : Add auto-generated docs for configuration variables
* `#760 <http://github.com/CoolProp/CoolProp/issues/760>`_ : Add `surface tension` to tabular backend
* `#759 <http://github.com/CoolProp/CoolProp/issues/759>`_ : Add comprehensive docs for REFPROP interface
* `#757 <http://github.com/CoolProp/CoolProp/issues/757>`_ : Cannot evaluate PT (or PH?) below p_triple
* `#756 <http://github.com/CoolProp/CoolProp/issues/756>`_ : HAPropsSI does not converge for T= 299.8 K
* `#754 <http://github.com/CoolProp/CoolProp/issues/754>`_ : Failure with sat derivative with QT and tables
* `#753 <http://github.com/CoolProp/CoolProp/issues/753>`_ : Relative humidity calculation error
* `#751 <http://github.com/CoolProp/CoolProp/issues/751>`_ : D-P is far slower than it should be
* `#750 <http://github.com/CoolProp/CoolProp/issues/750>`_ : Invalid index to calc_first_saturation_deriv in TabularBackends
* `#747 <http://github.com/CoolProp/CoolProp/issues/747>`_ : Plotting example on coolprop.org does not work - potentially related to issue #351
* `#746 <http://github.com/CoolProp/CoolProp/issues/746>`_ : Implement viscosity models for HFO (ECS?)
* `#745 <http://github.com/CoolProp/CoolProp/issues/745>`_ : Undocumented high level interface for saturation derivatives
* `#742 <http://github.com/CoolProp/CoolProp/issues/742>`_ : Expedite the D+Y flash routines
* `#741 <http://github.com/CoolProp/CoolProp/issues/741>`_ : Expedite the single-phase T+Y flash routines
* `#740 <http://github.com/CoolProp/CoolProp/issues/740>`_ : HapropsSI("T", "B", 299.15, "R", 0, "P", 101325) lead to an error
* `#739 <http://github.com/CoolProp/CoolProp/issues/739>`_ : Quality-related updates with tabular backend
* `#738 <http://github.com/CoolProp/CoolProp/issues/738>`_ : TTSE ranges
* `#737 <http://github.com/CoolProp/CoolProp/issues/737>`_ : Missing bib entry IAPWS-SurfaceTension-1994
* `#735 <http://github.com/CoolProp/CoolProp/issues/735>`_ : phase is wrong for water at STP
* `#734 <http://github.com/CoolProp/CoolProp/issues/734>`_ : F is missing from mixture interaction parameters on the web
* `#733 <http://github.com/CoolProp/CoolProp/issues/733>`_ : Typo in excess term in mixture docs
* `#731 <http://github.com/CoolProp/CoolProp/issues/731>`_ : Add EOS for Novec 649 from McLinden
* `#730 <http://github.com/CoolProp/CoolProp/issues/730>`_ : Merge references from paper about CoolProp into main bib file
* `#727 <http://github.com/CoolProp/CoolProp/issues/727>`_ : HapropsSI("T", "B", 299.15, "R", 0, "P", 101325) lead to an error
* `#726 <http://github.com/CoolProp/CoolProp/issues/726>`_ : Improve caching of derivative terms when using mixtures
* `#725 <http://github.com/CoolProp/CoolProp/issues/725>`_ : Implement dipole moment

5.1.1
-----

New features:

* A wrapper for the R language
* Tabular integration with tables from REFPROP only for now
* The Python wrapper is now also available on binstar: https://binstar.org/CoolProp/coolprop
* ... and a lot of little bugfixes (see issues)

Issues Closed:

* `#724 <http://github.com/CoolProp/CoolProp/issues/724>`_ : Gibbs not working as output (mass or molar)
* `#722 <http://github.com/CoolProp/CoolProp/issues/722>`_ : Predefined mixtures crash python
* `#721 <http://github.com/CoolProp/CoolProp/issues/721>`_ : v5.1.1
* `#714 <http://github.com/CoolProp/CoolProp/issues/714>`_ : Possible error in isobaric thermal expansion coefficient
* `#713 <http://github.com/CoolProp/CoolProp/issues/713>`_ : Bicubic backend and first_saturation_deriv
* `#712 <http://github.com/CoolProp/CoolProp/issues/712>`_ : Expose saturation derivatives from PropsSI [wishlist]
* `#708 <http://github.com/CoolProp/CoolProp/issues/708>`_ : CoolPropsetup.m needs to be installed
* `#707 <http://github.com/CoolProp/CoolProp/issues/707>`_ : conda builds
* `#703 <http://github.com/CoolProp/CoolProp/issues/703>`_ : 2/ HapropsSI ( "T" , "B" , ValueB, "W" , 0 , "P" , 101325) lead to an error
* `#702 <http://github.com/CoolProp/CoolProp/issues/702>`_ : 1 : HapropsSI ( "T" , "H" , ValueH, "W" , 0 , "P" , 101325) lead to an error
* `#700 <http://github.com/CoolProp/CoolProp/issues/700>`_ : If git is not found, still compile properly
* `#699 <http://github.com/CoolProp/CoolProp/issues/699>`_ : Fugacity using Python wrapper
* `#697 <http://github.com/CoolProp/CoolProp/issues/697>`_ : Get State (old-style) class working with predefined mixtures
* `#696 <http://github.com/CoolProp/CoolProp/issues/696>`_ : cp0 broken for tabular backends
* `#695 <http://github.com/CoolProp/CoolProp/issues/695>`_ : Problem with reference state
* `#691 <http://github.com/CoolProp/CoolProp/issues/691>`_ : variable names for second_partial_deriv
* `#688 <http://github.com/CoolProp/CoolProp/issues/688>`_ : PropsSI in saturation region
* `#685 <http://github.com/CoolProp/CoolProp/issues/685>`_ : Problem with Hazard output
* `#684 <http://github.com/CoolProp/CoolProp/issues/684>`_ : some problem and questions for calc in Excel
* `#681 <http://github.com/CoolProp/CoolProp/issues/681>`_ : Mix call failure after release update
* `#680 <http://github.com/CoolProp/CoolProp/issues/680>`_ : Tabular backend data range too small for (P,H) inputs and R245fa
* `#675 <http://github.com/CoolProp/CoolProp/issues/675>`_ : Get consistency plots working with Tabular backends
* `#674 <http://github.com/CoolProp/CoolProp/issues/674>`_ : QT inputs do not work for Tabular backends
* `#673 <http://github.com/CoolProp/CoolProp/issues/673>`_ : Mass-based saturation derivatives not supported
* `#672 <http://github.com/CoolProp/CoolProp/issues/672>`_ : Tabular methods returns hmolar for smolar for saturation
* `#671 <http://github.com/CoolProp/CoolProp/issues/671>`_ : MATLAB on OSX cannot load REFPROP
* `#670 <http://github.com/CoolProp/CoolProp/issues/670>`_ : Low-Level interfacing with MATLAB
* `#668 <http://github.com/CoolProp/CoolProp/issues/668>`_ : R wrapper
* `#664 <http://github.com/CoolProp/CoolProp/issues/664>`_ : Re-enable triple point for REFPROP backend for mixtures
* `#663 <http://github.com/CoolProp/CoolProp/issues/663>`_ : Vapor mass quality = 1 generates error for pseudo-pures
* `#662 <http://github.com/CoolProp/CoolProp/issues/662>`_ : Write function to determine phase after an update with PT and a guess for rho
* `#661 <http://github.com/CoolProp/CoolProp/issues/661>`_ : Predefined mixtures not working properly with Tabular backends
* `#660 <http://github.com/CoolProp/CoolProp/issues/660>`_ : T,X and PS, PD, PU not working with BICUBIC, but does with TTSE
* `#659 <http://github.com/CoolProp/CoolProp/issues/659>`_ : Add "PIP" as parameter
* `#658 <http://github.com/CoolProp/CoolProp/issues/658>`_ : Implement PIP for REFPROP
* `#657 <http://github.com/CoolProp/CoolProp/issues/657>`_ : Describe how to call REFPROP
* `#654 <http://github.com/CoolProp/CoolProp/issues/654>`_ : Add ability to calculate Ideal curves
* `#653 <http://github.com/CoolProp/CoolProp/issues/653>`_ : Implement update_with_guesses for P,T for REFPROP backend
* `#652 <http://github.com/CoolProp/CoolProp/issues/652>`_ : Implement solver for "true" critical point using REFPROP
* `#650 <http://github.com/CoolProp/CoolProp/issues/650>`_ : MATLAB examples not on website
* `#648 <http://github.com/CoolProp/CoolProp/issues/648>`_ : Link to examples broken
* `#647 <http://github.com/CoolProp/CoolProp/issues/647>`_ : Implement the new REFPROP header file and make necessary changes
* `#646 <http://github.com/CoolProp/CoolProp/issues/646>`_ : Add B,C virial coefficients for REFPROP backend
* `#645 <http://github.com/CoolProp/CoolProp/issues/645>`_ : PQ_INPUTS don't work with TTSE backend
* `#644 <http://github.com/CoolProp/CoolProp/issues/644>`_ : Get first_two_phase_deriv working with Tabular backends
* `#641 <http://github.com/CoolProp/CoolProp/issues/641>`_ : Install psyrc file
* `#640 <http://github.com/CoolProp/CoolProp/issues/640>`_ : Expose saturation_ancillary equation through python
* `#639 <http://github.com/CoolProp/CoolProp/issues/639>`_ : Incorrect error when non two-phase inputs to two-phase deriv
* `#638 <http://github.com/CoolProp/CoolProp/issues/638>`_ : Heavy Water Viscosity Unavailable
* `#636 <http://github.com/CoolProp/CoolProp/issues/636>`_ : Error surface tension in CoolProp v5.1.0
* `#635 <http://github.com/CoolProp/CoolProp/issues/635>`_ : Implement first_saturation_deriv for TTSE/BICUBIC
* `#631 <http://github.com/CoolProp/CoolProp/issues/631>`_ : Methane conductivity
* `#630 <http://github.com/CoolProp/CoolProp/issues/630>`_ : Make HS use DH rather than PH
* `#629 <http://github.com/CoolProp/CoolProp/issues/629>`_ : Handle PT inputs around saturation in a better way with BICUBIC
* `#628 <http://github.com/CoolProp/CoolProp/issues/628>`_ : Dry air enthalpy
* `#627 <http://github.com/CoolProp/CoolProp/issues/627>`_ : Test that H and S are the same for all the state points
* `#626 <http://github.com/CoolProp/CoolProp/issues/626>`_ : Improve docs for low-level interface
* `#622 <http://github.com/CoolProp/CoolProp/issues/622>`_ : TTSE fails around saturated liquid
* `#617 <http://github.com/CoolProp/CoolProp/issues/617>`_ : Block Tabular backend use with PropsSI somehow

5.1.0
-----

New features:

* Tabular interpolation using TTSE or Bicubic interpolation (http://www.coolprop.org/coolprop/Tabular.html)
* Equation of state for heavy water
* Added IF97 backend for industrial formulation for properties of pure water
* Lots of little bugfixes (see issues)

Issues Closed:

* `#624 <http://github.com/CoolProp/CoolProp/issues/624>`_ : Stability in two-phase region
* `#621 <http://github.com/CoolProp/CoolProp/issues/621>`_ : TTSE Input Param (Water)
* `#620 <http://github.com/CoolProp/CoolProp/issues/620>`_ : TTSE Problem (Water)
* `#618 <http://github.com/CoolProp/CoolProp/issues/618>`_ : H,S not working for pseudo-pure
* `#615 <http://github.com/CoolProp/CoolProp/issues/615>`_ : Ammonia T-P saturation calculation deviation
* `#614 <http://github.com/CoolProp/CoolProp/issues/614>`_ : Typos in parameter descriptions.
* `#612 <http://github.com/CoolProp/CoolProp/issues/612>`_ : Added missing cell "Input/Output" for enthalpy row.
* `#611 <http://github.com/CoolProp/CoolProp/issues/611>`_ : Splined Output Doubt
* `#609 <http://github.com/CoolProp/CoolProp/issues/609>`_ : Some Windows builds fail (error removing non-existent directory)
* `#608 <http://github.com/CoolProp/CoolProp/issues/608>`_ : MinGW builds fail
* `#605 <http://github.com/CoolProp/CoolProp/issues/605>`_ : CMake changes
* `#602 <http://github.com/CoolProp/CoolProp/issues/602>`_ : TTSE fails for two-phase H,P with heavy water
* `#601 <http://github.com/CoolProp/CoolProp/issues/601>`_ : Benzene conductivity bibtex is wrong
* `#599 <http://github.com/CoolProp/CoolProp/issues/599>`_ : Something is messed up with water properties
* `#595 <http://github.com/CoolProp/CoolProp/issues/595>`_ : add DOIs to bibliography
* `#591 <http://github.com/CoolProp/CoolProp/issues/591>`_ : Request for extension: table of quantities in the documentation for HAPropsSI like for PropsSI
* `#588 <http://github.com/CoolProp/CoolProp/issues/588>`_ : matplotlib and numpy should not be explicit dependencies
* `#586 <http://github.com/CoolProp/CoolProp/issues/586>`_ : HAProps humidity ratio calculation issue
* `#585 <http://github.com/CoolProp/CoolProp/issues/585>`_ : HAProps at low humidity ratio
* `#584 <http://github.com/CoolProp/CoolProp/issues/584>`_ : [Tabular] pure fluid AbstractState returns the wrong mole fractions
* `#583 <http://github.com/CoolProp/CoolProp/issues/583>`_ : Development docs only available on dreamhosters
* `#579 <http://github.com/CoolProp/CoolProp/issues/579>`_ : Issue with Excel Wrapper for Coolprop for OS X Excel 2011
* `#578 <http://github.com/CoolProp/CoolProp/issues/578>`_ : Update examples to show how to call TTSE and BICUBIC backends
* `#577 <http://github.com/CoolProp/CoolProp/issues/577>`_ : Unicode characters in bibtex not being escaped properly
* `#575 <http://github.com/CoolProp/CoolProp/issues/575>`_ : Phase envelopes should be able to be constructed for pure fluids too
* `#574 <http://github.com/CoolProp/CoolProp/issues/574>`_ : Methane (and pentane) transport properties
* `#573 <http://github.com/CoolProp/CoolProp/issues/573>`_ : Bug in derivatives from Matlab
* `#570 <http://github.com/CoolProp/CoolProp/issues/570>`_ : Implement EOS for heavy water
* `#569 <http://github.com/CoolProp/CoolProp/issues/569>`_ : REFPROP SPLNval for rhomolar_vap wrong
* `#568 <http://github.com/CoolProp/CoolProp/issues/568>`_ : Reference of state not working for Refprop backend
* `#567 <http://github.com/CoolProp/CoolProp/issues/567>`_ : Add IF97 Backend
* `#566 <http://github.com/CoolProp/CoolProp/issues/566>`_ : Retrieve phase envelopes from REFPROP using SPLNVAL function
* `#564 <http://github.com/CoolProp/CoolProp/issues/564>`_ : Molecular Formulas as Trivial Property
* `#562 <http://github.com/CoolProp/CoolProp/issues/562>`_ : Add docs about how to set the reference state
* `#556 <http://github.com/CoolProp/CoolProp/issues/556>`_ : [Tabular] Saturation curves for mixtures
* `#555 <http://github.com/CoolProp/CoolProp/issues/555>`_ : [Tabular] Re-enable the PHI0dll function for REFPROP
* `#552 <http://github.com/CoolProp/CoolProp/issues/552>`_ : IsFluidType function
* `#549 <http://github.com/CoolProp/CoolProp/issues/549>`_ : Implement up to 4th order derivatives of all Helmholtz terms (except SAFT)
* `#548 <http://github.com/CoolProp/CoolProp/issues/548>`_ : Problem with HAPropsSI
* `#546 <http://github.com/CoolProp/CoolProp/issues/546>`_ : Small speed enhancement for Julia wrapper
* `#541 <http://github.com/CoolProp/CoolProp/issues/541>`_ : Update CoolProp.jl
* `#540 <http://github.com/CoolProp/CoolProp/issues/540>`_ : Update CoolProp.jl
* `#539 <http://github.com/CoolProp/CoolProp/issues/539>`_ : Add SATTP to REFPROP wrapper
* `#537 <http://github.com/CoolProp/CoolProp/issues/537>`_ : [Tabular] rebuild tables if limits (especially enthalpies) have shifted
* `#536 <http://github.com/CoolProp/CoolProp/issues/536>`_ : Add low level interface to Julia wrapper as discussed in #534 + Fixes #497
* `#535 <http://github.com/CoolProp/CoolProp/issues/535>`_ : When using high-level wrapper of low-level interface, errors don't bubble properly
* `#534 <http://github.com/CoolProp/CoolProp/issues/534>`_ : Add error handling to Julia's wrapper
* `#532 <http://github.com/CoolProp/CoolProp/issues/532>`_ : More Coverity cleanups
* `#530 <http://github.com/CoolProp/CoolProp/issues/530>`_ : When reference state is changed, reducing/critical and hs_anchor states need to be changed
* `#529 <http://github.com/CoolProp/CoolProp/issues/529>`_ : First bunch of Coverity Scan static analysis warning fixes
* `#528 <http://github.com/CoolProp/CoolProp/issues/528>`_ : PQ Flash Failure for CO2+Water
* `#527 <http://github.com/CoolProp/CoolProp/issues/527>`_ : Silence all output to screen when building phase envelopes
* `#526 <http://github.com/CoolProp/CoolProp/issues/526>`_ : When building phase envelopes, stop when the composition is almost pure
* `#524 <http://github.com/CoolProp/CoolProp/issues/524>`_ : set_reference_state does not create expected output
* `#523 <http://github.com/CoolProp/CoolProp/issues/523>`_ : error: thermal conductivity R32:  _phase is unknown
* `#522 <http://github.com/CoolProp/CoolProp/issues/522>`_ : [Tabular] Implement solver when one of the inputs is not a native input
* `#521 <http://github.com/CoolProp/CoolProp/issues/521>`_ : [Tabular] Fix derivatives, and c_p
* `#520 <http://github.com/CoolProp/CoolProp/issues/520>`_ : [Tabular] Fix transport properties
* `#519 <http://github.com/CoolProp/CoolProp/issues/519>`_ : [Tabular] Fix cells close to the saturation curves
* `#518 <http://github.com/CoolProp/CoolProp/issues/518>`_ : Tabular methods implemented
* `#517 <http://github.com/CoolProp/CoolProp/issues/517>`_ : Isobaric expansion coefficient is not implemented
* `#516 <http://github.com/CoolProp/CoolProp/issues/516>`_ : [Tabular] Actually zip up the tables using zlib
* `#515 <http://github.com/CoolProp/CoolProp/issues/515>`_ : Kill off the CRT deprecate warning (#512)
* `#513 <http://github.com/CoolProp/CoolProp/issues/513>`_ : Primitive structures simplification attempt 2
* `#512 <http://github.com/CoolProp/CoolProp/issues/512>`_ : Kill off the CRT deprecate warning
* `#511 <http://github.com/CoolProp/CoolProp/issues/511>`_ : Python version should be 5.1.0dev, not just 5.1.0
* `#508 <http://github.com/CoolProp/CoolProp/issues/508>`_ : Add a ways of using the shared_ptr directly through shared library
* `#507 <http://github.com/CoolProp/CoolProp/issues/507>`_ : Add possibility to disable a backend at compile-time
* `#506 <http://github.com/CoolProp/CoolProp/issues/506>`_ : [Tabular] Add docs for TTSE and bicubic usage
* `#497 <http://github.com/CoolProp/CoolProp/issues/497>`_ : Julia and C++ Low Level Interface for faster Computation
* `#490 <http://github.com/CoolProp/CoolProp/issues/490>`_ : Add partial pressure of water as an output in HAPropsSI
* `#481 <http://github.com/CoolProp/CoolProp/issues/481>`_ : A bug is found when pressure approximates Critical Pressure for Air
* `#455 <http://github.com/CoolProp/CoolProp/issues/455>`_ : HS Inputs in PropsSI function working in two-phase region?
* `#297 <http://github.com/CoolProp/CoolProp/issues/297>`_ : Call matlab script from command line, with no window, catching errors, and never going interactive
* `#296 <http://github.com/CoolProp/CoolProp/issues/296>`_ : Update examples for v5
* `#262 <http://github.com/CoolProp/CoolProp/issues/262>`_ : Re-implement tabular methods
* `#43 <http://github.com/CoolProp/CoolProp/issues/43>`_ : [Tabular] Warn about tabular folder size

5.0.8
-----

New features:

* Added a Smath Studio native wrapper (thanks to Mike Kaganski for all his help)
* Lots of little cleanups to the code (thanks to Mike Kaganski)

Issues Closed:

* `#510 <http://github.com/CoolProp/CoolProp/issues/510>`_ : const, ref and iterator optimization
* `#509 <http://github.com/CoolProp/CoolProp/issues/509>`_ : Exceptions restructured
* `#505 <http://github.com/CoolProp/CoolProp/issues/505>`_ : AbstractState in python should implement phase() function
* `#504 <http://github.com/CoolProp/CoolProp/issues/504>`_ : More ref args
* `#503 <http://github.com/CoolProp/CoolProp/issues/503>`_ : Add compressibility factor for humid air
* `#502 <http://github.com/CoolProp/CoolProp/issues/502>`_ : thread_local broken on OSX
* `#501 <http://github.com/CoolProp/CoolProp/issues/501>`_ : thread_local: one more (hopefully portable) attempt
* `#500 <http://github.com/CoolProp/CoolProp/issues/500>`_ : Fix directory size calculations
* `#499 <http://github.com/CoolProp/CoolProp/issues/499>`_ : Longdouble remap
* `#498 <http://github.com/CoolProp/CoolProp/issues/498>`_ : HAProp - Conductivity & Viscosity
* `#496 <http://github.com/CoolProp/CoolProp/issues/496>`_ : Implement checking of directory size
* `#495 <http://github.com/CoolProp/CoolProp/issues/495>`_ : CoolPropDbl
* `#493 <http://github.com/CoolProp/CoolProp/issues/493>`_ : Avoid copying of parameters; some fixes for _HAPropsSI_inputs
* `#492 <http://github.com/CoolProp/CoolProp/issues/492>`_ : Add docs for Low-Level Interface
* `#488 <http://github.com/CoolProp/CoolProp/issues/488>`_ : Some more static analyser warning fixes
* `#487 <http://github.com/CoolProp/CoolProp/issues/487>`_ : Cannot use REFPROP to get reducing state variables
* `#485 <http://github.com/CoolProp/CoolProp/issues/485>`_ : Rewrite HAPropsSI to call _HAPropsSI
* `#484 <http://github.com/CoolProp/CoolProp/issues/484>`_ : Kill off all warnings in 64-bit compilation
* `#483 <http://github.com/CoolProp/CoolProp/issues/483>`_ : Problems noted by VS2013 static analysis
* `#479 <http://github.com/CoolProp/CoolProp/issues/479>`_ : RelativeHumidity simplification
* `#478 <http://github.com/CoolProp/CoolProp/issues/478>`_ : Julia 0.3 wrapper
* `#476 <http://github.com/CoolProp/CoolProp/issues/476>`_ : buildbot failure messages don't have the correct URL
* `#473 <http://github.com/CoolProp/CoolProp/issues/473>`_ : Wrapper for Julia 0.3
* `#472 <http://github.com/CoolProp/CoolProp/issues/472>`_ : Fix potential buffer overflow with get_parameter_information_string
* `#471 <http://github.com/CoolProp/CoolProp/issues/471>`_ : Document which inputs are possible in Props1SI
* `#470 <http://github.com/CoolProp/CoolProp/issues/470>`_ : Consider evaluating water at Tdb,p for transport properties in humid air
* `#469 <http://github.com/CoolProp/CoolProp/issues/469>`_ : Initialize fluids in HAProps_Aux
* `#468 <http://github.com/CoolProp/CoolProp/issues/468>`_ : Sanitize internal code in HAPropsSI
* `#467 <http://github.com/CoolProp/CoolProp/issues/467>`_ : Cp in HAPropsSI cannot be calculated in 5.0.7
* `#466 <http://github.com/CoolProp/CoolProp/issues/466>`_ : Prandtl number cannot be returned directly


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