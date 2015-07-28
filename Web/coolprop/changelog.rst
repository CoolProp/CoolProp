Changelog for CoolProp
======================

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