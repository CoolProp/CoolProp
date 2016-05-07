.. _mathcadprime:

********************************************************
CoolProp Wrapper for MathCAD Prime 3.x or later (64-bit)
********************************************************

| Copyright Scott Polak and Ian Bell, 2013
| Updated by Jeff Henning, 2016

Pre-compiled binaries
=====================

Pre-compiled binaries can be downloaded from :sfdownloads:`MathCADPrime`.  Development binaries coming from the buildbot server can be found at :sfnightly:`MathCADPrime`.

To Use
------

* Copy CoolPropMathcadWrapper.dll file to C:\\Program Files\\PTC\\Mathcad Prime 3.1\\Custom Functions

* Open the CoolPropFluidProperties.xmcd file in MathCAD, all CoolProp functions should evaluate properly. If not, press <Ctrl>-F9 to force recalculation of the entire workbook.

User-compiled binaries
======================

Common Requirements
-------------------

* Compilation of the MathCAD 15 wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

* You will need to have Visual Studio 2008 or later installed.  You will need the professional version of Visual Studio C++, or at least Visual Studio Express 2013 or later, as Mathcad Prime libraries are 64-bit and require the 64-bit compiler.

To Build
--------

* Recursively clone the CoolProp library to a local repository.::

    git clone https://github.com/CoolProp/CoolProp --recursive

* Change directory (cd) to the CoolProp directory you just created::

    cd CoolProp

* Go to the top level CoolProp directory and make a build directory (something like \buildprime)::

    mkdir buildprime
    cd buildprime

* Build the makefile using CMake (adjust root string for correct version of Prime)::

    cmake .. -DCOOLPROP_PRIME_MODULE=ON 
             -DCOOLPROP_PRIME_ROOT="C:/Program Files/PTC/Mathcad Prime 3.1" 
             -G "Visual Studio 10 2010 Win64" 
             -DCMAKE_VERBOSE_MAKEFILE=ON

    ( *Note: Mathcad Prime is 64-bit, so the 'Win64' option is necessary in the Visual Studio string.* )         
             
* Make the dynamic library (DLL)::

    cmake --build . --config Release

To Use
------

* Copy CoolProp\\buildprime\\Release\\CoolPropMathcadWrapper.dll file to C:\\Program Files\\PTC\\Mathcad Prime 3.1\\Custom Functions

* Open the CoolPropFluidProperties.xmcd file in MathCAD, all CoolProp functions should evaluate properly. If not, press <Ctrl>-F9 to force recalculation of the entire workbook.