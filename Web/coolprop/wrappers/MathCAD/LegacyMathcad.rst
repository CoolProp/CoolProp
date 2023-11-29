.. _mathcad15:

***********************************************************
CoolProp Wrapper for Legacy Mathcad (32-bit) - Discontinued
***********************************************************

| By Scott Polak and Ian Bell, 2013. 
| Updated by Jeff Henning, 2016.

Pre-compiled binaries (up to CoolProp v6.4.1)
=============================================

Pre-compiled binaries can be downloaded from :sfdownloads:`MathCAD15`.  These binaries are no longer generated as of CoolProp version 6.4.2 and there are no nightly builds.  If you have a perpetual Legacy Mathcad license, you can attempt to build your own User-Compiled Binary using the instructions further below.

To Use
------

* Copy CoolPropMathcadWrapper.dll to C:\\Program Files (x86)\\Mathcad\\Mathcad 15\\userefi 
    
* Copy CoolProp_EN.xml to C:\\Program Files (x86)\\Mathcad\\Mathcad 15\\doc\\funcdoc 
    
* Open the CoolPropFluidProperties.xmcdz file in Mathcad, all CoolProp functions should evaluate properly. If not, press <Ctrl>-F9 to force recalculation of the entire workbook.

* CoolProp functions can be inserted from the Mathcad Insert Functions panel under the function category: CoolProp.  Input parameters and a brief description of each function will be shown.


User-compiled binaries
======================

Common Requirements
-------------------

* Compilation of the Legacy Mathcad wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

* You will need to have Microsoft Visual Studio 2015 or later installed (Express or Community version is fine).

To Build
--------

* Recursively clone the CoolProp library to a local repository.::

    git clone https://github.com/CoolProp/CoolProp --recursive

* Change directory (cd) to the CoolProp directory you just created::

    cd CoolProp

* Go to the top level CoolProp directory and make a build directory** (something like \build15)::

    mkdir build15 
    cd build15

* Build the makefile using CMake (Note: Mathcad 15 is 32-bit!)::

    cmake .. -DCOOLPROP_MATHCAD15_MODULE=ON 
             -DCOOLPROP_MATHCAD15_ROOT="C:/Program Files (x86)/Mathcad/Mathcad 15"  
             -G "Visual Studio 10 2010" 
             -DCMAKE_VERBOSE_MAKEFILE=ON 

* Make the dynamic library (DLL)::

    cmake --build . --config Release

To Use
------

* Copy CoolProp\\build15\\Release\\CoolPropMathcadWrapper.dll to C:\\Program Files (x86)\\Mathcad\\Mathcad 15\\userefi 
    
* Copy CoolProp\\wrapper\\Mathcad\\CoolProp_EN.xml to C:\\Program Files (x86)\\Mathcad\\Mathcad 15\\doc\\funcdoc 
    
* Open the CoolPropFluidProperties.xmcdz file in Legacy Mathcad, all CoolProp functions should evaluate properly. If not, press <Ctrl>-F9 to force recalculation of the entire workbook.

* CoolProp functions can be inserted from the Mathcad Insert Functions panel under the function category: CoolProp.  Input parameters and a brief description of each function will be shown.