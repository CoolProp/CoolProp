CoolProp Wrapper for MathCAD Prime 3.x or later (64-bit)
==========================================================

| Copyright Scott Polak and Ian Bell, 2013
| Updated by Jeff Henning, 2016

Prerequisites
==============

* You will need to have Visual Studio 2008 or later installed.  You will need the professional version of Visual Studio C++, or at least Visual Studio Express 2013 or later, as Mathcad Prime libraries are 64-bit and require the 64-bit compiler.

* You will need CMake version 2.8.12 or later from https://cmake.org/download/

* You will need to install Git-SCM for Windows.  You can install this from https://git-for-windows.github.io

* You will need Anaconda/Miniconda Python, which you can get from https://store.continuum.io/cshop/anaconda
	

To Build
========

* **Recursively clone the CoolProp library to a local repository.**::

	git clone https://github.com/CoolProp/CoolProp --recursive

* **Change directory (cd) to the CoolProp directory you just created.**::

	cd CoolProp

* **Go to the top level CoolProp directory and make a build directory** (something like \buildprime)::

    mkdir buildprime
    cd buildprime

* **Build the makefile using CMake (adjust root string for correct version of Prime)**::

    cmake .. -DCOOLPROP_PRIME_MODULE=ON 
             -DCOOLPROP_PRIME_ROOT="C:/Program Files/PTC/Mathcad Prime 3.1" 
             -G "Visual Studio 10 2010 Win64" 
             -DCMAKE_VERBOSE_MAKEFILE=ON

	( *Note: Mathcad Prime is 64-bit, so the 'Win64' option is necessary in the Visual Studio string.* )		 
			 
* **Make the dynamic library (DLL)**::

    cmake --build . --config Release


To Use
======

* Copy CoolProp\\buildprime\\Release\\CoolPropMathcadWrapper.dll file to C:\\Program Files\\PTC\\Mathcad Prime 3.1\\Custom Functions

* Open the CoolPropFluidProperties.xmcd file in MathCAD, all CoolProp functions should evaluate properly. If not, press <Ctrl>-F9 to force recalculation of the entire workbook.

