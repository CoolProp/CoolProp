CoolProp Wrapper for MathCAD Prime 3.x or later
==================================================

| Copyright Scott Polak and Ian Bell, 2013
| Updated by Jeff Henning, 2016

Prerequisites
==============

* You will need to have Visual Studio 2008 installed.  Alternatively newer versions of Microsoft Visual Studio C++ should be fine.  Unfortunately, you will need the professional version of Visual Studio C++ as Mathcad Prime libraries are 64-bit and require the 64-bit compiler; the Express version of Visual Studio will not work.

* You will need CMake version 2.8.12 or later from https://cmake.org/download/

* You will need to install Git-SCM for Windows.  You can install this from https://git-for-windows.github.io

* You will need Anaconda/Miniconda Python, which you can get from https://store.continuum.io/cshop/anaconda
	

To Build
========

* **Go to the top level CoolProp directory and make a build directory** (something like \buildprime)::

	mkdir buildprime
	cd buildprime

* **Build the makefile using CMake (adjust root string for correct version of Prime)**::

	cmake .. -DCOOLPROP_PRIME_MODULE=ON 
	         -DCOOLPROP_PRIME_ROOT="C:/Program Files/PTC/Mathcad Prime 3.1" 
			 -G "Visual Studio 10 2010 Win64" 
			 -DCMAKE_VERBOSE_MAKEFILE=ON
	
* **Make the static library**::

	cmake --build . --config Release


To Use
======

* Copy CoolProp\\buildprime\\Release\\CoolPropMathcadWrapper.dll file to C:\\Program Files\\PTC\\Mathcad Prime 3.1\Custom Functions

* Open the CoolPropFluidProperties.xmcd file in MathCAD, all CoolProp functions should evaluate properly. If not, press <Ctrl>-F9 to force recalculation of the entire workbook.

