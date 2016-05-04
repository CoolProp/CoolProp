CoolProp Wrapper for MathCAD 15
==================================

| Copyright Scott Polak and Ian Bell, 2013. 
| Updated by Jeff Henning, 2016.

Prerequisites
==============

* You will need to have Visual Studio 2008 installed (Express version is fine).  Alternatively newer versions of Microsoft Visual Studio C++ should be fine.

* You will need CMake version 2.8.12 or later from https://cmake.org/download/

* You will need to install Git-SCM for Windows.  You can install this from https://git-for-windows.github.io

* You will need Anaconda/Miniconda Python, which you can get from https://store.continuum.io/cshop/anaconda


To Build
========

* **Go to the top level CoolProp directory and make a build directory** (something like \build15)::

	mkdir build15 
	cd build15

* **Build the makefile using CMake**::

    cmake .. -DCOOLPROP_MATHCAD15_MODULE=ON 
             -DCOOLPROP_MATHCAD15_ROOT="C:/Program Files (x86)/Mathcad/Mathcad 15"  
             -G "Visual Studio 10 2010" 
             -DCMAKE_VERBOSE_MAKEFILE=ON 

* **Make the static library**::

	cmake --build . --config Release


To Use
======

* Copy CoolProp\\build15\\Release\\CoolPropMathcadWrapper.dll to C:\\Program Files (x86)\\Mathcad\\Mathcad 15\\userefi 
	
* Copy CoolProp\\wrapper\\Mathcad\\CoolProp_EN.xml to C:\\Program Files (x86)\\Mathcad\\Mathcad 15\\doc\\funcdoc 
	
* Open the CoolPropFluidProperties.xmcd file in MathCAD, all CoolProp functions should evaluate properly. If not, press <Ctrl>-F9 to force recalculation of the entire workbook.

* CoolProp functions can be inserted from the Mathcad Insert Functions panel under the function category: CoolProp.  Input parameters and a brief description of each function will be shown.

