CoolProp Wrapper for Mathcad 15 (32-bit)
==========================================

| Copyright Scott Polak and Ian Bell, 2013. 
| Updated by Jeff Henning, 2016. 

 ..  
 
 > NOTE: Legacy Mathcad (through version 15) is no longer supported.  PTC has discontinued distribution and licensing of Legacy Mathcad.  Users should now be using only the top level of this wrapper to compile and develop the Mathcad Prime wrapper.  Users who may still have a perpetual license may still be able to compile this wrapper for Legacy Mathcad, but it is not being actively developed or tested and binaries on SourForge ceased with CoolProp 6.4.1. 
 
There are two ways to get the CoolProp wrapper add-in for Legacy Mathcad.  Download the pre-compiled binary from SourceForge or Compile your own.

Precompiled Binaries
===================
The binary DLL file (up through v6.4.1) can be downloaded and placed in the Mathcad installation directory. 
 
1. Goto the [SourceForge CoolProp](https://sourceforge.net/projects/coolprop/files/) web page. 
 
2. Browse to the desired version (6.4.1 or older). 
 
3. Under the MathCAD15 folder, download all three files 

   * Place the ``CoolPropMathcadWrapper.DLL`` file in the Mathcad 15 installation directory under the ``usrefi`` directory. 

   * Copy the ``CoolProp_EN.xml`` file to the ``doc\\funcdoc`` folder. 

4. Restart Mathcad 15. The add-in, ``CoolPropMathcadWrapper.DLL`` will automatically load the from the ``usrefi`` directory. 

5. Open the ``CoolPropFluidProperties.xmcdz`` file in Legacy Mathcad for examples on how to use the functions. 

6. The CoolProp functions will be added to the Insert Functions menu under the CoolProp category for easy insertion into your workbooks. 


Build Your Own
==============


Prerequisites
--------------

* You will need to have Microsoft Visual Studio 2008 or later installed (Express version is fine).

* You will need CMake version 2.8.12 or later from https://cmake.org/download/

* You will need to install Git-SCM for Windows.  You can install this from https://git-for-windows.github.io

* You will need Anaconda/Miniconda Python, which you can get from https://store.continuum.io/cshop/anaconda


To Build
--------

* **Recursively clone the CoolProp library to a local repository.**::

	git clone https://github.com/CoolProp/CoolProp --recursive

* **Change directory (cd) to the CoolProp directory you just created.**::

	cd CoolProp

* **Go to the top level CoolProp directory and make a build directory** (something like \build15)::

	mkdir build15 
	cd build15

* **Build the makefile using CMake** (Note: Mathcad 15 is 32-bit)::

    cmake .. -DCOOLPROP_MATHCAD15_MODULE=ON 
             -DCOOLPROP_MATHCAD15_ROOT="C:/Program Files (x86)/Mathcad/Mathcad 15"  
             -G "Visual Studio 10 2010" 
             -DCMAKE_VERBOSE_MAKEFILE=ON 

* **Make the dynamic library (DLL)**::

	cmake --build . --config Release


To Use
------

* Copy CoolProp\\build15\\Release\\CoolPropMathcadWrapper.dll to C:\\Program Files (x86)\\Mathcad\\Mathcad 15\\userefi 
	
* Copy CoolProp\\wrapper\\Mathcad\\CoolProp_EN.xml to C:\\Program Files (x86)\\Mathcad\\Mathcad 15\\doc\\funcdoc 
	
* Open the CoolPropFluidProperties.xmcd file in MathCAD, all CoolProp functions should evaluate properly. If not, press <Ctrl>-F9 to force recalculation of the entire workbook.

* CoolProp functions can be inserted from the Mathcad Insert Functions panel under the function category: CoolProp.  Input parameters and a brief description of each function will be shown.

