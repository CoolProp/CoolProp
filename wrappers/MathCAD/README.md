CoolProp Wrapper for PTC Mathcad Prime 7.0 or later (64-bit)
==========================================================

| Copyright Scott Polak and Ian Bell, 2013
| Updated by Jeff Henning, 2016

 ..  

There are two ways to get the CoolProp wrapper add-in for Mathcad Prime; download the pre-compiled binary from SourceForge or Compile your own. 


Precompiled Binaries
====================
The binary DLL file (including nightly builds and release versions) can be downloaded and placed in the Mathcad Prime installation directory. 

1. Goto the [SourceForge CoolProp](https://sourceforge.net/projects/coolprop/files/) web page. 

2. Browse to the desired version. 

3. Under the MathCADPrime folder, download all three files 

   * The README.txt files contains these instructions. 

   * Place the ``CoolPropMathcadWrapper.DLL`` file in the Mathcad Prime installation directory under the ``Custom Functions`` directory. 

   * Copy the ``CoolProp_EN.xml`` file to the ``Custom Functions\docs`` folder, creating it if it doesn't exist (this file is optional and for features in development). 

4. Restart Mathcad Prime. The add-in, ``CoolPropMathcadWrapper.DLL`` will automatically load the from the ``Custom Functions`` directory. 
5. Open the ``CoolPropFluidProperties.mcdx`` file (*also found in the Prime folder of this repository*) in Mathcad Prime for function usage examples. 
6. The CoolProp functions will _**not**_ be added to the internal Functions menu as this functionality has not been implemented in PTC Mathcad Prime as of version 9.0.0.0.  Users should refer to the usage ``.mcdx`` file referenced above and/or the companion PDF of that file. 


Build Your Own
==============

Prerequisites
-------------

* You will need to have Visual Studio 2015 or later installed.  You will need the professional version of Visual Studio C++, or at least Visual Studio Express 2015 or later, as Mathcad Prime libraries are 64-bit and require the 64-bit compiler.  The community edition of Visual Studio will also work.

* You will need CMake version 2.8.12 or later from https://cmake.org/download/.

* You will need to install Git-SCM for Windows.  You can install this from https://git-for-windows.github.io or https://www.git-scm.com.

* You will need Anaconda/Miniconda Python, which you can get from https://store.continuum.io/cshop/anaconda.  This is needed because the CMake process runs some python scripts that build the fluids and mixtures libraries.
	

To Build
--------

* **Recursively clone the CoolProp library to a local repository.**::

	git clone https://github.com/CoolProp/CoolProp --recursive

* **Change directory (cd) to the CoolProp directory you just created.**::

	cd CoolProp

* **In the top level CoolProp directory and make a build directory** (something like \buildprime)::

    mkdir buildprime
    cd buildprime

* **Build the makefile using CMake (adjust root string for correct version of Prime)**::

    cmake .. -DCOOLPROP_PRIME_MODULE=ON 
             -DCOOLPROP_PRIME_ROOT="C:/Program Files/PTC/Mathcad Prime x.0.0.0" 
             -G "Visual Studio 14 2015 Win64" 
             -DCMAKE_VERBOSE_MAKEFILE=ON

	( Note: Replace 'x.0.0.0' with your latest version.                                              )		 
	( Note: Replace 'Visual Studio 14 2015 Win64'  with your current Visual Studio version.          )		 
	( Note: Mathcad Prime is 64-bit, so the 'Win64' option is necessary in the Visual Studio string. )		 
			 
* **Make the dynamic library (DLL)**::

    cmake --build . --config Release

    (*Alternatively, open the VS solution created in the build directory with Visual Studio, set the config to Release and x64, and compile the DLL under the Build menu.*)

To Use
------

* Copy CoolProp\\buildprime\\Release\\CoolPropMathcadWrapper.dll file to C:\\Program Files\\PTC\\Mathcad Prime ``x.0.0.0``\\Custom Functions (replace ``x.0.0.0`` with your current version)

* Open the CoolPropFluidProperties.xmcd file in MathCAD, all CoolProp functions should evaluate properly.  If not, press <Ctrl>-F9 to force recalculation of the entire workbook.

