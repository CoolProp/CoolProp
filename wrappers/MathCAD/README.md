CoolProp Wrapper for PTC Mathcad Prime 7.0 or later (64-bit)
==========================================================

| Copyright Scott Polak and Ian Bell, 2013  
| Updated by Jeff Henning, 2016    

There are two ways to get the CoolProp wrapper add-in for Mathcad Prime; download the pre-compiled binary from SourceForge or Compile your own. 


Precompiled Binaries
====================
The binary DLL file (including nightly builds and release versions) can be downloaded and placed in the Mathcad Prime installation directory. 

1. Goto the [SourceForge CoolProp](https://sourceforge.net/projects/coolprop/files/) web page. 

2. Browse to the desired version. 

3. Under the MathCADPrime folder, download all three files 

   * The README.txt files contains these instructions. 

   * Place the ``CoolPropMathcadWrapper.DLL`` file in the Mathcad Prime installation directory under the ``Custom Functions`` directory. 

   * (Optional) Install the Custom Functions Add-in, [CustFunc](https://github.com/henningjp/CustFunc) and copy the associated ``CoolProp_EN.xml`` file to the ``Custom Functions\docs`` folder, creating it if it doesn't exist (this file is optional and for features in development). The [CustFunc](https://github.com/henningjp/CustFunc) add-in provides a pop-up (when pressing ``<F3>``) that gives brief descriptions of each implemented function and its input parameters and facilitates inserting these functions into the worksheet.

4. Restart Mathcad Prime. The add-in, ``CoolPropMathcadWrapper.DLL`` will automatically load the from the ``Custom Functions`` directory. 

5. Open the ``CoolPropFluidProperties.mcdx`` file (*also found in the Mathcad folder of this repository*) in Mathcad Prime for function usage examples. Press `<Ctrl>-<F5>` to force recalculation of the entire workbook. 

6. The CoolProp functions will _**not**_ be added to the internal Functions menu as this functionality has not been implemented in PTC Mathcad Prime as of version 11.0.0.0.  Users should refer to the usage ``.mcdx`` file referenced above and/or the companion PDF of that file.  Pop-up insert functionality can be implemented using the additional, optional [CustFunc](https://github.com/henningjp/CustFunc) add-in. 


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

             -G "Visual Studio 17 2022" -A x64"  

             -DCMAKE_VERBOSE_MAKEFILE=ON

	> NOTE: Replace `x.0.0.0` with your latest version. Example: `10.0.0.0`                                              	

	> NOTE: Replace `Visual Studio 17 2022` with your current Visual Studio compiler version.  
    
    > NOTE: Visual Studio versions earlier than Visual Studio 2017 will not recognize the bitness flag (`-A`) and it must be combined into the geneorator `-G` string with the `Win64` suffix, such as:  
     `-G "Visual Studio 14 2015 Win64"`         		 
			 
* **Make the dynamic library (DLL)**::

    cmake --build . --config Release

    (*Alternatively, open the VS solution created in the build directory with Visual Studio, set the config to Release and x64, and compile the DLL under the Build menu.*)

To Use
------

* Copy CoolProp\\buildprime\\Release\\CoolPropMathcadWrapper.dll file to C:\\Program Files\\PTC\\Mathcad Prime ``x.0.0.0``\\Custom Functions (replace ``x.0.0.0`` with your current version)

* Open the `CoolPropFluidProperties.mcdx` file in MathCAD and press `<Ctrl>-<F5>` to force recalculation of the entire workbook.  This file provides usage examples for the implemented CoolProp Functions.

* (Optional) Install the Custom Functions Add-in, [CustFunc](https://github.com/henningjp/CustFunc) and copy the associated ``CoolProp_EN.xml`` file to the ``Custom Functions\docs`` folder, creating it if it doesn't exist (this file is optional and for features in development). The [CustFunc](https://github.com/henningjp/CustFunc) add-in provides a pop-up (when pressing ``<F3>``) that gives brief descriptions of each implemented function and its input parameters and facilitates inserting these functions into the worksheet.