.. _mathcadprime:

********************************************************
CoolProp Wrapper for Mathcad Prime 7.0 or later (64-bit)
********************************************************

| Copyright Scott Polak and Ian Bell, 2013
| Updated by Jeff Henning, 2016

Pre-compiled binaries
=====================

Pre-compiled release binaries can be downloaded from :sfdownloads:`MathcadPrime`.  The latest development binaries coming from the nightly builds can be found at :sfnightly:`MathcadPrime`.

To Use
------

* Copy CoolPropMathcadWrapper.dll file to C:\\Program Files\\PTC\\Mathcad Prime x.0.0.0\\Custom Functions; where x.0.0.0 is the Mathcad Prime version being used (7 or later).
 
* Open the CoolPropFluidProperties.mcdx file in Mathcad Prime. Press `<Ctrl>-F5`` to force recalculation of the entire workbook. 
 
* (Optional) Use the Custom Functions add-in at `CustFunc <https://github.com/henningjp/CustFunc>`_ to provide a Custom Function pop-up (by pressing `<F3>`) that gives descriptions of each available function and its inputs and allows easy insertion of the CoolProp functions into the worksheet.

User-compiled binaries
======================

Common Requirements
-------------------

* Compilation of the Mathcad Prime wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

* You will need to have Visual Studio 2015 or later installed.  You will need the professional or community version of Visual Studio C++, or at least Visual Studio Express 2015 or later, as Mathcad Prime libraries are 64-bit and require the 64-bit compiler.

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
             -DCOOLPROP_PRIME_ROOT="C:/Program Files/PTC/Mathcad Prime 10.0.0.0"  
             -G "Visual Studio 17 2022" -A x64  
             -DCMAKE_VERBOSE_MAKEFILE=ON  
  
.. note::   
   Mathcad Varsion: Use your version of Mathcad Prime in place of **"10.0.0.0"** (with `7.0.0.0` or later)

.. note::
   Compiler Bitness: Mathcad Prime is 64-bit, so the `-A x64` option is necessary in the Visual Studio string for VS2017 or later.  

.. note::
   Older Compilers: Visual Studio versions earlier than 2017 will specify compiler bitness using the format **-G "Visual Studio 14 2015 Win64"** for the compiler string and will not use the `-A` option.
             
* Make the dynamic library (DLL)::

    cmake --build . --config Release

To Use
------

* Copy `CoolProp\\buildprime\\Release\\CoolPropMathcadWrapper.dll` file to `C:\\Program Files\\PTC\\Mathcad Prime 7.0.0.0\\Custom Functions` (or your current version) 
 
* Open the `CoolPropFluidProperties.mcdx` file in Mathcad Prime and Press `<Ctrl>-F5` to force recalculation of the entire workbook, all CoolProp functions should evaluate properly. 
 
* (Optional) Use the Custom Functions add-in at `CustFunc <https://github.com/henningjp/CustFunc>`_  to provide a Custom Function pop-up (by pressing `<F3>`) that gives descriptions of each available function and its inputs and allows easy insertion of the CoolProp functions into the worksheet.
