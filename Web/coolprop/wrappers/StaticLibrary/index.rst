.. _static_library:

**************
Static Library
**************

Static libraries can be used to compile all of CoolProp into one compilation unit, and then link that with user code.  This can be advantageous as CoolProp only needs to be compiled once and fast compilation of the user-defined code is then possible.

When writing your own C++ code, it is advised to compile CoolProp to a static library and then link CoolProp and your own code

Pre-compiled Binaries
=====================
Pre-compiled release binaries can be downloaded from :sfdownloads:`static_library`.  Development binaries coming from the buildbot server can be found at :sfnightly:`static_library`.  These static libraries are only useful if the compiler used to make the static library agrees with the static library that will be used to build your other code.  So best to follow the below instructions to build your own static library.

.. warning::

    Make sure that the compiler used to build the static library is the same used to build your application using the static library.  Otherwise you will get linking errors.

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the static library requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Compile
-------

You can build the static library using::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir -p build && cd build
    # Build the makefile using CMake
    cmake .. -DCOOLPROP_STATIC_LIBRARY=ON
    # Make the static library
    cmake --build .
    
.. note::

    If you use mingw port of gcc, you should add the generator name to the first call of cmake so that it reads something like::
    
        cmake ../.. -G "MinGW Makefiles" -DCOOLPROP_STATIC_LIBRARY=ON

Usage
-----

On linux and OSX, you can use the compiled static library in your application by doing something like this (starting in the directory ``build`` relative to the root of the source as in the above compilation step)::

    g++ -ldl -L. -I../../include main.cpp -lCoolProp
    
.. warning::
    
    In gcc and mingw ports of gcc, make sure that the `-lCoolProp` is the last argument in the line, otherwise you will certainly get linking errors.  See also: http://www.mingw.org/wiki/specify_the_libraries_for_the_linker_to_use

where main.cpp could have the contents for instance of::

    #include "CoolProp.h"
    #include <iostream>

    int main()
    {
        std::cout << CoolProp::PropsSI("T","P",101325,"Q",0,"Water") << std::endl;
        return 1;
    }
    