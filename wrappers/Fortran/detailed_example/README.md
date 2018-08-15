FORTRAN example with module
===========================

To build static library (windows)
---------------------------------

Start in root folder of repo

    mkdir build
    cd build
    mkdir gccstatic && cd gccstatic
    cmake ../.. -G "MinGW Makefiles" -DCOOLPROP_STATIC_LIBRARY=ON
    cmake --build .

This will generate the file libCoolProp.a which is a GCC static library that can be linked with GCC/GFORTRAN code

Make sure that the macro EXTERNC is defined - this will give extern "C" decorations for all the CoolPropLib.h functions.

Copy this .a file into the directory with the coolprop FORTRAN example

To build FORTRAN example
------------------------

    gfortran -c -Wall cpinterface.f90 cool_fortran_bind.f90
    gfortran -o main *.o libCoolProp.a -lstdc++
    main
    
.. warning::

    You MUST(!!!) put the -lstdc++ standard library *after* the cool_fortran_bind.f90  Same thing if you compile the fortran to object file, static library must always be at the end.

Compiling on Windows
--------------------

At the moment, the most reliable mixed compilation seems to be using the mingw-provided gfortran/gcc combination from mingw-get.  These are the versions used as of June 20, 2014::

    >gfortran --version
    GNU Fortran (GCC) 4.8.1
    Copyright (C) 2013 Free Software Foundation, Inc.

    >gcc --version
    gcc (GCC) 4.8.1
    Copyright (C) 2013 Free Software Foundation, Inc.