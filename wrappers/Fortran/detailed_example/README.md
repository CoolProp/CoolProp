FORTRAN example without module
==============================

To build static library
-----------------------

Start in root folder of repo

    mkdir build
    cd build
    mkdir gccstatic
    cd gccstatic
    cmake ../.. -G "MinGW Makefiles" -DCOOLPROP_STATIC_LIBRARY=ON
    cmake --build .

This will generate the file libCoolProp.a which is a GCC static library that can be linked with GCC/GFORTRAN code

Make sure that the macro COOLPROP_LIB is defined and that the macro CONVENTION=__cdecl is set.

Copy this .a file into the directory with the coolprop FORTRAN example

To build FORTRAN example
------------------------

    gfortran -c cool_fortran_bind.f90
    gfortran libCoolProp.a cool_fortran_bind.f90 -o main -lstdc++
    main