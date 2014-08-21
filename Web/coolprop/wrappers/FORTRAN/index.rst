.. _FORTRAN:

***************
FORTRAN Wrapper
***************

FORTRAN 95/2003 using ISO_C_BINDING
==================================

Build static library (windows + MinGW)
--------------------------------------

Start in root folder of recursively-cloned CoolProp repo::

    mkdir build && cd build
    mkdir gccstatic && cd gccstatic
    cmake ../.. -G "MinGW Makefiles" -DCOOLPROP_STATIC_LIBRARY=ON
    cmake --build .

This will generate the file libCoolProp.a which is a GCC static library that can be linked with GCC/GFORTRAN code

Make sure that the macro EXTERNC is defined - this will give extern "C" decorations for all the CoolPropLib.h functions.

Copy this .a file into the directory with the coolprop FORTRAN example

The simple example file ``cool_fortran_bind.f90``:

.. code-block:: fortran

    !Example calculates density of saturated liquid propane at 300 K:
    program simple

        USE cpinterface
        
        implicit none

        !Initialize the variables used in the example
        double precision T, Q, dens1
        character(LEN=32) fluid, out1, n1, n2

        T = 300                  ! Temperature [K]
        Q = 0                    ! Quality [-]

        out1 = "D"//CHAR(0)      ! String with of the output Property
        n1  = "T"//CHAR(0)       ! String with of the input Property #1
        n2  = "Q"//CHAR(0)       ! String with of the input Property #2
        fluid    = "Propane"//CHAR(0)   ! String with the fluid name
          
        dens1 = PropsSI(out1, n1, T, n2, Q, fluid)

        Print *, dens1

    end program simple
    
with the interface file ``cpinterface.f90``:

.. code-block:: fortran

    MODULE CPINTERFACE
        INTERFACE
            FUNCTION PropsSI (output, name1, prop1, name2, prop2, fluidname) BIND(C, NAME='PropsSI')
                use iso_c_binding
                real(C_DOUBLE) :: PropsSI
                character(KIND=c_char) :: output(*)
                character(c_char) :: name1(*)
                real(C_DOUBLE), VALUE :: prop1
                character(c_char) :: name2(*)
                real(C_DOUBLE), VALUE :: prop2
                character(kind=c_char) :: fluidname(*)
                    
            END FUNCTION PropsSI
        END INTERFACE
    END MODULE CPINTERFACE

In order to link all the files together, do::

    gfortran -c -Wall cpinterface.f90 cool_fortran_bind.f90
    gfortran -o main *.o libCoolProp.a -lstdc++
    main
    
.. warning::

    You MUST(!!!) put the -lstdc++ standard libary *after* the cool_fortran_bind.f90  Same thing if you compile the fortran to object file, static library must always be at the end.

Compiling on Windows
--------------------

At the moment, the most reliable mixed compilation seems to be using the mingw-provided gfortran/gcc combination from mingw-get.  Theese are the versions used as of June 20, 2014::

    >gfortran --version
    GNU Fortran (GCC) 4.8.1
    Copyright (C) 2013 Free Software Foundation, Inc.

    >gcc --version
    gcc (GCC) 4.8.1
    Copyright (C) 2013 Free Software Foundation, Inc.
    

FORTRAN77
=========

Pre-Compiled Binaries
---------------------

* Download the appropriate shared library for your architecture from from :sfdownloads:`shared_library`, or the development versions from the buildbot server at :bbbinaries:`shared_library`. 

Run
---

Use the sample file ``example.f77`` given by:

.. code-block:: fortran

    double precision T, Q, D, h, s
    character(LEN=32) Ref,Output, Name1, Name2
    double precision outVal, Prop1, Prop2

    T = 285
    Q = 0
    D = 1250;

    Output = "P"//CHAR(0)
    Name1  = "T"//CHAR(0)
    Prop1  = T
    Name2  = "Q"//CHAR(0)
    Prop2  = Q
    Ref    = "R134a"//CHAR(0)
    outval = 9999999

    write(*,*) "Saturation pressure for R134a: "
    call propssi(Output, Name1, Prop1, Name2, Prop2, Ref, outVal)
    write(*,*) "Result was: ", outVal/1e2, " bar"
    write(*,*) "-----------------------------------------------"
    
    end program

Place the shared library and the sample file in the same directory.  Build and run the example.f77 file with::

    gfortran -g -o example example.f77 -L. -lCoolProp
    example

