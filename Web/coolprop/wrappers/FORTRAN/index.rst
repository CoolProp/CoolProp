.. _FORTRAN:

***************
FORTRAN Wrapper
***************

For FORTRAN, there are two fundamental choices.  Choose one option:

* [Recommended] For FORTRAN 95 and newer, compile a static library of CoolProp and link it with FORTRAN code following the instructions here: :ref:`F95 and newer <FORTRAN95>`
* For FORTRAN 77 and newer, call a shared library of CoolProp using the instructions here: :ref:`F77 and newer <FORTRAN77>`

Common Requirements
===================
Compilation of the Fortran wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Compilers
=========

On linux, you need gcc and gfortran, which are easy to install using your package manager.

On windows, the most reliable mixed compilation seems to be using the gfortran/gcc combination from `MinGW-w64 <http://sourceforge.net/projects/mingw-w64/files>`_, whose installer allows you to install different versions of GCC, typically up to the last one. Version 5.3.0 is the one used as of February 10, 2016::

    >C:\>gfortran --version
    GNU Fortran (x86_64-posix-seh-rev0, Built by MinGW-W64 project) 5.3.0
    Copyright (C) 2015 Free Software Foundation, Inc.

    GNU Fortran comes with NO WARRANTY, to the extent permitted by law.
    You may redistribute copies of GNU Fortran
    under the terms of the GNU General Public License.
    For more information about these matters, see the file named COPYING


    >C:\>g++ --version
    g++ (x86_64-posix-seh-rev0, Built by MinGW-W64 project) 5.3.0
    Copyright (C) 2015 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

.. warning::
    gfortran in GCC version 5.1 has a bug preventing to open external files (a segmentation error is prompted), which makes this   version almost useless for fortran users.  


On OSX, the default compiler that comes with XCode is clang, gcc and g++ at the command prompt are just aliases to clang.  See for instance::

    $ gcc --version
    Configured with: --prefix=/Applications/Xcode.app/Contents/Developer/usr --with-gxx-include-dir=/usr/include/c++/4.2.1
    Apple LLVM version 5.1 (clang-503.0.40) (based on LLVM 3.4svn)
    Target: x86_64-apple-darwin13.3.0
    Thread model: posix

The most reliable mixed compilation on linux seems to be with the true-gcc/g++ toolchain rather than clang+gfortran.  To install gcc/g++, install gcc/g++ using homebrew with something like::

    brew install gcc-4.9

where you can change ``4.9`` to the most up to date version.  Search `braumeister <http://braumeister.org/formula/gcc>`_ for the most recipe for gcc.  Looks like ``gcc`` also works.

.. _FORTRAN95:

FORTRAN 95/2003 using ISO_C_BINDING
===================================

Build static library of CoolProp
--------------------------------

On all platforms, start in root folder of recursively-cloned CoolProp repo.

A) On linux,  do::

    mkdir build && cd build
    mkdir gccstatic && cd gccstatic
    cmake ../..  -DCOOLPROP_STATIC_LIBRARY=ON -DCOOLPROP_EXTERNC_LIBRARY=ON -DCMAKE_VERBOSE_MAKEFILE=ON
    cmake --build .

B) On Windows, the call to CMake should be done using the MinGW generator, but otherwise the procedure is the same::

    mkdir build && cd build
    mkdir gccstatic && cd gccstatic
    cmake ../.. -G "MinGW Makefiles" -DCOOLPROP_STATIC_LIBRARY=ON -DCOOLPROP_EXTERNC_LIBRARY=ON -DCMAKE_VERBOSE_MAKEFILE=ON
    cmake --build .

C) On OSX, cmake must use the true, real, gcc/g++ compiler (not clang).  Thus you must do something like this to make sure that it finds the right (true) gcc/g++ (see above)::

    mkdir build && cd build
    mkdir gccstatic && cd gccstatic
    cmake ../.. -DCOOLPROP_STATIC_LIBRARY=ON -DCOOLPROP_EXTERNC_LIBRARY=ON -DCMAKE_C_COMPILER="/usr/local/bin/gcc-4.9" -DCMAKE_CXX_COMPILER="/usr/local/bin/g++-4.9" -DCMAKE_VERBOSE_MAKEFILE=ON
    cmake --build .

If you are using a different version of gcc, change the version number for g++ and gcc

This will generate the file libCoolProp.a which is a GCC static library that can be linked with GCC/GFORTRAN code.  Copy this .a file into the directory with the coolprop FORTRAN example ``cool_fortran_bind.f90``:

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
    gcc -o main *.o libCoolProp.a -lstdc++ -ldl -lgfortran -lm
    main

On windows, you can leave off the ``-ldl`` and also the ``-lm`` might not be required.

On OSX, you must do the linking stage with true gcc so that it finds the right standard library.  Or alternatively, provide the full path to the libstdc++ static library and link with gfortran with something like::

    gfortran -o main *.o libCoolProp.a /usr/lib/libstdc++.a -ldl

.. warning::

    You MUST(!!!) put the -lstdc++ standard library *after* libCoolProp.a.  Same thing if you compile the fortran to object file, static library must always be at the end.

.. _FORTRAN77:

FORTRAN77
=========

Pre-Compiled Binaries
---------------------

* Download the appropriate shared library for your architecture from from :sfdownloads:`shared_library`, or the development versions from the buildbot server at :sfnightly:`shared_library`. Or you can built it yourself given the instructions at :ref:`shared_library`.

Run
---

Use the sample file ``example.for`` given by:

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
    write(*,*) "Result was: ", outVal/1e5, " bar"
    write(*,*) "-----------------------------------------------"

    end program

Place the shared library and the sample file in the same directory.

On linux, build and run the example.for file with::

    gfortran -g -o example example.for -L. -lCoolProp
    LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH example

On windows, the current folder is always searched for DLL, so you can just do::

    gfortran -g -o example example.for -L. -lCoolProp
    example
