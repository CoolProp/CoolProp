Build and Run
=============

To build this very simple example, do something like this::

    gfortran -c fortmyfunc.f90
    gcc -c -std=c++11 myfunc.cpp
    gfortran fortmyfunc.o myfunc.o -o main -lstdc++
    main

The``-lstdc++`` is required to link in the c++ standard libraries

Running ``main`` should yield something like::

               1
       42.000000000000000
       0.0000000000000000
       1.0000000000000000
       2.0000000000000000
       3.0000000000000000
       4.0000000000000000
       5.0000000000000000
       6.0000000000000000
       7.0000000000000000
       8.0000000000000000
       9.0000000000000000
       10.000000000000000
       11.000000000000000
       12.000000000000000
       13.000000000000000
       14.000000000000000
       15.000000000000000
       16.000000000000000
       17.000000000000000
       18.000000000000000
       19.000000000000000
       
Compiling on Windows
====================

At the moment, the most reliable mixed compilation seems to be using the mingw-provided gfortran/gcc combination from mingw-get.  These are the versions used as of June 20, 2014::

    >gfortran --version
    GNU Fortran (GCC) 4.8.1
    Copyright (C) 2013 Free Software Foundation, Inc.

    >gcc --version
    gcc (GCC) 4.8.1
    Copyright (C) 2013 Free Software Foundation, Inc.
