.. _static_library:

**************
Static Library
**************

.. contents:: :depth: 2

Introduction
============

Static libraries can be used to compile all of CoolProp into one compilation unit, and then link that with user code.  This can be advantageous as CoolProp only needs to be compiled once and fast incremental compilation of the user-defined code is then possible.

When writing your own C++ code, it is advised to compile CoolProp to a static library and then link CoolProp and your own code.

CMake Integration
=================

If you are using CMake, the process is quite trivial to integrate the CoolProp static library into your build system.  For instance, if you have the folder layout::

    main
     |- CMakeLists.txt (For your project, see below)
     |- mycode.cpp
     |- externals
        |- CoolProp
            |- src
            |- include
            |- ...
            |- CMakeLists.txt
            |-

Then CMakeLists.txt might have the contents::

    # See also http://stackoverflow.com/a/18697099
    cmake_minimum_required (VERSION 2.8.11)
    project (main)
    set(COOLPROP_STATIC_LIBRARY true)
    add_subdirectory ("${CMAKE_SOURCE_DIR}/externals/CoolProp" CoolProp)
    add_executable (main "${CMAKE_SOURCE_DIR}/mycode.cpp")
    target_link_libraries (main CoolProp)

Pre-compiled Binaries
=====================
Pre-compiled release binaries can be downloaded from :sfdownloads:`static_library`.  Development binaries coming from the buildbot server can be found at :sfnightly:`static_library`.  These static libraries are ONLY useful if the compiler used to build the static library agrees with the static library that will be used to build your other code.  This cannot be guaranteed, and no effort is made to that effect.  Build your own static libraries following the instructions below.

.. warning::

    Make sure that the compiler used to build the static library is the same used to build your application using the static library.  Otherwise you will get linking errors.   For instance, on windows, a .lib file generated with Microsoft Visual Studio 10 is not compatible with Microsoft Visual Studio 9, or 11, and definitely not with mingw+gcc.

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

    If you use mingw port of gcc on windows, you should add the generator name to the first call of cmake so that it reads something like::
    
        cmake .. -G "MinGW Makefiles" -DCOOLPROP_STATIC_LIBRARY=ON
        
.. note::

    If you use Microsoft Visual Studio, you should tell cmake what exact version of visual studio you would like it to use, by doing something like::
    
        cmake .. -G "Visual Studio 12 2013 Win64" -DCOOLPROP_STATIC_LIBRARY=ON
        
    which is a 64-bit build for Microsoft Visual Studio 2013 (even express version) for instance.  You can get the full list of supported generators on your machine by doing `cmake --help`.
    
.. note::
    
    If you use gcc with libstdc++ (like on ubuntu) and want to build the debug library, you should add the proper cxx flags to link to the correct debug libstdc++ librariries::
    
        cmake .. -DCOOLPROP_DEBUG=ON -DCMAKE_CXX_FLAGS_DEBUG='-g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC'

Usage
-----

For all platforms we start with a simple example file here called main.cpp; the "Hello world" of CoolProp ::

    #include "CoolProp.h"
    #include <iostream>

    int main()
    {
        std::cout << CoolProp::PropsSI("T","P",101325,"Q",0,"Water") << std::endl;
        return 1;
    }

Linux and OSX
^^^^^^^^^^^^^

On linux and OSX, you can use the compiled static library in your application by doing something like this (starting in the directory ``build`` relative to the root of the source as in the above compilation step)::

    g++ -ldl -L. -I../../include main.cpp -lCoolProp

This will result in an executable which can be run by the user.

.. warning::
    
    In gcc and mingw ports of gcc, make sure that the `-lCoolProp` is the last argument in the line, otherwise you will certainly get linking errors.  See also: http://www.mingw.org/wiki/specify_the_libraries_for_the_linker_to_use .
    
Windows
^^^^^^^

On windows the two main compiler families are Visual Studio and MINGW+GCC.

**Mingw+gcc**: If you use mingw, follow the instructions like for linux and OSX, and leave off the ``-ldl`` argument to the compilation.

**Visual Studio**: 

a) Generate the static library following the command line instructions above, ensuring that you have selected the proper visual studio version as well as ``Win64`` in your generator if you would like a 64-bit static library
b) Create a new empty project in visual studio, change to 64-bit (x64) build type if you would like
c) Add the include directory of CoolProp to the list of include directories the ``C/C++->General`` tab in visual studio
d) Add the directory where the .lib file is to the list of library directories in the ``Linker->General`` tab of the properties
e) Add ``CoolProp.lib`` to the list of .lib files in the ``Linker->Input`` tab in visual studio
f) Add the above main.cpp file to your project
