.. _R:

************
R Wrapper
************

.. contents:: :depth: 2

Pre-compiled Binaries
=====================

Pre-compiled binaries can be downloaded from :sfdownloads:`R`.  Development binaries coming from the buildbot server can be found at :sfnightly:`R`.  Download the shared library for your platform, as well as the ``CoolProp.R`` module

Usage
-----
At the R console, run::

    dyn.load(paste("CoolProp", .Platform$dynlib.ext, sep=""))
    source("CoolProp.R")
    cacheMetaData(1)
    PropsSI("T","P",101325,"Q",0,"Water")
    
There is example code :ref:`at the end of this page <r_example>`

.. warning::

    For OSX users that wish to call REFPROP, you may be required to set the environmental variable ``DYLD_LIBRARY_PATH`` to the folder containing your REFPROP.dylib shared library, which is probably ``/opt/refprop``
    
.. warning::

    If you want to use ``Rscript`` rather than ``R``, you need to pass the argument ``--default-packages=methods`` to get it to load the necessary packages for calling SWIGG-ed code.  Or call ``library(methods)`` at the top of the file (before ``source("CoolProp.R")``).  See also http://stackoverflow.com/a/19468533\n

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the R wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Also, you need to install R

Windows
-------

You must install the MINGW compiler suite.  The `TDM suite <http://tdm-gcc.tdragon.net/>`_ is recommended.

Build
-----

Linux and OSX
^^^^^^^^^^^^^

Once the dependencies are installed, you can run the builder and tests using::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    mkdir -p  CoolProp/build && cd CoolProp/build
    # Build the makefile using CMake
    cmake .. -DCOOLPROP_R_MODULE=ON -DR_BIN="/usr/bin" -DCMAKE_BUILD_TYPE=Release
    # Make the R shared library
    cmake --build .

Windows (32-bit and 64-bit)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

You need to just slightly modify the building procedure::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir build && cd build
    # Build the makefile using CMake (Must use MINGW since R is built with MINGW)
    cmake .. -DCOOLPROP_R_MODULE=ON -G "MinGW Makefiles" -DR_BIN="c:\Program Files\R\R-3.2.1\bin\x64" -DCMAKE_BUILD_TYPE=Release
    # Make the R shared library
    cmake --build .

.. _r_example:

Example Code
============

.. literalinclude:: Example.R

Example Code Output
===================

.. literalinclude:: Example.out
