.. _Octave:

**************
Octave Wrapper
**************

.. contents:: :depth: 2

Pre-compiled Binaries
=====================
Pre-compiled release binaries can be downloaded from :sfdownloads:`Octave`.  Development binaries coming from the buildbot server can be found at :sfnightly:`Octave`.  Download the oct file appropriate to your system.

On Linux systems you can put the generated .oct file in
``/usr/share/octave/?octave.version.number?/m`` folder. You will need superuser
privileges to do this.

If you place .oct file somewhere outside octave path, you have to use
"addpath" function at beginning of your code.

Example: adding the folder that contains CoolProp.oct file to the Octave path::

    addpath('/home/?user_name?/Some_folder/CoolProp')
    
There is example code :ref:`at the end of this page <octave_example>`

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the Octave wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

You will also need:

* SWIG (see :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`)
* Octave (and development headers)

Linux
-----

For ubuntu and friends, you can install build dependencies using::

    sudo apt-get install swig octave liboctave-dev

OSX
---
For OSX, your best best is a binary installer (see http://wiki.octave.org/Octave_for_MacOS_X), alternatively, you can install from Homebrew, though as of July 6, 2014, this functionality was broken in OSX 10.9.  If you use the installer, you might want to add the octave binary folder onto the path.  To do so, add to the file .profile (or create it) in your home directory::

    export PATH="/usr/local/octave/3.8.0/bin:$PATH"

Windows
-------
For windows, the situation is ok, but not great.  Only the MinGW builds are supported, and not comfortably

1. Download a MinGW build from `Octave for windows <http://wiki.octave.org/Octave_for_Microsoft_Windows>`_.

2. Extract the zip file to somewhere on your computer without any spaces in the path (c:\\octave-x.x.x is a good choice)

3. Rename the sh.exe in the bin folder of your installation to _sh.exe

.. warning::
    MinGW has problems with the latest version of CoolProp.  This seems to be a GCC-related 
    issue and using a more up-to-date version of GCC helps.  Unfortunately, MinGW is stuck 
    at GCC 4.8.  You could try the `TDM-GCC distribution <http://tdm-gcc.tdragon.net>`_ 
    that comes with the latest GCC. This version seems to work fine.

Build
-----

Once the dependencies are installed, you can run the builder and tests using::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir -p build &&  cd build
    # Build the makefile using CMake
    cmake .. -DCOOLPROP_OCTAVE_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_CXX_STANDARD=11
    # Make the OCT files
    cmake --build .
    # Run the integration tests (optional)
    ctest --extra-verbose

On windows, you need to just slightly modify the building procedure::

    # The folder containing the folders bin, mingw, swigwin, include, etc...
    # is required
    set OCTAVE_ROOT=c:\octave-x.y.z
    set PATH=c:\MinGW;c:\octave-x.y.z\bin;c:\swigwin-x.y.z;%PATH%
    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir build
    # Move into that folder
    cd build
    # Build the makefile using CMake
    cmake .. -G "MinGW Makefiles" -DCOOLPROP_OCTAVE_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_CXX_STANDARD=11
    # Make the OCT files
    cmake --build .
    # Run the integration tests (optional)
    ctest --extra-verbose

.. _octave_example:

Example Code
============

.. literalinclude:: Example.m
   :language: octave

Example Code Output
===================

.. literalinclude:: Example.out
