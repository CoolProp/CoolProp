.. _Octave:

**************
Octave Wrapper
**************

Pre-compiled Binaries
=====================
Pre-compiled binaries can be downloaded from XXXXXXXXXXXXXX

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the Octave wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Octave Requirements
-------------------
* SWIG
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
Due to difficulties with interfacing CMake/SWIG/Visual Studio, the Visual Studio compiled versions of octave are not supported as of version 5.  The only windows port of Octave that is supported is the MinGW compiled version.

Build
-----

Once the dependencies are installed, you can run the builder and tests using::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir -p build/Octave
    # Move into that folder
    cd build/Octave
    # Build the makefile using CMake
    cmake ../.. -DCOOLPROP_OCTAVE_MODULE=ON -DBUILD_TESTING=ON
    # Make the OCT files (by default files will be generated in folder install_root/Octave relative to CMakeLists.txt file)
    make install
    # Run the integration tests
    ctest --extra-verbose

On windows, you need to just slightly modify the building procedure::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir build/Octave
    # Move into that folder
    cd build/Octave
    # Build the makefile using CMake
    cmake ../.. -G "MinGW Makefiles" -DCOOLPROP_OCTAVE_MODULE=ON -DBUILD_TESTING=ON
    # Make the OCT files (by default files will be generated in folder install_root/Octave relative to CMakeLists.txt file)
    make install
    # Run the integration tests
    ctest --extra-verbose

This building process on windows requires that make and the ``bin`` folder of Octave are on the ``PATH``.  The MinGW version of Octave conveniently also includes the MinGW C++ compiler used to build it, look inside the distribution.

Usage
=====
On Linux systems you can put generated .oct file in
``/usr/share/octave/?octave.version.number?/m`` folder. You will need superuser
privileges to do this.

If you place .oct file somewhere outside octave path, you have to use
"addpath" function at begining of your code.

Example: adding the folder that contains CoolProp.oct file to the Octave path::

    addpath('/home/?user_name?/Some_folder/CoolProp')
