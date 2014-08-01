.. _EES:

***********
EES Wrapper
***********

EES is an acausal solver that can be used to solve a wide range of technical problems.  It can be obtained from http://www.fchart.com/ees/.  Though EES has its own set of thermodynamic properties, CoolProp also implements a number of things that are not in EES (incompressibles, interpolation methods, etc.).

Pre-compiled Binaries
=====================
Pre-compiled binaries can be downloaded from XXXXXXXXXXXXXX

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the EES wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Build
-----

Once the dependencies are installed, you can run the builder and tests using::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir build/EES
    # Move into that folder
    cd build/EES
    # Build the makefile using CMake
    cmake ../.. -DCOOLPROP_EES_MODULE=ON -DBUILD_TESTING=ON
    # Make the DLF file and the installer (by default files will be generated in folder install_root/EES relative to CMakeLists.txt file)
    make install
    # Run the integration tests
    ctest --extra-verbose
