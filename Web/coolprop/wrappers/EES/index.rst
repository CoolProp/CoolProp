.. _EES:

***********
EES Wrapper
***********

EES is an acausal solver that can be used to solve a wide range of technical problems.  It can be obtained from http://www.fchart.com/ees/.  Though EES has its own set of thermodynamic properties, CoolProp also implements a number of things that are not in EES (incompressibles, interpolation methods, etc.).

Pre-compiled Binaries
=====================
Pre-compiled binaries can be downloaded from :sfdownloads:`EES` - follow the instructions there

User-Compiled Binaries
======================

Requirements
-------------------
Compilation of the EES wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Additionally, you must download `InnoSetup <http://www.jrsoftware.org/isinfo.php>`_ and add it to the system path.

Build
-----

Once the dependencies are installed, you can run the installer with::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir build/EES
    # Move into that folder
    cd build/EES
    # Build the makefile using CMake
    cmake ../.. -DCOOLPROP_EES_MODULE=ON
    # Make the DLF file and the installer (by default installer will be generated in folder install_root/EES relative to CMakeLists.txt file)
    cmake --build . --target install
    
Usage
-----
Open EES, you should see the external function COOLPROP_EES.  The functions ``PropsSI`` takes the same inputs as described in the high-level API for C-only inputs.  You can use something like
``fluid$='REFPROP-MIX:'||'R32'||'['||yy$||']'||and$||'R125'||'['||xx$||']'`` where yy and xx are mole fractions of R32 and R125 respectively to encode the string in EES.

The function ``PropsSIZ`` takes the normal inputs, but then also takes the mole fractions as an array rather than encoding them in the string.  The example file for EES demonstrates all of these types of inputs
