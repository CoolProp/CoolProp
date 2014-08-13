.. _MATLAB:

**************
MATLAB Wrapper
**************

Pre-compiled Binaries
=====================
Pre-compiled binaries can be downloaded from :sfdownloads:`MATLAB`.  Download the files appropriate to your installation of MATLAB.

Usage
-----

Place the mex files somewhere on the MATLAB path.

If you place mex file somewhere outside MATLAB path, you have to use
"addpath" function at begining of your code.

Example: adding the folder that contains CoolProp.mexw32 file to the Octave path::

    addpath('/home/USERNAME/Some_folder/CoolProp')

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the MATLAB wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`
    
Linux
-----

Install MATLAB using installer downloaded from www.mathworks.com.  As of version R2014a, only 64-bit MATLAB is available

OSX
---

Install MATLAB using installer downloaded from www.mathworks.com.  As of version R2014a, only 64-bit MATLAB is available

Windows
-------

Install MATLAB using installer downloaded from www.mathworks.com.  As of version R2014a, both of 32-bit and 64-bit MATLAB is available


SWIG+MATLAB
-----------
On Ubuntu, to cross-compile, add the package ``mingw-w64`` in the Synaptic Package manager.  This will install a cross-platform installation toolchain

Build
-----

Linux and OSX
^^^^^^^^^^^^^

Once the dependencies are installed, you can run the builder and tests using::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir build && cd build
    # Set an environmental variable that points to your MATLAB installation for use in CMake (adjust if needed)
    export MATLAB_ROOT=/usr/local/MATLAB/R2014a # or /Applications/MATLAB_R2014a.app
    # Build the makefile using CMake
    cmake .. -DCOOLPROP_MATLAB_MODULE=ON -DBUILD_TESTING=ON
    # Make the MEX files (by default files will be generated in folder install_root/MATLAB relative to CMakeLists.txt file)
    make install
    # Run the integration tests
    ctest --extra-verbose

Windows (32-bit and 64-bit)
^^^^^^^^^^^^^^^^^^^^^^^^^^^ 

You need to just slightly modify the building procedure::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir build && cd build
    # Set an environmental variable that points to your MATLAB installation for use in CMake (adjust if needed)
    set MATLAB_ROOT=c:\Program Files\MATLAB\R2014a
    # Build the makefile using CMake
    cmake .. -DCOOLPROP_MATLAB_MODULE=ON -DBUILD_TESTING=ON
    # Make the MEX files (by default files will be generated in folder install_root/MATLAB relative to CMakeLists.txt file)
    make install
    # Run the integration tests
    ctest --extra-verbose

