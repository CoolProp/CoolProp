.. _Java:

************
Java Wrapper
************

Pre-compiled Binaries
=====================
Pre-compiled binaries can be downloaded from :sfdownloads:`Java`, which come from the continuous development buildbot server at :ref:`http://www.coolprop.dreamhosters.com:8010/binaries`.

When the binaries are found, java should by executed () with 

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the Java wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`
    
Linux & OSX
-----------

1. Download a zip file of the Java Development Kit (JDK) for Java from `Oracle downloads <http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html>`_. If you are in a 32-bit system, download the 32-bit system, else download the 64-bit version.

2. Expand the zip file you downloaded

3. Add the ``bin`` folder of the JDK that you installed.  For instance, add the following to ~/.profile:
      
      export /path/to/java/SDK/bin:$PATH 
      
  to ~/.profile where the path ``/path/to/java/SDK/bin`` points to the absolute path appropriate to


Windows
-------

Install MATLAB using installer downloaded from www.mathworks.com.  As of version R2014a, both of 32-bit and 64-bit MATLAB is available

Build
-----

Linux and OSX
^^^^^^^^^^^^^

Once the dependencies are installed, you can run the builder and tests using::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    mkdir -p  CoolProp/build/Java && cd CoolProp/build/Java
    # Build the makefile using CMake
    cmake ../.. -DCOOLPROP_MATLAB_MODULE=ON -DBUILD_TESTING=ON
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
    mkdir build/MATLAB
    # Move into that folder
    cd build/MATLAB
    # Set an environmental variable that points to your MATLAB installation for use in CMake (adjust if needed)
    set MATLAB_ROOT=c:\Program Files\MATLAB\R2014a
    # Build the makefile using CMake
    cmake ../.. -DCOOLPROP_MATLAB_MODULE=ON -DBUILD_TESTING=ON
    # Make the MEX files (by default files will be generated in folder install_root/MATLAB relative to CMakeLists.txt file)
    make install
    # Run the integration tests
    ctest --extra-verbose

Usage
=====

Place the mex files somewhere on the MATLAB path.

If you place mex file somewhere outside MATLAB path, you have to use
"addpath" function at begining of your code.

Example: adding the folder that contains CoolProp.mexw32 file to the Octave path::

    addpath('/home/USERNAME/Some_folder/CoolProp')
