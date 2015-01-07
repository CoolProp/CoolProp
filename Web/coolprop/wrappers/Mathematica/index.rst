.. _Mathematica:

*******************
Mathematica Wrapper
*******************

Pre-compiled Binaries
=====================
Pre-compiled binaries can be downloaded from :sfdownloads:`Mathematica`, which come from :sfnightly:`the nightly snapshots <Mathematica>`.

Place the shared library for your platform (windows: CoolProp.dll, OSX: CoolProp.dylib, linux: CoolProp.so) in the directory that you can obtain from Mathematica::

    FileNameJoin[{$BaseDirectory, "SystemFiles", "LibraryResources", $SystemID}]

If this directory doesn't yet exist, create it. At the command prompt, you should be able to call ``FindLibrary["CoolProp"]``

There is a small example file ``example.nb`` that demonstrates how to call the functions in the DLL which can be downloaded from :sfdownloads:`Mathematica`

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the Mathematica wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Linux and OSX
^^^^^^^^^^^^^

Once the dependencies are installed, you can run the builder and tests using::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    mkdir -p CoolProp/build && cd CoolProp/build
    # Build the makefile using CMake
    cmake .. -DCOOLPROP_MATHEMATICA_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON
    # Make the shared library
    cmake --build .

Windows
^^^^^^^

You need to just slightly modify the building procedure

1. Checkout and preparation::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir build && cd build

2. Pick a toolchain (A or B)

    A: Building using Visual Studio::

        # Build the makefile using CMake
        cmake .. -DCOOLPROP_MATHEMATICA_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON -G "Visual Studio 10 2010 Win64"
        
    If you have a different version of Visual Studio installed, replace the generator name in the ``-G`` argument
        
    B: Building using MinGW::

        # Build the makefile using CMake
        cmake .. -DCOOLPROP_MATHEMATICA_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON -G "MinGW Makefiles"
    
3. Actually do the build::
    
    # Make the shared library
    cmake --build . --config Release
