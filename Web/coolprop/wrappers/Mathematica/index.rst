.. _Mathematica:

*******************
Mathematica Wrapper
*******************

Pre-compiled Binaries
=====================
Pre-compiled binaries are no longer uploaded to the  :sfdownloads:`latest Release` on SourceForge or :sfnightly:`the nightly snapshots <Mathematica>`.  This is because:
* Usage is one of the lowest of all the wrappers supported by CoolProp with few active developers
* The wrapper should really be built by the user for their specific version of Wolfram Mathematica on their system.  Major and minor updates to Mathematica are just too frequent to keep up with. 

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the Mathematica wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

A version of Wolfram Mathematica (v13.0 or higher) or Wolfram Engine

Linux and OSX
^^^^^^^^^^^^^

Once the dependencies are installed, you can run the builder and tests using::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp
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
    git clone https://github.com/CoolProp/CoolProp
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir build && cd build

2. Pick a toolchain (A or B)

    A: Building using Visual Studio::

        # Build the makefile using CMake
        cmake .. -DCOOLPROP_MATHEMATICA_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON -G "Visual Studio 17 2022" -A x64"
        
.. note::   
    If you have a different version of Visual Studio installed, replace the generator name in the ``-G`` argument.  Note that the DLL must be 64-bit by using the ``-A x64`` option at the end.  Prior to VS 2017, this "bitness" was embedded in the ``-G`` argument using the format **-G "Visual Studio 14 2015 Win64"** and will not use the ``-A`` option.
        
    B: Building using MinGW::

        # Build the makefile using CMake
        cmake .. -DCOOLPROP_MATHEMATICA_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON -G "MinGW Makefiles"
    
3. Actually do the build::
    
    # Make the shared library
    cmake --build . --config Release

Install and Verify
==================
Place the shared library for your platform (Windows: ``CoolProp.dll``, OSX: ``CoolProp.dylib``, Linux: ``CoolProp.so``) in the directory that you can obtain from Mathematica::

    FileNameJoin[{$BaseDirectory, "SystemFiles", "LibraryResources", $SystemID}]

If this directory doesn't yet exist, create it. At the command prompt, you should be able to call ``FindLibrary["CoolProp"]``

There is a small example file ``example.nb`` in the CoolProp repository under ``\CoolProp\wrappers\Mathematica`` that demonstrates how to call the functions in the installed DLL.

Future Work
===========
As of Wolfram Mathematica 14.0, the Wolfram Language has the ability to make API calls *directly to the CoolProp shared library*.  This can be done from a Wolfram Language ``paclet`` that obviates the need for compiling a C++ wrapper and can be handled solely through Mathematica's **Paclet Manager**.  This means that the wrapper can be built with automation scripts by any user for any version of Mathematica greater than v14.0.

We welcome any Wolfram Language users to take this on and contribute to the CoolProp project.
