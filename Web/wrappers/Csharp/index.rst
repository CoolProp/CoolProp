.. _Csharp:

**********
C# Wrapper
**********

Pre-compiled Binaries
=====================
Pre-compiled binaries can be downloaded from XXXXXXXXXXXXXX

User-Compiled Binaries
======================

Requirements
------------
* CMake
* git
* SWIG
* C#
* python (to build the headers for CoolProp)

Linux
-----

For ubuntu and friends, you will need to install Mono C# as well as the compiler (and other dependencies) using::

    sudo apt-get install cmake swig mono-mcs mono-runtime

Windows
-------
For windows, download the Visual Studio 2010 version of C# (other versions should be fine too)

Compile
-------

Once mono c# is installed, you can run the builder and tests using::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir -p build/Csharp
    # Move into that folder
    cd build/Csharp
    # Build the makefile using CMake
    cmake ../.. -DCOOLPROP_JAVA_MODULE=ON -DBUILD_TESTING=ON
    # Make the C# files (by default files will be generated in folder install_root/Csharp relative to CMakeLists.txt file)
    make install
    # Run the integration tests
    ctest --extra-verbose
