.. _Csharp:

**********
C# Wrapper
**********

Pre-compiled Binaries
=====================
Pre-compiled release binaries can be downloaded from :sfdownloads:`Csharp`.  Development binaries coming from the buildbot server can be found at :bbbinaries:`Csharp`.  Download the files appropriate to your system.

To Use
------

Copy all the .cs files to a location you want.  You will need to have a copy of some version of C#.

At the command prompt, run::

    call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
    csc *.cs /platform:x86
    call Example

where you might need to update the path to visual studio depending on your version installed.  

32-bit: Use `/platform:x86` to avoid PINVOKE errors!
64-bit: Use `/platform:x64` to avoid PINVOKE errors!

Alternatively, you can add all the .cs files to a visual studio project.

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the C# wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Octave Requirements
-------------------
* SWIG
* C#

Linux
-----

For ubuntu and friends, you will need to install Mono C# as well as the compiler (and other dependencies) using::

    sudo apt-get install swig mono-mcs mono-runtime

Windows
-------
For Windows, download the Visual Studio 2010 version of C# (other versions should probably be fine too)

Compile
-------

Once mono c# is installed, you can run the builder and tests using::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir -p build/Csharp && cd build
    # Build the makefile using CMake
    cmake .. -DCOOLPROP_CSHARP_MODULE=ON -DBUILD_TESTING=ON
    # Make the C# files (by default files will be generated in folder install_root/Csharp relative to CMakeLists.txt file)
    make install
    # Run the integration tests
    ctest --extra-verbose
