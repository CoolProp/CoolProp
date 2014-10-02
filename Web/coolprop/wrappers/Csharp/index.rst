.. _Csharp:

**********
C# Wrapper
**********

Pre-compiled Binaries
=====================
Pre-compiled release binaries can be downloaded from :sfdownloads:`Csharp`.  Development binaries coming from the buildbot server can be found at :bbbinaries:`Csharp`.

To Use
------

Copy all the platform-independent .cs files to a folder on your computer you want, here we call it ``platform-independent``.  Copy the DLL for your system architecture to the same location.  Copy the Example.cs file to the same location.  You will need to have a copy of some version of C#.

Windows
^^^^^^^

At the command prompt, run::

    call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
    csc Example.cs platform-independent/*.cs -platform:x64
    Example

where you might need to update the path to visual studio depending on your version installed.  Use `-platform:x86` to tell C# that your DLL is 32-bit if you are on 32-bit, or `-platform:x64` if you are on 64-bit.

Alternatively, you can add all the .cs files to a visual studio project.  

Linux/OSX
^^^^^^^^^

Same idea as windows, but command line is just a bit different::

    mcs Example.cs platform-independent/*.cs -platform:x64
    ./Example

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
    mkdir build && cd build
    # Build the makefile using CMake
    cmake .. -DCOOLPROP_CSHARP_MODULE=ON -DBUILD_TESTING=ON
    # Make the C# files (by default files will be generated in folder install_root/Csharp relative to CMakeLists.txt file)
    make install
    # Run the integration tests
    ctest --extra-verbose
