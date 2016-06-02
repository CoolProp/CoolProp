.. _Csharp:

**********
C# Wrapper
**********

.. contents:: :depth: 2

Pre-compiled Binaries
=====================

To Use
------

Pre-compiled binaries can be downloaded from :sfdownloads:`Csharp`.  Development binaries coming from the buildbot server can be found at :sfnightly:`Csharp`.

Download the ``platform-independent.7z`` file and expand it to a folder called ``platform-independent`` using 7-zip.  Download the special C# shared library for your system architecture to the same location from either :sfdownloads:`Csharp` (release) or :sfnightly:`Csharp` (development).  Copy the Example.cs file to the same location.  You will need to have a copy of some version of C#.

When you are finished, you should have a folder layout something like ::

    main
     |- CoolProp.dll
     |- Example.cs
     |- platform-independent
        |- AbstractState.cs
        |- Configuration.cs
        |- ...
        
There is example code :ref:`at the end of this page <csharp_example>`

Windows
^^^^^^^

At the command prompt, run::

    call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
    csc Example.cs platform-independent/*.cs -platform:x64
    Example

where you might need to update the path to visual studio depending on your version installed.  Use `-platform:x86` to tell C# that your DLL is 32-bit if you are on 32-bit, or `-platform:x64` if you are on 64-bit.

Alternatively, you can add all the .cs files to a visual studio project.  If you do that, add the DLL to the project as well, right-click on the DLL, and select the option to copy it to the output directory.

Linux/OSX
^^^^^^^^^

Same idea as windows, but command line is just a bit different::

    mcs Example.cs platform-independent/*.cs -platform:x64
    ./Example
    
Use `-platform:x86` to tell C# that your shared library is 32-bit if you are on 32-bit, or `-platform:x64` if you are on a 64-bit platform.

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the C# wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Additionally, you will need:
* SWIG (see :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`)
* C#

OSX
---

For OSX, to install the necessary tools using homebrew, you can do::
    
    homebrew install mono

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
    cmake --build .
    # Run the integration tests (optional)
    ctest --extra-verbose

If you want to change the package that CoolProp resides in, you can do so by changing the cmake call to read::

    cmake .. -DCOOLPROP_CSHARP_MODULE=ON -DBUILD_TESTING=ON -DCOOLPROP_SWIG_OPTIONS="-namespace package.name"

where ``package.name`` is replaced with the desired name    
    
.. _csharp_example:

Example Code
============

.. literalinclude:: Example.cs
   :language: csharp

Example Code Output
===================

.. literalinclude:: Example.out
