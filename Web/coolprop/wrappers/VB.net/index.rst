.. _VBdotnet:

**************
VB.net Wrapper
**************

Pre-compiled Binaries
=====================

Information
-----------

As C# and VB.net share the .NET framework, it is possible to call the C# wrapper code from VB.net, thus making it possible to access the entire code of CoolProp from VB.net.  The overall process is not too painful.

Pre-compiled binaries of the C# wrapper can be downloaded from :sfdownloads:`Csharp`.  Development binaries coming from the buildbot server can be found at :sfnightly:`Csharp`.  Download the ``platform-independent.7z`` file and expand it to a folder called ``platform-independent`` using 7-zip.  Download the special C# shared library for your system architecture to the same location from either :sfdownloads:`Csharp` (release) or :sfnightly:`Csharp` (development). 

When you are finished, you should have a folder layout something like ::

    main
     |- CoolProp.dll
     |- Example.vb
     |- platform-independent
        |- AbstractState.cs
        |- Configuration.cs
        |- ...

Open Visual Studio 2012 (or any other version of Visual Studio).  Even the express version works.

Create a VB console application.  Add a project that is a C# Class library to the solution.  Add all the platform-independent .cs files to the class library project.  Add the Example.vb file to the VB console application project.  Add the CoolProp.dll file as an existing file to the VB console project.  Right click on the DLL, select the property "Copy to Output Directory" and make sure it is set to "Copy always".  Right click on the console project, "Add Reference...", expand Solution, then projects. Select the class library project.  This will copy the class library dll to the output directory and allow Intellisense to properly parse the C# code.

Caveats
-------

The architecture of the solution/projects should match that of the CoolProp.dll file.  So if you are using the 64-bit DLL, make sure the architecture of the console application project is set to ``x64`` rather than ``Any CPU``

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the C# wrapper for VB.net requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Build on 64-bit windows
^^^^^^^^^^^^^^^^^^^^^^^

The build procedure:

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir build && cd build
    # Build the makefile using CMake
    cmake .. -DCOOLPROP_CSHARP_MODULE=ON -DCOOLPROP_SWIG_OPTIONS="-namespace CoolProp" -G "Visual Studio 11 2012 Win64"
    # Make the C# shared library
    cmake --build . --config Release