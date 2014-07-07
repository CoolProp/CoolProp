.. _wrappers:

******************
Available Wrappers
******************

CoolProp at its core is a C++ library, but it can be of interest to use this code base from other programming environments.  For that reason, wrappers have been constructed for most of the programming languages of technical interest to allow users to seamlessly interface CoolProp and existing codebases.  Downloads and instructions for each wrapper are included in the page for the wrapper.

=======================                    ===========================  =======================================
Language                                   Operating Systems            Notes
=======================                    ===========================  =======================================
:ref:`Octave <Octave>`                     linux, OSX, win              Wrapper is SWIG based
:ref:`C# <Csharp>`                         linux, OSX, win              Wrapper is SWIG based
Java                                       linux, OSX, win              Wrapper is SWIG based
Scilab                                     linux, OSX, win              Wrapper is SWIG based (experimental)
MATLAB                                     linux, OSX, win              Wrapper is SWIG based (VERY experimental)
Python                                     linux, OSX, win
Modelica                                   linux, OSX, win
Javascript                                 linux, OSX, win              Also works in all internet browsers
Maple                  
Mathematica                      
FORTRAN                                    linux, OSX, win
EES                                        windows only
Microsoft Excel                            windows only
=======================                    ===========================  =======================================

.. _wrapper_common_prereqs:

Common Wrapper Prerequisites
============================

On all platforms for which CoolProp is supported, the compilation of one of the wrappers requires a few common prerequisites, described here. They are:

* git (to interface with the CoolProp repository at https://github.com/CoolProp/CoolProp)
* CMake (platform-independent software to generate makefiles)
* C++ compiler (see below)

Windows
-------
On Windows, download the newest binary installer for CMake from `CMake downloads <http://www.cmake.org/cmake/resources/software.html>`_.  Run the installer.  Check that at the command prompt you can do::

    C:\Users\XXXX>cmake -version
    cmake version 2.8.12.2
    
For git, your best best is the installer from http://msysgit.github.io/.  Check that at the command prompt you can do something like::

    C:\Users\XXXX>git --version
    git version 1.9.4.msysgit.0

For the C++ compiler, the options are a bit more complicated.  There are multiple (binary incompatible) versions of Visual Studio, as well as G++ ports for windows (MinGW).  Unless you are compiling the python wrappers, you can compile with MinGW, so you should obtain the `MinGW installer <http://sourceforge.net/projects/mingw/files/Installer/mingw-get-setup.exe/download>`_ and run it.  You should install all the packages available, and you must install to a path without spaces. ``C:\MinGW`` is recommended as an installation path.  

If you are compiling for Python 2.7, you can install Visual Studio 2008 Express from `VS2008Express installer <http://go.microsoft.com/?linkid=7729279>`_.

If you are compiling for Python 3.x, you can install Visual Studio 2010 Express from `VS2010Express installer <http://www.visualstudio.com/en-us/downloads#d-2010-express>`_.

All three compilers should co-exist happily on the path, so you should be fine installing all three, but they are rather sizeable installs.

Linux
-----    
On debian based linux distributions (ubuntu, etc.), you can simply do::

    sudo apt-get install cmake git g++
    
although ``git`` is probably already packaged with your operating system; ``g++`` probably isn't

OSX
---
OSX should come with a c++ compiler (clang), for git and cmake your best bet is `Homebrew <http://brew.sh/>`_.  With Homebrew installed, you can just do::

    sudo brew install cmake git

.. toctree::
    :hidden:

    Octave/index.rst
    Csharp/index.rst