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

.. wrapper_common_prereqs

Common Wrapper Prerequisites
============================

On all platforms for which CoolProp is supported, the compilation of one of the wrappers requires a few common prerequisites, described here. They are:

* git (to interface with the CoolProp repository at https://github.com/CoolProp/CoolProp)
* CMake (platform-independent software to generate makefiles)
* C++ compiler (see below)

Windows
-------
On Windows, download the newest binary installer for CMake from http://www.cmake.org/cmake/resources/software.html.  Run the installer.  Check that at the command prompt you can do::

    C:\Users\XXXX>cmake -version
    cmake version 2.8.12.2

Linux
-----    
On debian based linux distributions (ubuntu, etc.), you can simply do::

    sudo apt-get install cmake git g++
    
although ``git`` is probably already packaged with your operating system; ``g++`` probably isn't

OSX
---
OSX should come with a c++ compiler (clang), for git and cmake your best bet is `Homebrew <http://brew.sh/>`_.  With Homebrew installed, you can just do::

    brew install cmake git

.. toctree::
    :hidden:

    Octave/index.rst
    Csharp/index.rst