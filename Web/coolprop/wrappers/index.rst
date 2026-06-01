.. _wrappers:

******************
Available Wrappers
******************

CoolProp at its core is a C++ library, but it can be of interest to use this code base from other programming environments.  For that reason, wrappers have been constructed for most of the programming languages of technical interest to allow users to seamlessly interface CoolProp and existing codebases.

There are also installer packages available on the :ref:`page on installation packages <Installers>`.

Downloads and instructions for each wrapper are included in the page for the wrapper given in the table below.

======================================================= ===========================  =====================================================
Target                                                  Operating Systems            Notes
======================================================= ===========================  =====================================================
:ref:`Static library <static_library>`                  linux, OSX, win
:ref:`Shared library (DLL) <shared_library>`            linux, OSX, win              Included in the Windows :ref:`installer <Installers>`
:ref:`Python <Python>`                                  linux, OSX, win              Wrapper is Cython based
:ref:`Octave <Octave>`                                  linux, OSX, win              Wrapper is SWIG based
:ref:`C# <Csharp>`                                      linux, OSX, win              Wrapper is SWIG based
:ref:`VB.net <VBdotNet>`                                Windows only                 Wrapper is SWIG based
:ref:`MATLAB <MATLAB>`                                  linux, OSX, win              Wrapper is SWIG based
:ref:`Java <Java>`                                      linux, OSX, win              Wrapper is SWIG based
:ref:`R <R>`                                            linux, OSX, win              Wrapper is SWIG based
:ref:`Scilab <Scilab>`                                  linux, OSX, win              Wrapper is SWIG based (experimental)
:ref:`Julia <Julia>`                                    linux, OSX, win              
`Modelica <https://github.com/modelica/ExternalMedia>`_ linux, OSX, win
:ref:`PHP <PHP>`                                        linux, OSX, win              Mostly used on linux
:ref:`Javascript <Javascript>`                          cross-platform               Works in all internet browsers
:ref:`Labview <Labview>`                                Windows only
:ref:`Maple <Maple>`                                    linux, OSX, win              CoolProp is included in Maple 2016
:ref:`Mathcad <MathCAD>`                                Windows only
:ref:`SMath Studio <SMath>`                             linux, OSX, win
:ref:`Mathematica <Mathematica>`
:ref:`FORTRAN <FORTRAN>`                                linux, OSX, win
:ref:`EES <EES>`                                        Windows only                 Included in the Windows :ref:`installer <Installers>`
:ref:`Microsoft Excel <Excel>`                          Windows only                 Included in the Windows :ref:`installer <Installers>`
:ref:`LibreOffice <LibreOffice>`                        Windows, linux
:ref:`Delphi & Lazarus <Delphi>`                        linux, OSX, win
:ref:`iOS (iPhone) <ios>`                       
:ref:`Android <Android>`                       
======================================================= ===========================  =====================================================

.. _wrapper_common_prereqs:

Common Wrapper Prerequisites
============================

On all platforms for which CoolProp is supported, the compilation of one of the wrappers requires a few common prerequisites, described here. They are:

* git (to interface with the CoolProp repository at https://github.com/CoolProp/CoolProp)
* python (to generate the header files and convert them to binary file)
* CMake (platform-independent software to generate makefiles)
* C++ compiler
* 7-zip

Windows
-------
On Windows, download the newest binary installer for CMake from `CMake downloads <https://www.cmake.org/cmake/resources/software.html>`_.  Run the installer.  Check that at the command prompt you can do something like::

    C:\Users\XXXX>cmake -version
    cmake version 2.8.12.2

For git, your best best is the installer from https://msysgit.github.io/.  Check that at the command prompt you can do something like::

    C:\Users\XXXX>git --version
    git version 1.9.4.msysgit.0

For 7-zip, download the installer from https://www.7-zip.org/ .  Check that at the command prompt you can do something like::

    C:\Users\XXXX>7z

    7-Zip [64] 9.20  Copyright (c) 1999-2010 Igor Pavlov  2010-11-18

    Usage: 7z <command> [<switches>...] <archive_name> [<file_names>...]
           [<@listfiles...>]

For python, you should be using `Anaconda/Miniconda <https://www.anaconda.com/download>`_ for your python installation.  Or, you can just install `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_, which is sufficient.  You will also need to install the ``pip`` and ``six`` packages, which can be achieved by the following at a command prompt: ``conda install pip six``.

For the C++ compiler, the options are a bit more complicated.  There are multiple (binary incompatible) versions of Visual Studio, as well as G++ ports for windows (MinGW).  Unless you are compiling the python wrappers, you can compile with MinGW.  The modern way to obtain a MinGW-w64 toolchain on Windows is via `MSYS2 <https://www.msys2.org/>`_; after installing it, install the toolchain with ``pacman -S mingw-w64-x86_64-gcc`` and add the ``mingw64\bin`` directory to your PATH.  Install to a path without spaces.

If you want to build 64-bit extensions, install the free Community edition of Visual Studio (the "Desktop development with C++" workload).  Python/Cython is built with a particular MSVC toolset, so the Visual Studio version should match your Python version as shown in the table below.  Since Visual Studio 2015 the C runtime (UCRT) is stable, so VS2015/2017/2019/2022 are mutually compatible for modern Python versions.

+---------------------------+-------------+-------------------------+
| Visual Studio             | Visual C++  | Python Versions         |
+===========================+=============+=========================+
| 2008 Professional         | 9.0         | 2.6, 2.7, 3.0, 3.1, 3.2 |
+---------------------------+-------------+-------------------------+
| 2010 Professional         | 10.0        | 3.3, 3.4                |
+---------------------------+-------------+-------------------------+
| 2015/2017/2019/2022       | 14.x        | 3.5 and newer           |
| Community                 |             |                         |
+---------------------------+-------------+-------------------------+

Otherwise, for wrappers other than Python, you can select a Visual Studio version freely.

The older compilers can co-exist on the path with the newer ones.  Downloads of VS2008 and 2010 are getting harder to find, so this may influence your Python version of choice.  See this `WindowsCompilers <https://wiki.python.org/moin/WindowsCompilers>`_ page at wiki.python.org for more info on the Windows C++ compilers needed for the various Python/Cython versions and where to download them (most for free).

Linux
-----
On debian based linux distributions (ubuntu, etc.), you can simply do::

    sudo apt-get install cmake git g++ p7zip libpython3-dev

although ``git`` is probably already packaged with your operating system; ``g++`` probably isn't.  Python is (probably) included in your distribution, but the headers aren't.  For python, you need the ``six`` package, a ``pip install six`` should do it.

OSX
---
OSX should come with a c++ compiler (clang), for git and cmake your best bet is `Homebrew <https://brew.sh/>`_.  With Homebrew installed, you can just do::

    brew install cmake git p7zip
    
OSX includes a python version, but you should be using `Anaconda/Miniconda <https://www.anaconda.com/download>`_ for your python installation.  Or, you can just install `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_, which is sufficient. For python, you need the ``six`` package, a ``pip install six`` should do it.

If you have never done any command-line compilation before on OSX, chances are that you do not have the utilities needed. Thus you need to first install Xcode: see the description on the page https://guide.macports.org/#installing.xcode . After installing, you need to accept the license by running the following command in the Terminal::

    xcodebuild -license

and explicitly typing "agree" before closing. Then you can use the compiler from the command line.

.. toctree::
    :hidden:

    Android/index.rst
    Octave/index.rst
    Csharp/index.rst
    MATLAB/index.rst
    MathCAD/index.rst
    FORTRAN/index.rst
    PHP/index.rst
    EES/index.rst
    iOS/index.rst
    Java/index.rst
    Javascript/index.rst
    Julia/index.rst
    Labview/index.rst
    Python/index.rst
    LibreOffice/index.rst
    Excel/index.rst
    Maple/index.rst
    Mathematica/index.rst
    Scilab/index.rst
    SMath/index.rst
    StaticLibrary/index.rst
    SharedLibrary/index.rst
    DelphiLazarus/index.rst
    VB.net/index.rst
    R/index.rst
    Installers/index.rst