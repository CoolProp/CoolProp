.. _MSYS2:

***********
MSYS2 Setup
***********

.. image:: https://www.msys2.org/logo.png
   :alt: MSYS2 Logo
   :width: 75px
   :align: left

This is not a wrapper in itself, but guidance on using the MSYS2 toolchain to compile libraries and wrappers on Windows.

Background
==========

Prior versions of CoolProp supported **MINGW** as a Windows port of gcc.  However, the original **MINGW** has been frozen since 2013 and is no longer compatible.  Mingw-x64 (*with 32 and 64-bit compilers*) was forked form MinGW in 2007 to add 64-bit support and coverage of newer Windows APIs absent from the original project.  `MSYS2 <https://www.msys2.org/>`__ followed in the 2010s, building on Mingw-x64 while adding the ``pacman`` package manager and a modernized distribution model. There is also a `VS Code extension <https://code.visualstudio.com/docs/cpp/config-mingw#_prerequisites>`_ facilitating use of **VS Code** as the IDE to run the g++ compiler and GDB debugger.  `MSYS2 <https://www.msys2.org/>`__ is the best place to download this set of GNU compiler tools as well as help on

* Installing and setting up the MSYS2/UCRT64 environment
* Terminals for command line operations in the UCRT64 environment
* Compilers and toolchains, including the most recent gcc version
* A UCRT64 specific version of CMake and other compiler tools
* Ninja CMake generator for auto-selection of toolchain based on current terminal environment


Installation
============

Install **MSYS2** from `MSYS2.org <https://www.msys2.org/>`_.

* Requires 64-bit Windows 10 (1809+) or newer.  Support for Windows 8.1 was dropped in February, 2026.
* Install to ``C:\msys64``` only.  This will make life much easier!
* When given a choice, choose the **UCRT64** environment!

  * The **UCRT** C standard library is newer and also used by MS Visual Studio by default.  It has better compatibility with MSVC, both at build and runtime, and is for modern C++ development.
  * The **MSYS** environment uses the cygwin C Library and is really only for package management.
  * There is a **MINGW64** environment, but it is being phased out of **MSYS2**.
  * The **CLANG64** or **CLANGARM64** environments can be used (especially if compiling on ARM architecture), but have not been tested with CoolProp (*yet*).

In addition to installation, read through the `MSYS2 Documentation <https://www.msys2.org/docs/what-is-msys2/>`_.  There is a lot of information there that will help with environment and tools setup.

Adding Compiler and Tools
=========================

Instead of MinGW or Ming-w64, install the **UCRT64** compiler toolchain::

    pacman -S --needed base_devel mingw-w64-ucrt-x86_64-toolchain

This will install:

* g++
* gcc
* ld
* standard C++ libraries
* modern Windows CRT (UCRT) libraries

.. note::
   **MinGW** is no longer supported and **Ming-w64** is being phased out of MSYS2 in favor of **UCRT64**.

Tools Needed for the CoolProp Build Process
======================================================

Other than the compilers, there are a few additional tools that will be needed.  These can be installed separately with the pacman package manager, or all at once by listing them all on the packman package list.

**Install UCRT64 CMake**

MSYS2 has its own version of cmake that is **UCRT64** aware.  Install it with::

   pacman -S --needed mingw-w64-ucrt-x86_64-cmake


**Install UCRT64 Git**

This is practically identical to the standalone "Git for Windows", but it is neatly bundled and pre-configured by the MSYS2 team to seamlessly integrate with your UCRT64 pathing and shell settings.  Install with::

    pacman -S --needed mingw-w64-ucrt-x86_64-git

.. note::
   MSYS2 has its own **git** package primarily because standard Windows programs and Unix/Linux programs handle file paths, lines endings, and system processes fundamentally differently.


**Install UCRT64 Ninja**

MSYS2 **Ninja** is a compact, ultra-fast, low-level build system optimized strictly for execution speed. On `MSYS2 <https://www.msys2.org/>`__.  Install it with::

    pacman -S --needed mingw-w64-ucrt-x86_64-ninja


.. note::
   Ninja is actually the default backend for the MSYS2 version of CMake. If you do not pass a -G flag to CMake, it automatically looks for a Ninja installation.

Verify the Installation
=======================

Open an ``MSYS2/MSYS2 UCRT64`` terminal from the Windows start menu and type::

  which gcc g++ cmake git

all found paths should be in ``/ucrt64/bin/``. This means they are all in the ``UCRT64`` environment path.

Keeping MSYS2 Up-to-Date
========================

MSYS2 uses a global update process, rather than a package update, in order to ensure resolution of package dependencies.  Periodically, the MSYS2 environment should be updated with the command::

  pacman -Syu

UCRT64 Shell Environment
========================

**Always, always, always** build CoolProp in a ``UCRT64`` terminal (see `MSYS2 Terminals <https://www.msys2.org/docs/terminals/>`_ on setting this up).  Inside a UCRT64 Shell, the normal Windows environment variable and paths are completely replaced with a UCRT64 environment.  This means that running git, cmake, or g++ will use the MSYS2 UCRT64 versions of these codes instead of the Windows native or Cygwin versions if they are installed on the machine.

Pro-Tip:
========

1. If using VS Code as your IDE/Text Editor, add the `MSYS2 recommended json <https://www.msys2.org/docs/terminals/>`_ so that the ``MSYS2 UCRT64`` profile is available when launching a terminal window.
2. Add the path to VS Code to your .bashrc file so that it is available within a UCRT64 Shell and can be launched by typing ``code`` at the shell prompt.
