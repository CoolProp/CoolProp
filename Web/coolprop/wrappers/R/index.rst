.. _R:

************
R Wrapper
************

.. contents:: :depth: 2

Pre-compiled Binaries
=====================
Pre-compiled R binaries are no longer supplied in CoolProp's SourceForge repository.  The last release that included them is 6.4.1 (August 2020), which has files for 64-bit Windows and 64-bit macOS only, never for Linux.  There are simply too many OS, hardware, and R version combinations to keep current.  The R wrapper must be built from source (see below).

.. warning::

    The old 6.4.1 files are still online at https://sourceforge.net/projects/coolprop/files/CoolProp/6.4.1/R/ as the last known binaries.  They were built against R releases from 2020 and may not load in a current R.

Usage
-----
The shared library (``CoolPropR.so`` on Linux, ``CoolPropR.dll`` on Windows, ``CoolPropR.so`` on macOS) and the ``CoolProp.R`` file are placed in the same folder (the build folder, or the folder they were installed to).  From that folder, at the R console, run::

    dyn.load(paste("CoolPropR", .Platform$dynlib.ext, sep=""))
    source("CoolProp.R")
    cacheMetaData(1)
    PropsSI("T","P",101325,"Q",0,"Water")

There is example code :ref:`at the end of this page <r_example>`

.. warning::

    For OSX users that wish to call REFPROP, you may be required to set the environmental variable ``DYLD_LIBRARY_PATH`` to the folder containing your REFPROP.dylib shared library, which is probably ``/opt/refprop``
    
.. note::

    Old versions of ``Rscript`` did not load the ``methods`` package, which SWIG-ed code needs, so ``--default-packages=methods`` or ``library(methods)`` was required.  Current R loads it by default (tested with R 4.5.2 on Ubuntu 26.04, where the steps above work in both ``R`` and ``Rscript``).  The example code still calls ``library(methods)`` for older R versions.  See also https://stackoverflow.com/a/19468533

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the R wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Also, R and its development headers need to be installed, along with SWIG.

Linux
-----

On **Ubuntu** and friends, the build dependencies are installed this way::

    sudo apt-get install r-base-dev swig cmake build-essential git

The ``r-base-dev`` package provides the R headers (``R.h``) and pulls in ``r-base``.  These steps were tested on Ubuntu 26.04 LTS with R 4.5.2, SWIG 4.4.0, CMake 4.2, and GCC 15.2.

.. note::

    The CMake variable ``R_BIN`` is required.  It is the folder that contains the ``R`` executable, ``/usr/bin`` for a distro-packaged R.  CMake stops with an error if it is missing.

Build
-----

Linux
^^^^^

Once the dependencies are installed, the wrapper is built and tested with::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp
    # Move into the folder just created
    mkdir -p  CoolProp/build && cd CoolProp/build
    # Build the makefile using CMake
    cmake .. -DCOOLPROP_R_MODULE=ON -DR_BIN="/usr/bin" -DCMAKE_BUILD_TYPE=Release
    # Make the R shared library (-j runs the compilers in parallel)
    cmake --build . -j8
    # Optional: run the example script through Rscript
    ctest

The build folder then holds ``CoolPropR.so``, ``CoolProp.R``, and ``Example.R``.  ``cmake --install .`` copies them to ``R/`` and ``R/Linux_64bit`` under the install prefix.

.. note::

    The ``ctest`` example calls REFPROP in its last blocks.  Without REFPROP installed those blocks print a "REFPROP example skipped" message and the rest of the example still runs.

OSX
^^^

.. attention::

    The OSX instructions below have not been re-verified recently (no macOS machine was available for testing).  An update PR is welcome.

The commands are the same as for Linux, with ``R_BIN`` pointing to the folder that contains the ``R`` executable of the R installation in use.

Windows
-------

The R wrapper is built with the **MSYS2 UCRT64** toolchain.  Current R for Windows (4.2 and later, from https://cran.r-project.org/bin/windows/base/) is itself built for the UCRT runtime, so the wrapper has to be built with a UCRT toolchain too.  The older TDM-GCC and generic MinGW-w64 toolchains (which used the old ``msvcrt`` runtime) are not recommended.  The Visual Studio generators are not supported (CMake stops with an error).

Install `MSYS2 <https://www.msys2.org/>`_, then, from an *MSYS2 UCRT64* shell, install the tools::

    pacman -S mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-cmake mingw-w64-ucrt-x86_64-ninja mingw-w64-ucrt-x86_64-swig mingw-w64-ucrt-x86_64-python git

``R_BIN`` is the folder of CRAN's R that holds ``R.exe`` and ``R.dll`` (use forward slashes), here ``C:/Program Files/R/R-4.6.2/bin/x64``.  Still in the UCRT64 shell::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp
    cd CoolProp
    # Make a build folder
    mkdir build && cd build
    # Build the makefile using CMake
    cmake .. -G Ninja -DCOOLPROP_R_MODULE=ON -DR_BIN="C:/Program Files/R/R-4.6.2/bin/x64" -DCMAKE_BUILD_TYPE=Release
    # Make the R shared library
    cmake --build .

The shared library is ``CoolPropR.dll``, found in the ``Release`` sub-folder of the build folder (next to ``CoolProp.R`` one level up; copy it next to ``CoolProp.R`` before using the Usage steps above).  It only imports Windows system DLLs and ``R.dll``, so it loads in a plain R without the MSYS2 folders on the ``PATH``.

These steps were tested with MSYS2 GCC 16.2, CMake 4.4, SWIG 4.5.1, and CRAN R 4.6.2 (``R.version$crt`` is ``ucrt``), which printed ``373.1243`` for ``PropsSI("T","P",101325,"Q",0,"Water")``.  Building against the R from MSYS2's own ``mingw-w64-ucrt-x86_64-r`` package has not been tested.

.. _r_example:

Example Code
============

.. literalinclude:: Example.R

Example Code Output
===================

.. literalinclude:: Example.out
