
.. _shared_library:

********************
Shared Library (DLL)
********************

General Information
===================

Shared libraries are compiled code that can be accessed by other programs.  On windows, they are DLL files, on other platforms .os (linux) or .dylib (OSX)

There are a few things that need to be considered when determining what shared library you should build/use:

* `Calling convention <https://en.wikipedia.org/wiki/Calling_convention>`_: ``__stdcall`` or ``__cdecl`` - only a consideration on 32-bit windows
* Architecture: 32-bit or 64-bit
* Compiler: Visual Studio, Mingw, GCC, clang

Pre-Compiled Binaries
======================

Download the appropriate shared library for your architecture from :sfdownloads:`shared_library`, or the nightly development versions at :sfnightly:`shared_library`.

Users of the Borland compiler might find these instructions useful for using Visual Studio compiled shared libraries with Borland compiler: https://bcbjournal.org/articles/vol4/0012/Using_Visual_C_DLLs_with_CBuilder.htm

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of a shared library requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Windows
-------
On windows, the greatest amount of complexity is experienced.

Your compiler options are:

* `MSYS2 <https://www.msys2.org/>`_ UCRT64 (formerly MinGW, a windows port of GCC)
* Visual Studio, multiple versions can be installed in parallel

By default, cmake will use your most up to date version of MS Visual Studio it finds.

.. warning::
    The old MinGW GCC may have problems building the latest version of CoolProp.  Use a more recent toolchain (for instance via `MSYS2 <https://www.msys2.org/>`_) that ships a modern GCC.

Your calling convention options are:

* 32-bit __stdcall (``-DCOOLPROP_SHARED_LIBRARY=ON -DCOOLPROP_STDCALL_LIBRARY=ON``)
* 32-bit __cdecl (``-DCOOLPROP_SHARED_LIBRARY=ON -DCOOLPROP_CDECL_LIBRARY=ON``)
* 64-bit (``-DCOOLPROP_SHARED_LIBRARY=ON``)

You can select the compiler in the call to cmake below.

.. note::
   If using `MSYS2/UCRT64 <https://www.msys2.org/>`_, perform these steps in a UCRT64 terminal.  If using Visual Studio, you can use a command window (cmd) or a PowerShell terminal.

1. Check out sources and go into the build directory::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir build && cd build

2. Generate the build file.  Here is where it gets a little complicated.

   A. If you use `MSYS2/UCRT64 <https://www.msys2.org/>`_, these are your options:

      For 64-bit DLL (default)::

        cmake .. -DCOOLPROP_SHARED_LIBRARY=ON -G "Ninja" -DCMAKE_BUILD_TYPE=Release

      For 32-bit __stdcall DLL::

        cmake .. -DCOOLPROP_SHARED_LIBRARY=ON -DCOOLPROP_STDCALL_LIBRARY=ON -G "Ninja" -DCMAKE_BUILD_TYPE=Release

      For 32-bit __cdecl DLL::

        cmake .. -DCOOLPROP_SHARED_LIBRARY=ON -DCOOLPROP_CDECL_LIBRARY=ON -G "Ninja" -DCMAKE_BUILD_TYPE=Release

      You can cross-compile by forcing a non-native bitness by using the additional flags ``-DFORCE_BITNESS_32=ON`` and ``-DFORCE_BITNESS_64=ON``.

   .. note::
      The commands above should be executed *within* an **MSYS2/UCRT64** environment window (see `MSYS2 <https://www.msys2.org/>`_ for setup).  This window will remove windows system paths to native Windows tools and, instead, rely on the UCRT64 toolchain and your ``.bashrc`` file. **UCRT64** versions of ``gcc``, ``CMake`` and ``Ninja`` should be installed as part of this tool chain using the ``pacman`` package manager.  The ``Ninja`` generator simplifies the process by choosing the MingW, Ming-w64, or ucrt64 gcc/g++ compiler, depending on the environment window in which the commands are executed.

   B. If using Visual Studio, replace the Visual Studio version with the version that you are using (see Note below).  Running the command ``cmake`` at the command prompt will tell you which generators are supported

      For 64-bit DLL (Watch out for the 64-bit flag with Win64)::

        cmake .. -DCOOLPROP_SHARED_LIBRARY=ON -G "Visual Studio 17 2022" -A x64

      For 32-bit __stdcall DLL::

        cmake .. -DCOOLPROP_SHARED_LIBRARY=ON -DCOOLPROP_STDCALL_LIBRARY=ON -G "Visual Studio 17 2022" -A x64

      For 32-bit __cdecl DLL::

        cmake .. -DCOOLPROP_SHARED_LIBRARY=ON -DCOOLPROP_CDECL_LIBRARY=ON -G "Visual Studio 17 2022" -A x64

   .. note::
      With Visual Studio 2019 or later, the bitness is no longer embedded in the ``-G`` parameter, but specified with an additional ``-A`` parameter followed by either ``Win32`` or ``x64``.  As of Visual Studio 2017, 64-bit is the default.

      Example: for 64-bit,

        -G "Visual Studio 17 2022" -A x64

      or for 32-bit,

        -G "Visual Studio 17 2022" -A Win32


   Since the bitness for Visual Studio is already specified with the ``-A`` switch (*or embedded in the ``-G`` string for older VS versions*), passing an additional flag to force a certain bitness will cause an error and the process will terminate.

3. Do the build::

    cmake --build . --config Release

   .. note::
      Alternatively, if using MSYS2 UCRT64, you can simply issue the comman: ``ninja``
|

Linux & OSX
-----------

On linux and OSX there is no calling convention to worry about, only options are 32-bit and 64-bit compilation. Also here you can force cross-compilation using ``-DFORCE_BITNESS_32=ON`` and ``-DFORCE_BITNESS_64=ON``.

For 32-bit compilation::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir build && cd build
    # Generate builder
    cmake .. -DCOOLPROP_SHARED_LIBRARY=ON -DFORCE_BITNESS_32=ON -DCMAKE_BUILD_TYPE=Release
    # Build
    cmake --build .

For 64-bit compilation::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir build && cd build
    # Generate builder (defaults to 64-bit on 64-bit machine)
    cmake .. -DCOOLPROP_SHARED_LIBRARY=ON -DCMAKE_BUILD_TYPE=Release
    # Build
    cmake --build .

.. note::

    In order to get more verbose debugging of the build process, see build flags actually used, etc., you can add the flag ``-DCMAKE_VERBOSE_MAKEFILE=ON`` to the generator step above

On Linux, installation could be done by::

    # Change "32" to match your system bitness
    sudo cp libCoolProp.so /usr/local/lib/libCoolProp.so.32.:version:
    pushd /usr/local/lib
    sudo ln -sf libCoolProp.so.32.:version: libCoolProp.so.5
    sudo ln -sf libCoolProp.so.5 libCoolProp.so
    popd


Using
=====

Windows
-------

Here is a small example for calling the shared library from C on windows, as contributed by Philipp Rollmann, Guentner::

    #include "windows.h"
    #include "stdio.h"
    int main(){
        // Define DLL functions
        typedef double (WINAPI *Props1SI)(char Refrigerant[20], char PropertyToReturn[20]);
        typedef double (WINAPI *PropsSI)(char PropertyToReturn[20], char InputProperty1[20], double InputValue1, char InputProperty2[20], double InputValue2, char Refrigerant[20]);

        // addresses
        Props1SI Props1SIAddress;
        PropsSI PropsSIAddress;
        double result1, result2;

        // load DLL; change this path as needed
        HINSTANCE CoolPropDll = LoadLibraryA("C:\\CoolProp\\CoolProp.dll");

        if (CoolPropDll)
        {
             // addresses
             Props1SIAddress = (Props1SI) GetProcAddress(CoolPropDll, "_Props1SI@8");
             PropsSIAddress = (PropsSI) GetProcAddress(CoolPropDll, "_PropsSI@32");
             // call function
             if (Props1SIAddress && PropsSIAddress)
             {
                    result1 = (*Props1SIAddress) ("R410A", "Tcrit");
                    printf("R410A Tcrit: %g\n", result1);
                    result2 = (*PropsSIAddress) ("Dmass", "T", 298.15, "P", 101325, "R410A");
                    printf("R410A density: %g\n", result2);
             }

             // unload DLL
             FreeLibrary(CoolPropDll);
        }
        else{
            printf("Could not load CoolProp DLL.");
        }
    }

Here is another snippet of using the shared library in windows when (for your application), you MUST use a Visual Studio 32-bit stdcall dll of CoolProp for compatibility with other tools::

    // This is to get all the function prototypes from the header
    #define EXPORT_CODE extern "C" __declspec(dllimport)
    #define CONVENTION __stdcall
    #include "CoolProp/CoolPropLib.h"
    #undef EXPORT_CODE
    #undef CONVENTION

    #include <iostream>
    //---------------------------------------------------------------------------
    int main(int argc, char* argv[])
    {
        std::cout << PropsSI("T","P",101325,"Q",0,"Water") << std::endl;
        return 1;
    }


Linux
-----

Based on the discussion on `GitHub <https://github.com/CoolProp/CoolProp/issues/1600>`_, you can use the following steps to link against the CoolProp libraries::

    cat <<EOF > main.cpp
    #include <iostream>
    #include "CoolProp/CoolPropLib.h"
    int main(){
        double T{293.15};
        double P{1e5};
        //double res{0.0};
        double res{PropsSI("D", "T", T, "P", P, "Water")};
        std::cout << "Density: " << res << std::endl;
    }

    EOF
    g++ -std=c++11 -Wall -O2 -o main -DCOOLPROP_LIB -I../include main.cpp libCoolProp.so -ldl
    LD_LIBRARY_PATH=`pwd` ./main


