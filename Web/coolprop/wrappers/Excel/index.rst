
.. _Excel:

*************
Excel Wrapper
*************

Pre-compiled Binaries for windows
=================================
Pre-compiled release binaries can be downloaded from :sfdownloads:`MicrosoftExcel`.  Development binaries coming from the buildbot server can be found at :bbbinaries:`MicrosoftExcel`.

Download all the files. The basic protocol is:

Part 1:
-------
1. Copy the files CoolProp.dll and CoolProp_x64.dll to c:\\CoolProp folder. Technically you only need the DLL that matches your system architecture (CoolProp.dll = 32-bit, CoolProp_x64.dll = 64-bit), but it can't hurt to copy both if you don't know which system architecture version you have.

Part 2:
-------
1. Open Excel
2. Go to the menu File-->Options-->Add-Ins
3. At the bottom, select Manage: Excel Add-ins, then click the Go.. button
4. Click the browse button
5. Browse to the file CoolProp.xlam you downloaded, select it
6. Make sure the CoolProp Add-in is selected.
7. Open the file TestExcel.xlsx and try to re-evaluate one of the cells - they should work now

Pre-compiled Binaries for OSX
=============================

Part 1:
-------
Download pre-compiled release binaries for OSX from :sfdownloads:`shared_library/Darwin/32bit__cdecl_calling_convention/`.  Development binaries coming from the buildbot server can be found at :bbbinaries:`shared_library/Darwin/32bit__cdecl_calling_convention/`. Place the downloaded file called libCoolProp.dylib in the ${HOME}/lib folder (make this folder if needed).

Part 2:
-------
Download the xlam from :sfdownloads:`MicrosoftExcel` or the development version from :bbbinaries:`MicrosoftExcel`.

1. Open Excel 2011 for Mac
2. Go to the menu Tools-->Add-Ins
3. Click the "Select..." button
4. Browse to the file CoolProp.xlam you downloaded, select it
5. Make sure the CoolProp Add-in is selected.
6. Add this code to a cell - it should work::

    =PropsSI("T","P",101325,"Q",0,"Water")

User-compiled Binaries
======================

Common Requirements
-------------------
Compilation of the Excel wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Build
-----

The instructions here are for a 64-bit windows system that will compile both 64-bit and 32-bit versions of the DLL::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder for the 32-bit DLL
    mkdir build/32bit__stdcall && cd build/32bit__stdcall
    # Build the MSVC project using CMake
    cmake ../.. -G "Visual Studio 10" -DCOOLPROP_32BIT_STDCALL_SHARED_LIBRARY=ON
    # Make the shared library
    cmake --build . -C Release
    cd ../..
    # Make a build folder for the 64-bit DLL
    mkdir build/64bit && cd build/64bit
    # Build the MSVC project using CMake
    cmake ../.. -G "Visual Studio 10 Win64" -DCOOLPROP_64BIT_SHARED_LIBRARY=ON
    # Make the shared library
    cmake --build . -C Release
    cd ../..
    # Copy the generated DLL
    copy build\32bit__stdcall\CoolProp.dll c:\CoolProp
    copy build\64bit\CoolProp.dll c:\CoolProp