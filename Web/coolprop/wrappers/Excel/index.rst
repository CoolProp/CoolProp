
.. _Excel:

*************
Excel Wrapper
*************

Pre-compiled Binaries for windows
=================================
Pre-compiled release binaries can be downloaded from :sfdownloads:`MicrosoftExcel`.  Development binaries coming from the buildbot server can be found at :sfnightly:`MicrosoftExcel`.

Download all the files. The basic protocol is:

Part 1:
-------
1. Copy the files CoolProp.dll and CoolProp_x64.dll to c:\\CoolProp folder. Technically you only need the DLL that matches your system architecture (CoolProp.dll = 32-bit, CoolProp_x64.dll = 64-bit), but it can't hurt to copy both if you don't know which system architecture version you have.

Part 2:
-------
1. Open Excel
2. Determine the trusted directory for addins by going to File->Options...->Trust Center->Trust Center Settings...->Trusted Locations . Note the add-ins directory. Put your CoolProp.xlam in that directory
3. Go to the menu File-->Options-->Add-Ins
4. At the bottom, select Manage: Excel Add-ins, then click the Go.. button
5. Click the browse button
6. Browse to the file CoolProp.xlam you downloaded, select it (if necessary, CoolProp might already be visible)
7. Make sure the CoolProp Add-in is selected.
8. Open the file TestExcel.xlsx and try to re-evaluate one of the cells - they should work now

If you are trying to over-write the CoolProp xlam, you should:

a) In Excel, determine the location of the CoolProp xlam that it wants to load by going to Addins menu
b) Move the xlam somewhere else temporarily
c) In Excel, uncheck the CoolProp add-in - it might have already warned about the addin not being in the right place
d) Restart Excel. Check that the addin has been removed, or at least the path is not the old path to CoolProp
e) Determine the trusted directory for addins by going to File->Options...->Trust Center->Trust Center Settings...->Trusted Locations . Note the add-ins directory. Put your CoolProp.xlam in that directory
f) In Excel, add the add-in for CoolProp again by checking it and browsing to the correct location
g) Restart Excel, add-in should still be there.

.. note::

    If you are having problems with the path to the XLAM in your worksheet appearing as the complete path to the xlam (but not the correct path), you might need to go into ``Data->Update Links`` in Excel, select CoolProp.xlam, and select "Change Source..." and select the location of your xlam file.  That should update all the links.
    
Pre-compiled Binaries for OSX
=============================

Part 1:
-------
Download pre-compiled release binaries for OSX from :sfdownloads:`shared_library/Darwin/32bit/`.  Development binaries coming from the buildbot server can be found at :sfnightly:`shared_library/Darwin/32bit/`. Place the downloaded file called libCoolProp.dylib in the ${HOME}/lib folder (make this folder if needed).

Part 2:
-------
Download the xlam from :sfdownloads:`MicrosoftExcel` or the development version from :sfnightly:`MicrosoftExcel`.

1. Open Excel 2011 for Mac
2. Go to the menu Tools-->Add-Ins
3. Click the "Select..." button
4. Browse to the file CoolProp.xlam you downloaded, select it
5. Make sure the CoolProp Add-in is selected.
6. Add this code to a cell - it should work::

    =PropsSI("T","P",101325,"Q",0,"Water")

If it doesn't work and you get error number 53, it might be because you have a 64-bit .dylib file and you want a 32-bit .dylib file.  For instance when you run the ``file`` command on your .dylib, you should see something like:

    $ file libCoolProp.dylib
    libCoolProp.dylib: Mach-O dynamically linked shared library i386

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
    cmake ../.. -G "Visual Studio 10" -DCOOLPROP_SHARED_LIBRARY=ON -DCOOLPROP_STDCALL_LIBRARY=ON
    # Make the shared library
    cmake --build . --config Release
    cd ../..
    # Make a build folder for the 64-bit DLL
    mkdir build/64bit && cd build/64bit
    # Build the MSVC project using CMake
    cmake ../.. -G "Visual Studio 10 Win64" -DCOOLPROP_SHARED_LIBRARY=ON 
    # Make the shared library
    cmake --build . --config Release
    cd ../..
    # Copy the generated DLL
    copy build\32bit__stdcall\CoolProp.dll c:\CoolProp
    copy build\64bit\CoolProp.dll c:\CoolProp
