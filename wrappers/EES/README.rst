CoolProp for EES

Ian Bell, January 2013

Run the installer to install the CoolProp wrapper for EES.

You can see the docs at `<CoolProp.htm>`_

There is an example file included in the c:\\EES32\\Userlib\\CoolProp_EES folder when installed
(c:\\EES64\\Userlib64\\CoolProp_EES for the 64-bit EES)

Bitness
-------
EES cannot load an external library of the wrong bitness.  The 32-bit program
reads COOLPROP_EES.dlf and CoolProp.LIB from Userlib, the 64-bit program
EES64.exe reads COOLPROP_EES.dlf64 and CoolProp.LIB64 from Userlib64.  Both are
built, published and offered by the installer, see `<EES64.md>`_ for the details
and for what still needs a manual check in EES.

Debug
-----
1. Create a new DLL project in Visual Studio (or use the CMake project).
2. Add main.cpp and all cpp files in ROOT/CoolProp folder into the project
3. Set the output file to c:\\EES32\\Userlib\\COOLPROP_EES\\COOLPROP_EES.dlf
   (c:\\EES64\\Userlib64\\COOLPROP_EES\\COOLPROP_EES.dlf64 for the 64-bit build)
4. Set the include folder to ROOT/CoolProp
5. Set c:\\EES32\\ees as the command to run in Debugging (c:\\EES64\\ees64 for
   the 64-bit build)
6. Build
7. Load your failing EES file
