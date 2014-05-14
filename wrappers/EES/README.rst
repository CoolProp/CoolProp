CoolProp for EES

Ian Bell, January 2013

Run the installer to install the CoolProp wrapper for EES.

You can see the docs at `<CoolProp.htm>`_

There is an example file included in the c:\\EES\\Userlib\\CoolProp_EES folder when installed


Debug
-----
1. Create a new DLL project in Visual Studio 2010 (or similar).
2. Add main.cpp and all cpp files in ROOT/CoolProp folder into the project
3. Set the output file to c:\EES32\Userlib\COOLPROP_EES\COOLPROP_EES.dlf
4. Set the include folder to ROOT/CoolProp
5. Set c:\EES32\ees as the command to run in Debugging
6. Build
7. Load your failing EES file