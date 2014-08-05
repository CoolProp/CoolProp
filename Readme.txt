CoolProp is a thermophysical property database and wrappers for a selection of programming environments

It offers similar functionality to REFPROP, but CoolProp is open-source and free, with flexible licensing terms

It was originally developed by Ian Bell, currently a post-doc at the University of Liege, in Liege, Belgium.

To see whats new, go to http://coolprop.sourceforge.net/changelog.html

Developers
==========

Light debugging
---------------
1. Install COOLPROP_EES
2. Append '$DEBUG' to the fluid name
3. Open the log.txt and log_stdout.txt files in c:\ees32\userlib\COOLPROP_EES

Hardcore debugging
------------------
To use a debug DLL, do (from root of repo)

mkdir build/EES
cd build/EES
cmake ../.. -G "Visual Studio 10 2010" -DCOOLPROP_EES_MODULE=ON

open the visual studio project, for the COOLPROP_EES project:

1. Change the output directory to C:\ees32\userlib\COOLPROP_EES (this is where the DLF will go)
2. Under debugging, set the command to c:\ees32\ees
3. Set a breakpoint somewhere that it will get hit (in the COOLPROP_EES function for instance)
4. Run the project, it will build and start EES, open your code or call some inputs for EES
5. Debugger should stop at your breakpoint



