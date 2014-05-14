Building
========

Requirements
------------
SWIG (http://www.swig.org/)
A compiler (the Visual Studio 2010 C++ Express version should be fine)

To Build
--------
Run the script build_x64.bat - adjust the paths if necessary to the include folders for your java installation

If on 32-bit windows, run the build_win32.bat file

Each script will put the DLL in the corresponding folder (win32 for 32-bit, x64 for 64-bit)

Running
=======
At the console, run::

    javac *.java
    java runme
    
which should output::

    702.820647167934
    
Hiccups
=======
If the bin folder of the installation for java is not on the path, you may need to add it.

