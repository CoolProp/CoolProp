To Use
======

Get the dll folder from Mathematica by doing::

    FileNameJoin[{$BaseDirectory, "SystemFiles", "LibraryResources", $SystemID}]
    
If this path doesn't exist, make it.

Put the DLL in this folder. (or any other folder in ``$LibraryPath``)

Open the file example.nb to see how to use.

To Build
========

Get "WolframLibrary.h" from C:\\Program Files\\Wolfram Research\\Mathematica\\9.0\\SystemFiles\\IncludeFiles\\C (or similary for your computer) and copy it to the source folder wrappers/Mathematica

Run the file BuildDLL.bat (requires visual studio 2010 with 64-bit support).  Some paths have been hardcoded and might need to be changed to work properly.
