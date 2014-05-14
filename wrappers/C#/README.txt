Ian Bell, Ph.D. (ian.h.bell@gmail.com)
November 2013

Info
----
These are the files needed for Csharp applications on Windows. A similar build process should be used for 
non-Windows application though no knowledge is available for non-Windows applications

To Use
------

Copy all the .cs files to a location you want.  You will need to have a copy of some version of C#.

At the command prompt, run::

    call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
    csc *.cs /platform:x86
    call Example

where you might need to update the path to visual studio depending on your version installed.  The `/platform:x86` is important to ensure that it uses 32-bit calling conventions to avoid PINVOKE errors!

Alternatively, you can add all the .cs files to a visual studio project.

Build
-----
The build batch file BuildCsharpDLL.bat should be run to generate the C# files. This requires SWIG!

Make sure you keep it in "x86" rather than "Any CPU" architecture so that it uses 32-bit calling conventions.  Otherwise you will get PINVOKE errors!!!
