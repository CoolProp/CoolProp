Building
====================

There are a number of ways to build the library:

Directly with Make
-------------------

There is a non-windows makefile in the wrappers/SharedLibrary that can be used to make the library.


With CMAKE
-----------

The platform independent cmake program can be used. In using cmake, CoolProp uses the standard procedure:

1) Make a build directory and cd to build
2) cmake ..
3) make
4) make 

Although step 1 and 2 does depend on the OS, step 3 and 4 will depend on how the program is being built.
If on windows the cmake gui can be used which will perform steps 1 and 2, if using visual studio then
the build process will need to be run twice.

Make needs to be called twice, the first make step will dynamically generate a number of files from the
JSON fluid definitions - the second make run will actually generate the program.

Testing
-------

CMake generates a target for testing. You can build the test executable with `make testRunner`.

Building wrappers using CMake
-----------------------------

To build the 32-bit __cdecl DLL for instance, from the root directory execute

```
mkdir build
cd build
mkdir 32bitcdecl
cd 32bitcdecl
cmake ..\.. -G "Visual Studio 10" -DCOOLPROP_32BIT_CDECL_SHARED_LIBRARY=ON
cmake --build . --config Release --target INSTALL
```

which will build the DLL and put it in the bin/shared_library/Windows/32bit/__cdecl_calling_convention folder.  Alternatively you could build using MINGW GCC using

```
mkdir build
cd build
mkdir 32bitcdeclgcc
cd 32bitcdeclgcc
cmake ..\.. -G "MinGW Makefiles" -DCOOLPROP_32BIT_CDECL_SHARED_LIBRARY=ON
cmake --build . --config Release --target install
```

which will also be the same protocol on linux