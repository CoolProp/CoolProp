
.. _mathcad:

***************
MathCAD Wrapper
***************

Pre-Compiled Binaries
=====================

To be continued.


User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the static library requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Compile
-------

To build::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir -p build && cd build
    # Build the makefile using CMake
    cmake .. -DCOOLPROP_MATHCAD_MODULE=ON -DCOOLPROP_MATHCAD_ROOT="C:/Program Files/PTC/Mathcad Prime 3.0" -G "Visual Studio 10 2010 Win64" -DCMAKE_VERBOSE_MAKEFILE=ON
    # Make the static library
    cmake --build . --config Release
    
The .dll will be in the Release folder, copy it to ``C:/Program Files/PTC/Mathcad Prime 3.0/Custom Functions``.