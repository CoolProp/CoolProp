.. _EES:

***********
EES Wrapper
***********

EES is an acausal solver that can be used to solve a wide range of technical problems.  It can be obtained from https://www.fchart.com/ees/.  Though EES has its own set of thermodynamic properties, CoolProp also implements a number of things that are not in EES (incompressibles, interpolation methods, etc.).

Users
=====

32-bit and 64-bit EES
---------------------
EES comes in a 32-bit and a 64-bit flavour, and the two cannot share an external library.  The 32-bit program loads ``COOLPROP_EES.dlf`` and ``CoolProp.LIB`` from its ``Userlib`` folder (by default ``c:\EES32\Userlib``), the 64-bit program ``EES64.exe`` loads ``COOLPROP_EES.dlf64`` and ``CoolProp.LIB64`` from ``Userlib64`` (by default ``c:\EES64\Userlib64``).  Both are built and published, in the ``EES/EES`` and ``EES/EES64`` folders of the downloads, and the Windows installer offers one task per flavour.

Automated Installation
----------------------
As of October 2016, the EES wrapper can be installed from the Windows package as described on the :ref:`page on installation packages <Installers>`. Please refer to the documentation there for issues related to the installation process.

Manual Install
--------------
Pre-compiled binaries can be downloaded from :sfdownloads:`EES` - follow the instructions there.  Or you can download an installer from the development preview at :sfnightly:`EES`.  The ``EES`` folder of the downloads holds both flavours, the 32-bit files in ``EES`` and the 64-bit files in ``EES64``.  To install by hand, copy the four files of the matching flavour into ``c:\EES32\Userlib\COOLPROP_EES`` or ``c:\EES64\Userlib64\COOLPROP_EES`` and restart EES.

Usage
-----
Open EES, you should see the external function COOLPROP_EES.  The Function Information dialog shows an example call for it.  The function ``PropsSI`` takes the same inputs as described in the :ref:`High-Level API <high_level_api>`.  You can use something like::

    xx = string$(0.5)
    yy = string$(0.5)
    fluid$='HEOS::Methane['||yy||']&Ethane['||xx||']'

which is a 50/50 molar blend of methane and ethane.

The function ``PropsSIZ`` takes the normal inputs, but then also takes the mole fractions as an array rather than encoding them in the string.  The example file for EES demonstrates all of these types of inputs

Errors and units
----------------
A call that CoolProp cannot evaluate stops the calculation and shows the CoolProp error message, rather than returning zero and letting the solve continue with that number.  A warning is shown without stopping the calculation.

EES skips its unit check for a ``COOLPROP_EES`` call: the units of the arguments depend on the property keys inside the fluid string, which EES does not pass when it asks an external function for units.  The unit system itself is still checked by the library file, and each function has its own: ``PropsSI`` and ``PropsSIZ`` require K, Pa, J and mass, and the deprecated ``coolprop`` requires K, kPa, kJ and mass.  The other deprecated function, ``coolpropsi``, does not work at all: it calls the external function with an undefined string variable, so use ``PropsSI``.

Debugging
---------
1. Install CoolProp EES wrapper
2. Append ``'$DEBUG'`` to the fluid name
3. Open the log.txt and log_stdout.txt files to see the error.  The wrapper opens them by name, so they are written to the working directory of the EES process, not to the user library folder.

Developers
==========

Requirements
------------
Compilation of the EES wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Additionally, you must download `InnoSetup <https://www.jrsoftware.org/isinfo.php>`_ and add it to the system path.

Build
-----

Once the dependencies are installed, you can run the installer with::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp
    # Make a build folder
    mkdir CoolProp/build && cd CoolProp/build
    # Build the makefile using CMake, 32-bit for the 32-bit EES
    cmake .. -G "Visual Studio 17 2022" -A Win32 -DCOOLPROP_EES_MODULE=ON
    # Make the DLF file
    cmake --build . --target COOLPROP_EES --config Release

The 64-bit library is built from the same sources, only the architecture changes::

    cd .. && mkdir build64 && cd build64
    cmake .. -G "Visual Studio 17 2022" -A x64 -DCOOLPROP_EES_MODULE=ON
    cmake --build . --target COOLPROP_EES --config Release

This creates ``COOLPROP_EES.dlf64`` next to ``CoolProp.LIB64``.  Both bitnesses are built and packaged automatically by the ``COOLPROP_WINDOWS_PACKAGE_INSTALLER`` target, which is what the nightly and the tagged releases run.  ``wrappers/EES/DEVELOPER.md`` describes the interface, the build and the packaging for both flavours.

Low-level debugging
-------------------
To make and use a debug DLL, do (from root of repo)::

    mkdir build && cd build
    cmake .. -G "Visual Studio 10 2010" -DCOOLPROP_EES_MODULE=ON

This will make a Visual Studio Project called CoolProp.sln defaulting to 32-bit build.  Open the visual studio project, for the COOLPROP_EES project:

1. Change the output directory to C:\\EES32\\Userlib\\COOLPROP_EES (this is where the DLF will go)
2. Under debugging, set the command to c:\\EES32\\ees.  You can also set the arguments to the file that you want EES to open
3. Set a breakpoint somewhere that it will get hit (in the COOLPROP_EES function for instance)
4. Right-click on the COOLPROP_EES project, select "Set as StartUp Project"
5. Run the project, it will build and start EES, open your code or call some inputs for EES
6. Debugger should stop at your breakpoint


