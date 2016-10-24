.. _EES:

***********
EES Wrapper
***********

EES is an acausal solver that can be used to solve a wide range of technical problems.  It can be obtained from http://www.fchart.com/ees/.  Though EES has its own set of thermodynamic properties, CoolProp also implements a number of things that are not in EES (incompressibles, interpolation methods, etc.).

Users
=====

Automated Installation
----------------------
As of October 2016, the EES wrapper can be installed from the Windows package as described on the :ref:`page on installation packages <Installers>`. Please refer to the documentation there for issues related to the installation process.

Manual Install
--------------
Pre-compiled binaries can be downloaded from :sfdownloads:`EES` - follow the instructions there.  Of you can download an installer from the development preview at :sfnightly:`EES`.

Usage
-----
Open EES, you should see the external function COOLPROP_EES.  The function ``PropsSI`` takes the same inputs as described in the :ref:`High-Level API <high_level_api>`.  You can use something like::

    xx = string$(0.5)
    yy = string$(0.5)
    fluid$='HEOS::Methane['||yy||']&Ethane['||xx||']'

which is a 50/50 molar blend of methane and ethane.

The function ``PropsSIZ`` takes the normal inputs, but then also takes the mole fractions as an array rather than encoding them in the string.  The example file for EES demonstrates all of these types of inputs

Debugging
---------
1. Install CoolProp EES wrapper
2. Append ``'$DEBUG'`` to the fluid name
3. Open the log.txt and log_stdout.txt files in c:\\ees32\\userlib\\COOLPROP_EES to see the error.

Developers
==========

Requirements
------------
Compilation of the EES wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Additionally, you must download `InnoSetup <http://www.jrsoftware.org/isinfo.php>`_ and add it to the system path.

Build
-----

Once the dependencies are installed, you can run the installer with::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Make a build folder
    mkdir CoolProp/build && cd CoolProp/build
    # Build the makefile using CMake
    cmake .. -DCOOLPROP_EES_MODULE=ON
    # Make the DLF file
    cmake --build . --target COOLPROP_EES

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


