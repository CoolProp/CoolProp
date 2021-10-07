.. _Csharp:

**********
C# Wrapper
**********

.. contents:: :depth: 2

NuGet packages (3-party wrappers)
=================================

SharpFluids
-----------

This C# NuGet package uses CoolProp to perform all Fluid Properties lookups. It combines the speed from the low-level lookup with a units of measurement system packed into a easy-to-use system. If you are new to using CoolProp, this is a good place to start.

How to start 

- Create a new C# Console App(.NET Framework) project in Visual studio
- Right click your new project and press 'Manage NuGet Packages'
- Go to 'Browse' and search for 'SharpFluids' and press 'Install'
- Add this to the top of your code :: 

    using SharpFluids;
    using UnitsNet;

- Add this example code in your main ::

    Fluid Water = new Fluid(FluidList.Water);
    Water.UpdatePT(Pressure.FromBars(1.013), Temperature.FromDegreesCelsius(13));
    Console.WriteLine("Density of water at 13°C: " + Water.RHO);
    Console.ReadLine();

- Press 'F5' or 'Start' to check that it is working

- If you have problems or questions, `Find SharpFluids at Github <https://github.com/MadsKirkFoged/SharpFluids>`_.

SharpProp
---------

It is a simple, full-featured, lightweight CoolProp wrapper for C#. SharpProp gets published on `NuGet <https://www.nuget.org/packages/SharpProp/>`_.

All CoolProp features are included: thermophysical properties of pure fluids, mixtures and humid air (all in *SI units*).
Also you can easily convert the results to JSON, add new properties or inputs for lookups, and more.

Examples
^^^^^^^^

Don't forget to add ``using SharpProp;`` at the top of the code.

To calculate the specific heat of saturated water vapour at *101325 Pa*: ::

    var waterVapour = new Fluid(FluidsList.Water);
    waterVapour.Update(Input.Pressure(101325), Input.Quality(1));
    Console.WriteLine(waterVapour.SpecificHeat); // 2079.937085633241

To calculate the dynamic viscosity of propylene glycol aqueous solution with *60 %* mass fraction at *101325 Pa* and *253.15 K*: ::

    var propyleneGlycol = new Fluid(FluidsList.MPG, 0.6);
    propyleneGlycol.Update(Input.Pressure(101325), Input.Temperature(253.15));
    Console.WriteLine(propyleneGlycol.DynamicViscosity); // 0.13907391053938847

To calculate the density of ethanol aqueous solution (with ethanol *40 %* mass fraction) at *200 kPa* and *277.15 K*: ::

    var mixture = new Mixture(new List<FluidsList> {FluidsList.Water, FluidsList.Ethanol}, new List<double> {0.6, 0.4});
    mixture.Update(Input.Pressure(200e3), Input.Temperature(277.15));
    Console.WriteLine(mixture.Density); // 883.3922771627759

To calculate the wet bulb temperature of humid air at *99 kPa*, *303.15 K* and *50 %* relative humidity: ::

    var humidAir = new HumidAir();
    humidAir.Update(InputHumidAir.Pressure(99e3), InputHumidAir.Temperature(303.15),
        InputHumidAir.RelativeHumidity(0.5));
    // or use:
    // var humidAir = HumidAir.WithState(InputHumidAir.Pressure(99e3), InputHumidAir.Temperature(303.15),
    //     InputHumidAir.RelativeHumidity(0.5));
    Console.WriteLine(humidAir.WetBulbTemperature); // 295.0965785590792

For any questions or more examples, `see SharpProp on GitHub <https://github.com/portyanikhin/SharpProp>`_.

Pre-compiled Binaries
=====================

To Use
------

Pre-compiled binaries can be downloaded from :sfdownloads:`Csharp`.  Development binaries coming from the buildbot server can be found at :sfnightly:`Csharp`.

Download the ``platform-independent.7z`` file and expand it to a folder called ``platform-independent`` using 7-zip.  Download the special C# shared library for your system architecture to the same location from either :sfdownloads:`Csharp` (release) or :sfnightly:`Csharp` (development).  Copy the Example.cs file to the same location.  You will need to have a copy of some version of C#.

When you are finished, you should have a folder layout something like ::

    main
     |- CoolProp.dll
     |- Example.cs
     |- platform-independent
        |- AbstractState.cs
        |- Configuration.cs
        |- ...
        
There is example code :ref:`at the end of this page <csharp_example>`

Windows
^^^^^^^

At the command prompt, run::

    call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
    csc Example.cs platform-independent/*.cs -platform:x64
    Example

where you might need to update the path to visual studio depending on your version installed.  Use `-platform:x86` to tell C# that your DLL is 32-bit if you are on 32-bit, or `-platform:x64` if you are on 64-bit.

Alternatively, you can add all the .cs files to a visual studio project.  If you do that, add the DLL to the project as well, right-click on the DLL, and select the option to copy it to the output directory.

Linux/OSX
^^^^^^^^^

Same idea as windows, but command line is just a bit different::

    mcs Example.cs platform-independent/*.cs -platform:x64
    ./Example
    
Use `-platform:x86` to tell C# that your shared library is 32-bit if you are on 32-bit, or `-platform:x64` if you are on a 64-bit platform.

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the C# wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Additionally, you will need:
* SWIG (see :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`)
* C#

OSX
---

For OSX, to install the necessary tools using homebrew, you can do::
    
    homebrew install mono

Linux
-----

For ubuntu and friends, you will need to install Mono C# as well as the compiler (and other dependencies) using::

    sudo apt-get install swig mono-mcs mono-runtime

Windows
-------
For Windows, download the Visual Studio 2010 version of C# (other versions should probably be fine too)

Compile
-------

Once mono c# is installed, you can run the builder and tests using::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir build && cd build
    # Build the makefile using CMake
    cmake .. -DCOOLPROP_CSHARP_MODULE=ON -DBUILD_TESTING=ON
    # Make the C# files (by default files will be generated in folder install_root/Csharp relative to CMakeLists.txt file)
    cmake --build .
    # Run the integration tests (optional)
    ctest --extra-verbose

If you want to change the package that CoolProp resides in, you can do so by changing the cmake call to read::

    cmake .. -DCOOLPROP_CSHARP_MODULE=ON -DBUILD_TESTING=ON -DCOOLPROP_SWIG_OPTIONS="-namespace package.name"

where ``package.name`` is replaced with the desired name    
    
.. _csharp_example:

Example Code
============

.. literalinclude:: Example.cs
   :language: csharp

Example Code Output
===================

.. literalinclude:: Example.out
