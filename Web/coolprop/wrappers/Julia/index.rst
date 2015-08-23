.. _Julia:

*************
Julia Wrapper
*************

There are two fundamental options for interfacing Julia and CoolProp. Choose one option:

* Call the shared library directly using only Julia code (very fast!)
* Get access to all the python interface using ``PyCall`` package

Option A: Calling shared library directly
=========================================

.. warning::

    For now, only 64-bit architectures are supported.
    
Download a precompiled shared library appropriate to the computer you are using from :sfdownloads:`sourceforge <shared_library>` or the development version from :sfnightly:`the nightly snapshots <shared_library>`.  

NB: For Linux users, the file libCoolProp.so must be renamed into CoolProp.so and copied into Julia's library folder (e.g. /usr/lib/x86_64-linux-gnu/julia/).
If you compile/download the library often, you may prefer to make a symbolic link from your downloaded/compiled library to Julia's library folder::

    sudo ln -s libCoolProp.so /usr/lib/x86_64-linux-gnu/julia/CoolProp.so
    
If you want to use the library in an other directory, once renamed (or after creating a `CoolProp.so` symbolic link) you can set the system shared library search path to the current directory before running Julia::

    export LD_LIBRARY_PATH=${PWD}

Or directly ask Julia to search in the current directory::

    push!(DL_LOAD_PATH,".")

Note that you can replace `${PWD}` (present working directory) by any path you want, and similarly for the `.` of the `"."` in Julia.

Download the module ``CoolProp.jl`` from :sfdownloads:`sourceforge <Julia>` or the development version from :sfnightly:`the nightly snapshots <Julia>` and place in the same folder as the shared library (in `~/.julia/$version/CoolProp/src/` for Linux users).

Usage
-----
At the console, do something like this in the folder that contains CoolProp.jl and the shared library:

High level interface::

    julia> import CoolProp
    
    julia> CoolProp.PropsSI("T","P",101325.0,"Q",0.0,"Water")
    373.1242958476879

Low level interface::

    julia> handle = CoolProp.AbstractState_factory("HEOS", "Water")
    0
    
    julia> CoolProp.AbstractState_update(handle,CoolProp.get_input_pair_index("PT_INPUTS"),101325, 300)
    
    julia> CoolProp.AbstractState_keyed_output(handle,CoolProp.get_param_index("C"))
    4180.635776569655

    julia> CoolProp.AbstractState_free(handle)
    
Option B: Using PyCall package in Julia
=======================================

You must have python installed on your computer, and have CoolProp installed in your python installation.

Open julia in console and check for installed packages::
    
    # 2 required packages:
    # - PyCall            0.4.8
    # - PyPlot            1.2.9    
    julia> Pkg.status()

If required packages are not installed, use the following commands::

    julia> Pkg.add("PyCall")
    julia> Pkg.add("PyPlot")
    julia> Pkg.update()

Usage::

    # Enable PyCall package:
    julia> using PyCall

    # Import CoolProp module:
    julia> @pyimport CoolProp.CoolProp as CP

    # Call some CoolProp properties:
    julia> Tbp = CP.PropsSI("T","P",101325.0,"Q",0.0,"Water")

User-Compiled Binaries
======================

Build the 64-bit shared library for your architecture following the instructions at :ref:`shared_library`.
