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

NB: For Linux users, the file libCoolProp.so must be renamed into CoolProp.so and copied into Julia's library folder (e.g. /usr/lib/x86_64-linux-gnu/julia/)

Download the module ``CoolProp.jl`` from :sfdownloads:`sourceforge <Julia>` or the development version from :sfnightly:`the nightly snapshots <Julia>` and place in the same folder as the shared library.
The wrapper should be valid for Julia 0.4 and above. For Julia 0.3 and lower, use the one in the 0.3 folder.

Usage
-----
At the console, do something like this in the folder that contains CoolProp.jl and the shared library::

    julia> import CoolProp
    
    julia> CoolProp.PropsSI("T","P",101325.0,"Q",0.0,"Water")
    373.1242958476879
    
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
