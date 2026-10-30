.. _Julia:

*************
Julia Wrapper
*************

.. note::

    The Julia wrapper is an official Julia package, registered in the Julia General registry.  It is maintained outside of this repository at `CoolProp/CoolProp.jl <https://github.com/CoolProp/CoolProp.jl>`_; use that repository to file issues and pull requests.  See also the `CoolProp page on JuliaHub <https://platform.juliahub.com/ui/Packages/General/CoolProp>`_.

.. warning::

    The version of CoolProp inside the Julia package may lag behind the latest CoolProp release.  The compiled library is built and registered separately by the Julia community (see :ref:`julia_binaries`), and the wrapper must then allow the new library version.  See :ref:`julia_check_version` to find out which version you have.

There are two fundamental options for interfacing Julia and CoolProp. Choose one option:

A. Install the registered ``CoolProp`` package.

   - Pure Julia code calling a compiled library that is downloaded for you automatically.
   - Very fast.
   - Nothing to download, compile, or copy by hand.

B. Use the ``PyCall`` package to call the Python interface.

   - Needs your own Python installation with CoolProp.
   - Gives access to everything in the Python interface, and to whichever CoolProp version your Python has (including development builds).

Option A: The CoolProp Julia package
====================================

Step 1: Get Julia
-----------------

Julia 1.3 or later is required.  The package has been tested with recent Julia versions.  Skip this step if Julia is already installed, or if you will use JuliaHub (below).

Windows
^^^^^^^

Install Julia from the Microsoft Store, using either the `MSIX App Installer <https://install.julialang.org/Julia.appinstaller>`_ or this command in a command prompt:

.. code-block:: bat

    winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore

A GUI installer is also available from the `Julia downloads page <https://julialang.org/downloads/>`_.  Start Julia from the Start menu, or type ``julia`` in a terminal.

Linux and macOS
^^^^^^^^^^^^^^^

Run this in a terminal.  It installs the latest stable Julia and the ``juliaup`` tool for managing Julia versions:

.. code-block:: bash

    curl -fsSL https://install.julialang.org | sh

Then start Julia by typing ``julia``.

JuliaHub (in the cloud)
^^^^^^^^^^^^^^^^^^^^^^^

`JuliaHub <https://juliahub.com>`_ runs Julia in a web browser, with nothing to install on your computer.  JuliaHub has a free tier with a monthly allowance of compute hours; see the JuliaHub site for its current terms.

1. Sign in at `platform.juliahub.com <https://platform.juliahub.com/ui/Home>`_.
2. On the **Home** menu, under **Applications**, launch **Julia IDE**, choose the machine size and a time limit, and start it.
3. When the instance is ready, click **Connect**.  A VS Code environment opens in your browser.
4. Press Ctrl+Shift+P and run **Julia: Start REPL**, then continue with Step 2 in that REPL.

The menu names may change as JuliaHub updates its site; see the `JuliaHub documentation <https://help.juliahub.com/juliahub/dev/tutorials/julia_ide/>`_ for current instructions.

**Alternatively**, put the commands in a Pluto Notebook:

1. Sign in as indicated above.
2. On the **Home** menu, under **Applications**, launch **Pluto.jl**.
3. When the instance is ready, click **Connect**.  A Pluto.jl notebook will open in your browser and you can begin creating notebook cells.

   - Enter markdown text in cells to document the notebook calculations.
   - Enter the julia commands to make CoolProp calls.
   - Only one command allowed per cell.

Step 2: Add the CoolProp package
--------------------------------

In the Julia REPL, add the package (this is done once):

.. code-block:: julia

    using Pkg
    Pkg.add("CoolProp")

This works the same on Windows, Linux, and JuliaHub.  The compiled CoolProp library is downloaded automatically, and no separate shared library needs to be downloaded, renamed, or placed on a search path.  The first load takes a short time while the package precompiles.

If you want unit-aware calls (see below), also add Unitful:

.. code-block:: julia

    Pkg.add("Unitful")

Step 3: Use it
--------------

The functions are described in :ref:`high_level_api`.  Load the package and call them directly.

High level interface, the boiling point of water at atmospheric pressure:

.. code-block:: julia

    julia> using CoolProp

    julia> PropsSI("T", "P", 101325.0, "Q", 0.0, "Water")
    373.12429584766636

Humid air, the enthalpy at 20 °C (293.15 K), 50 % relative humidity, and 101325 Pa:

.. code-block:: julia

    julia> HAPropsSI("H", "T", 293.15, "R", 0.5, "P", 101325.0)
    38622.838923867734

Unit-aware calls, the saturation pressure of water at 100 °C.  The result carries its unit:

.. code-block:: julia

    julia> using CoolProp, Unitful

    julia> PropsSI("P", "T", 100u"°C", "Q", 0, "Water")
    101417.99666001559 Pa

Low level interface, the specific heat of water at 101325 Pa and 300 K:

.. code-block:: julia

    julia> handle = AbstractState_factory("HEOS", "Water")
    0

    julia> AbstractState_update(handle, get_input_pair_index("PT_INPUTS"), 101325, 300)

    julia> AbstractState_keyed_output(handle, get_param_index("C"))
    4180.635776555601

    julia> AbstractState_free(handle)

The last digits of results may differ slightly between computers and CoolProp versions.

.. _julia_check_version:

Which CoolProp version do I have?
---------------------------------

Ask the library itself, and look at the package versions:

.. code-block:: julia

    julia> using CoolProp, Pkg

    julia> get_global_param_string("version")
    "7.2.0"

    julia> Pkg.status("CoolProp_jll"; mode=Pkg.PKGMODE_MANIFEST)
    Status `~/.julia/environments/v1.13/Manifest.toml`
      [692b3bcd] JLLWrappers v1.8.0
    ⌅ [3351c21f] CoolProp_jll v7.2.0+0
      [56f22d72] Artifacts v1.11.0
      [8f399da3] Libdl v1.11.0
    Info Packages marked with ⌅ have new versions available but compatibility constraints restrict them from upgrading. To see why use `status --outdated -m`

.. note::

    The ``mode`` option is needed because ``CoolProp_jll`` is installed as a dependency of ``CoolProp`` rather than added directly; without it, ``Pkg.status`` reports "No Matches".  The path on the first line depends on your Julia version and environment.

.. note::

    When a newer ``CoolProp_jll`` exists but cannot be installed because the ``CoolProp`` package does not yet allow it, ``Pkg.status`` marks it with a ⌅ symbol (Julia 1.8 and newer; older versions of Julia do not show this marker, so compare the version shown with the latest ``CoolProp_jll`` release in the `CoolProp_jll repository <https://github.com/JuliaBinaryWrappers/CoolProp_jll.jl>`_ instead).  ``Pkg.update()`` picks up newer versions as they become allowed.

.. _julia_binaries:

Where the compiled library comes from
-------------------------------------

The ``CoolProp`` package does not contain the CoolProp library itself.  It depends on the ``CoolProp_jll`` package, which holds prebuilt binaries.  Those binaries are built from the CoolProp source release by the community build system `Yggdrasil <https://github.com/JuliaPackaging/Yggdrasil/tree/master/C/CoolProp>`_ and published as `CoolProp_jll <https://github.com/JuliaBinaryWrappers/CoolProp_jll.jl>`_.  Binaries are provided for Windows, Linux (glibc and musl), macOS, and FreeBSD, on both 64-bit and several other processor types.

Julia chooses the binary that matches your computer.  The package has no setting to use your own shared library instead.  To use a development or custom build of CoolProp from Julia, use Option B below with a Python that has that build installed.

Option B: Using the PyCall package in Julia
===========================================

You must have Python installed on your computer, and have CoolProp installed in that Python installation.  Use the full path to the Python executable that you will set as ``PYTHON`` (below), so that CoolProp is installed into the same Python that PyCall uses.  The paths shown are only examples; use your own.

On Windows:

.. code-block:: bat

    C:\Users\me\AppData\Local\Python\python.exe -m pip install CoolProp

On Linux or macOS:

.. code-block:: bash

    /usr/bin/python3 -m pip install CoolProp

.. note::

    If the path contains spaces (for example ``C:\Program Files\...``), put it in quotes on the command line.

    On many Linux distributions and on Homebrew, the system Python is "externally managed" and ``pip`` refuses to install into it with an ``externally-managed-environment`` error.  The ``/usr/bin/python3`` that comes with macOS is also an old stub.  In these cases, create a virtual environment (``python3 -m venv ~/coolprop-venv``), install CoolProp with ``~/coolprop-venv/bin/python -m pip install CoolProp``, and set ``PYTHON`` to ``~/coolprop-venv/bin/python`` (written out as a full path).

Tell Julia which Python to use, then add the PyCall package.  Set ``PYTHON`` to the full path of that same Python (the path below is only an example):

.. code-block:: julia

    julia> ENV["PYTHON"] = raw"C:\Users\me\AppData\Local\Python\python.exe"

    julia> using Pkg

    julia> Pkg.add("PyCall")

    julia> Pkg.build("PyCall")

On Linux or macOS, use a path such as ``/usr/bin/python3``.  If PyCall was already installed with a different Python, the ``Pkg.build`` step switches it to the new one.

Usage:

.. code-block:: julia

    julia> using PyCall  # Enable PyCall package

    julia> CP = pyimport("CoolProp.CoolProp")  # Import CoolProp module

    julia> Tbp = CP.PropsSI("T", "P", 101325.0, "Q", 0.0, "Water")  # Call some CoolProp properties
    373.12429584766636

    julia> CP.get_global_param_string("version")  # Which CoolProp version does Python have?
    "8.0.0"

The version reported is the one installed in your Python, which can be newer than the one in the Option A package.
