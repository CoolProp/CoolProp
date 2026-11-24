.. _octave_rhel:

****************
RHEL Octave Info
****************

.. contents:: :depth: 2

The Octave packages in RHEL 9's standard and EPEL repositories are typically several releases behind upstream.  Flatpak is the fastest way to get the current GNU Octave (11.3.0 when this page was written) on RHEL 9 without dependency conflicts, using the official `GNU Octave package on Flathub <https://flathub.org/en/apps/org.octave.Octave>`_.  The steps below were tested on RHEL 9.5.

This page covers only what differs from the generic Linux instructions on the :ref:`Octave Wrapper <Octave>` page: installing Octave and SWIG, and building inside the Flatpak sandbox.

Why Flatpak Works Well Here
===========================

The Flathub package bundles the compilers (``gcc``, ``g++``, ``gfortran``, ``make``), CMake, Python, and Octave's own development headers inside its runtime.  Because the wrapper is built on this toolchain, the requirement that ``CoolProp.oct`` be built with the exact same toolchain as Octave (see the :ref:`Octave Wrapper <Octave>` page) is met automatically.  Add-on packages installed with ``pkg install`` are stored under ``~/.var/app/org.octave.Octave/`` and survive Flatpak updates.

There are a few limits to keep in mind:

* **SWIG is not included**, but it is easily added (see below).
* **Host libraries are invisible.**  The sandbox cannot see ``-devel`` headers or libraries installed on the RHEL host.  This does not matter for CoolProp, which is self-contained, but it does for add-ons that need system libraries such as MPI or HDF5.
* **Host folders outside the home directory are invisible.**  The sandbox sees the home directory by default.  Anything else (``/opt``, ``/srv``, ...) has to be shared explicitly::

    flatpak override --user --filesystem=/path/to/dir org.octave.Octave

Building Octave itself from source on RHEL 9 is only worth the effort when an add-on must link against specialized host libraries or an HPC/MPI stack.

Installing Octave
=================

RHEL 9 ships with Flatpak, and Flathub was already enabled on the RHEL 9 machine used for testing.  Check to see if Flathub is enabled with::

    flatpak remotes

If ``flathub`` is not in the list, it can be added once with::

    flatpak remote-add --if-not-exists flathub https://dl.flathub.org/repo/flathub.flatpakrepo

(This command comes from Flathub's own setup instructions, but it was not needed or run during testing.)

Then install Octave from the command line::

    flatpak install -y flathub org.octave.Octave

Alternatively (**Easiest Route**), ``GNOME Software`` can be used.  If Flathub was enabled above close and reopen it (or run ``killall gnome-software``) so the new remote is picked up.  Then search for "GNU Octave", and click *Install*.

.. note::
    GNOME Software installs Flatpaks system-wide rather than per-user.  The ``flatpak run`` and ``flatpak override --user`` commands on this page work either way.

Check the installed version with::

    flatpak run org.octave.Octave --version

Installing SWIG
===============

SWIG is needed inside the Flatpak environment for the CoolProp wrapper build but doesn't have to be built there.  SWIG only generates wrapper source code, so it does not have to match Octave's toolchain and can be built outside the Flatpak sandbox.  It is best built and installed under the home directory, since the sandbox can already see that.

1. Install the build prerequisites on the host (once)::

    sudo dnf install -y gcc-c++ make pcre2-devel

2. Download the SWIG source archive (``swig-<version>.tar.gz``) from `SourceForge <https://sourceforge.net/projects/swig/files/swig/>`_ with a browser.  Pick the release folder (e.g., 4.5.1) and then the ``.tar.gz`` file.

   .. note::
        SourceForge blocks command-line downloads (``curl``, ``wget``), so the file has to be downloaded by hand through the web interface.

3. Unpack, configure, build, and install it under the home directory::

    mkdir -p ~/swig-build
    mv ~/Downloads/swig-<version>.tar.gz ~/swig-build
    cd ~/swig-build
    tar -xf swig-<version>.tar.gz
    cd swig-<version>
    ./configure --prefix=$HOME/swig-build/install
    make -j$(nproc)
    make install

4. Check that it runs::

    ~/swig-build/install/bin/swig -version

Building CoolProp Inside the Flatpak Shell
==========================================

Building with the Flatpak runtime's own ``cmake``, ``g++``, and ``octave-config`` avoids any compiler or ``glibc`` mismatch with the RHEL host, and is the way this page was tested.

.. warning::
    Do not run the build from the host's own ``cmake``.  CoolProp's CMake finds Octave through ``octave-config`` and then looks for Octave's headers and libraries on the file system it can see, so a host build either fails to find Octave or picks up a different host Octave (such as the older ``dnf`` package).  The host compiler (GCC 11 on RHEL 9) is also older than the Flatpak runtime's (GCC 13/14), which risks C++ ABI mismatches.

1. Clone CoolProp (from a regular RHEL terminal) into the home directory, for example ``~/CoolProp``, so the sandbox can see it.  A source folder elsewhere on the host first needs the ``flatpak override --filesystem=...`` command shown above.

2. Open a shell inside the Octave Flatpak::

    flatpak run --command=bash org.octave.Octave

3. Configure, build, and (optionally) install, as in the generic Linux instructions, adding the extra ``-DSWIG_EXECUTABLE`` flag so CMake finds the SWIG built above::

    # move to the CoolProp local repo folder
    cd ~/CoolProp
    # create a build directory and change to that directory
    mkdir -p buildOct && cd buildOct
    cmake .. -DCOOLPROP_OCTAVE_MODULE=ON -DSWIG_EXECUTABLE=$HOME/swig-build/install/bin/swig
    cmake --build . --parallel $(nproc)
    # (OPTIONAL) Install, and register it with ~/.octaverc for auto-loading
    cmake --install .

   Inside the shell, ``cmake``, ``g++``, ``python3``, ``mkoctfile``, and ``octave-config`` are all on the ``PATH`` and resolve the Flatpak's own Octave headers and libraries.  
   
   .. Tip::
      Instead of the ``-DSWIG_EXECUTABLE`` flag, ``export PATH="$HOME/swig-build/install/bin:$PATH"`` in the Flatpak shell before running ``cmake`` works as well.

Running the Example
-------------------

Flatpak sees the home directory, so the generated example is run the same way as in the generic instructions, from the same Flatpak shell.  CMake generates ``Example.m`` in the build folder, so it can be run from there without installing::

    cd ~/CoolProp/buildOct
    octave Example.m

If the optional ``cmake --install .`` step was run, a copy is also placed under the install prefix, and the path to it can be given from any folder, for example (with the default install prefix)::

    octave ~/CoolProp/install_root/Octave/Example.m

.. _octave_rhel_refprop:

REFPROP Under Flatpak
=====================

The Flatpak sandbox does not inherit the host shell's environment, so the environment variables that CoolProp needs to find REFPROP (see :ref:`the notes on the Octave wrapper page <octave_calls_notes>`) must be given to Flatpak itself.  (A shell opened with ``flatpak run --command=bash`` also reads ``~/.bashrc``, so exports placed there apply inside that shell only.)  The overrides below are persistent and apply to every way of starting this Octave.  For all users of the machine::

    sudo flatpak override \
        --env=COOLPROP_REFPROP_ROOT=/path/to/refprop \
        --env=RPPREFIX=/path/to/refprop \
        org.octave.Octave

or for the current user only::

    flatpak override --user \
        --env=COOLPROP_REFPROP_ROOT=/path/to/refprop \
        --env=RPPREFIX=/path/to/refprop \
        org.octave.Octave

The REFPROP folder has to be visible to the sandbox: anywhere under the home directory works, and any other location needs the ``--filesystem=`` override shown earlier.

.. raw:: html

    <div style="width: 100%; height: 2px; background-color: #72727254; margin-top: 2rem; margin-bottom: 0px; display: block; clear: both;"></div>
    <p style="margin-top: 4px; margin-bottom: 2rem; font-style: italic;">
        Page Last Updated and verified: 11-22-2026
    </p>
