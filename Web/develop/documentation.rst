.. _developer_documentation:

*************
Documentation
*************

General
-------

The CoolProp documentation at https://coolprop.org is built by running a collection of `reStructuredText <https://docutils.sourceforge.net/rst.html>`_ (.rst) files through Sphinx, which converts them to html.  This is done so that pages can be easily edited/created in .rst format and iPython examples of running CoolProp can be embedded and executed dynamically using the latest version of CoolProp, displaying the example output in the resulting html documentation.  The build runs on GitHub Actions (see ``.github/workflows/docs_docker-run.yml``) in the ``ghcr.io/coolprop/coolprop_docs_02_builder`` container; the documentation is published to GitHub Pages on tagged releases.  The instructions below serve to document how the automated builds are done as well as provide instruction for generating local documentation, useful for creating, editing, viewing, and debugging new contributor content.

Build Sphinx documentation
--------------------------

0. On a Mac system, you need a working Latex distribution and some more command line tools that can be installed 
   with ``brew install jpeg libpng libtiff libtool imagemagick``. Additionally, you should add the local string 
   to your bash environment in ``/Users/username/.bash_profile``::

    export LC_ALL=en_US.UTF-8
    export LANG=en_US.UTF-8
    

1. Check out the sources in the CoolProp/Web folder::

    git clone https://github.com/CoolProp/CoolProp
    cd CoolProp/Web

2. Create the build environment.  The authoritative list of dependencies is
   ``dev/docker/conda_environment.yml`` (the same file the CI container is built
   from).  The simplest way to reproduce it is with conda::

    conda env create -f dev/docker/conda_environment.yml
    conda activate docs   # the environment is named "docs"

   You will also need a working install of the CoolProp Python wheel built from
   the same checkout (``pip wheel .`` then ``pip install`` the resulting wheel),
   so that the embedded examples run against the matching version of CoolProp.

3. Run setup script, ``CoolProp/Web/scripts/__init__.py`` to generate dynamic content::

    cd scripts
    python __init__.py
    cd ..

4. To build the documentation, go into the CoolProp/Web folder and run::

    make html

5. Move the generated docs in ``_build`` to wherever you want


Build Doxygen documentation
---------------------------

All the configuration is done in the ``Doxyfile`` file.

1. If you don't have doxygen, install it from `doxygen downloads <https://www.doxygen.nl/download.html>`_.  Or on Linux a ``sudo apt-get install doxygen`` should do it.

2. Simply fire up a shell in the root of the repo, and type::

    doxygen Doxyfile


Building the documentation on Linux
-----------------------------------

1. Make sure you have what you need::

    sudo apt-get install build-essential imagemagick git gfortran cmake doxygen pandoc python3-dev python3-pip
    sudo apt-get install libblas-dev liblapack-dev # numpy / scipy
    sudo apt-get install libpng-dev tk libfreetype6-dev # matplotlib

2. Follow the instructions above to create the conda environment.


Building the documentation on Windows
-------------------------------------

Building on Windows works the same way as on Linux/macOS: create the conda
environment from ``dev/docker/conda_environment.yml`` (see above), install a
matching CoolProp Python wheel, run the setup scripts, and then ``make html``.
The ``conda`` package manager (Miniconda or Anaconda) is the simplest way to get
a clean, MKL-backed scientific Python stack on Windows.

**Prerequisites**

1. Install GNUwin32 Make at https://gnuwin32.sourceforge.net/packages/make.htm
   and add ``C:\Program Files (x86)\GnuWin32\bin`` to the ``PATH`` environment
   variable, so that ``make`` is available from any command window.

2. Install `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ (or
   Anaconda) for 64-bit Windows, then create the build environment::

    conda env create -f dev/docker/conda_environment.yml
    conda activate docs

3. If you want to generate the wrapper examples in the documentation, you will
   need the software packages that they use (e.g. Octave, REFPROP, etc.).  If
   they aren't installed, those examples are skipped with warnings that can be
   ignored.

**Run Setup Scripts**

The setup scripts generate dynamic content for the documentation.  Open a
command prompt, ``conda activate docs``, ``cd`` to ``CoolProp\Web\scripts`` and
run::

    python __init__.py

(On Linux/macOS the same entry point is used; see ``Web/scripts/__init__.py``.)
These only need to be run once, and some may emit warnings depending on which
optional packages you have installed.

**Build the Documentation**

From the ``CoolProp\Web`` directory, with the ``docs`` environment active::

    make html

This creates a ``CoolProp\Web\_build`` directory containing the html pages; open
``_build\html\index.html`` in any browser.  ``make`` skips parts of the build
that have not changed.

To remove the build, use ``make clean``; to force a completely fresh build use
``make clean html``.
