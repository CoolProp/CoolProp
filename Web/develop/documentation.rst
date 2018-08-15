.. _developer_documentation:

*************
Documentation
*************

General
-------

The CoolProp documentation at http://www.Coolprop.org ( and http://www.CoolProp.org/dev ) is built dynamically (nightly, in the case of the dev pages) by running a collection of `reStructuredText <http://docutils.sourceforge.net/rst.html>`_ (.rst) files through Sphinx, which converts them to html.  This is done so that pages can be easily edited/created in .rst format and iPython examples of running CoolProp can be embedded and executed dynamically using the latest version of CoolProp, displaying the example output in the resulting html documentation.  The instructions below serve to document how the automated builds are done as well as provide instruction for generating local documentation, useful for creating, editing, viewing, and debugging new contributor content.

Build Sphinx documentation
--------------------------

0. On a Mac system, you need a working Latex distribution and some more command line tools that can be installed 
   with ``brew install jpeg libpng libtiff libtool imagemagick``. Additionally, you should add the local string 
   to your bash environment in ``/Users/username/.bash_profile``::

    export LC_ALL=en_US.UTF-8
    export LANG=en_US.UTF-8
    

1. Check out the sources in the CoolProp/Web folder::

    git clone https://github.com/CoolProp/CoolProp --recursive
    cd CoolProp/Web

2. Make a virtualenv::

    virtualenv --distribute ~/env/py27
    source ~/env/py27/bin/activate # Turn on this virtual env, should see (py27) in your command shell next to the prompt to tell you that environment is active

3. Then install prerequisites into this virtualenv::

    pip install --upgrade setuptools cython
    pip install numpy scipy matplotlib python-dateutil pytz pandas
    pip install sphinx ipython sphinxcontrib-bibtex sphinxcontrib-doxylink sphinxcontrib-napoleon cloud-sptheme
    pip install wxpython


4. Run setup script, `CoolProp\Web\scripts\__init__.py` to generate dynamic content::

    cd scripts
    python __init__.py
    cd ..

5. To build the documentation, go into the CoolProp/Web folder and run::

    make html

6. Move the generated docs in ``_build`` to wherever you want


Build Doxygen documentation
---------------------------

All the configuration is done in the ``Doxyfile`` file.

1. If you don't have doxygen, install it from `doxygen downloads <http://www.stack.nl/~dimitri/doxygen/download.html>`.  Or on Linux a ``sudo apt-get install doxygen`` should do it.

2. Simply fire up a shell in the root of the repo, and type::

    doxygen Doxyfile


Building the documentation on Linux
-----------------------------------

1. Make sure you have what you need::

    sudo aptitude install build-essential imagemagick git gfortran cmake doxygen ipython python-pip python-virtualenv
    sudo aptitude install libatlas-base-dev libatlas3-base # numpy
    sudo aptitude install gcc gfortran python-dev libblas-dev liblapack-dev libblas3 liblapack3 # scipy
    sudo aptitude install python-dev libpng-dev tk libfreetype6-dev # matplotlib
    sudo aptitude install libxml2-dev libxslt1-dev libxslt1.1 python-all-dev # pandas

2. Follow the instructions above to create the virtual environment.


Building the documentation on Windows
-------------------------------------

Use these instructions to build the CoolProp documentation (html) on 64-bit windows systems.  This should only be needed if you need to reproduce the site at http://www.CoolProp.org or http://www.CoolProp.org/dev on a local machine.  This is extremely useful (if not necessary) when contributing changes to the documentation to see how changes to the .rst files will look in the html documentation.  There are probably many ways to do this that will work, but the method below makes use of the **Conda** environments to create clean virtual environments with the requisite packages needed to build most of the documentation.  This is not meant to be a comprehensive build, as some of the scripts don't run well under windows Anaconda (e.g. pdfnup, rsync, etc.).  It will however assemble the bulk of the documentation for those who need to contribute or edit doc pages.

**Prerequisites**

1. Install GNUwin32 Make at http://gnuwin32.sourceforge.net/packages/make.htm and add ``C:\Program Files (x86)\GnuWin32\bin`` to the PATH environment variable.  This will ensure that `make` is available from any command window.  


2. Install Anaconda using the 64-bit Windows Installer with the following installation options:

 - Use the ``Per User`` installation (vice installing for all users).  This will allow the installation of python packages without having to open a window with administrative privileges.  Windows (8, 8.1, and 10) does a "lock-down" of the ``Program Files`` directories such that they can only be modified with administrative privileges.

 - Choose Python version 2.7.x or 3.6.x

 - Add Anaconda to the path.  This will ensure that the Anaconda/python commands (like Conda) are available from any command window.  

 
3. If you want to generate the wrapper examples in the documentation, you will need to install the software packages that they use (e.g. Octave, REFPROP, etc.).  If they aren't installed, these examples will be skipped, but will generate some warnings and errors in the steps below.  These can be ignored.  
	
.. note::

    The 64-bit Windows version of Anaconda will automatically install 64-bit versions of packages like **numpy** and **scipy** by default; something that is not guaranteed in other python environments for Windows.  Also, Anaconda is packaged with MKL-powered binary versions of some of the most popular numerical/scientific Python libraries for improved performance; including **numpy** and **scipy**.  This is a requisite for CoolProp Python scripts, and another reason to use Anaconda (or Miniconda) with CoolProp.  

.. note:: 

    For help in using Conda commands, this `Navigator Cheat Sheet <https://docs.continuum.io/_downloads/Anaconda_CheatSheet.pdf>`_ and `Conda Cheat Sheet <http://conda.pydata.org/docs/_downloads/conda-cheatsheet.pdf>`_ can be very useful.
   
**Setup a Virtual Python Environment**

Set up a virtual python environment and name it something like CP27 (that's what is used in the examples below) or CP36 if you are using Python 3.6.  This virtual environment will contain all the modules needed to build the CoolProp documentation.  Setting up a virtual environment is a very simple thing to do in the Anaconda Navigator (Graphical Interface), but you can also set up the environment using the following commands in a command line window::

    conda create --name CP27 python=2.7
    activate CP27 # Will cause command prompt to be prefixed with (CP27)
	
.. note::

   Any instructions below that take place on the Windows are assumed to be in the virtual environment created here.  For example, when opening a Windows command prompt (cmd), *activate CP27* (or CP36, or whatever you named your virtual environment) first before issing any other commands.

**Install Required Python Modules**

Most of the following modules can be installed from the Anaconda Navigator, and it is much simpler to do so.  Even if the module is already in the Anaconda environment, it is best to update the packages listed, making sure to get the latest Windows 64-bit versions.  However, some of the requisite modules are not in the Conda library and have to be installed from the command line using `pip install`.  These will be designated as such.  The versions of each module that have been tested for generating the docs on Windows are shown in the table below. 

+-------------------------+-------------+-------------+
| Modules (Conda)         | Python 2.7  | Python 3.6  |
+=========================+=============+=============+
| setuptools              | 27.2.0      | 36.4.0      |
+-------------------------+-------------+-------------+
| cython                  | 0.25.2      | 0.26        |
+-------------------------+-------------+-------------+
| numpy                   | 1.12.1      | 1.13.1      |
+-------------------------+-------------+-------------+
| scipy                   | 0.19.0      | 0.19.1      |
+-------------------------+-------------+-------------+
| matplotlib              | 2.0.2       | 2.0.2       |
+-------------------------+-------------+-------------+
| python-dateutil         | 2.6.0       | 2.6.1       |
+-------------------------+-------------+-------------+
| pytz                    | 2017.2      | 2017.2      |
+-------------------------+-------------+-------------+
| pandas                  | 0.20.1      | 0.20.3      |
+-------------------------+-------------+-------------+
| sphinx                  | 1.5.6       | 1.6.3       |
+-------------------------+-------------+-------------+
| ipython                 | 5.3.0       | 6.1.0       |
+-------------------------+-------------+-------------+
| wxpython                | 3.0         |     n/a     |
+-------------------------+-------------+-------------+

The following modules will need to be installed using *pip install*.

+-------------------------+-------------+-------------+
| Modules (pip)           | Python 2.7  | Python 3.6  |
+=========================+=============+=============+
| sphinxcontrib-bibtex    | 0.3.5       | 0.3.5       |
+-------------------------+-------------+-------------+
| sphinxcontrib-doxylink  | 1.3         | 1.3         |
+-------------------------+-------------+-------------+
| sphinxcontrib-napoleon  | 0.6.1       | 0.6.1       |
+-------------------------+-------------+-------------+
| cloud-sptheme           | 1.9.4       | 1.9.4       |
+-------------------------+-------------+-------------+
| wxPython                |     n/a     | 4.0.0b1     |
+-------------------------+-------------+-------------+


.. note::

   Errors may occur during the install of some of the above Python modules or while making the html docs.  Make sure that all of the above Python modules are the most recent versions using pip install --upgrade.  There are known issues with the default sphinx version installed with Anaconda on Windows.
   
If updating the dev version of the documentation (i.e. the latest dev source files are loaded), the dev version of the CoolProp module will need to be installed.  This ensures that any example code will run properly with the latest CoolProp functions.::

    pip install -vvv --pre
     --trusted-host www.coolprop.dreamhosters.com
     --find-links http://www.coolprop.dreamhosters.com/binaries/Python/ -U 
     --force-reinstall CoolProp


**Run Setup Scripts**

There are a number of setup scripts that have to be run to generate dynamic content for the web documentation.  To run them, open a windows command prompt (cmd) and cd to the CoolProp\\Web\\scripts directory.  Then run the following scripts::

    python ..\..\dev\scripts\examples\win64run.py
    python coolprop.tabular.speed.py
    python fluid_properties.phase_envelope.py
    python fluid_properties.PurePseudoPure.py
    python fluid_properties.Mixtures.py
    python coolprop.parametric_table.py
    python coolprop.configuration.py
    # The next three scripts can take a while to run
    python fluid_properties.Consistency.py  # Many errors/warnings generated; this is normal
    python logo_2014.py
    python fluid_properties.REFPROPcomparison.py  # Only if you have REFPROP

.. note::

   These scripts are normally run by the BuildBot using the Python 2.x initialization script, `CoolProp\\Web\\scripts\\__init__.py <https://github.com/CoolProp/CoolProp/blob/master/Web/scripts/__init__.py#L107>`_ (lines 107-122).  This script could be run to help with automation, however, there are linux and OSX shell scripts included that will not run on Windows.  Also, these scripts only need to be run once, and many may generate errors and warning messages that will be useful in debugging your python environment.  Once the dynamic content from these scripts has been generated, you're ready to build the documentation.
   
**Build the Documentation**

There is a Makefile that will build the entire site.  This can take a while, especially the first time.  Open a Windows command prompt, activate your virtual environment, cd to the *CoolProp\\Web* directory and type::

    make html

This will create a new *CoolProp\\Web\\_build* directory that will contain the html pages.  If the build already exists, make will skip the parts of the documentation that have not changed.  You can look at built documentation pages by opening *CoolProp\\Web\\_build\\html\\index.html* in any web browser, or by double clicking on this file.

To remove the docs, use::

    make clean
	
Or to build a completely fresh version of the documentation (this will take longer) use::

    make clean html
	
.. note::

   Warnings and errors will be generated, depending on what packages you have installed.  If pages or sections are missing from your build, check that you have the prerequisites installed, have all the up-to-date python modules, and have run the setup scripts successfully.