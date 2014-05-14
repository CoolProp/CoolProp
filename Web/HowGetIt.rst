Downloading CoolProp
====================

How to get it?
--------------

Option 1 (easiest)
^^^^^^^^^^^^^^^^^^

Head to https://sourceforge.net/projects/coolprop/files/CoolProp and download the most recent version.  Each language has instructions on what you should do.  All the files are already compiled and should work out of the box.

Option 1a (for Python users)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Nightly build installers are also available at https://sourceforge.net/projects/coolprop/files/CoolProp/Nightly for a limited subset of python configurations and are updated every night to be current with the main developer's personal codebase.

.. warning::

    Nightly build may break your code, give wacky results, or otherwise. Use at your own risk.

Option 2 (for Python users)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

CoolProp is now on `PyPI <http://pypi.python.org/pypi/CoolProp>`_.  If you already have `cython <http://www.cython.org>`_ (version > 0.17) installed and your default compiler is already configured (see below), you can just do::

    easy_install CoolProp
    
Or if you already have CoolProp installed, you can upgrade it with::

    easy_install -U CoolProp
    
Or using pip::

    pip install CoolProp
    
Option 3 (developers and the courageous)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CoolProp is an open-source project, and is actively looking for developers.  The project is hosted in a git repository on github at::

    https://github.com/ibell/coolprop
    
and you can check out the sources by doing::

    svn checkout https://github.com/ibell/coolprop coolprop-code

or if you want to just browse the repository, you can go to https://github.com/ibell/coolprop.

Compiler Configuration
----------------------
If you are on OSX or linux/unix, you probably don't have to do anything at all since python will just use the most recent version of gcc.

If you are a windows user and you have installed Visual Studio 2008 (even the `free express version <http://www.microsoft.com/visualstudio/en-us/products/2008-editions/express>`_ works) python will default to this compiler and everything should go just fine.  Make sure you do not install the 2010 version since python 2.x versions are compiled with Visual Studio 2008 compiler.  Yes I know that is annoying.

The only thing that is a bit tricky if if you have not installed Visual Studio 2008 and instead want to use the MINGW compiler (a windows version of the gcc compiler).  It can be installed using the `python(x,y) <http://www.pythonxy.com>`_ distribution for instance.  In that case, if you are running a command line build of CoolProp, you need to do something like::

    python setup.py build --compiler=mingw32 install
    
Or if you want to use ``easy_install``, you need to create a distutils configuration file called ``distutils.cfg`` located in ``c:\\Python27\\lib\\distutils\\distutils.cfg`` with the contents::

    [build]
    compiler = mingw32

Uninstall
---------
If you don't want CoolProp anymore, just delete the CoolProp folder in the Lib/site-packages folder for your distribution, as well as the CoolProp .egg file in Lib/site-packages
