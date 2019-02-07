.. _MATLAB:

**************
MATLAB Wrapper
**************

.. contents:: :depth: 2

The only interface for MATLAB that is supported as of version 6.2 is via Python (you can also call the DLL directly, see below). The maintenance hurdles of supporting the old SWIG MATLAB wrapper proved to be too difficult to surmount.

Via Python
==========

You will need to acquire a Python interpreter, the easiest method to do so if you do not already have Python installed on your computer is to download the installer from python.org: `link to downloads <https://www.python.org/downloads/>`_. There is also information on the Mathworks website. The Python wrapper only uses methods that are in the Python standard library, so the standard installation of Python would work fine. If you would like to have a more full-featured Python installation, you can install a full-fledged Python installation from Anaconda: `Anaconda download <https://www.anaconda.com/download/>`_.

In your MATLAB shell, you can inquire about what Python version MATLAB intends to use with a command like::

    >> pyversion

           version: '3.6'
        executable: 'D:\Anaconda\python.exe'
           library: 'D:\Anaconda\python36.dll'
              home: 'D:\Anaconda'
          isloaded: 0

Good. It found Python, and has not loaded it yet. You are ready!

If you have multiple copies of Python on your computer already, then you can tell MATLAB which one you want it to use by passing the absolute path to the python executable to pyversion. For instance::


    >> pyversion d:\Anaconda\envs\py36\python.exe
    >> pyversion

           version: '3.6'
        executable: 'd:\Anaconda\envs\py36\python.exe'
           library: 'd:\Anaconda\envs\py36\python36.dll'
              home: 'd:\Anaconda\envs\py36'
          isloaded: 0

Finally, you need to install CoolProp into your given copy of python. This one-liner calls the ``pip`` program of Python to install the CoolProp package from the PYPI package index. Watch out for the spaces in the arguments, they are important!::

    >> [v,e] = pyversion; system([e,' -m pip install --user -U CoolProp'])

Use
---

Then you can calculate the normal boiling point temperature of water::

    >> py.CoolProp.CoolProp.PropsSI('T','P',101325,'Q',0,'Water')

    ans = 

      373.1243

Similar approaches are possible for the low-level interface. Addition of docs documenting how to use the low-level interface in MATLAB would be welcome.  Also, more advanced wrappers for MATLAB are available, and are stored on github: https://github.com/CoolProp/CoolProp/tree/master/wrappers/MATLAB

.. _low_level_high_level_matlab: 

Calling Low-Level interface through DLL
=======================================

.. literalinclude:: LowLevelHighLevelMATLAB.m
   :language: matlab
   