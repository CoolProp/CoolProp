
.. _Installers:

******************
Installer Packages
******************

.. warning:: We only have an installer for Windows and it is still under development. Please file an issue if you experience any problems with it.

Packaged binaries can be downloaded from :sfdownloads:`Installers`.  Development binaries coming from the continuous integration server can be found at :sfnightly:`Installers`.

*****************
Windows Installer
*****************

Shared Libraries
================

The installer for Windows includes the shared libraries (dll) for 32bit (stdcall and cdecl calling conventions) and 64bit.  The libraries are installed into the user's ``AppData/CoolProp`` folder and the is added to the user's search path for libraries and executables via the environment variable ``%PATH%``. 

Microsoft Excel
===============

The add-in for Microsoft Excel can also be installed from this package and it automatically determines which of the two available add-ins is needed.  There is ``CoolProp.xla`` for Microsoft Office versions prior to version 2007 and ``CoolProp.xlam`` for more recent Office installations. The add-ins are installed into the user's add-in directory and they are activated by default.  Please note that the Excel add-in includes the installation of some of the shared libraries.  Upon installation, an example file is placed on the user's desktop.  This file demonstrates some of the features of the Excel wrapper and can be moved or deleted freely, it will not be removed by the uninstaller. Have a look at the dedicated :ref:`page <Excel>` for more information.

Engineering Equation Solver
===========================

A custom library for the Engineering Equation Solver (EES) is also part of this package.  The installer puts the files into ``C:\EES32\Userlib\COOLPROP_EES`` and you have to move them manually if you installed EES to another location.  There is also an example file for EES that will be placed on the user's desktop.  This file demonstrates some of the features of the Excel wrapper and can be moved or deleted freely, it will not be removed by the uninstaller.  Have a look at the dedicated :ref:`page <EES>` for more information.
