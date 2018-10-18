
.. _LibreOffice:

*******************
LibreOffice Wrapper
*******************

General Information
-------------------

The wrapper for LibreOffice uses the CoolProp python package on Linux and Windows systems. The template spreadsheet file :download:`LibreOffice_coolprop.ods` includes a python script and Basic macros to create the cell functions ``=PropsSI()`` and ``=PhaseSI()`` that wraps the appropriate functions from the underlying CoolProp python package.

In LibreOffice the macro security level must be changed to medium to allow the unsigned macros to run. At the medium security level a confirmation is required before executing macros from an opened file.

The Basic macros can be directly changed with the LibreOffice Basic editor. To change the included python script (e.g. to add more functions), unzip the ods-file, edit the file ``Scripts/python/coolprop.py``, compress the changed contents to a zip archive and change the file extension to ods.


Linux
-----

1. On systems that split the LibreOffice package, install the necessary python script provider. On Ubuntu this can be done by::

    sudo apt-get install libreoffice-script-provider-python

2. If the python package manager pip is not already installed on the system, then install it with the package manager, e.g.::

    # on Ubuntu LibreOffice uses python 3
    sudo apt-get install python3-pip

3. Install CoolProp python package with pip:: 

    # system wide installation
    sudo pip3 install coolprop

    # install in user directory
    pip3 install coolprop --user

4. Open the example file :download:`LibreOffice_coolprop.ods`


Windows
-------

1. Install LibreOffice from official installer (that already includes a bundled python)

2. Check the python version bundled with Libreoffice from the install directory::
 
    # e.g. c:\programs\LibreOffice\program\python-core-3.5.5 --> python 3.5

3. Download coolprop wheel file for that python version from::

    # choose win32 or amd64 depending on the used LibreOffice installer
    # e.g. CoolProp-6.1.0-cp35-cp35m-win_amd64.whl for python 3.5 and 64bit LibreOffice
    https://pypi.org/project/CoolProp/#files

4. Unzip that wheel file into the LibreOffice python site-packages directory::

    # e.g. c:\programs\LibreOffice\program\python-core-3.5.5\lib\site-packages

5. Open the example file :download:`LibreOffice_coolprop.ods`

