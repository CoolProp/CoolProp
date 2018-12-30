
.. _LibreOffice:

*******************
LibreOffice Wrapper
*******************

General Information
-------------------

This is a LibreOffice extension, that makes the CoolProp High-Level Interface functions available in LibreOffice Calc. All functions and their parameters contains descriptions, that are available in the function wizard.

The LibreOffice extension itself is platform-independent. In the background the extension uses the Python CoolProp package for all calculations. Thus, the CoolProp Python package must be installed for the Python interpreter that is used by LibreOffice.

The extension contains a helper function to automate the installation of the Python package. This function installs the package inside the LibreOffice extension directory, e.g.::

    # on Windows
    C:\Users\user\AppData\Roaming\LibreOffice\4\user\uno_packages\cache\uno_packages\lu376486is.tmp_\CoolProp.oxt\pythonpath

    # on Linux
    /home/user/.config/libreoffice/4/user/uno_packages/cache/uno_packages/lu18931y72pq5.tmp_/CoolProp.oxt/pythonpath

    # on macOS
    /Users/user/Library/Application Support/LibreOffice/4/user/uno_packages/cache/uno_packages/lu104274oq0.tmp_/CoolProp.oxt/pythonpath

Installation
------------

1. Download the CoolProp.oxt Extension for LibreOffice (don't rename the file) and the example spreadsheet file ``TestLibreOffice.ods`` from :sfdownloads:`LibreOffice`

2. On Linux systems that split the LibreOffice package, install the necessary python script provider. On Ubuntu this can be done by::

    sudo apt-get install libreoffice-script-provider-python

3. Install the CoolProp Extension by double-clicking the oxt-file (install only for user)

4. Test the installation with the ``TestLibreOffice.ods`` file

5. If the CoolProp functions don't work, install the CoolProp Python package:

   a. In LibreOffice go to the menu Tools | Options | CoolProp | Installation
   b. Click the „Install CoolProp“ button
   c. The helper function checks if the CoolProp package is already installed. Otherwise it tries to download and install the Python package.

6. Open the file ``TestLibreOffice.ods``. All CoolProp formulas should be working.


Manual installation of the CoolProp Python package
--------------------------------------------------

Alternatively, you can also install the CoolProp Python package by yourself. If your LibreOffice installation uses the system Python (e.g. on Linux) then you can easily install the CoolProp Python package with pip::

    # system wide installation
    sudo pip3 install coolprop

    # install in user directory
    pip3 install coolprop --user


Another options is to download the CoolProp Python package at pypi.org. First check which Python version is used be LibreOffice. On Windows LibreOffice contains a bundled Python interpreter::
 
    # e.g. c:\programs\LibreOffice\program\python-core-3.5.5 --> python 3.5


Download the CoolProp wheel file for that python version from and unzip it to the ``pythonpath`` folder in the extension directory (please see above)::

    # choose win32 or amd64 depending on your Python interpreter
    # e.g. CoolProp-6.2.1-cp35-cp35m-win_amd64.whl for python 3.5 and 64bit LibreOffice
    https://pypi.org/project/CoolProp/#files
