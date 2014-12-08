
.. _LibreOffice:

*******************
LibreOffice Wrapper
*******************

Windows using DLL (shared library)
----------------------------------

A very simple example can be found here: :download:`TestLibreOffice.ods`.  You will need to download the release 32-bit __stdcall shared library for windows from :sfdownloads:`sourceforge <shared_library/Windows/32bit__stdcall_calling_convention>` or from :sfnightly:`the development snapshots <shared_library/Windows/32bit__stdcall_calling_convention>`.  Place the downloaded DLL in c:\\CoolProp

Linux using python
------------------

1. You will need two files - :download:`coolprop.py` and :download:`TestLibreOfficePy.ods`

2. At the command prompt do something like::

    # Make a folder for the script to reside in
    mkdir -p ~/.config/libreoffice/4/user/Scripts/python/
    # Copy the dowloaded .py file to the folder you made (adjust download path as necessary)
    cp ~/Downloads/coolprop.py ~/.config/libreoffice/4/user/Scripts/python/

    # Install python script provider (uses python 3)
    sudo apt-get install libreoffice-script-provider-python
    # pyton scripting in LibreOffice uses python3, make sure we populate the right python site-packages folder
    sudo apt-get install python3-pip
    sudo pip3 install Cython
    # Compile CoolProp v5+ manually
    # you need to have numpy and matplotlib already - you can do sudo apt-get install python3-matplotlib
    git clone http://github.com/CoolProp/CoolProp --recursive
    cd CoolProp/wrappers/Python
    sudo python3 setup.py install

3. Open example in LibreOffice