
.. _Excel:

*************
Excel Wrapper
*************

.. contents:: :depth: 2

Pre-compiled Binaries for windows
=================================

Automated Installation
----------------------
As of October 2016, the Excel wrapper for Microsoft Windows can be installed from the Windows package as described on the :ref:`page on installation packages <Installers>`. Please refer to the documentation there for issues related to the installation process on Windows.

.. [#] If you get security warnings after you installed the CoolProp add-in, you can try adding the installation folder to your trusted directories: Got to ``File->Options->Trust Center->Trust Center Settings...->Trusted Locations`` and add the folder where your ``xlam``-file lives. This is most likely ``C:\Users\username\AppData\Roaming\Microsoft\AddIns\``.


Manual Installation
-------------------
For a manual installation,   

1.  Get the shared libraries at :sfdownloads:`shared_library/Windows/32bit__stdcall` and :sfdownloads:`shared_library/Windows/64bit`, put the DLLs into a folder of your choice and rename them to `CoolProp_stdcall.dll` and `CoolProp_x64.dll`. Make sure to add that folder to your path.  Technically you only need the DLL that matches your system architecture (`CoolProp_stdcall.dll` = 32-bit, `CoolProp_x64.dll` = 64-bit), but it can’t hurt to copy both if you don’t know which system architecture version you have.  The Excel macro will select the correct one and use it.
2.  Download the xla and xlam files from :sfdownloads:`MicrosoftExcel` and activate the add-in from Excel as described below. and copy the files to a convenient accessible location.
3.  The **TestExcel.xlsx** from :sfdownloads:`MicrosoftExcel` can be copied to a working directory in ``My Documents``.


1.  Open Excel
2.  Go to the menu ``File–>Options–>Add-Ins``
3.  At the bottom, select ``Manage: Excel Add-ins``, then click the ``Go...`` button
4.  Click the ``Browse`` button on the Add-in Manager panel
5.  Browse to the file **CoolProp.xlam** you downloaded and select it
6.  Make sure the CoolProp Add-in is selected (box checked) and close the Add-in Manager.
7.  Open the file **TestExcel.xlsx** and try to re-evaluate one of the cells; the CoolProp formulas should all be working now. (To recalculate the entire worksheet, press ``Ctrl``-``Alt``-``Shift``-``F9`` ) [#]_

.. [#] **Alternate DLL Location** - Some environments, lock down the folders included in the binary search path for normal users for security reasons.  If this is the case, you will need to put the DLL files in an alternate location (possibly on a shared network location for all users).  Follow the instructions below:

  1. Place the CoolProp DLL files in the alternate location
  2. Place the CoolProp xlam file in a writable location and open it.
  3. You will get an Excel error, ``File not found - CoolProp_stdcall.dll``.  Clicking **Ok** on the error will only clear the first of many.  Press and hold the **``Enter``** key until all of the errors clear.
  4. Make sure that the Developer ribbon is visible on the main menu.  If not
  
     - Go to **File | Options** on the main menu and select Cusomize Ribbon
     - Make sure that the Developer main tab is checked (ON)
     
  5. Open the Visual Basic editor and use **Edit | Replace** to replace all occurrences of **CoolProp_stdcall.dll** with the full path to the alternate location for your DLL files, making sure to press the save button (disk image) or **File | Save** before exiting the VBA editor.
  6. Save the CoolProp.xlam file.


.. [#] If you are having problems with the path to the XLAM in your worksheet appearing as the complete path to the xlam (but not the correct path), you might need to go into ``Data->Update Links`` in Excel, select CoolProp.xlam, and select ``Change Source...`` and select the location of your xlam file.  That should update all the links.

 One possible cause of this problem is the security feature of Windows. You can avoid it by manually giving the permission for ordinary access to the xlam file as follows:

 1. Right-click the xlam file and open Properties.
 2. There might be a security message "This file came from another computer and ...".
 3. Check the **Unblock** checkbox (or button other than Windows 10) and click **Apply** at the bottom.
 4. After confirming the security message has disappeared, click **OK** to exit from Properties.
    
Pre-compiled Binaries for OSX
=============================

.. warning:: 

  There are now both 32-bit and 64-bit versions of Microsoft Excel on OSX.  You need to make sure that your bitness of the shared library for CoolProp (and perhaps REFPROP) match that of Excel.  

Part 1:
-------

There are several ways to determine the bitness of your Excel version.  The easiest is to open a terminal, and do something like::

    Ians-Mac-mini:~ ian$ file /Applications/Microsoft\ Excel.app/Contents/MacOS/Microsoft\ Excel 
    /Applications/Microsoft Excel.app/Contents/MacOS/Microsoft Excel: Mach-O 64-bit executable x86_64

Or you can go into Excel->About Excel.  If version is greater than 15.24, you are running a 64-bit version of Excel.

Part 2:
-------
We need to convince Microsoft Excel to load our shared library, and it seems the only place it is willing to look for shared libraries is in the folder ``/Users/${USER}/Library/Group Containers/UBF8T346G9.Office``, where ``${USER}`` should be replaced with your user name.  This is because Excel is now sandboxed.

Following http://apple.stackexchange.com/a/106814, save these contents as the file ``~/Library/LaunchAgents/my.startup.plist`` (obviously replace ``ihb`` with the appropriate user name)::

    <?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
    <plist version="1.0">
    <dict>
    <key>Label</key>
    <string>my.startup</string>
    <key>ProgramArguments</key>
    <array>
      <string>sh</string>
      <string>-c</string>
      <string>launchctl setenv DYLD_LIBRARY_PATH "/Users/ihb/Library/Group Containers/UBF8T346G9.Office"</string>
    </array>
    <key>RunAtLoad</key>
    <true/>
    </dict>
    </plist>

This ``.plist`` will be run as soon as the computer starts, and will set the ``DYLD_LIBRARY_PATH`` environmental variable, and Microsoft Excel will then read this variable, and be willing to load your shared library

Make sure to log out and log back in to have this ``.plist`` take effect.

Part 2a (optional):
-------------------
If you want to add additional paths to the terminal, you can add a line like this to your ``~/.bash_profile`` for instance to append paths to the ``DYLD_LIBRARY_PATH`` variable. It calls ``launchctl`` to extract the ``DYLD_LIBRARY_PATH`` environment variable and prepends ``/another/path`` to it::

    export DYLD_LIBRARY_PATH="/another/path:`launchctl getenv DYLD_LIBRARY_PATH`"

Part 3:
-------

Download the xlam from :sfdownloads:`MicrosoftExcel` or the development version from :sfnightly:`MicrosoftExcel`.

Place XLAM file in ``/Users/${USER}/Library/Group Containers/UBF8T346G9.Office``, where ``${USER}`` should be replaced with your user name

Follow the below instructions depending on the version of Excel you have.  If you can't figure out what version of Excel you have, it's fine to have both 32-bit and 64-bit versions of the .dylib sitting next to each other.

32-bit
^^^^^^

Download pre-compiled release binaries for OSX from :sfdownloads:`shared_library/Darwin/32bit/`.  Development binaries coming from the buildbot server can be found at :sfnightly:`shared_library/Darwin/32bit/`. 

Place the downloaded file ``libCoolProp.dylib`` in the folder ``/Users/${USER}/Library/Group Containers/UBF8T346G9.Office`` too, but RENAME it to ``libCoolProp_32bit.dylib`` (this is to ensure that there is no name clash with the standard 64-bit shared library).

64-bit
^^^^^^

Download pre-compiled release binaries for OSX from :sfdownloads:`shared_library/Darwin/64bit/`.  Development binaries coming from the buildbot server can be found at :sfnightly:`shared_library/Darwin/64bit/`. 

Place the downloaded file ``libCoolProp.dylib`` in the folder ``/Users/${USER}/Library/Group Containers/UBF8T346G9.Office``.

Part 4:
-------

Open Excel, go to ``Tools/Add-ins...``. In browse, go to the folder listed above with the ``BF8T346G9.Office`` in it. Select CoolProp.xlam.

Part 4b:
-------
Go to Tools/Macro/Visual_Basic_Editor and open Module 1 in CoolProp.xlam.  Replace all references to “libCoolProp.dylib” with references to "/Users/${USER}/Library/Group Containers/UBF8T346G9.Office/libCoolProp.dylib”, again changing ${USER} to your user name.  Save and close the Visual Basic Editor.

Part 5:
-------
Add this to a cell::

    =PropsSI("T","P",101325,"Q",0,"Water")

make sure you get something like 373.1242958 K.

Debugging
---------

* If it doesn't work and you get error number 53, it might be because you have a 64-bit .dylib file and you want a 32-bit .dylib file.  For instance when you run the ``file`` command on your .dylib, you should see something like::

    $ file libCoolProp_32bit.dylib
    libCoolProp.dylib: Mach-O dynamically linked shared library i386

  the ``i386`` is the important bit, that indicates that the shared library is 32-bit.

REFPROP support on OSX
======================

You can also call REFPROP through the Excel wrapper of CoolProp, but it requires a few tweaks to work properly

1. The refprop dylib (with the correct bitness!), as well as the ``fluids`` and ``mixtures`` folders of REFPROP should be placed in the folder ``refprop`` inside ``/Users/${USER}/Library/Group Containers/UBF8T346G9.Office``.  Make sure the shared library is called ``librefprop.dylib``.
2. An environment variable called ``COOLPROP_REFPROP_PATH`` should be set to the folder ``/Users/${USER}/Library/Group Containers/UBF8T346G9.Office/refprop`` (see above about how to do that in a ``.plist`` file).  The CoolProp xlam, on loading, will query this environment variable to determine which path to use for REFPROP.  It seems from my testing that this path MUST be a subfolder of ``/Users/${USER}/Library/Group Containers/UBF8T346G9.Office`` due to the sandboxing.

User-compiled Binaries
======================

Common Requirements
-------------------
Compilation of the Excel wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Build
-----

The instructions here are for a 64-bit windows system that will compile both 64-bit and 32-bit versions of the DLL::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder for the 32-bit DLL
    mkdir build/32bit__stdcall && cd build/32bit__stdcall
    # Build the MSVC project using CMake
    cmake ../.. -G "Visual Studio 10" -DCOOLPROP_SHARED_LIBRARY=ON -DCOOLPROP_STDCALL_LIBRARY=ON
    # Make the shared library
    cmake --build . --config Release
    cd ../..
    # Make a build folder for the 64-bit DLL
    mkdir build/64bit && cd build/64bit
    # Build the MSVC project using CMake
    cmake ../.. -G "Visual Studio 10 Win64" -DCOOLPROP_SHARED_LIBRARY=ON 
    # Make the shared library
    cmake --build . --config Release
    cd ../..
    # Copy the generated DLL
    copy build\32bit__stdcall\CoolProp.dll c:\CoolProp
    copy build\64bit\CoolProp.dll c:\CoolProp
