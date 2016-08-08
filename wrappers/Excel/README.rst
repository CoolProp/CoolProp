Excel Wrapper
=============

Pre-compiled Binaries for windows
---------------------------------

Pre-compiled release binaries can be downloaded from `MicrosoftExcel <http://sourceforge.net/projects/coolprop/files/CoolProp/6.0.0/MicrosoftExcel>`_ .  Development binaries coming from the buildbot server can be found at `MicrosoftExcel <http://sourceforge.net/projects/coolprop/files/CoolProp/nightly/MicrosoftExcel>`_.

Download all the files. The basic protocol is:

**Part 1:**

1.  Copy the files **CoolProp.dll** and **CoolProp_x64.dll** to **C:\\CoolProp** [#]_ folder. Technically you only need the DLL that matches your system architecture (`CoolProp.dll` = 32-bit, `CoolProp_x64.dll` = 64-bit), but it can’t hurt to copy both if you don’t know which system architecture version you have.  The excel macro will select the correct one and use it.
2.  Copy the file **CoolProp.xlam** to a convenient accessible location (it could go in the C:\\CoolProp folder as well).
3.  **CoolProp.xlsx** can be copied to a working directory in ``My Documents``.

**Part 2:**

1.  Open Excel
2.  Go to the menu ``File–>Options–>Add-Ins``
3.  At the bottom, select ``Manage: Excel Add-ins``, then click the ``Go...`` button
4.  Click the ``Browse`` button on the Add-in Manager panel
5.  Browse to the file **CoolProp.xlam** you downloaded and select it
6.  Make sure the CoolProp Add-in is selected (box checked) and close the Add-in Manager.
7.  Open the file **TestExcel.xlsx** and try to re-evaluate one of the cells; the CoolProp formulas should all be working now. (To recalculate the entire worksheet, press ``Ctrl``-``Alt``-``Shift``-``F9`` ) [#]_

.. [#] **Alertnate DLL Location** - Some environments, lock down the C:\\ drive from normal users for security reasons.  If this is the case, you will need to put the DLL files in an alternate location (possibly on a shared network location for all users).  Follow the instructions below:

  1. Place the CoolProp DLL files in the alternate location
  2. Place the CoolProp xlam file in a writable location and open it.
  3. You will get an Excel error, ``File not found - C:\CoolProp\CoolProp.dll``.  Clicking **Ok** on the error will only clear the first of many.  Press and hold the **``Enter``** key until all of the errors clear.
  4. Make sure that the Developer ribbon is visible on the main menu.  If not
  
     - Go to **File | Options** on the main menu and select Cusomize Ribbon
     - Make sure that the Developer main tab is checked (ON)
     
  5. Open the Visual Basic editor and use **Edit | Replace** to replace all occurances of **C:\\CoolProp** with the path to the alternate location for your CoolProp.dll files, making sure to press the save button (disk image) or **File | Save** before exiting the VBA editor.
  6. Save the CoolProp.xlam file.


.. [#] If you are having problems with the path to the XLAM in your worksheet appearing as the complete path to the xlam (but not the correct path), you might need to go into ``Data->Update Links`` in Excel, select CoolProp.xlam, and select ``Change Source...`` and select the location of your xlam file.  That should update all the links.


Pre-compiled Binaries for OSX
-----------------------------

**Part 1:**

Download pre-compiled release binaries for OSX from `shared_library/Darwin/32bit/ <http://sourceforge.net/projects/coolprop/files/CoolProp/6.0.0/shared_library/Darwin/32bit/>`_.  Development binaries coming from the buildbot server can be found at `shared_library/Darwin/32bit/ <http://sourceforge.net/projects/coolprop/files/CoolProp/nightly/shared_library/Darwin/32bit/>`_. Place the downloaded file called libCoolProp.dylib in the ${HOME}/lib folder (make this folder if needed).

**Part 2:**

Download the xlam from `MicrosoftExcel <http://sourceforge.net/projects/coolprop/files/CoolProp/6.0.0/MicrosoftExcel>`_ or the development version from `MicrosoftExcel <http://sourceforge.net/projects/coolprop/files/CoolProp/nightly/MicrosoftExcel>`_.

1.  Open Excel 2011 for Mac

2.  Go to the menu ``Tools–>Add-Ins``

3.  Click the ``Select...`` button

4.  Browse to the file CoolProp.xlam you downloaded, select it

5.  Make sure the CoolProp Add-in is selected.

6.  Add this code to a cell - it should work: ::

    =PropsSI("T","P",101325,"Q",0,"Water")


If it doesn’t work and you get error number 53, it might be because you have a 64-bit .dylib file and you want a 32-bit .dylib file.  For instance when you run the ``file`` command on your .dylib, you should see something like: ::

    $ file libCoolProp.dylib
    libCoolProp.dylib: Mach-O dynamically linked shared library i386


User-compiled Binaries
------------------------

**Common Requirements**

Compilation of the Excel wrapper requires a few `common wrapper pre-requisites <http://www.coolprop.org/coolprop/wrappers/index.html#wrapper-common-prereqs>`_


**Build (Windows)**

MS Excel requires the CoolProp shared library, or Dynamic Link Library (DLL) on Windows, and will use either the 64-bit version or the 32-bit, ``__stdcall`` version.  The instructions here are for a 64-bit windows system that will use Microsoft Visual Studio 2010 to compile *both* the 64-bit and 32-bit versions of the DLL.  

.. code-block:: bash

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
  copy build\32bit__stdcall\CoolProp.dll c:\CoolProp\CoolProp.dll
  copy build\64bit\CoolProp.dll c:\CoolProp\CoolProp_x64.dll

The above script should be adjusted for your specific compiler, replacing "Visual Studio 10" with your compiler name/release.

**Build (OSX)**

On OSX there is no calling convention to worry about, and only the 32-bit compilation is needed.  You can force 32-bit compilation using -DFORCE_BITNESS_32=ON.  These instructions will compile both the 32-bit library needed for Excel on OSX.

.. code-block:: bash

  # Check out the sources for CoolProp
  git clone https://github.com/CoolProp/CoolProp --recursive
  # Move into the folder you just created
  cd CoolProp
  # Make a build folder
  mkdir build && cd build
  # Generate builder
  cmake .. -DCOOLPROP_SHARED_LIBRARY=ON -DFORCE_BITNESS_32=ON -DCMAKE_BUILD_TYPE=Release
  # Build
  cmake --build .
  cd ..
  # Copy the generated DLL
  cp build/libCoolProp.dylib ${HOME}/lib

Usage
------------------------

The following CoolProp funcitons are implemented as Excel Functions (see CoolProp.xlsx for examples)


- **get_global_param_string("param")** 
- **PropsSI("OutputName","Input1Name",Value1,"Input2Name",Value2,"FluidString")** 
- **Props1SI("FluidString","OutputName")** 
- **PhaseSI("Input1Name",Value1,"Input2Name",Value2,"FluidString")** 
- **HAPropsSI("OutputName","Input1Name",Value1,"Input2Name",Value2,"Input3Name",Value3)** 

  
and utility routine

  
- **MixtureString(Names,Fractions)** 

  
where


"param" =      
    "version", "gitrevision", "fluids_list","parameter_list", or "predefined_mixtures"

"OutputName" =
    double quoted name of the fluid property to output
   
"Input1Name" =
    double quoted name of the first input state property
 
Value1 =  
    the value, in SI units, of the first input state property

"Input2Name" =  
    double quoted name of the second input state property

Value2 =  
    the value, in SI units, of the second input state property

"Input3Name" =  
    double quoted name of the third input state property; typically humidity ratio, "R"

Value3 =   
    the value, in SI units, of the third input state property
  
"FluidString" =
    double quoted name of the fluid to evaluate

Names =
    range of cells containing fluid names (as text)

Fractions =
    range of cells containing mole fractions of the corresponding fluids in the **Names** range


| See `CoolProp.org <http://www.coolprop.org>`_ for more details on usage.
| See **TestExcel.xlsx** for examples of usage of these functions.
