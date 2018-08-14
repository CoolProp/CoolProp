.. _SMath:

*******************
SMathStudio Wrapper
*******************

Pre-compiled Binaries
=====================

To Use
------

Pre-compiled binaries can be downloaded from :sfdownloads:`SMath`.  Development binaries coming from the buildbot server can be found at :sfnightly:`SMath`.

Extract the files in the .7z file using 7-zip.  Run the install.bat file in the zip file to copy the files to the right location.  It will make the installed version of CoolProp the default CoolProp in SMath.

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the SMath Studio wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Native Wrapper
--------------

0. Check out coolprop::

    git clone https://github.com/CoolProp/CoolProp --recursive
    
1. Create a build directory ``build``::

    mkdir build

2. In the ``build`` directory, run the command::

    cmake .. -DCOOLPROP_SMATH_MODULE=ON
    
  This will inject the version number in the appropriate format into a few template files using CMake
  
3. Open the solution file in ``build`` directory.

4. Make sure the mode is set to Release (not Debug).  Build the project, the generated and copied files will be in ``build/wrappers/SMath/coolprop_wrapper/bin/Release``.

5. From the :sfdownloads:`download page <shared_library>`, download the 64-bit DLL ``CoolProp.dll`` file from ``Windows/64bit`` folder and rename to ``CoolProp_x64.dll`` and place with the files in ``build/wrappers/SMath/coolprop_wrapper/bin/Release``.  Download the 32-bit ``__stdcall`` DLL ``CoolProp.dll`` from ``Windows/32bit__stdcall_calling_convention`` and rename to ``CoolProp_x86.dll`` and place with the files in ``build/wrappers/SMath/coolprop_wrapper/bin/Release``.

6. Run the ``build_zip.bat`` file that is in the ``build/wrappers/SMath`` folder.  It will create a zip file with the needed files. 

7. To install, unzip the ``coolprop_wrapper.7z`` zip file and run the ``install.bat`` script that was in the zip file

Old Method
==========

There are two ways to link CoolProp and SMathStudio :
Remark: I used the first one, the second one isn't quite clear to me.

1) The easiest way is to download and install the portable SMath Studio distribution with Maxima, CoolProp and many other plugins :
http://smath.info/wiki/%28S%28p1gq0245uny440554qdivzik%29%29/SMath%20with%20Plugins.ashx

2) To link CoolProp to SMathStudio_Desktop follow the identical instructions as given below for MathCAD 15 :
http://en.smath.info/forum/yaf_postst2568_CoolProp.aspx#post13310

Notes :
- In October 2014 SMath Studio with CoolProp worked on Windows only. 
- The portable SMath Studio distribution is provided on the SMath Studio site, however so far it's indicated as unofficial (October 2014) 

More information on SMath Studio is given under :
http://en.wikipedia.org/wiki/SMath_Studio
http://en.smath.info/