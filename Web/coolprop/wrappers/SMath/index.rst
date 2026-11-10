.. _SMath:

*******************
SMathStudio Wrapper
*******************

Pre-compiled Binaries
=====================

.. Note::
    The SMath wrapper files moved to the SMath open-source SVN repository around 2018 and the *official* source is no longer tracked in the `CoolProp/CoolProp repository <https://github.com/CoolProp/CoolProp/tree/master/wrappers/SMath>`__. Although the stale wrapper code is still there, it does not compile with current versions of CoolProp.  Information on the current *official* wrapper / Plugin can be found at `SMath > CoolProp Wrapper <https://smath.info/en-US/view/ca92ef03-c7da-4888-98ad-528482733e2f/summary>`__.  Instructions for using the *official* CoolProp Plugin for SMath are shown below.

To Use
------

Pre-compiled binaries are loaded through the SMath Plugin Manager in SMath Studio itself (Tested in v1.5.0.9678).  Simply,

1. Open SMath Studio
2. From the main menu, select **Tools | Plugins**
3. Change the Source location (top right corner of panel) from **Local storage** to **Online Gallery**
4. Type **CoolProp** in the search field at the bottom of the panel and select the CoolProp version from the gallery window
5. Click the **Install** button underneath the gallery

A list of the CoolProp function implementation is available on the official `SMath > CoolProp Wrapper > Features <https://smath.info/en-US/view/ca92ef03-c7da-4888-98ad-528482733e2f/features>`__ page.  All functions mirror the high-level interface calls prefixed with *CoolProp_*.

**Caveat:** Prefixing the high-level names does not always produce a valid SMath function. For example,
- Use `CoolProp_Props` for `PropsSI`
- Use `CoolProp_Props1` for `Props1SI`
- Use `CoolProp_HAProps` for `HAPropsSI`

**Wrapper questions and Issues** should be posted in the `SMath Forum (CoolProp Wrapper Topic) <https://smath.com/en-US/forum/topic/gcKGsU/CoolProp-Wrapper>`__. (Requires an SMath login to post)

.. Warning::
    This *official* version is still based on CoolProp 6.4.1 as evidenced by the calling::
        
        CoolProp_get_global_param_string("version") = 6.4.1

    Although the release date for the wrapper is more recent.

.. seealso::
    From the `SMath Forum (CoolProp Wrapper Topic) <https://smath.com/en-US/forum/topic/gcKGsU/CoolProp-Wrapper>`__ Post (#83) there is an *unoficial* version of the the Plugin (source and distribution binary) that supports **CoolProp v8.0.0** located on GitHub at `defencedog/smath-devel <https://github.com/defencedog/smath-devel/tree/main/plugin/CoolPropWrapper_src_win>`__.  The README in this repository also includes information on running the Plugin and building the wrapper on other platforms if needed.
    
    The installable zip file is located in the `dist directory <https://github.com/defencedog/smath-devel/tree/main/plugin/CoolPropWrapper_src_win/dist>`__ and supports Windows (32 and 64-bit) and  aarch64 (Linux).  Download the compressed zip file to any location on your machine (e.g., Downloads folder).  Use the Plugin Manager (as above), but point to the downloaded zip file in **Local storage**  using the **Install from file..** button (*SMath v1.5 or later*) to browse to the zip file location.

User-Compiled Binaries (Legacy Instructions)
============================================

.. Attention::
    While the legacy wrapper files still exist on GitHub in the CoolProp/CoolProp repository, these do no currently compile with more recent versions of CoolProp (7.x and onward).  Legacy instructions remain below for reference should the code ever be tracked again in the CoolProp repository in the future. 

Common Requirements
-------------------

Compilation of the SMath Studio wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Native Wrapper
--------------

1. Check out coolprop::

    git clone https://github.com/CoolProp/CoolProp
    
2. Create a build directory ``build``::

    mkdir build

3. In the ``build`` directory, run the command::

    cmake .. -DCOOLPROP_SMATH_MODULE=ON
    
  This will inject the version number in the appropriate format into a few template files using CMake
  
4. Open the solution file in ``build`` directory.

5. Make sure the mode is set to Release (not Debug).  Build the project, the generated and copied files will be in ``build/wrappers/SMath/coolprop_wrapper/bin/Release``.

6. From the :sfdownloads:`download page <shared_library>`, download the 64-bit DLL ``CoolProp.dll`` file from ``Windows/64bit`` folder and rename to ``CoolProp_x64.dll`` and place with the files in ``build/wrappers/SMath/coolprop_wrapper/bin/Release``.  Download the 32-bit ``__stdcall`` DLL ``CoolProp.dll`` from ``Windows/32bit__stdcall_calling_convention`` and rename to ``CoolProp_x86.dll`` and place with the files in ``build/wrappers/SMath/coolprop_wrapper/bin/Release``.

7. Run the ``build_zip.bat`` file that is in the ``build/wrappers/SMath`` folder.  It will create a zip file with the needed files. 

8. To install, unzip the ``coolprop_wrapper.7z`` zip file and run the ``install.bat`` script that was in the zip file


More information on SMath Studio is given under: 

- https://en.wikipedia.org/wiki/SMath_Studio  
- https://smath.com/en-US/

.. raw:: html

    <div style="width: 100%; height: 2px; background-color: #72727254; margin-top: 2rem; margin-bottom: 0px; display: block; clear: both;"></div>
    <p style="margin-top: 4px; margin-bottom: 2rem; font-style: italic;">
        Page Last Updated: 10-10-2026
    </p>