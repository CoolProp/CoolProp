An example of CoolProp for Visual Basic .NET
============================================

Pierre-Yves Manach and Ian Bell, 2013

To Use
------
1. Put the Windows CoolProp Shared Library (DLL) from `SourceForge <https://sourceforge.net/projects/coolprop/files/CoolProp>`_ in ``C:\CoolProp`` (or adjust the paths in coolprop.vb)  

   This file comes in two flavors 32-bit and 64-bit, but you will want the Windows 64-bit version.  You can get this file two ways:  

   - Download and run the Windows Installer (v6.1 or later) from `SourceForge <https://sourceforge.net/projects/coolprop/files/CoolProp>`_.  This will place **all** versions of files into the user's profile under the ``AppData\Roaming\CoolProp`` folder.  The VB.NET wrapper will use the  Windows 64-bit version, ``CoolProp.dll``, which you can relocate to ``C:\CoolProp`` if needed.  

   - Download ``CoolProp.dll`` directly from the ``shared_library/Windows/64-bit`` directory on `SourceForge <https://sourceforge.net/projects/coolprop/files/CoolProp>`_.  


2. Open the solution (*WpfApplication1.sln*) in Visual Studio 2012 or later.  


3. Run!