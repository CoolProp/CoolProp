
.. _Maple:

*************
Maple Wrapper
*************

A very simple example can be found here: :download:`sample_file.mw`.  You will need to download the release 64-bit shared library for windows from :sfdownloads:`sourceforge <shared_library/Windows/64bit>` or the development version of the shared library from :sfnightly:`the snapshot server <shared_library/Windows/64bit>`.

Place the downloaded DLL in c:\\CoolProp and rename to CoolProp_x64.dll

On other platforms, the protocol should be more or less the same.  Download a DLL appropriate to the version of Maple you are using from :sfdownloads:`sourceforge <shared_library>` or the development version from :sfnightly:`the nightly snapshots <shared_library>`.  You will probably need to change the path to the shared library in the Maple file, but aside from that, it should work fine.  File an issue at `github issues <https://github.com/CoolProp/CoolProp/issues>`_ if that is not the case.

A more in-depth example has been put together by the Maple Enginers: :download:`Analysis of a Refrigeration Cycle with CoolProp.mw`