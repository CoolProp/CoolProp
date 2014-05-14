Steps
-----

1. Download a WinXP 32-bit version virtual box image
2. Install Winxp 32-bit virtual box image, install 7zip, tortoisegit, git
3. Take the updated_vxworks63gccdist.ZIP from (http://www.ni.com/white-paper/5694/en/) and extract to c:\\gccdist so that you have folders c:\\gccdist\\docs, c:\\gccdist\\supplemental, c:\\gccdist\\WindRiver
4. Check out coolprop sources
5. Open a console and cd to CoolProp sources
6. cd to wrappers/Labview/vxWorks
7. run build.bat (It can take a really long time (several minutes) for some reason for Mixtures.cpp, be patient)
8. Upload the generated file CoolProp.out to ni-rt/system folder using the National Instruments Measurement & Automation Explorer.  Right click on the unit, then file transfer. Also see http://www.ni.com/white-paper/3365/en/

Debugging
---------
1. Go to the IP address of the cRIO in your browser (I found IE worked better than Chrome for this)
2. Create a new session
3. Go to the console (You will need to have this enabled on your cRIO)
4. Try to run the VI on the cRIO, it will try to load the module and print out the errors that occur (if any) when loading the module

Fixing Problems
---------------
Make variables static if at all possible
Do not enable Catch (it uses language features that are incompatible with vxWorks cRIO)

To be determined: What do do about _Dtest, _Nan, _Inf ??

See Also
--------
https://decibel.ni.com/content/docs/DOC-13537