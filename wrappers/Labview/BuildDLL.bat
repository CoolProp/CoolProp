REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"
REM ******* compile all the sources ***************

cl /MP3 /O2 /Oi /GL /I "..\..\CoolProp" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_USRDLL" /D "LABVIEW_RT_EXPORTS" /D "COOLPROP_LIB" /D "CONVENTION=__cdecl" /D "_WINDLL" /FD /EHsc /MT /Gy /W3 /c /Zi /TP ..\..\CoolProp\*.cpp

link /OUT:"CoolProp.dll" /INCREMENTAL:NO /NOLOGO /DLL /MANIFEST:NO /DEBUG /SUBSYSTEM:WINDOWS /OPT:REF /OPT:ICF /LTCG /DYNAMICBASE /NXCOMPAT /MACHINE:X86 /ERRORREPORT:PROMPT kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib *.obj

dumpbin /EXPORTS CoolProp.dll > exports.txt
erase *.obj
erase *.pdb
erase *.idb
erase *.lib
erase *.exp
