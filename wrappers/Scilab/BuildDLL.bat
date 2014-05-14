
REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"

REM ******* compile all the sources ***************
cl /MP3 /O2 /Oi /GL /I "..\..\CoolProp" /D "NDEBUG" /D "_WINDOWS" /D "_USRDLL" /D "CONVENTION=__cdecl"  /D "EXPORT_CODE=extern \"C\" __declspec(dllexport)" /D "_WINDLL" /FD /EHsc /MT /Gy /W3 /c /Zi /TP scilab_functions.cpp
cl /MP3 /O2 /Oi /GL /I "..\..\CoolProp" /D "NDEBUG" /D "_WINDOWS" /D "_USRDLL" /D "CONVENTION=__cdecl"  /D "EXPORT_CODE=extern \"C\" __declspec(dllexport)" /D "_WINDLL" /FD /EHsc /MT /Gy /W3 /c /Zi /TP ..\..\CoolProp\*.cpp

REM ******** link into DLL *****************
link /OUT:"CoolProp.dll" /INCREMENTAL:NO /NOLOGO /DLL /MANIFEST:NO /DEBUG /SUBSYSTEM:WINDOWS /OPT:REF /OPT:ICF /LTCG /DYNAMICBASE /NXCOMPAT /ERRORREPORT:PROMPT  *.obj

dumpbin /EXPORTS CoolProp.dll > exports.txt
erase *.obj
erase *.pdb
erase *.idb
erase *.lib
erase *.exp

REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat" amd64

REM ******* compile all the sources ***************
cl /MP3 /O2 /Oi /GL /I "..\..\CoolProp" /D "NDEBUG" /D "_WINDOWS" /D "_USRDLL" /D "COOLPROP_LIB" /D "_WINDLL" /FD /EHsc /MT /Gy /W3 /c /Zi /TP scilab_functions.cpp
REM /D "CONVENTION=__cdecl" 
cl /MP3 /O2 /Oi /GL /I "..\..\CoolProp" /D "NDEBUG" /D "_WINDOWS" /D "_USRDLL" /D "COOLPROP_LIB" /D "_WINDLL" /FD /EHsc /MT /Gy /W3 /c /Zi /TP ..\..\CoolProp\*.cpp

REM ******** link into DLL *****************
link /OUT:"CoolProp_x64.dll" /INCREMENTAL:NO /NOLOGO /DLL /MANIFEST:NO /DEBUG /SUBSYSTEM:WINDOWS /OPT:REF /OPT:ICF /LTCG /DYNAMICBASE /NXCOMPAT /ERRORREPORT:PROMPT  *.obj

dumpbin /EXPORTS CoolProp_x64.dll > exports_x64.txt
erase *.obj
erase *.pdb
erase *.idb
erase *.lib
erase *.exp
