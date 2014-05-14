REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
call "C:\Program Files\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

REM ******* compile all the sources ***************
cl /c /MP3 /O2 /EHsc /DNDEBUG /I../../include /DCOOLPROP_LIB ../../src/*.cpp
cl /c /MP3 /O2 /EHsc /DNDEBUG /I../../include /DCOOLPROP_LIB ../../src/Backends/*.cpp
cl /c /MP3 /O2 /EHsc /DNDEBUG /I../../include /DCOOLPROP_LIB ../../src/Fluids/*.cpp
cl /c /MP3 /O2 /EHsc /DNDEBUG /I../../include /DCOOLPROP_LIB ../../src/Tests/*.cpp

link /DLL *.obj /OUT:CoolProp.dll
dumpbin /EXPORTS CoolProp.dll > exports.txt

erase *.obj
erase *.exp
erase *.lib