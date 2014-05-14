REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
call "C:\Program Files\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

REM ******* compile all the sources ***************
cl /c /MP3 /I../../CoolProp /EHsc /DCOOLPROP_LIB /D"CONVENTION=__cdecl" ../../CoolProp/*.cpp

mkdir dll
link /DLL *.obj /OUT:dll\CoolProp.dll
copy dll\CoolProp.dll .
dumpbin /EXPORTS CoolProp.dll > exports.txt
erase *.obj
erase *.exp