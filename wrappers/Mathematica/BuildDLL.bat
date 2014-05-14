REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat" amd64
call "C:\Program Files\Microsoft Visual Studio 10.0\VC\vcvarsall.bat" amd64

copy "C:\Program Files\Wolfram Research\Mathematica\9.0\SystemFiles\IncludeFiles\C\WolframLibrary.h" .

REM ******* compile all the sources ***************
cl /c /MP3 /I../../CoolProp /EHsc CoolPropMathematica.cpp
cl /c /MP3 /I../../CoolProp /EHsc ../../CoolProp/*.cpp

link /DLL *.obj /OUT:CoolProp.dll

REM folder obtained from FileNameJoin[{$BaseDirectory, "SystemFiles", "LibraryResources", $SystemID}] in Mathematica
copy CoolProp.dll "C:\ProgramData\Mathematica\SystemFiles\LibraryResources\Windows-x86-64\"

dumpbin /EXPORTS CoolProp.dll > exports.txt
erase *.obj
erase *.exp
erase *.lib
erase WolframLibrary.h