
@echo on
erase CoolProp.dll
REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
call "C:\Program Files\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

REM ******* compile all the sources ***************
swig -java -c++ -outcurrentdir -I../../include ../../src/CoolProp.i
cl /c /MP3 /I../../include /I"C:\Program Files\Java\jdk1.7.0_40\include" /I"C:\Program Files\Java\jdk1.7.0_40\include\win32" /EHsc *.cxx
cl /c /MP3 /I../../include /I"C:\Program Files\Java\jdk1.7.0_40\include" /I"C:\Program Files\Java\jdk1.7.0_40\include\win32" /EHsc ../../src/*.cpp
cl /c /MP3 /I../../include /I"C:\Program Files\Java\jdk1.7.0_40\include" /I"C:\Program Files\Java\jdk1.7.0_40\include\win32" /EHsc ../../src/Backends/*.cpp
cl /c /MP3 /I../../include /I"C:\Program Files\Java\jdk1.7.0_40\include" /I"C:\Program Files\Java\jdk1.7.0_40\include\win32" /EHsc ../../src/Fluids/*.cpp
cl /c /MP3 /I../../include /I"C:\Program Files\Java\jdk1.7.0_40\include" /I"C:\Program Files\Java\jdk1.7.0_40\include\win32" /EHsc ../../src/Tests/*.cpp

link /DLL *.obj /OUT:CoolProp.dll
dumpbin /EXPORTS CoolProp.dll > exports.txt
mkdir win32
move CoolProp.dll win32
erase *.obj
erase *.exp
erase *.lib