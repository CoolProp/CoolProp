REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat" amd64
call "C:\Program Files\Microsoft Visual Studio 10.0\VC\vcvarsall.bat" amd64

REM ******* compile all the sources ***************
cl /c  /MP3 /I../../CoolProp /EHsc /DCOOLPROP_LIB ../../CoolProp/*.cpp

link /DLL CoolProp.obj *.obj /OUT:CoolProp_x64.dll
copy CoolProp_x64.dll c:\CoolPropx64.dll

dumpbin /EXPORTS CoolProp_x64.dll > exports_x64.txt
erase *.obj
erase *.exp
erase *.lib
