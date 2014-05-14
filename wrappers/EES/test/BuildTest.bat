cd ..
call BuildDLL
copy CoolProp.dll test
cd test

REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
call "C:\Program Files\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

REM ******* compile all the sources from CoolProp***************
cl /EHsc /I../../../CoolProp testdll.cpp ../../../CoolProp/CoolPropTools.cpp

erase *.obj

testdll.exe