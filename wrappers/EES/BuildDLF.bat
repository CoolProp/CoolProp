REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
call "C:\Program Files\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

REM ******* compile all the sources ***************
cl /c /MP3 /I../../CoolProp /EHsc ../../CoolProp/*.cpp

REM ******* compile the wrapper ***************
cl /c /EHsc /I../../CoolProp main.cpp
link /DLL *.obj /OUT:COOLPROP_EES.dlf

erase *.obj