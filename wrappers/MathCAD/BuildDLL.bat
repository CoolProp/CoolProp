REM ******** call the Microsoft Visual studio (these are for VS2008) ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
call "C:\Program Files\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

REM ******* compile all the sources ***************
cl /c /MP3 /I../../CoolProp /EHsc /DEXTERNC ../../CoolProp/*.cpp
cl /c /I../../CoolProp /EHsc /DEXTERNC /DCOOLPROPMATHCADWRAPPER_EXPORTS CoolPropMathcad.cpp

link /DLL *.obj mcaduser.lib /OUT:CoolPropMathcadWrapper.dll /ENTRY:"DllEntryPoint"
dumpbin /EXPORTS CoolPropMathcadWrapper.dll > exports.txt
erase *.obj
REM ~ erase *.exp
REM ~ erase *.lib