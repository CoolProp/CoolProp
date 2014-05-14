REM ******** call the Microsoft Visual studio (these are for VS2008) ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat" x64
call "C:\Program Files\Microsoft Visual Studio 10.0\VC\vcvarsall.bat" x64

REM ******* compile all the sources ***************
cl /c /I../../../CoolProp /EHsc /DEXTERNC /DCOOLPROPMATHCADWRAPPER_EXPORTS "CoolProp 4.0 Beta Mathcad Prime 3.0 Wrapper.cpp"
cl /c /I../../../CoolProp /EHsc /DEXTERNC ../../../CoolProp/*.cpp

link /DLL *.obj mcaduser.lib /OUT:CoolPropMathcadWrapper.dll
dumpbin /EXPORTS CoolPropMathcadWrapper.dll > exports.txt
erase *.obj