
@echo off

call:defineEnv 
REM Change to "call:defineEnv x86" or "call:defineEnv amd64" according to your needs

cl /c /MP /I../../CoolProp /EHsc /DCOOLPROP_LIB ../../CoolProp/*.cpp

link /DLL *.obj /OUT:CoolProp.dll
lib CoolProp.obj *.obj /OUT:CoolProp.lib

dumpbin /EXPORTS CoolProp.dll > exportsDLL.txt
dumpbin /HEADERS CoolProp.lib > exportsLIB.txt

erase *.obj
erase *.exp

goto:eof

rem ******** define some general functions ************
:defineEnv    - set the variables, accepts one argument
set stdpaths="C:\Program Files (x86)\Microsoft Visual Studio ","C:\Program Files\Microsoft Visual Studio "
rem this order assures that the latest version is used...
set versions="8.0","9.0","10.0","11.0","12.0"
set relPaths="\VC\vcvarsall.bat"
set filename=""
for %%i in (%stdpaths%) do (
  for %%j in (%versions%) do (
    for %%k in (%relPaths%) do (
      call:loadScript "%%~i%%~j%%~k" %%~j %~1
    )
  )
)
goto:eof

:loadScript
rem echo "%~1" "%~2"
if exist "%~1" (
  echo .
  echo Found Visual Studio v. %~2
  echo Calling "%~1" %~3 
  call "%~1" %~3 
  echo .
) 
REM else (
REM   echo Could not find "%~1"
REM )
goto:eof