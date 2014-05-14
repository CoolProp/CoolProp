REM ******** set the variables ************
REM call both to ensure that one works
call "C:\Program Files\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"

REM ******* compile all the sources from CoolProp ***************
cl /c /Ox /fp:fast /I../../../CoolProp /EHsc Example.cpp
if %errorlevel% neq 0 exit /b %errorlevel%
cl /c /Ox /MP3 /fp:fast /I../../../CoolProp /EHsc ../../../CoolProp/*.cpp
if %errorlevel% neq 0 exit /b %errorlevel%

link *.obj /OUT:Example.exe
if %errorlevel% neq 0 exit /b %errorlevel%
erase *.obj

call Example > Output.txt
if %errorlevel% neq 0 exit /b %errorlevel%
erase Example.exe