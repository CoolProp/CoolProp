call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"

swig -octave -c++ -outcurrentdir -o CoolProp_wrap.cpp ../../../CoolProp/CoolProp.i
    if %errorlevel% neq 0 exit /b %errorlevel%

REM : %%~nf gives just the file name, no path or extension
REM : %%f gives the full path and extension
for %%f in (..\..\..\CoolProp\*.cpp) do mkoctfile -v -c -I..\..\CoolProp -o %%~nf.o %%f
   if %errorlevel% neq 0 exit /b %errorlevel%
mkoctfile -v -c -I..\..\..\CoolProp -o CoolProp_wrap.o CoolProp_wrap.cpp
   if %errorlevel% neq 0 exit /b %errorlevel%
mkoctfile -v -o CoolProp.oct *.o
   if %errorlevel% neq 0 exit /b %errorlevel%
erase *.o
erase CoolProp.lib
erase CoolProp.exp
erase CoolProp_wrap.cpp
octave Example.m > Output.txt
    if %errorlevel% neq 0 exit /b %errorlevel%