
mkdir 3.6.1
mkdir 3.6.2
mkdir 3.6.4

swig -octave -c++ -outcurrentdir -o CoolProp_wrap.cpp  -I../../include ../../src/CoolProp.i

REM ~ call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"
REM ~ REM : %%~nf gives just the file name, no path or extension
REM ~ REM : %%f gives the full path and extension
REM ~ mkoctfile -v -c -I..\..\CoolProp -o CoolProp_wrap.o CoolProp_wrap.cpp
REM ~ if %errorlevel% neq 0 exit /b %errorlevel%
REM ~ for %%f in (..\..\CoolProp\*.cpp) do mkoctfile -v -c -I..\..\CoolProp -o %%~nf.o %%f
REM ~ if %errorlevel% neq 0 exit /b %errorlevel%
REM ~ mkoctfile -v -o 3.6.1\CoolProp.oct *.o
REM ~ if %errorlevel% neq 0 exit /b %errorlevel%
REM ~ erase *.o

call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
for %%f in (..\..\src\*.cpp) do c:\octave3_6_2\bin\mkoctfile -v -c -I..\..\include -o %%~nf.o %%f
for %%f in (..\..\src\Backends\REFPROP\*.cpp) do c:\octave3_6_2\bin\mkoctfile -v -c -I..\..\include -o %%~nf.o %%f
for %%f in (..\..\src\Backends\Helmholtz\*.cpp) do c:\octave3_6_2\bin\mkoctfile -v -c -I..\..\include -o %%~nf.o %%f
for %%f in (..\..\src\Backends\Helmholtz\Fluids\*.cpp) do c:\octave3_6_2\bin\mkoctfile -v -c -I..\..\include -o %%~nf.o %%f
for %%f in (..\..\src\Backends\Incompressible\*.cpp) do c:\octave3_6_2\bin\mkoctfile -v -c -I..\..\include -o %%~nf.o %%f
for %%f in (..\..\src\Tests\*.cpp) do c:\octave3_6_2\bin\mkoctfile -v -c -I..\..\include -o %%~nf.o %%f

c:\octave3_6_2\bin\mkoctfile -v -c -I..\..\include -o CoolProp_wrap.o CoolProp_wrap.cpp
c:\octave3_6_2\bin\mkoctfile -v -o 3.6.2\CoolProp.oct *.o
REM ~ erase *.o

REM ~ call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
REM ~ for %%f in (..\..\CoolProp\*.cpp) do c:\octave3_6_4\bin\mkoctfile -v -c -I..\..\CoolProp -o %%~nf.o %%f
REM ~ c:\octave3_6_4\bin\mkoctfile -v -c -I..\..\CoolProp -o CoolProp_wrap.o CoolProp_wrap.cpp
REM ~ c:\octave3_6_4\bin\mkoctfile -v -o 3.6.4\CoolProp.oct *.o
REM ~ erase *.o