
@echo on
copy ..\..\..\wrappers\Java\Example.java Example.java

REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
call "C:\Program Files\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

REM ******* compile all the sources ***************
swig -java -c++ -outcurrentdir ../../../CoolProp/CoolProp.i
if %errorlevel% neq 0 exit /b %errorlevel%
cl /c /I../../../CoolProp /I"C:\Program Files\Java\jdk1.7.0_40\include" /I"C:\Program Files\Java\jdk1.7.0_40\include\win32" /EHsc *.cxx
if %errorlevel% neq 0 exit /b %errorlevel%
cl /c /I../../../CoolProp /I"C:\Program Files\Java\jdk1.7.0_40\include" /I"C:\Program Files\Java\jdk1.7.0_40\include\win32" /EHsc ../../../CoolProp/*.cpp
if %errorlevel% neq 0 exit /b %errorlevel%
link /DLL *.obj /OUT:CoolProp.dll
if %errorlevel% neq 0 exit /b %errorlevel%
erase *.obj
erase *.exp
erase *.lib

javac *.java
if %errorlevel% neq 0 exit /b %errorlevel%
java Example > Output.txt
if %errorlevel% neq 0 exit /b %errorlevel%

erase *.java
erase *.class
erase CoolProp.dll
erase CoolProp_wrap.cxx
copy ..\..\..\wrappers\Java\Example.java Example.java