REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat" amd64
call "C:\Program Files\Microsoft Visual Studio 10.0\VC\vcvarsall.bat" amd64

REM ******* compile all the sources ***************
cl -Gz mult.c -LD -link

dumpbin /EXPORTS mult.dll > exports.txt
erase *.obj
erase *.exp
erase *.lib