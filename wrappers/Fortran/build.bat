@echo off
REM Some Windows commands to compile the Fortran interface for F77
setlocal 
set comPath=C:\MinGW32-xy\bin
set FC=%comPath%\gfortran
set CC=%comPath%\gcc
REM set comPath=C:\MinGW32\bin
REM set FC=%comPath%\g77
REM set CC=%comPath%\gcc
%FC% -g -c examplef.for
%CC% -g -c coolpropIF.c 
%FC% -g -o example examplef.o coolpropIF.o -L../SharedLibrary -lCoolProp
endlocal 