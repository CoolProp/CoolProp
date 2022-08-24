@echo off
REM Some Windows commands to compile the Fortran interface for F77
copy C:\Users\Belli\Documents\Code\coolprop-v5-new\build\32bitcdecl\Debug\Coolprop.dll .
C:\MinGW\bin\gfortran -g -o example examplef.for -L. -lCoolProp
example