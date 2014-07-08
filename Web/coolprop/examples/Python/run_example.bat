echo off
copy ..\..\..\wrappers\Python\examples\Example.py
if %errorlevel% neq 0 exit /b %errorlevel%
python Example.py > Output.txt
if %errorlevel% neq 0 exit /b %errorlevel%