@echo off

:: path to test
SET testPath="."

:: launch batch test
"%ProgramFiles(x86)%\SMath Studio\SMathStudio_Desktop.exe" -silent -t %testPath%
::"%ProgramFiles%\SMath Studio\SMathStudio_Desktop.exe" -silent -t %testPath%
::"C:\Program Files (x86)\SMath Studio\SMathStudio_Desktop.exe" -silent -t %testPath%

ECHO.
ECHO Tested files in: %testPath%
ECHO.

:: prevents auto-close
pause