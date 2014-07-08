@echo off
REM ~ make latex
REM ~ cd _build/latex
REM ~ pdflatex CoolPropdoc.tex
REM ~ pdflatex CoolPropdoc.tex
REM ~ pdflatex CoolPropdoc.tex
REM ~ pdflatex CoolPropdoc.tex
REM ~ pdflatex CoolPropdoc.tex
REM ~ copy /Y CoolPropdoc.pdf ..\..\_static\
REM ~ cd ..\..

sphinx-apidoc -T -f -o apidoc ../CoolProp
mingw32-make html