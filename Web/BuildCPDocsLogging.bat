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

rem sphinx-apidoc -f -o apidoc ../CoolProp
make html 2>&1 | wtee log.txt
