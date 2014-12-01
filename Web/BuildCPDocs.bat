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

sphinx-apidoc -T -f -e -o apidoc C:\\Miniconda\\lib\\site-packages\\coolprop-5.0.0-py2.7-win-amd64.egg\\CoolProp
mingw32-make html_release