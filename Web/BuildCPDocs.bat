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
python -V
python -c "import IPython.sphinxext.ipython_console_highlighting as I; print I.__file__"
mingw32-make html