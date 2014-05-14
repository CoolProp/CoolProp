Generating the documentation using DOXYGEN
==================================================
Christoph Richter, ch.richter@tu-bs.de, 2007-01-19

1. Get the latest version of doxygen that can be found at
   http://www.stack.nl/~dimitri/doxygen/
   and install it on your system

2. Open a command line shell and go to
   \path-to-external-media-library\Projects\Sources

3. Type in the following command
   doxygen Doxyfile
   This will create the html documentation in
   \path-to-external-media-library\Documentation\html
   and the LaTeX documentation in
   \path-to-external-media-library\Documentation\latex

Note that the first chapter of the documentation can be found in "documentation.h".