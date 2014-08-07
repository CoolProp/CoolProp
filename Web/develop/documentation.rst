.. _developer_documentation:

*************
Documentation
*************

Build Sphinx documentation
--------------------------

1. Check out the sources in the CoolProp/Web folder

    git clone https://github.com/CoolProp/CoolProp --recursive

Make a virtualenv - helps to keep sane ! 

    virtualenv ~/env/py27
    source ~/env/py27/bin/activate # Turn on this virtual env, should see (py27) in your command shell next to the prompt to tell you that environment is active

Then install prerequisites into this virtualenv:
  
    pip install Cython
    pip install CoolProp sphinx sphinxcontrib-doxylink cloud-sptheme ipython matplotlib numpy scipy

To build the documentation, go into the CoolProp/Web folder and run:

    make html
  
Build Doxygen documentation
---------------------------

All the configuration is done in the ``Doxyfile`` file.

1. If you don't have doxygen, install it from `doxygen downloads <http://www.stack.nl/~dimitri/doxygen/download.html>`.  Or on linux a ``sudo apt-get install doxygen`` should do it.

2. Simply fire up a shell in the root of the repo, and type 

  doxygen Doxyfile
  
