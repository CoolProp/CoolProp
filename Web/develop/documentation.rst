.. _developer_documentation:

*************
Documentation
*************

Build Sphinx documentation
--------------------------

0. On a Mac system, you need a working Latex distribution and some more command line tools that can be installed 
   with ``brew install jpeg libpng libtiff libtool imagemagick``. Additonally, you should add the local string 
   to your bash environment in ``/Users/username/.bash_profile``::

    export LC_ALL=en_US.UTF-8
    export LANG=en_US.UTF-8
    

1. Check out the sources in the CoolProp/Web folder::

    git clone https://github.com/CoolProp/CoolProp --recursive

2. Make a virtualenv::

    virtualenv --distribute ~/env/py27
    source ~/env/py27/bin/activate # Turn on this virtual env, should see (py27) in your command shell next to the prompt to tell you that environment is active

3. Then install prerequisites into this virtualenv::

    pip install --upgrade setuptools cython
    pip install cython numpy scipy matplotlib python-dateutil pytz pandas
    pip install sphinx ipython sphinxcontrib-bibtex sphinxcontrib-doxylink sphinxcontrib-napoleon cloud-sptheme
    pip install wxpython


4. To build the documentation, go into the CoolProp/Web folder and run::

    make html

5. Move the generated docs in ``_build`` to wherever you want


Build Doxygen documentation
---------------------------

All the configuration is done in the ``Doxyfile`` file.

1. If you don't have doxygen, install it from `doxygen downloads <http://www.stack.nl/~dimitri/doxygen/download.html>`.  Or on linux a ``sudo apt-get install doxygen`` should do it.

2. Simply fire up a shell in the root of the repo, and type::

    doxygen Doxyfile


Building the documentation on Linux
-----------------------------------

1. Make sure you have what you need::

    sudo aptitude install build-essential imagemagick git gfortran cmake doxygen ipython python-pip python-virtualenv
    sudo aptitude install libatlas-base-dev libatlas3-base # numpy
    sudo aptitude install gcc gfortran python-dev libblas-dev liblapack-dev libblas3 liblapack3 # scipy
    sudo aptitude install python-dev libpng-dev tk libfreetype6-dev # matplotlib
    sudo aptitude install libxml2-dev libxslt1-dev libxslt1.1 python-all-dev # pandas

2. Follow the instructions above to create the virtual environment.