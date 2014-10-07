.. _developer_documentation:

*************
Documentation
*************

Build Sphinx documentation
--------------------------

1. Check out the sources in the CoolProp/Web folder::

    git clone https://github.com/CoolProp/CoolProp --recursive

2. Make a virtualenv::

    virtualenv ~/env/py27
    source ~/env/py27/bin/activate # Turn on this virtual env, should see (py27) in your command shell next to the prompt to tell you that environment is active

3. Then install prerequisites into this virtualenv::
  
    pip install --upgrade Cython
    pip install numpy scipy matplotlib
    pip install sphinx ipython sphinxcontrib-bibtex sphinxcontrib-doxylink sphinxcontrib-napoleon cloud-sptheme


4. To build the documentation, go into the CoolProp/Web folder and run::

    make html
    
5. Move the generated docs in ``_build`` to wherever you want

6. In case you try to set up a new buildbot (from http://docs.buildbot.net/current/tutorial/firstrun.html)::

    pip install sqlalchemy==0.7.10 buildbot-slave
    buildslave create-slave your-slave coolprop.dreamhosters.com:port your-slave pass
  
Build Doxygen documentation
---------------------------

All the configuration is done in the ``Doxyfile`` file.

1. If you don't have doxygen, install it from `doxygen downloads <http://www.stack.nl/~dimitri/doxygen/download.html>`.  Or on linux a ``sudo apt-get install doxygen`` should do it.

2. Simply fire up a shell in the root of the repo, and type::

    doxygen Doxyfile
  

Creating a documentation slave 
------------------------------

1. Make sure you have what you need. For Linux::

    sudo aptitude install build-essential gfortran python-matplotlib python-pip python-dev cmake ipython git 
    sudo aptitude install libblas-dev liblapack-dev ngrep doxygen 
    sudo pip install --upgrade Cython
    sudo pip install numpy scipy matplotlib
    sudo pip install --upgrade sphinx ipython sphinxcontrib-bibtex sphinxcontrib-doxylink sphinxcontrib-napoleon cloud-sptheme
    
2. Create a new buildbot slave (from http://docs.buildbot.net/current/tutorial/firstrun.html)::

    sudo pip install sqlalchemy==0.7.10 buildbot-slave
    buildslave create-slave your-slave coolprop.dreamhosters.com:port your-slave password
    buildslave start your-slave
    
3. Go to the homepage http://coolprop.dreamhosters.com:8010 and watch your slave work.
  