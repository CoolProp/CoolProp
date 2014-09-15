.. _Java:

************
Java Wrapper
************

Pre-compiled Binaries
=====================
Pre-compiled binaries can be downloaded from :sfdownloads:`Java`, which come from  :bbbinaries:`the continuous development buildbot server <Java>`.

Usage
-----
At the console, run::

    javac *.java
    java Example

User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the Java wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`
    
Linux & OSX
-----------

1. Download a zip file of the Java Development Kit (JDK) for Java from `Oracle downloads <http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html>`_. If you are in a 32-bit system, download the 32-bit system, else download the 64-bit version.

2. Expand the zip file you downloaded

3. Add the ``bin`` folder of the JDK that you installed.  For instance, add this::
      
    export /path/to/java/SDK/bin:$PATH 
      
  to ``~/.profile`` where the path ``/path/to/java/SDK/bin`` points to the absolute path for the ``bin`` folder of your Java installation.

Windows
-------

Download and run the JDK installer from `Oracle downloads <http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html>`_. If you are in a 32-bit system, download the 32-bit system, else download the 64-bit version.

Build
-----

Linux and OSX
^^^^^^^^^^^^^

Once the dependencies are installed, you can run the builder and tests using::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    mkdir -p  CoolProp/build && cd CoolProp/build
    # Build the makefile using CMake
    cmake .. -DCOOLPROP_JAVA_MODULE=ON -DBUILD_TESTING=ON
    # Make the java files
    make install
    # Run the integration tests
    ctest --extra-verbose

Windows (32-bit and 64-bit)
^^^^^^^^^^^^^^^^^^^^^^^^^^^ 

You need to just slightly modify the building procedure::

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    cd CoolProp
    # Make a build folder
    mkdir build && cd build
    # Build the makefile using CMake
    cmake .. -DCOOLPROP_JAVA_MODULE=ON -DBUILD_TESTING=ON
    # Make the Java shared library
    make install
    # Run the integration tests
    ctest --extra-verbose
