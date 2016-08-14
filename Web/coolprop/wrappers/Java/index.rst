.. _Java:

************
Java Wrapper
************

.. contents:: :depth: 2

Pre-compiled Binaries
=====================

Pre-compiled binaries can be downloaded from :sfdownloads:`Java`.  Development binaries coming from the buildbot server can be found at :sfnightly:`Java`.

Download the ``platform-independent.7z`` file and expand it to a folder called ``platform-independent`` using 7-zip.  Download the special Java shared library for your system architecture to the same location from either :sfdownloads:`Java` (release) or :sfnightly:`Java` (development).  Copy the Example.java file to the same location.  You will need to have a copy of some version of java.

When you are finished, you should have a folder layout something like ::

    main
     |- CoolProp.dll
     |- Example.java
     |- platform-independent
        |- AbstractState.java
        |- Configuration.java
        |- ...

Usage
-----
At the console, run::

    javac *.java
    java Example
    
There is example code :ref:`at the end of this page <java_example>`

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
    cmake .. -DCOOLPROP_JAVA_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
    # Make the java files
    cmake --build .
    # Run the integration tests
    ctest --extra-verbose

If you want to change the package that CoolProp resides in, you can do so by changing the cmake call to read::

    cmake .. -DCOOLPROP_JAVA_MODULE=ON -DBUILD_TESTING=ON -DCOOLPROP_SWIG_OPTIONS="-package package.name"

where ``package.name`` is replaced with the desired name

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
    cmake --build . --config Release (drop the --config Release if you use MinGW)
    # Run the integration tests
    ctest --extra-verbose

If you want to change the package that CoolProp resides in, you can do so by changing the cmake call to read::

    cmake .. -DCOOLPROP_JAVA_MODULE=ON -DBUILD_TESTING=ON -DCOOLPROP_SWIG_OPTIONS="-package package.name"

where ``package.name`` is replaced with the desired name

.. _java_example:

Example Code
============

.. literalinclude:: Example.java
   :language: java

Example Code Output
===================

.. literalinclude:: Example.out