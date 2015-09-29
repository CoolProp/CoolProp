.. _Android:

************
Android Wrapper
************

.. contents:: :depth: 2


User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the Android wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Linux & OSX
-----------

* Install NDK from here: https://developer.android.com/ndk/downloads/index.html

* Do the build with something like:

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Move into the folder you just created
    mkdir -p CoolProp/build && cd CoolProp/build
    # Construct the makefile using CMake
    cmake .. -DCOOLPROP_ANDROID_WRAPPER=ON -DNDK_PATH=~/Downloads/android-ndk-r10e
    # Now actually do the build
    cmake --build .