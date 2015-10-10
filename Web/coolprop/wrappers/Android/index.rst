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

* To Build 

    # Check out the sources for CoolProp
    git clone https://github.com/CoolProp/CoolProp --recursive
    # Open the CMakeLists.txt file in the CoolProp directory for editing
    # Within the COOLPROP_ANDROID_MODULE modify the command:  set(ANDROID_PACKAGE_NAME "CoolProp") to reflect your package name
    # Make and move into a build folder
    mkdir -p CoolProp/build && cd CoolProp/build
    # If a target architecture other than the default (armeabi) is desired, the ndk-build can be modified by editing wrappers/Android/Android.mk.template.  For details see https://developer.android.com/ndk/guides/index.html
    # Construct the makefile using CMake
    cmake .. -DCOOLPROP_ANDROID_MODULE=ON -DNDK_PATH=~/Downloads/android-ndk-r10e (change path based on your installation)
    # Now actually do the build
    cmake --build .


* To Incorporate into an Android Studio Project
    # Copy the java files from the package directory (i.e build/com/example/myprogram) to the pacakge directory of the android proejct
    # Copy the the build/libs/armeabi folder containing the .so file to the app/src/main/jniLibs folder of the Android Project (you may have reate the jniLibs folder).
    # In the main activity of the Android Project add:

		static {
			System.loadLibrary("CoolProp");
		}

	# Now CoolProp functions can be called as described in the java wrapper documentation.
