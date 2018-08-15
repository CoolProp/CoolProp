.. _Android:

***************
Android Wrapper
***************

.. contents:: :depth: 2


User-Compiled Binaries
======================

Common Requirements
-------------------
Compilation of the Android wrapper requires a few :ref:`common wrapper pre-requisites <wrapper_common_prereqs>`

Furthermore, you will need to install swig, which can be obtained from http://www.swig.org/download.html for windows users (add the unzipped directory to your PATH variable), or from your package manager (apt, homebrew, etc.) on other platforms.

Windows
-------

* Follow the instructions for Linux & OSX, but the first call to cmake should read::

    cmake .. -DCOOLPROP_ANDROID_MODULE=ON -G "MinGW Makefiles" -DNDK_PATH=c:/Downloads/android-ndk-r10e (change path based on your installation)

* This change (telling CMake to use the MinGW generator with ``-G "MinGW Makefiles"``) is needed because by default on windows it tries to use the most-recent installed version of Visual Studio, which conflicts with the Android SDK.

.. warning::

    As of Aug 2016, the version 12 of the NDK does not compile CoolProp correctly.  You must use 10e for some reason.  See also https://github.com/CoolProp/CoolProp/issues/1178

Linux & OSX
-----------

* Install NDK from here: https://developer.android.com/ndk/downloads/index.html

* To Build
    - Check out the sources for CoolProp::

        git clone https://github.com/CoolProp/CoolProp --recursive

    - Open the CMakeLists.txt file in the CoolProp directory for editing
    - Within the COOLPROP_ANDROID_MODULE modify the command:  ``set(ANDROID_PACKAGE_NAME "CoolProp")`` to reflect your package name
    - Make and move into a build folder::

        mkdir -p CoolProp/build && cd CoolProp/build

    - If a target architecture other than the default (armeabi) is desired, the ndk-build can be modified by editing wrappers/Android/Android.mk.template.  For details see https://developer.android.com/ndk/guides/index.html
    - Construct the makefile using CMake::

        cmake .. -DCOOLPROP_ANDROID_MODULE=ON -DNDK_PATH=~/Downloads/android-ndk-r10e (change path based on your installation)

    - Now actually do the build::

        cmake --build .

* To Incorporate into an Android Studio Project
    - Copy the java files from the package directory (i.e build/com/example/myprogram) to the package directory of the android project
    - Copy the the build/libs/armeabi folder containing the .so file to the app/src/main/jniLibs folder of the Android Project (you may have create the jniLibs folder).
    - In the main activity of the Android Project add the code below so the CoolProp functions can be called as described in the java wrapper documentation::

        static {
            System.loadLibrary("CoolProp");
            }


Example Code
======================

* mainactivity.java

.. code:: java

    package com.example.coolprop.androidexample;

    import android.support.v7.app.ActionBarActivity;
    import android.os.Bundle;
    import android.view.Menu;
    import android.view.MenuItem;
    import android.widget.TextView;


    public class MainActivity extends ActionBarActivity {

        static {
            System.loadLibrary("CoolProp");
        }

        @Override
        protected void onCreate(Bundle savedInstanceState) {
            super.onCreate(savedInstanceState);
            setContentView(R.layout.activity_main);

            TextView CoolPropsDisplay = (TextView) findViewById(R.id.CoolPropsDisplay);
            CoolPropsDisplay.setText(Double.toString(CoolProp.PropsSI("T", "P", 101300, "Q", 0, "Water")));
        }

        @Override
        public boolean onCreateOptionsMenu(Menu menu) {
            // Inflate the menu; this adds items to the action bar if it is present.
            getMenuInflater().inflate(R.menu.menu_main, menu);
            return true;
        }

        @Override
        public boolean onOptionsItemSelected(MenuItem item) {
            // Handle action bar item clicks here. The action bar will
            // automatically handle clicks on the Home/Up button, so long
            // as you specify a parent activity in AndroidManifest.xml.
            int id = item.getItemId();

            //noinspection SimplifiableIfStatement
            if (id == R.id.action_settings) {
                return true;
            }

            return super.onOptionsItemSelected(item);
        }
    }

* activity_main.xml

.. code:: xml

    <RelativeLayout xmlns:android="http://schemas.android.com/apk/res/android"
        xmlns:tools="http://schemas.android.com/tools" android:layout_width="match_parent"
        android:layout_height="match_parent" android:paddingLeft="@dimen/activity_horizontal_margin"
        android:paddingRight="@dimen/activity_horizontal_margin"
        android:paddingTop="@dimen/activity_vertical_margin"
        android:paddingBottom="@dimen/activity_vertical_margin" tools:context=".MainActivity">

        <LinearLayout
            android:layout_width="fill_parent"
            android:layout_height="fill_parent"
            android:orientation="vertical">

            <TextView
                android:layout_width="fill_parent"
                android:layout_height="fill_parent"
                android:textSize="40dp"
                android:gravity="center"
                android:id="@+id/CoolPropsDisplay"
                android:layout_weight="1"
                android:textAlignment="gravity" />

        </LinearLayout>

    </RelativeLayout>
