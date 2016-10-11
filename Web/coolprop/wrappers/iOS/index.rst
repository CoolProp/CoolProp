.. _ios:

************
iOS Tutorial
************

This tutorial is based on a tutorial put together by `Babak Samareh <mailto:babak.samareh@gmail.com>`_.


1. Pull the latest version of CoolProp from github::

    git clone https://github.com/CoolProp/CoolProp --recursive
    
2. Make a build folder and generate the headers::

    cd CoolProp && mkdir build && cd build
    cmake .. -G Xcode -DCOOLPROP_STATIC_LIBRARY=ON -DCOOLPROP_IOS_TARGET=ON
    
3. Open the generated CoolProp.xcodeproj file in xCode.

4. Select the project root then under TARGETS select the CoolProp (Cocoa Touch Static Library) and open the Build Settings tab. You should end up with something like this:

    .. image:: xcode.png
        :height: 283px

5. Setup your project.

    Under Architectures: Make sure that Build Active Architecture Only is set to NO
    
    Under Deployment: Targeted Device Family: 1,2

6. Change your scheme from Debug to Release. Build the project. 

.. warning:: 

    If you are building while on simulator, you end up with a static library .a file with 32/64 bit compatibility. If you want to build an actual arm static library, you have to compile on an actual device. 

7. Connect your phone, select it from the drop down list and build again. Now under products you should be able to see this:

    .. image:: xcode2.png
        :height: 216px

    That libCoolProp.a is the file you are looking for. It contains all the required ARM editions. Right click on it and click reveal in finder.
