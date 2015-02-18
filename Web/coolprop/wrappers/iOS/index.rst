.. _ios:

************
iOS Tutorial
************

This tutorial is based on a tutorial put together by `Babak Samareh <mailto:babak.samareh@gmail.com>`_.


1. Pull the latest version of CoolProp from github::

    git clone https://github.com/CoolProp/CoolProp --recursive
    
2. Generate the headers::

    cd CoolProp/dev && python generate_headers.py
    
3. Open XCode and create a new project, for the template, under iOS, choose Cocoa Touch Static Library and save the project somewhere.  Now there are two files generated with your project name, one .h and one .m. Get rid of them.

4. Right click on your project name and select "add files". Go to the folder where you have compiled CoolProp and add /src and /include folders. You should end up with something like this:

    .. image:: xcode.png
        :height: 200px

    Don't forget to remove the .i and .cxx files, otherwise you will get errors.
    
5. Set up your project.

    Under PROJECT -> Build Settings change the following::
    
        Under Architectures: Build Active Architecture Only : No
        Under Deployment: Targeted Device Family: iPhone/iPad
    
    You also have to sign your code, so under Code Signing -> Code Signing identity, select the relevant profile.

    Under TARGETS -> Build Settings::
    
        Under Deployment: Targeted Device Family: iPhone/iPad
        Under Linking: Other Linker Flags: -ObjC -all
        Under Search Paths: Header Search Paths: add the location for /externals and /Include, both of them recursive. 
    
    It should look like this:

    .. image:: xcode2.png
        :height: 200px

6. Change your scheme from Debug to Release. Build the project. 

.. warning:: 

    If you are building while on simulator, you end up with a static library .a file with 32/64 bit compatibility. If you want to build an actual arm static library, you have to compile on an actual device. 

7. Connect your phone, select it from the drop down list and build again. Now under products you should be able to see this:

    .. image:: xcode3.png
        :height: 200px

    That libCoolPropStaticLib.a is the file you are looking for. It contains all the required ARM editions. Right click on it and click reveal in finder.