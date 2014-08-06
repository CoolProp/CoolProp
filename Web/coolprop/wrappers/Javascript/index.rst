.. _Javascript:

******************
Javascript Wrapper
******************


Users
=====

Precompiled binaries
--------------------

* Download the precompiled binaries from :sfdownloads:`Javascript`

* Load your js file into your website, following the structure of `the example here <https://github.com/CoolProp/CoolProp/blob/master/wrappers/Javascript/index.html>`_

Developers
==========

On linux, but binaries generated are cross-platform, follow the instructions from `emscripten.org <http://kripken.github.io/emscripten-site/docs/getting_started/downloads.html>`_ - download the portable emscripten SDK `emsdk` for linux.

0. First download node.js, clang++ and llvm using 
    
    sudo apt-get install nodejs clang++ llvm
    
1. Expand the SDK zip file linked above
2. At the console in the folder that contains the file emsdk run the commands 

    emsdk update # This will fetch the list of things to download
    
    emsdk install latest # This will download and install the full toolchain
    
3. Go enjoy a nice walk or a cup of coffee - it could be a while
    
4. Set the environmental variable ``EMSCRIPTEN`` to point to the root of the toolchain.  For instance in your /etc/rc.local file you can add the line

    export EMSCRIPTEN=/path/to/folder/containing/emsdk
    
5. 

