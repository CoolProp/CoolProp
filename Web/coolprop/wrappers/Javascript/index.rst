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

On linux, but binaries generated are cross-platform.  We are following the instructions from `emscripten.org <http://kripken.github.io/emscripten-site/docs/getting_started/downloads.html>`_ - download the portable emscripten SDK `emsdk` for linux.

1. First download node.js, clang++ and llvm using 
    
    sudo apt-get install nodejs clang++ llvm
    
2. Expand the SDK zip file linked above

3. At the console in the folder that contains the file emsdk run the commands 

    emsdk update # This will fetch the list of things to download
    
    emsdk install latest # This will download and install the full toolchain

4. Go enjoy a nice walk or a cup of coffee - it will be a while

5. Activate the SDK just compiled

    emsdk activate latest # This will make the file ~/.emscripten with the paths to most of the binaries compiled in SDK
    
6. Modify the file ``~/.emscripten`` to make NODE_JS path equal to `which nodejs`

7. Check out coolprop

    git clone https://github.com/CoolProp/CoolProp
    
8. Folder creating

    mkdir -p build/JS && cd build/JS
    
9. cmake ../.. 

