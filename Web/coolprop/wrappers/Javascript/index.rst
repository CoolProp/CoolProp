.. _Javascript:

******************
Javascript Wrapper
******************


Pre-Compiled Binaries
=====================

* Download the precompiled binaries from :sfdownloads:`Javascript`

* Load your js file into your website, following the structure of `the example here <https://github.com/CoolProp/CoolProp/blob/master/wrappers/Javascript/index.html>`_, which is also included at the above download link

User-Compiled Binaries
======================

Linux
-----
We are following the instructions from `emscripten.org <http://kripken.github.io/emscripten-site/docs/getting_started/downloads.html>`_ - download the portable emscripten SDK `emsdk` for linux.

1. First download node.js, clang++ and llvm using::
    
    sudo apt-get install nodejs clang++ llvm
    
2. Expand the SDK zip file linked above

3. At the console in the folder that contains the file emsdk run the commands::

    emsdk update # This will fetch the list of things to download
    
    emsdk install latest # This will download and install the full toolchain

4. Go enjoy a nice walk or a cup of coffee - it will be a while

5. Activate the SDK just compiled::

    emsdk activate latest # This will make the file ~/.emscripten with the paths to most of the binaries compiled in SDK
    
6. Modify the file ``~/.emscripten`` to make NODE_JS path equal to ``which nodejs`` if it doesn't already

7. Make an environmental variable (in ~/.profile) ``export EMSCRIPTEN=/path/to/emsdk`` that points to the folder that contains ``emc++``, ``emcc``, etc.

8. Check out coolprop::

    git clone https://github.com/CoolProp/CoolProp --recursive
    
9. Folder creating::

    mkdir -p build && cd build
    
10. Build the Javascript module::

     cmake .. -DCOOLPROP_JAVASCRIPT_MODULE=ON -DCMAKE_TOOLCHAIN_FILE=${EMSCRIPTEN}/cmake/Platform/Emscripten.cmake

Windows
-------
1. Download the `EMSDK installer <http://kripken.github.io/emscripten-site/docs/getting_started/downloads.html>`_, run the web download installer, that will install everything, and get you ready.

2. In the ``%HOME%/.emscripten`` file, make sure that there is only one entry for NODE_JS and it points to the right place.  Mine looks like::

    import os
    SPIDERMONKEY_ENGINE = ''
    LLVM_ROOT='C:/Program Files/Emscripten/clang/e1.21.0_64bit'
    NODE_JS='C:/Program Files/Emscripten/node/0.10.17_64bit/node.exe'
    PYTHON='C:/Program Files/Emscripten/python/2.7.5.3_64bit/python.exe'
    EMSDK_GIT='C:/Program Files/Emscripten/git/1.8.3'
    EMSCRIPTEN_ROOT='C:/Program Files/Emscripten/emscripten/1.21.0'
    V8_ENGINE = ''
    TEMP_DIR = 'c:/users/belli/appdata/local/temp'
    COMPILER_ENGINE = NODE_JS
    JS_ENGINES = [NODE_JS]

3. Open an Emscripten command Prompt (Start->Emscripten->Emscripten command Prompt)

4. Make an environmental variable ``EMSCRIPTEN`` that points to the folder that contains ``emc++``, ``emcc``, etc.

5. Navigate to the root of the source

6. Build the build folder::

    mkdir build && cd build

7. Build the Javascript module::

    cmake ../.. -DCOOLPROP_JAVASCRIPT_MODULE=ON -DCMAKE_TOOLCHAIN_FILE=%EMSCRIPTEN%/cmake/Platform/Emscripten.cmake
