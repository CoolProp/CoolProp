.. _Javascript:

******************
Javascript Wrapper
******************


Pre-Compiled Binaries
=====================

* Download the precompiled binaries from :sfdownloads:`Javascript`, or the development versions from the buildbot server at :sfnightly:`Javascript`. Remember to get both files, ``coolprop.js`` and ``coolprop.js.mem``.

* You can load your js file into your website, following the structure of `the example here <https://github.com/CoolProp/CoolProp/blob/master/Web/coolprop/wrappers/Javascript/index.html>`_, which is also included at the above download. 

* Alternatively, you can link to our server directly to make sure that you always have the latest version of CoolProp. To do so, include the address ``<script src="http://www.coolprop.sourceforge.net/jscript/coolprop.js"></script>`` in your HTML header instead of the relative path ``<script src="coolprop.js"></script>``.

* A live demo of the Javascript library in action can also be found `online <http://www.coolprop.sourceforge.net/jscript/index.html>`_.

* There is a bug in emscripten causing problems when a ``Release`` build of CMake is used.  Switching to ``RelWithDebInfo`` config seems to solve it, just delete the generated .wast file.

Serving the JS
==============

The .wasm file can cause some problems for hosting emscripten-compiled JS files.  When locally hosting, you need to ensure that the .wasm file extension is served with MIME type ``application/wasm``.  A simple webserver config in Python 3 for testing CoolProp (with proper hosting of the WASM) could read::

    import http.server
    from http.server import BaseHTTPRequestHandler, HTTPServer

    port=8000
    print("Running on port %d" % port)

    http.server.SimpleHTTPRequestHandler.extensions_map['.wasm'] = 'application/wasm'
    httpd = HTTPServer(('localhost', port), http.server.SimpleHTTPRequestHandler)

    print("Mapping \".wasm\" to \"%s\"" %
    http.server.SimpleHTTPRequestHandler.extensions_map['.wasm'])
    httpd.serve_forever()

On Apache, you need to `setup the server appropriately <https://emscripten.org/docs/compiling/WebAssembly.html?highlight=apache#web-server-setup>`_.

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

     cmake .. -DCOOLPROP_JAVASCRIPT_MODULE=ON -DCMAKE_TOOLCHAIN_FILE=${EMSCRIPTEN}/cmake/Modules/Platform/Emscripten.cmake

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

    cmake ../.. -DCOOLPROP_JAVASCRIPT_MODULE=ON -DCMAKE_TOOLCHAIN_FILE=%EMSCRIPTEN%/cmake/Modules/Platform/Emscripten.cmake
