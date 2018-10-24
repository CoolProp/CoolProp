#!/bin/bash
mkdir -p /src/Javascript
cd /src/Javascript
cmake .. -DCOOLPROP_JAVASCRIPT_MODULE=ON -DCMAKE_TOOLCHAIN_FILE=${EMSCRIPTEN}/cmake/Modules/Platform/Emscripten.cmake
cmake --build . --target install
exit 0
