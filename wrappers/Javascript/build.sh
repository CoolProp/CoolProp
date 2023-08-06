#!/bin/bash

mkdir /bldbind
cd /bldbind
cmake /src -DCMAKE_TOOLCHAIN_FILE=/emsdk/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake -DCOOLPROP_JAVASCRIPT_MODULE=ON
cmake --build . --target coolprop

cp coolprop.* /src/wrappers/Javascript

exit 0
