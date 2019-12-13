#!/bin/bash
mkdir /cpbuild || 0
cd /cpbuild
cmake /coolprop -DCOOLPROP_CLANG_ADDRESS_SANITIZER=ON -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_C_COMPILER=/usr/bin/clang-9 -DCMAKE_CXX_COMPILER=/usr/bin/clang++-9
ASAN_OPTIONS=verbosity=1 cmake --build .