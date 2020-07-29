#!/bin/bash
mkdir /cpbuild || 0
cd /cpbuild
ls /usr/bin/clang*
cmake /coolprop -DCOOLPROP_CLANG_ADDRESS_SANITIZER=ON -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_C_COMPILER=/usr/bin/clang-10 -DCMAKE_CXX_COMPILER=/usr/bin/clang++-10 -DPYTHON_EXECUTABLE=/usr/bin/python3
ASAN_OPTIONS=verbosity=1:atexit=1 cmake --build .
