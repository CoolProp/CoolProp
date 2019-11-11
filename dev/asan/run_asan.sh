#!/bin/bash
ls /
mkdir /cpbuild && cd /cpbuild
cmake /coolprop -DCOOLPROP_CLANG_ADDRESS_SANITIZER=ON -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_C_COMPILER=/usr/bin/clang-9 -DCMAKE_CXX_COMPILER=/usr/bin/clang++-9
#DYLD_LIBRARY_PATH=/usr/lib/clang/3.9.0/lib/darwin/ ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer  ASAN_OPTIONS=verbosity=1 cmake --build .
LSAN_OPTIONS=verbosity=1:log_threads=1 cmake --build .
