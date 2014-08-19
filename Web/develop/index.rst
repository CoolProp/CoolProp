
**************************
Information for Developers
**************************

.. toctree::

  cmake.rst
  buildbot.rst
  documentation.rst
  
Address Sanitizer
-----------------
1. Clang >3.6, build from source on linux following instructions
2. Check out CoolProp using git:

       git clone https://github.com/CoolProp/CoolProp --recursive

3. cd CoolProp
4. mkdir build/asan && cd build/asan
5. sudo ln -s /path/to/your/clang/bin/clang /usr/bin/clang36
5. CC=clang36 CXX=clang36 cmake ../.. -DCOOLPROP_CLANG_ADDRESS_SANITIZER
6. ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.5 ASAN_OPTIONS=verbosity=1 ./CatchTestRunner

The verbosity=1 is to make sure that ASAN is actually running
