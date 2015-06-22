.. _develop:

**************************
Information for Developers
**************************

.. toctree::

  code.rst
  backends.rst
  cmake.rst
  buildbot.rst
  documentation.rst
  release.rst
  
Address Sanitizer
-----------------

Address sanitizer is a module of the clang compiler that can help to pinpoint several memory problems, like addressing memory that is out of range.  

The instructions here explain how to get address sanitizer working for CoolProp for testing purposes.  

The easiest solution is to use OSX and download the binaries for LLVM+clang from http://llvm.org/releases/download.html.  You will need to expand the file with something like::

    tar -xJf clang+llvm-3.5-x86_64-apple-darwin10.9.tar.xz

1. Check out CoolProp using git::

    git clone https://github.com/CoolProp/CoolProp --recursive

2. Move into folder::

    cd CoolProp && mkdir -p build/asan && cd build/asan
    
3. Set environmental variable to the root of the clang installation::

    export CLANG_ROOT=/Users/Ian/Downloads/clang+llvm-3.5.0-macosx-apple-darwin

4. Run the cmake call with the special clang version with asan support::

    cmake ../.. -DCOOLPROP_CLANG_ADDRESS_SANITIZER=ON -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_C_COMPILER=${CLANG_ROOT}/bin/clang -DCMAKE_CXX_COMPILER=${CLANG_ROOT}/bin/clang++

5. Build:: 
    
    cmake --build .

6. Execute::

    DYLD_LIBRARY_PATH=${CLANG_ROOT}/lib/clang/3.5.0/lib/darwin/ ASAN_SYMBOLIZER_PATH=${CLANG_ROOT}/bin/llvm-symbolizer  ASAN_OPTIONS=verbosity=1 ./CatchTestRunner

The ``verbosity=1`` is to make sure that ASAN is actually running
