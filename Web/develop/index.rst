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
  
Address Sanitizer
-----------------

Address sanitizer is a module of the clang compiler that can help to pinpoint several memory problems, like addressing memory that is out of range.  

The instructions here explain how to get address sanitizer working for CoolProp for testing purposes.

1. You need Clang >3.6, build from source on linux following instructions

2. Check out CoolProp using git::

    git clone https://github.com/CoolProp/CoolProp --recursive

3. ``cd CoolProp``

4. ``mkdir build/asan && cd build/asan``

5. ``sudo ln -s /path/to/your/clang/bin/clang /usr/bin/clang36``

6. ``CC=clang36 CXX=clang36 cmake ../.. -DCOOLPROP_CLANG_ADDRESS_SANITIZER``

7. ``ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.5 ASAN_OPTIONS=verbosity=1 ./CatchTestRunner``

The ``verbosity=1`` is to make sure that ASAN is actually running
