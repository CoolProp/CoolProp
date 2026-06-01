
***********************
Diagnostics and Testing
***********************

Continuous Integration
----------------------

CoolProp is tested on GitHub Actions; the workflow definitions live in
``.github/workflows/``.  Pushes and pull requests trigger the test matrix
automatically.

Address Sanitizer
-----------------

Address sanitizer (ASan) is an instrumentation mode of the clang and gcc
compilers that can help to pinpoint several memory problems, like addressing
memory that is out of range.

Modern clang and gcc ship ASan as part of the standard toolchain, so no special
compiler download is needed — any reasonably recent ``clang``/``clang++`` (or
``gcc``/``g++``) will do.

1. Check out CoolProp using git::

    git clone https://github.com/CoolProp/CoolProp

2. Make a build folder::

    cd CoolProp && mkdir -p build/asan && cd build/asan

3. Run the cmake call, enabling the address sanitizer build::

    cmake ../.. -DCOOLPROP_CLANG_ADDRESS_SANITIZER=ON -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++

4. Build::

    cmake --build .

5. Execute, asking ASan to be verbose so you can confirm it is active::

    ASAN_OPTIONS=verbosity=1 ./CatchTestRunner

The ``verbosity=1`` is to make sure that ASan is actually running.
