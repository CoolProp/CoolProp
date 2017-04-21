
***********************
Diagnostics and Testing
***********************

Travis Builds
-------------


Coverity Scan
-------------

1. Download coverity and expand the zip file

2. Make up a little batch script to build the coverity analysis.  Obviously change the paths to suit your machine.

   .. note:: I had to set the path to the g++ compiler, and it didn't like visual studio, so I used mingw.  On linux, this isn't much of an issue as it just uses g++
    
   Here it is::

       set PATH_TO_COV=C:\Users\Belli\Downloads\cov-analysis-win64-7.6.0\cov-analysis-win64-7.6.0\bin\cov-build
       set PATH=C:\TDM-GCC-64\bin;%PATH%
       cmake ../.. -DCOOLPROP_SHARED_LIBRARY=ON -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
       %PATH_TO_COV% --dir cov-int cmake --build . --config Release
    
3. Zip up the results

4. Upload to coverity scan

.. note::

    In April 2017, the upload to Coverity Scan has been automated and all you have to do is to merge and push new code to the branch coverity_scan, Travis CI will do the rest.


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
