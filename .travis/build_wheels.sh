#!/bin/bash

SETUP_PY_ARGS="$1"

# https://github.com/pypa/python-manylinux-demo/blob/master/travis/build-wheels.sh
set -e -x

# Get the directory containing this script
# see http://stackoverflow.com/a/246128/1360263
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Install a system package required by our library
#yum install -y atlas-devel

#yum install -y cmake

#if [ "$SETUP_PY_ARGS" = *"32" ]; then
#    CMAKE_URL="https://cmake.org/files/v3.6/cmake-3.6.3-Linux-i386.tar.gz"
#else
#    CMAKE_URL="https://cmake.org/files/v3.7/cmake-3.7.2-Linux-x86_64.tar.gz"
#fi
#mkdir cmake && wget --no-check-certificate --quiet -O - ${CMAKE_URL} | tar --strip-components=1 -xz -C cmake
#export PATH=${DIR}/cmake/bin:${PATH}

mkdir -p /io/wheelhouse_tmp
mkdir -p /io/wheelhouse

OLD_PATH=${PATH}

# Compile wheels
for PYBIN in /opt/python/*/bin; do  
    PYV_MAJOR=`"${PYBIN}/python" -c "import sys;print(list(sys.version_info[:2])[0])";`
    PYV_MINOR=`"${PYBIN}/python" -c "import sys;print(list(sys.version_info[:2])[1])";`
    echo Detected Python ${PYV_MAJOR}.${PYV_MINOR}
    if [ "${PYV_MAJOR}" -le "2" -a "${PYV_MINOR}" -lt "7" ]; then
        continue
    fi
    export PATH="${PYBIN}:$OLD_PATH"
    #ls -lh "${PYBIN}"
    pip install cython wheel
    #"${PYBIN}/pip" install scikit-build cmake
    pushd /io/wrappers/Python
    python setup.py bdist_wheel ${SETUP_PY_ARGS}
    cp dist/*.whl /io/wheelhouse_tmp/
    popd
    #deactivate
    #"${PYBIN}/pip" install cython wheel
    #"${PYBIN}/pip" wheel /io/wrappers/Python --wheel-dir /io/wheelhouse_tmp/ --build-options ${SETUP_PY_ARGS}
    #"${PYBIN}/pip" wheel /io/wrappers/Python -w /io/wheelhouse_tmp/
done

export PATH="$OLD_PATH"

# Bundle external shared libraries into the wheels
for whl in /io/wheelhouse_tmp/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done

## Install packages and test
#for PYBIN in /opt/python/*/bin/; do
#    "${PYBIN}/pip" install python-manylinux-demo --no-index -f /io/wheelhouse
#    (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
#done
