#!/bin/bash

SETUP_PY_ARGS="$1"

# https://github.com/pypa/python-manylinux-demo/blob/master/travis/build-wheels.sh
set -e -x

# Install a system package required by our library
#yum install -y atlas-devel

yum install -y cmake

mkdir /io/wheelhouse_tmp
mkdir /io/wheelhouse

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    #ls -lh "${PYBIN}"
    "${PYBIN}/pip" install cython wheel
    pushd /io/wrappers/Python
    "${PYBIN}/python" setup.py bdist_wheel ${SETUP_PY_ARGS}
    cp dist/*.whl /io/wheelhouse_tmp/
    popd
    deactivate
    #"${PYBIN}/pip" install cython wheel
    #"${PYBIN}/pip" wheel /io/wrappers/Python --wheel-dir /io/wheelhouse_tmp/ --build-options ${SETUP_PY_ARGS}
    #"${PYBIN}/pip" wheel /io/wrappers/Python -w /io/wheelhouse_tmp/
done

# Bundle external shared libraries into the wheels
for whl in /io/wheelhouse_tmp/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done

## Install packages and test
#for PYBIN in /opt/python/*/bin/; do
#    "${PYBIN}/pip" install python-manylinux-demo --no-index -f /io/wheelhouse
#    (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
#done
