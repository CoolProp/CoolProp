#!/bin/bash
# https://github.com/pypa/python-manylinux-demo/blob/master/travis/build-wheels.sh
set -e -x

# Install a system package required by our library
#yum install -y atlas-devel

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install cython wheel
    "${PYBIN}/pip" wheel /io/wrappers/Python -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done

## Install packages and test
#for PYBIN in /opt/python/*/bin/; do
#    "${PYBIN}/pip" install python-manylinux-demo --no-index -f /io/wheelhouse
#    (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
#done
