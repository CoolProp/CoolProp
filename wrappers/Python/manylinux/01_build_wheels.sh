#!/bin/bash

SETUP_PY_ARGS="$1"

# Stop on errors
set -ex 

# Get the directory containing this script
# see http://stackoverflow.com/a/246128/1360263
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# # Build the .tar.gz with the sources
# pushd ${DIR}/../pypi
# source /py27/bin/activate
# python prepare_pypi.py --dist-dir=`pwd`/../dist
# deactivate
# popd 

pushd ${DIR}/..

for PYBIN in /py*; do
    source ${PYBIN}/bin/activate
    c++ --version
    python setup.py bdist_wheel ${SETUP_PY_ARGS}
    deactivate
done

# Bundle external shared libraries into the wheels
for whl in dist/*.whl; do
    # auditwheel comes in base image
    auditwheel repair $whl -w ../../install_root/Python
done

popd 

exit 0