#!/bin/bash

# Stop on errors
set -ex

# Get the directory containing this script
# see http://stackoverflow.com/a/246128/1360263
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Move to repository root (two levels up from wrappers/Python/manylinux)
pushd ${DIR}/../../..

for PYBIN in /py*; do
    source ${PYBIN}/bin/activate
    c++ --version
    # Use pip wheel with scikit-build-core
    python -m pip wheel --no-deps --no-build-isolation -w dist .
    deactivate
done

# Bundle external shared libraries into the wheels
for whl in dist/*.whl; do
    # auditwheel comes in base image
    auditwheel repair $whl -w wrappers/Python/install_root
done

popd

exit 0