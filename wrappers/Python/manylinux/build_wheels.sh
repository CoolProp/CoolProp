#!/bin/bash

# Get the directory containing this script
# see http://stackoverflow.com/a/246128/1360263
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Into the pypi directory
cd ${DIR}/../pypi
# Build the .tar.gz with the sources
source /py27/bin/activate
python prepare_pypi.py --dist-dir=`pwd`/../dist
deactivate

for PYBIN in /py*; do
    source ${PYBIN}/bin/activate
    pip wheel `ls ../dist/*.tar.gz`
    deactivate
done