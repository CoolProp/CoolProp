#!/bin/bash
set -v

# This script is intended to be run INSIDE the docker container, and
# if all goes to plan, the _build/html folder should contain the contents
# of the build

# Turn on our conda environment
source activate docs

# Build/Install CoolProp and check
cd /coolprop/wrappers/Python
python setup.py bdist_wheel --dist-dir dist cmake=default,64
pip install -vvv --force-reinstall --ignore-installed --upgrade --no-index `ls dist/CoolProp*.whl`
rm -rf dist
cd /coolprop
python -c "import CoolProp; print(CoolProp.__gitrevision__)"
python -c "import CoolProp; print(CoolProp.__file__)"

# Run the slow stuff
cd /coolprop/Web/scripts
python -u __init__.py False

# Doxygen
cd /coolprop
doxygen --version && doxygen Doxyfile

# api documentation
cd /coolprop/Web
sphinx-apidoc -T -f -e -o apidoc ../wrappers/Python/CoolProp

# All the rest of the docs
cd /coolprop/Web
make html