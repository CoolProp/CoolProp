#!/bin/bash
set -v
set -e 

# This script is intended to be run INSIDE the docker container, and
# if all goes to plan, the _build/html folder should contain the contents
# of the build

# Turn on our conda environment
source activate docs

# Try to install dependencies on the fly, or rely on the existing environment
#conda install six numpy cython matplotlib requests jinja2 pyyaml

# Build/Install CoolProp and check
cd /coolprop/wrappers/Python
python setup.py bdist_wheel --dist-dir dist cmake=default,64
pip install -vvv --force-reinstall --ignore-installed --upgrade --no-index `ls dist/CoolProp*.whl`
rm -rf dist
cd /coolprop
python -c "import CoolProp; print(CoolProp.__gitrevision__)"
python -c "import CoolProp; print(CoolProp.__file__)"

# Run the slow stuff, if needed, or demanded
cd /coolprop/Web/scripts
python -u __init__.py $1

# Doxygen
cd /coolprop
doxygen --version && doxygen Doxyfile

# api documentation
cd /coolprop/Web
sphinx-apidoc -T -f -e -o apidoc ../wrappers/Python/CoolProp

# All the rest of the docs
cd /coolprop/Web
make html