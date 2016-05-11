#!/bin/bash

# Stop on errors
set -ex

# Fire up docker
docker-machine restart default

# Set up the docker environment
eval "$(docker-machine env default)"

# Build the docker image with the virtual environments, cython, etc.
docker build -t coolprop/manylinux .
# Go up three directories
pushd ../../..
# Run the build script
docker run --rm -v `pwd`:/io coolprop/manylinux /io/wrappers/Python/manylinux/build_wheels.sh