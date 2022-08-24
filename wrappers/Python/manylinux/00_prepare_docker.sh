#!/bin/bash

BITNESS="$1"
if [[ "$BITNESS" == "64" ]] 
then
  DOCKER_IMG_NAME="coolprop/manylinux"
  SETUP_PY_ARGS="cmake=default,64"
elif [[ "$BITNESS" == "32" ]] 
then
  DOCKER_IMG_NAME="coolprop/manylinux32"
  SETUP_PY_ARGS="cmake=default,32"
else
  echo "Received an unknown argument, aborting."
  exit 1
fi 

DOCKER_MACHINE_TAG="v2.1.0"

# Stop on errors
set -ex

# Only start docker-machine if needed
# https://gist.github.com/ladyrassilon/31fa6877a9c24ae4e8f0

# Check if docker-machine is available as command
if hash docker-machine 2>/dev/null
then
  docker_running=$(docker-machine ls | grep default)
else
  echo "No docker-machine available."
  docker_running=""
fi

# If docker-machine is available, chack if we have 
# to start virtual machine that the containers run in.
if [[ "$docker_running" == *"Stopped"* ]] 
then
  docker-machine start default
  eval "$(docker-machine env default)"
  env | grep "DOCKER_HOST"
elif [[ "$docker_running" == *"Running"* ]]
then
  eval "$(docker-machine env default)"
  env | grep "DOCKER_HOST"
else 
  echo "Direct access to docker required, assuming native installation."
fi

# Get the directory containing this script
# see http://stackoverflow.com/a/246128/1360263
CUR_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Go up three directories
pushd ${CUR_DIR}/../../..

# Run the build script
chmod +x ${CUR_DIR}/01_build_wheels.sh
# Rewmove the manylinux container if it didn't complete properly last time (see also https://stackoverflow.com/a/38225298/1360263)
docker stop manylinux || true && docker rm manylinux || true
# docker run --rm -v `pwd`:/io ${DOCKER_IMG_NAME}:${DOCKER_MACHINE_TAG} /io/wrappers/Python/manylinux/01_build_wheels.sh ${SETUP_PY_ARGS}
docker run -itd --name manylinux ${DOCKER_IMG_NAME}:${DOCKER_MACHINE_TAG} bash
docker cp . manylinux:/io
docker exec manylinux /io/wrappers/Python/manylinux/01_build_wheels.sh ${SETUP_PY_ARGS}
docker cp manylinux:/io/install_root install_root
docker stop manylinux 
docker rm manylinux

popd 

exit 0
