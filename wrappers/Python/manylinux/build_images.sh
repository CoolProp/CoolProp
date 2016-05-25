#!/bin/bash

BITNESS="$1"
if [[ "$BITNESS" == "64" ]] 
then
  DOCKER_MACHINE_NAME="coolprop/manylinux"
  DOCKER_MACHINE_ARGS="cmake=default,64"
elif [[ "$BITNESS" == "32" ]] 
then
  DOCKER_MACHINE_NAME="coolprop/manylinux32"
  DOCKER_MACHINE_ARGS="cmake=default,32"
else
  DOCKER_MACHINE_NAME="coolprop/manylinux"
  DOCKER_MACHINE_ARGS="cmake=default,64"
fi 

DOCKER_MACHINE_TAG="v1.4.2"

# Stop on errors
set -ex

# Only start docker-machine if needed
# https://gist.github.com/ladyrassilon/31fa6877a9c24ae4e8f0
if hash docker-machine 2>/dev/null
then
  docker_running=$(docker-machine ls | grep default)
else
  echo "No docker-machine available."
  docker_running=""
fi

if [[ "$docker_running" == *"Stopped"* ]] 
then
  docker-machine start default
  eval "$(docker-machine env default)"
  env | grep "DOCKER_HOST"
elif [[ "$docker_running" == *"Running"* ]]
then
    eval "$(docker-machine env default)"
fi

# Go up three directories
pushd ../../..

# Run the build script
docker run --rm -v `pwd`:/io ${DOCKER_MACHINE_NAME}:${DOCKER_MACHINE_TAG} /io/wrappers/Python/manylinux/build_wheels.sh ${DOCKER_MACHINE_ARGS}