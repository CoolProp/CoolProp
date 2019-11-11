FROM ubuntu:18.04

RUN apt-get -y -m update && apt-get install -y cmake python3-numpy g++ gfortran git zip python3-six lcov nano python3-pip wget lsb-core software-properties-common

RUN bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"

ENV PATH "$PATH:/usr/bin"
