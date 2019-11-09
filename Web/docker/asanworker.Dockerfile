FROM ubuntu:18.04

RUN apt-get -y -m update && apt-get install -y cmake python3-numpy g++ gfortran git zip python3-six lcov nano python3-pip wget lsb-core software-properties-common

RUN bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"

# This ADD block forces a build (invalidates the cache) if the git repo contents have changed, otherwise leaves it untouched.
# See https://stackoverflow.com/a/39278224
ADD https://api.github.com/repos/CoolProp/CoolProp/git/refs/heads/master CoolProp.json
RUN git clone --recursive https://github.com/CoolProp/CoolProp /CoolProp

WORKDIR /CoolProp
ENV PATH "$PATH:/usr/bin"
CMD mkdir build && cd build && cmake .. -DCOOLPROP_CLANG_ADDRESS_SANITIZER=ON -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_C_COMPILER=/usr/bin/clang-9 -DCMAKE_CXX_COMPILER=${CLANG_ROOT}/bin/clang++-9
