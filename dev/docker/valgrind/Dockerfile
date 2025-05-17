FROM ubuntu:24.04

RUN apt-get -y -m update &&  DEBIAN_FRONTEND=noninteractive apt-get install -y python3 python3-numpy python3-pip cmake g++ git zip valgrind

# This ADD block forces a build (invalidates the cache) if the git repo contents have changed, otherwise leaves it untouched.
# See https://stackoverflow.com/a/39278224
ADD https://api.github.com/repos/coolprop/coolprop/git/refs/heads/master CPcmake-version.json
RUN git clone --recursive --shallow-submodules https://github.com/coolprop/coolprop

WORKDIR /coolprop/dev
RUN python3 generate_headers.py

# Build the catch tests
WORKDIR /coolprop
RUN mkdir build && cd build && ls -l && \
    cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCOOLPROP_CATCH_MODULE=ON && \
    cmake --build . 

COPY run_all.py /coolprop/build

WORKDIR /coolprop/build
CMD python3 -u run_all.py > /output/log.txt 2>&1