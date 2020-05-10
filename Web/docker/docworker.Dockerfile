FROM continuumio/miniconda3

RUN mkdir /usr/share/man/man1/

RUN apt-get -y -m update && \
    apt-get install -y \
        g++ make cmake swig doxygen p7zip-full \
        mono-mcs \
        octave liboctave-dev \
        r-base-dev \
        default-jre default-jdk \
        texlive-extra-utils \
        imagemagick rsync

        
ADD conda_environment.yml /environment.yml
RUN conda env create -f /environment.yml

# This ADD block forces a build (invalidates the cache) if the git repo contents have changed, otherwise leaves it untouched.
# See https://stackoverflow.com/a/39278224
ADD https://api.github.com/repos/usnistgov/REFPROP-cmake/git/refs/heads/master RPcmake-version.json
RUN git clone --recursive https://github.com/usnistgov/REFPROP-cmake /REFPROP-cmake

ADD REFPROP_sources /REFPROP_sources

WORKDIR /REFPROP-cmake
SHELL ["/bin/bash", "-c"] # https://github.com/moby/moby/issues/7281#issuecomment-389440503
RUN source activate docs && \
    python -c "import numpy; print(numpy.__file__)" && \
    mkdir build && \
    cd build && \
    cmake .. -DREFPROP_FORTRAN_PATH=/REFPROP_sources/FORTRAN && \
    cmake --build . && \
    mkdir -p /opt/refprop && \
    cp librefprop.so /opt/refprop && \
    cp -r /REFPROP_sources/FLUIDS /opt/refprop && \
    cp -r /REFPROP_sources/MIXTURES /opt/refprop 

RUN python -m pip install pybtex
