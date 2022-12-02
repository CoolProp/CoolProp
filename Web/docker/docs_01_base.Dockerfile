FROM continuumio/miniconda3

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
RUN mkdir -p /opt
