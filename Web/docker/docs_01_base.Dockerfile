FROM continuumio/miniconda3

RUN apt-get -y -m update && \
    apt-get install -y \
        g++ make cmake swig doxygen p7zip-full \
        mono-mcs \
        octave liboctave-dev \
        r-base-dev \
        default-jre default-jdk \
        texlive-extra-utils pdfjam \
        imagemagick rsync && \
    apt-get autoclean

ADD conda_environment.yml /environment.yml
RUN conda env create -f /environment.yml && conda clean --all --yes
RUN mkdir -p /opt
