# Set the base image
FROM coolprop/slavebase

# File Author / Maintainer
MAINTAINER Jorrit Wronski (jowr@ipu.dk)

# Just to make sure
USER root
WORKDIR /

# Install system tools for compilation of Python extensions numpy, scipy and pandas
RUN apt-get update -qq && \
    apt-get install -y build-essential gcc python-dev \
        libhdf5-serial-dev libnetcdf-dev liblapack-dev libatlas-dev libblas-dev && \
    apt-get clean

# Run as user
USER buildbot
WORKDIR /home/buildbot

# Install miniconda Python distribution
RUN curl -o miniconda.sh http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && \
    chmod +x miniconda.sh && ./miniconda.sh -b && rm miniconda.sh && \
    echo "export PATH=/home/buildbot/miniconda3/bin:/home/buildbot/miniconda/bin:$PATH" >> /home/buildbot/.bash_profile
env PATH /home/buildbot/miniconda3/bin:/home/buildbot/miniconda/bin:$PATH

# Add the miniconda environments for the Python builds
RUN \
conda create -n CoolProp27 python=2.7 cython pip jinja2 pyyaml pycrypto numpy scipy matplotlib pandas wheel && \
conda create -n CoolProp33 python=3.3 cython pip jinja2 pyyaml pycrypto numpy scipy matplotlib pandas wheel && \
conda create -n CoolProp34 python=3.4 cython pip jinja2 pyyaml pycrypto numpy scipy matplotlib pandas wheel && \
conda clean -yilts
