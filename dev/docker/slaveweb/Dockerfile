# Set the base image
FROM coolprop/slavebase

# File Author / Maintainer
MAINTAINER Jorrit Wronski (jowr@ipu.dk)

# Just to make sure
USER root
WORKDIR /

# Install system tools for compilation of Python extensions numpy, scipy and pandas
RUN apt-get update -qq && \
    apt-get install -y build-essential git gcc gfortran cmake python-dev \
        libhdf5-serial-dev libnetcdf-dev liblapack-dev libatlas-dev libblas-dev \
        libxml2-dev libxslt1-dev libxslt1.1 doxygen libtool \
        libatlas3-base libblas3 liblapack3 libatlas-dev libblas-dev liblapack-dev \
        libpng-dev tk libfreetype6-dev && \
    apt-get clean

# Run as user
USER buildbot
WORKDIR /home/buildbot

# Install miniconda Python distribution
RUN curl -o miniconda.sh http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && \
    chmod +x miniconda.sh && ./miniconda.sh -b && rm miniconda.sh && \
    echo "export PATH=/home/buildbot/miniconda3/bin:/home/buildbot/miniconda/bin:$PATH" >> /home/buildbot/.bash_profile
env PATH /home/buildbot/miniconda3/bin:/home/buildbot/miniconda/bin:$PATH

# Add the additional packages to the conda environment
RUN \
conda create -n CoolPropWeb python=2.7 setuptools cython numpy scipy matplotlib pandas ipython wxpython && \
conda clean -yilts

# Add some pip packages that only exist there
RUN \
source activate CoolPropWeb && \
pip install sphinx cloud-sptheme sphinxcontrib-bibtex sphinxcontrib-doxylink sphinxcontrib-napoleon && \
source deactivate 
