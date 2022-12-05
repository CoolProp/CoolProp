# This is a dockerfile for building the docker container that is
# the foundation for the documentation builder. All components of
# this image are publicly available and it could be built using
# github actions or some other CI tool. However, it does not change
# frequently and building the image is expensive which is why it
# is not part the frequent CI runs.
# 
# Normally, a CI workflow should take care of executing the commands 
# to build the docker image. However, you can also use an access 
# token to manually build the new image and push it to github.
# 
# $ cat your_token | docker login ghcr.io -u USERNAME --password-stdin
# $ docker build --file docs_01_base.Dockerfile --tag ghcr.io/coolprop/coolprop_docs_01_base:dev .
# $ docker push ghcr.io/coolprop/coolprop_docs_01_base:dev

FROM continuumio/miniconda3

RUN apt-get -y -m update && \
    apt-get install -y \
        g++ make cmake swig doxygen p7zip-full \
        mono-mcs \
        octave liboctave-dev \
        r-base-dev \
        default-jre default-jdk \
        texlive-extra-utils \
        imagemagick rsync && \
    apt-get autoclean
	
# Allow ImageMagick to invoke Ghostscript
RUN sed -i '/disable ghostscript format types/,+6d' /etc/ImageMagick-6/policy.xml

ADD conda_environment.yml /environment.yml
RUN conda update -n base -c defaults conda && conda env create -f /environment.yml && conda clean --all --yes
RUN mkdir -p /opt
