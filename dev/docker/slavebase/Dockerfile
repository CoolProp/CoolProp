# Set the base image to Debian
FROM debian:stable

# File Author / Maintainer
MAINTAINER Jorrit Wronski (jowr@ipu.dk)

# Just to make sure
USER root
WORKDIR /

# Install system tools for compilation
RUN apt-get update -qq && \
    apt-get install -y build-essential gfortran gcc cmake git && \
    apt-get install -y bash rsync curl wget && \
    apt-get install -y python-dev python-pip && \
    apt-get clean

# Install the buildbot slave software and add the appropriate user
RUN pip install buildbot-slave
RUN groupadd -r buildbot && \
    useradd -r -d /home/buildbot -m -s /bin/bash -g buildbot buildbot # && mkdir /buildslave && chown buildbot:buildbot /buildslave

# Make sure bash is the standard shell
RUN rm /bin/sh && ln -sf /bin/bash /bin/sh
USER buildbot
WORKDIR /home/buildbot

# Add some standard hosts
RUN \
mkdir .ssh && \
ssh-keyscan github.com >> .ssh/known_hosts && \
ssh-keyscan bitbucket.org >> .ssh/known_hosts && \
ssh-keyscan coolprop.org >> .ssh/known_hosts && \
ssh-keyscan coolprop.dreamhosters.com >> .ssh/known_hosts && \
cat .ssh/known_hosts

# Add the entrypoint that creates a repository on first launch
COPY ./Dockerfile.slave.entrypoint.sh /
ENTRYPOINT ["/Dockerfile.slave.entrypoint.sh"]
