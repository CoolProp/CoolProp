# Set the base image to Debian
FROM debian:stable
# File Author / Maintainer
MAINTAINER {{ author }} ({{ email }})
#
#
RUN apt-get update -qq && \
    apt-get install -y{% for pkg in lin_dev_pkgs %} {{ pkg }}{% endfor %} && \
    apt-get clean
#
RUN pip install{% for pkg in pip_dev_pkgs %} {{ pkg }}{% endfor %}
#
RUN groupadd -r buildbot && useradd -r -g buildbot buildbot && \
    mkdir /buildslave && chown buildbot:buildbot /buildslave
#
USER buildbot
WORKDIR /buildslave
RUN buildslave create-slave . {{ masterhost }} {{ slavename }} {{ slavepassword }}
ENTRYPOINT ["/usr/local/bin/buildslave"]
CMD ["start", "--nodaemon"]
