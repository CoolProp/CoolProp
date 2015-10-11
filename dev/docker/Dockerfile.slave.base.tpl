# Set the base image to Debian
FROM debian:testing
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
RUN groupadd -r buildbot && \
    useradd -r -d /home/buildbot -m -s /bin/bash -g buildbot buildbot # && mkdir /buildslave && chown buildbot:buildbot /buildslave
#
RUN rm /bin/sh && ln -sf /bin/bash /bin/sh
USER buildbot
WORKDIR /home/buildbot
#
#RUN \{% for host in ssh_hosts %} 
#ssh-keyscan {{ host }} >> .ssh/known_hosts && \{% endfor %}
#cat .ssh/known_hosts
#
#RUN buildslave create-slave ${SLAVEDIR} ${MASTERHOST} ${SLAVENAME} ${SLAVEPASSWORD}
#ENTRYPOINT ["/usr/local/bin/buildslave"]
#CMD ["start", "--nodaemon"]
#
COPY ../Dockerfile.slave.entrypoint.sh /
ENTRYPOINT ["/Dockerfile.slave.entrypoint.sh"]
