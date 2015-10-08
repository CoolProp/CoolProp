# Set the base image to Debian
FROM coolprop:slavebase
# File Author / Maintainer
MAINTAINER {{ author }} ({{ email }})
#
#
USER buildbot
WORKDIR /home/buildbot
#
RUN curl -o miniconda.sh http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && \
    chmod +x miniconda.sh && ./miniconda.sh -b && rm miniconda.sh && \
    echo "export PATH=/home/buildbot/miniconda3/bin:/home/buildbot/miniconda/bin:$PATH" >> /home/buildbot/.bash_profile
#
env PATH /home/buildbot/miniconda3/bin:/home/buildbot/miniconda/bin:$PATH
#
RUN \{% for py in cnd_env_pyt %}
conda create -n {{ py }}{% for pkg in cnd_dev_pkgs %} {{ pkg }}{% endfor %} && \{% endfor %}
conda clean
#
RUN \{% for py in cnd_env %}
source activate {{ py }} && \
pip install{% for pkg in pip_add_pkgs %} {{ pkg }}{% endfor %} && \
source deactivate && \{% endfor %}
conda clean -yilts
