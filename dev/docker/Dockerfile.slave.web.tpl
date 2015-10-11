# Set the base image to Debian
FROM coolprop/slavepython
# File Author / Maintainer
MAINTAINER {{ author }} ({{ email }})
#
RUN apt-get update -qq && \
    apt-get install -y{% for pkg in lin_web_pkgs %} {{ pkg }}{% endfor %} && \
    apt-get clean
#
RUN \{% for py in cnd_env_pyt %}
conda install -n {{ py }}{% for pkg in cnd_web_pkgs %} {{ pkg }}{% endfor %} && \{% endfor %}
conda clean -yilts
#
RUN \{% for py in cnd_env %}
source activate {{ py }} && \
pip install{% for pkg in pip_web_pkgs %} {{ pkg }}{% endfor %} && \
source deactivate && \{% endfor %}
conda clean -yilts
