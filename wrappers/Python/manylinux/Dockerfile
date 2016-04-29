
FROM quay.io/pypa/manylinux1_x86_64

# Install Cmake
RUN wget --no-check-certificate http://www.cmake.org/files/v2.8/cmake-2.8.10.2.tar.gz && tar -zxvf cmake-2.8.10.2.tar.gz && cd cmake-2.8.10.2 && ./bootstrap && make && make install

RUN /opt/python/cp27-cp27m/bin/pip install virtualenv

# Create virtual environments
RUN /opt/python/cp27-cp27m/bin/virtualenv /py27 --python=/opt/python/cp27-cp27m/bin/python
RUN /opt/python/cp27-cp27m/bin/virtualenv /py27mu --python=/opt/python/cp27-cp27mu/bin/python
RUN /opt/python/cp27-cp27m/bin/virtualenv /py33 --python=/opt/python/cp33-cp33m/bin/python
RUN /opt/python/cp27-cp27m/bin/virtualenv /py34 --python=/opt/python/cp34-cp34m/bin/python
RUN /opt/python/cp27-cp27m/bin/virtualenv /py35 --python=/opt/python/cp35-cp35m/bin/python

# Install pip packages needed
RUN /py27/bin/pip install cython virtualenv auditwheel
RUN /py27mu/bin/pip install cython virtualenv
RUN /py33/bin/pip install cython virtualenv
RUN /py34/bin/pip install cython virtualenv
RUN /py35/bin/pip install cython virtualenv