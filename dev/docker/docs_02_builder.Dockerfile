# This is a dockerfile for building the docker container that can 
# create the documentation. It requires REFPROP and cannot be made 
# publicly available. 
# 
# Normally, a CI workflow should take care of executing the commands 
# to build the docker image. However, you can also use an access 
# token to manually build the new image and push it to github.
# 
# $ copy the REFPROP sources to the directory of this file.
# $ cat your_token | docker login ghcr.io -u USERNAME --password-stdin
# $ docker build --file docs_02_builder.Dockerfile --tag ghcr.io/coolprop/coolprop_docs_02_builder:dev .
# $ docker push ghcr.io/coolprop/coolprop_docs_02_builder:dev


# Use an intermediate container to build REFPROP
FROM ghcr.io/coolprop/coolprop_docs_01_base:dev as intermediate

# This ADD block forces a build (invalidates the cache) if the git repo contents have changed, otherwise leaves it untouched.
# See https://stackoverflow.com/a/39278224
ADD https://api.github.com/repos/usnistgov/REFPROP-cmake/git/refs/heads/master RPcmake-version.json
RUN git clone --recursive https://github.com/usnistgov/REFPROP-cmake /REFPROP-cmake

# Add the REFPROP source code to the repository, manage the build context accordingly 
ADD REFPROP_sources /REFPROP_sources

# Build the sources using the Fortran compiler
SHELL ["/bin/bash", "-c"] # https://github.com/moby/moby/issues/7281#issuecomment-389440503
RUN source activate docs && \
    python -c "import numpy; print(numpy.__file__)" && \
    cmake -B /REFPROP-build -S /REFPROP-cmake -DREFPROP_FORTRAN_PATH=/REFPROP_sources/FORTRAN && \
    cmake --build /REFPROP-build
	
# Install the REFPROP files
SHELL ["/bin/bash", "-c"] # https://github.com/moby/moby/issues/7281#issuecomment-389440503
RUN mkdir -p /opt/refprop && \
    cp /REFPROP-build/librefprop.so /opt/refprop && \
    cp -r /REFPROP_sources/FLUIDS /opt/refprop && \
    cp -r /REFPROP_sources/MIXTURES /opt/refprop
	
# Delete the sources to avoid distributing them
RUN rm -rf /REFPROP_sources


# Start with the second stage image
FROM ghcr.io/coolprop/coolprop_docs_01_base:dev
# Use the output of the earlier build
COPY --from=intermediate /opt/refprop /opt/refprop
