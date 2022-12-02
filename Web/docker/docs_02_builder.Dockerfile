# Use an intermediate container to build REFPROP
FROM ghcr.io/coolprop/coolprop_docs_01_base as intermediate

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
FROM ghcr.io/coolprop/coolprop_docs_01_base
# Use the output of the earlier build
COPY --from=intermediate /opt/refprop /opt/refprop
