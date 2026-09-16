# Ubuntu 26.04 LTS.
#
# The current Ubuntu long term support release, and the distribution most PSOPT
# users are on.  Everything PSOPT needs is packaged here, IPOPT included, so
# this image is also the one that says what a plain "apt-get install" build of
# PSOPT gets: which IPOPT version, which Eigen, which compiler.
#
# Only the dependency stanza belongs in this file.  The build, the tests, the
# install and the consumer check are containers/common/build_and_test.sh and
# are identical on every distribution in this matrix.
FROM ubuntu:26.04

ENV DEBIAN_FRONTEND=noninteractive

# coinor-libipopt-dev is the reason this list is short.  Where a distribution
# does not have it, containers/common/ensure_ipopt.sh builds IPOPT from source
# instead; the packages below are enough for either path, since the source
# build needs a Fortran compiler, LAPACK, patch, wget and unzip.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        build-essential \
        gfortran \
        cmake \
        git \
        pkg-config \
        python3 \
        patch \
        wget \
        unzip \
        file \
        libeigen3-dev \
        libblas-dev \
        liblapack-dev \
        gnuplot-nox \
        coinor-libipopt-dev \
    && rm -rf /var/lib/apt/lists/*

# Anything built from source here goes under /usr, and a source build puts its
# pkg-config files and libraries in /usr/lib even on distributions whose own
# libraries live in /usr/lib64.  Neither pkg-config nor the loader necessarily
# looks in both, so both are named once, here, and not in the shared script.
ENV PKG_CONFIG_PATH=/usr/lib/pkgconfig:/usr/lib64/pkgconfig
ENV LD_LIBRARY_PATH=/usr/lib:/usr/lib64

# Before the source tree is copied, so that editing PSOPT does not rebuild
# IPOPT.  Each script reports where its dependency came from and neither fails
# merely because the distribution does not package it.
COPY containers/common/ensure_eigen.sh containers/common/ensure_ipopt.sh /opt/
RUN /opt/ensure_eigen.sh && /opt/ensure_ipopt.sh

COPY . /src
RUN /src/containers/common/build_and_test.sh

CMD ["/bin/bash"]
