# Fedora 44.
#
# The fast-moving end of the Red Hat family, and the place where a new compiler
# reaches users first.  Fedora is usually a version or two of GCC ahead of the
# LTS distributions, so a construct that a newer compiler has become stricter
# about fails here before it fails anywhere else.  Fedora also packages IPOPT,
# as coin-or-Ipopt-devel, which makes this the Red Hat side counterpart of the
# Debian images.
#
# Only the dependency stanza belongs in this file.  The build, the tests, the
# install and the consumer check are containers/common/build_and_test.sh and
# are identical on every distribution in this matrix.
FROM fedora:44

# 'which' is here because configure scripts inside the IPOPT source build use
# it and a Fedora container does not have it by default.
RUN dnf install -y \
        gcc \
        gcc-c++ \
        gcc-gfortran \
        make \
        cmake \
        git \
        pkgconf-pkg-config \
        python3 \
        patch \
        wget \
        unzip \
        file \
        which \
        eigen3-devel \
        blas-devel \
        lapack-devel \
        gnuplot \
        meson \
        ninja-build \
        curl \
        coin-or-Ipopt-devel \
        MUMPS-devel \
    && dnf clean all

# Anything built from source here goes under /usr, and a source build puts its
# pkg-config files and libraries in /usr/lib even on distributions whose own
# libraries live in /usr/lib64.  Neither pkg-config nor the loader necessarily
# looks in both, so both are named once, here, and not in the shared script.
# On Fedora this is not hypothetical: pkgconf searches /usr/lib64/pkgconfig and
# a source-built IPOPT lands in /usr/lib/pkgconfig.
ENV PKG_CONFIG_PATH=/usr/lib/pkgconfig:/usr/lib64/pkgconfig
ENV LD_LIBRARY_PATH=/usr/lib:/usr/lib64

COPY containers/common/ensure_eigen.sh containers/common/ensure_ipopt.sh /opt/
RUN /opt/ensure_eigen.sh && /opt/ensure_ipopt.sh

COPY . /src
RUN /src/containers/common/build_and_test.sh

CMD ["/bin/bash"]
