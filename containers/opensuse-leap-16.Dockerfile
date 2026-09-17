# openSUSE Leap 16.0.
#
# Leap is built from the same sources as SUSE Linux Enterprise, so it stands in
# this matrix for the SUSE side of the enterprise world in the way that the UBI
# image stands for the Red Hat side.  Its package set is conservative and it
# does not carry IPOPT in the main repository, so this image ordinarily
# exercises the source path: distribution compiler and distribution Eigen,
# IPOPT and MUMPS built by coinbrew.  That combination is worth testing on its
# own account, since it is what a user on any distribution without an IPOPT
# package ends up with.
#
# Only the dependency stanza belongs in this file.  The build, the tests, the
# install and the consumer check are containers/common/build_and_test.sh and
# are identical on every distribution in this matrix.
FROM opensuse/leap:16.0

RUN zypper --non-interactive refresh && \
    zypper --non-interactive install --no-recommends \
        gcc \
        gcc-c++ \
        gcc-fortran \
        make \
        cmake \
        git \
        pkg-config \
        python3 \
        patch \
        wget \
        unzip \
        file \
        which \
        tar \
        gzip \
        eigen3-devel \
        blas-devel \
        lapack-devel \
        gnuplot \
        meson \
        ninja \
        curl \
    && zypper clean --all

# Anything built from source here goes under /usr, and a source build puts its
# pkg-config files and libraries in /usr/lib even though this distribution's own
# libraries live in /usr/lib64.  Neither pkg-config nor the loader looks in both
# by default, so both are named once, here, and not in the shared script.
ENV PKG_CONFIG_PATH=/usr/lib/pkgconfig:/usr/lib64/pkgconfig
ENV LD_LIBRARY_PATH=/usr/lib:/usr/lib64

COPY containers/common/ensure_eigen.sh containers/common/ensure_ipopt.sh /opt/
RUN /opt/ensure_eigen.sh && /opt/ensure_ipopt.sh

COPY . /src
RUN /src/containers/common/build_and_test.sh

CMD ["/bin/bash"]
