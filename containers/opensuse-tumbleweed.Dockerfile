# openSUSE Tumbleweed.
#
# The rolling release, and the opposite end of the same distribution from Leap:
# the newest GCC, the newest CMake, the newest Eigen, all changing weekly.  It
# is in the matrix as an early warning.  A change in a dependency that will
# reach the LTS distributions in a year or two reaches Tumbleweed within days,
# so a failure here is usually a report about the future rather than about
# openSUSE, and the right response to one is often to fix PSOPT before anybody
# else meets the problem.
#
# Because the image moves, a failure here is also the one most likely to be
# transient.  Read it against the previous run before treating it as a defect.
#
# Only the dependency stanza belongs in this file.  The build, the tests, the
# install and the consumer check are containers/common/build_and_test.sh and
# are identical on every distribution in this matrix.
FROM opensuse/tumbleweed

# findutils is named here and in no other image in this matrix; diffutils here
# and in Fedora.  The openSUSE base images carry neither, and what that cost is
# worth recording.  IPOPT's configure ran without cmp, diff or xargs,
# reported each as "command not found" eleven times and produced a working
# IPOPT anyway, and then build_qp_backends.sh could not find osqp-config.cmake
# -- because find was not there either, and every call to it in the shared
# scripts discards stderr so that an empty prefix reads as an empty answer.  An
# absent find therefore answered "not installed" for every file it was asked
# about, and the install log directly above the failure showed cmake writing
# the very file the check said did not exist.  containers/common/check_tools.sh
# now names the tools before anything uses one.
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
        findutils \
        diffutils \
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

COPY containers/common/check_tools.sh containers/common/ensure_eigen.sh containers/common/ensure_ipopt.sh /opt/
RUN /opt/check_tools.sh && /opt/ensure_eigen.sh && /opt/ensure_ipopt.sh

COPY . /src
RUN /src/containers/common/build_and_test.sh

CMD ["/bin/bash"]
