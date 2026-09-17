# Red Hat Enterprise Linux 10, through the Universal Base Image.
#
# NOT IN THE MATRIX, and kept anyway.
#
# .github/workflows/distros.yml does not build this image.  It was in the matrix
# for the first run and it was the only job that did not pass, and chasing it
# further was judged not to be worth the effort it was taking.  The file stays
# because it is a working starting point for anyone who needs PSOPT on RHEL and
# because throwing it away would discard what writing it established.  It is not
# a supported target and nothing checks it.  By hand, from the top of the
# repository:
#
#     docker build --progress=plain -f containers/rhel-10.Dockerfile .
#
# containers/README.md says more about why it is out.
#
# What makes this image hard is worth recording for whoever picks it up.
# redhat/ubi10 really is RHEL 10: the same glibc, the same compilers, the same
# runtime a user on a supported subscription has.  What it does not have is the
# same set of REPOSITORIES.  The UBI repositories are deliberately narrow, so
# packages an ordinary RHEL installation reaches through CodeReady Builder, and
# everything from EPEL, are absent until they are enabled, and enabling them
# from inside a container is the part that does not reliably work.
#
# Hence the shape below.  The CodeReady Builder and EPEL steps are allowed to
# fail: where they succeed, BLAS, LAPACK and Eigen come from packages and this
# tests PSOPT against RHEL's own libraries; where they do not, the ensure_
# scripts build what is missing from source and this tests PSOPT against a bare
# RHEL, which is what a user inside an air-gapped or subscription-limited estate
# actually has.  The log says which happened.
#
# And the claim to make if it ever does pass: nothing in PSOPT depends on what
# the wider repository set adds.  Not that a subscribed RHEL 10 behaves
# identically, which this cannot show.
#
# Only the dependency stanza belongs in this file.  The build, the tests, the
# install and the consumer check are containers/common/build_and_test.sh and are
# identical on every distribution that uses them.
FROM redhat/ubi10

# The part that must work.  If this fails there is no point continuing: there
# is no compiler.
RUN dnf install -y \
        gcc \
        gcc-c++ \
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
        tar \
        gzip \
        diffutils \
        findutils \
    && dnf clean all

# The part that is allowed to fail.  CodeReady Builder carries the -devel
# packages that are not in the base repository set, and EPEL carries Eigen.
# The repository identifier differs between UBI versions and between UBI and a
# subscribed RHEL, so both spellings are tried and neither is required.
RUN dnf config-manager --set-enabled ubi-10-codeready-builder-rpms || \
    dnf config-manager --set-enabled codeready-builder-for-rhel-10-x86_64-rpms || \
    echo "CodeReady Builder could not be enabled; continuing without it"
RUN dnf install -y \
      https://dl.fedoraproject.org/pub/epel/epel-release-latest-10.noarch.rpm \
      || echo "EPEL could not be installed; continuing without it"

# Wanted from packages if they can be had, and provided from source if they
# cannot.  EPEL carries coin-or-Ipopt-devel under the same name Fedora uses, so
# if the EPEL step above succeeded this image may get IPOPT from a package after
# all; the log will say which happened.  Installed one at a time so that one
# unavailable package does not take the rest of the list with it.
#
# gfortran is the one that matters most, because it is needed if IPOPT has to be
# built from source: MUMPS is Fortran.  gnuplot's absence costs nothing but a
# warning on stderr from the plotting routines.
RUN for p in gcc-gfortran blas-devel lapack-devel eigen3-devel \
             coin-or-Ipopt-devel gnuplot meson ninja-build curl ; do \
        dnf install -y "$p" || echo "optional package $p unavailable; continuing" ; \
    done ; dnf clean all

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
