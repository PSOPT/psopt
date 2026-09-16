# Arch Linux.
#
# Arch is in the matrix for two reasons.  It is the distribution the published
# PSOPT continuous integration image is built on, so a failure here breaks the
# ordinary build as well as this one; and it is a rolling release with an
# aggressive toolchain, which makes it the other early warning in this matrix
# alongside Tumbleweed.
#
# IPOPT is deliberately NOT taken from the AUR here, although coin-or-ipopt is
# there and the published CI image uses it.  The AUR build needs an unprivileged
# user, an AUR helper and a chain of makepkg steps whose outcome depends on what
# the AUR looked like that morning, and none of that is a property of PSOPT.
# containers/common/ensure_ipopt.sh builds IPOPT from source instead, which is
# reproducible and which is also what an Arch user who does not want an AUR
# helper does.  The AUR route keeps being exercised by the image published from
# the Dockerfile at the top of this repository.
#
# Only the dependency stanza belongs in this file.  The build, the tests, the
# install and the consumer check are containers/common/build_and_test.sh and
# are identical on every distribution in this matrix.
FROM archlinux:latest

# A rolling distribution has to be brought up to date before anything is
# installed: a partial upgrade on Arch is not a supported state and produces
# library version mismatches that look like PSOPT defects.
RUN pacman -Syu --noconfirm && \
    pacman -S --noconfirm --needed \
        base-devel \
        gcc-fortran \
        cmake \
        git \
        pkgconf \
        python \
        wget \
        unzip \
        eigen \
        blas \
        lapack \
        gnuplot \
    && pacman -Scc --noconfirm

# Anything built from source here goes under /usr.  Arch keeps its own
# libraries in /usr/lib, so this is a formality on this distribution, but the
# matrix is easier to read when every image says the same thing.
ENV PKG_CONFIG_PATH=/usr/lib/pkgconfig:/usr/lib64/pkgconfig
ENV LD_LIBRARY_PATH=/usr/lib:/usr/lib64

COPY containers/common/ensure_eigen.sh containers/common/ensure_ipopt.sh /opt/
RUN /opt/ensure_eigen.sh && /opt/ensure_ipopt.sh

COPY . /src
RUN /src/containers/common/build_and_test.sh

CMD ["/bin/bash"]
