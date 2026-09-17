# Manjaro.
#
# Manjaro is Arch with its package updates held back for a few weeks of testing,
# so it is not a redundant copy of the image above: at any given moment it has
# older versions of exactly the packages Arch has just moved.  That gap is the
# thing being tested.  When the Arch job fails and this one passes, the cause is
# a package that changed in the last fortnight, and the log of the two jobs side
# by side says which one.
#
# Everything said in arch.Dockerfile about not using the AUR applies here too.
#
# Only the dependency stanza belongs in this file.  The build, the tests, the
# install and the consumer check are containers/common/build_and_test.sh and
# are identical on every distribution in this matrix.
FROM manjarolinux/base:latest

# pacman 7 downloads packages as an unprivileged user inside a seccomp and
# landlock sandbox.  A container build cannot always install that filter, and
# pacman then stops before it has synchronised anything:
#
#   error: error restricting syscalls via seccomp: 22!
#   error: switching to sandbox user 'alpm' failed!
#   error: failed to synchronize all databases
#
# 22 is EINVAL: the filter was refused, not violated.  The sandbox limits what
# the download process can do to the machine it runs on, and here that machine
# is a container built for one test run and then discarded, so the alternative
# is not a safer build but no build at all.
#
# Turned off through pacman.conf and not through the --disable-sandbox flag,
# because an older pacman does not have that flag and would stop on it, whereas
# it merely warns about a directive it does not recognise.  The directive has to
# go inside [options]: appending it to the end of the file would place it in the
# last repository section, where it does nothing.
RUN sed -i '/^\[options\]/a DisableSandbox' /etc/pacman.conf

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
        meson \
        ninja \
        curl \
    && pacman -Scc --noconfirm

# Anything built from source here goes under /usr.  Arch keeps its own
# libraries in /usr/lib, so this is a formality on this distribution, but the
# matrix is easier to read when every image says the same thing.
ENV PKG_CONFIG_PATH=/usr/lib/pkgconfig:/usr/lib64/pkgconfig
ENV LD_LIBRARY_PATH=/usr/lib:/usr/lib64

COPY containers/common/check_tools.sh containers/common/ensure_eigen.sh containers/common/ensure_ipopt.sh /opt/
RUN /opt/check_tools.sh && /opt/ensure_eigen.sh && /opt/ensure_ipopt.sh

COPY . /src
RUN /src/containers/common/build_and_test.sh

CMD ["/bin/bash"]
