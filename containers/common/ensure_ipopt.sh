#!/usr/bin/env bash
#
# Make IPOPT available, whichever way this distribution allows.
#
# The Dockerfile has already tried the distribution's own package and is allowed
# to have failed doing so.  If pkg-config can see IPOPT now, nothing happens and
# the job records which version the distribution supplied.  If it cannot, IPOPT
# is built from source through coinbrew, which also brings ASL and MUMPS.
#
# Written this way on purpose.  Whether a distribution packages IPOPT is one of
# the questions this matrix exists to answer, and a job that merely dies when the
# package is absent answers it badly: it reports "failed" where the truth is
# "not packaged here, and everything works once it is built".  The two outcomes
# need to look different in the log, so they do.
#
set -euo pipefail

IPOPT_VERSION=${IPOPT_VERSION:-releases/3.14.19}
PREFIX=${PREFIX:-/usr}

if pkg-config --exists ipopt; then
    echo "IPOPT-SOURCE: distribution package"
    echo "IPOPT-VERSION: $(pkg-config --modversion ipopt)"
    exit 0
fi

echo "IPOPT-SOURCE: built from source (no usable distribution package found)"

# MUMPS is Fortran, so the source path needs a Fortran compiler that the
# package path does not.  Saying so here costs one line and saves reading a
# configure log to find out that the distribution's Fortran compiler is in a
# repository this image cannot reach.
if ! command -v "${FC:-gfortran}" >/dev/null 2>&1; then
    echo "No Fortran compiler (${FC:-gfortran}) is installed, and IPOPT has to be"
    echo "built from source on this distribution.  MUMPS is Fortran, so the build"
    echo "cannot proceed.  Install the distribution's gfortran package in the"
    echo "Dockerfile; on a narrow repository set it may be in an optional"
    echo "repository that the image has not enabled."
    exit 1
fi

echo "Building IPOPT ${IPOPT_VERSION} with coinbrew. This takes several minutes."

git clone --depth 1 https://github.com/coin-or/coinbrew /tmp/coinbrew
cd /tmp/coinbrew
./coinbrew fetch "Ipopt@${IPOPT_VERSION}" --no-prompt
./coinbrew build Ipopt --prefix="${PREFIX}" --no-prompt --verbosity=2
./coinbrew install Ipopt --no-prompt
cd /
rm -rf /tmp/coinbrew

ldconfig || true
export PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig:${PREFIX}/lib64/pkgconfig:${PKG_CONFIG_PATH:-}"

if ! pkg-config --exists ipopt; then
    echo "IPOPT was built but pkg-config still cannot see it."
    echo "PKG_CONFIG_PATH=${PKG_CONFIG_PATH}"
    find "${PREFIX}" -name 'ipopt.pc' 2>/dev/null || true
    exit 1
fi
echo "IPOPT-VERSION: $(pkg-config --modversion ipopt)"
