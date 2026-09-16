#!/usr/bin/env bash
#
# Make Eigen available, whichever way this distribution allows.
#
# The same argument as ensure_ipopt.sh, for the same reason.  Most of these
# distributions package Eigen and the Dockerfile has already installed it; a few
# do not reach it from the repositories a container image is given, Red Hat's
# UBI being the case this was written for.  A job that dies because one header
# collection is not in a narrow repository set reports "PSOPT fails on RHEL",
# which is not what happened.
#
# find_package(Eigen3 NO_MODULE) needs Eigen3Config.cmake, not merely the
# headers, so the source install is done through Eigen's own CMake install
# target and not by copying a directory.
#
set -euo pipefail

EIGEN_VERSION=${EIGEN_VERSION:-3.4.0}
PREFIX=${PREFIX:-/usr}

# Ask the question exactly as PSOPT's CMakeLists.txt asks it, by configuring a
# project that asks it.  Testing for eigen3.pc instead would be weaker -- a
# distribution can ship the pkg-config file and the CMake package separately,
# and it is the CMake package that decides whether PSOPT configures -- and
# `cmake --find-package` is a deprecated mode that is not present everywhere.
have_eigen() {
    local probe=/tmp/eigen_probe
    rm -rf "${probe}"
    mkdir -p "${probe}"
    cat > "${probe}/CMakeLists.txt" <<'EOF'
cmake_minimum_required(VERSION 3.16)
project(eigen_probe NONE)
find_package(Eigen3 REQUIRED NO_MODULE)
EOF
    cmake -S "${probe}" -B "${probe}/build" >/dev/null 2>&1
}

if have_eigen; then
    echo "EIGEN-SOURCE: distribution package"
    echo "EIGEN-VERSION: $(pkg-config --modversion eigen3 2>/dev/null || echo unknown)"
    exit 0
fi

echo "EIGEN-SOURCE: built from source (no usable distribution package found)"
git clone --depth 1 --branch "${EIGEN_VERSION}" \
    https://gitlab.com/libeigen/eigen.git /tmp/eigen
cmake -S /tmp/eigen -B /tmp/eigen/build -DCMAKE_INSTALL_PREFIX="${PREFIX}"
cmake --install /tmp/eigen/build
rm -rf /tmp/eigen

if ! have_eigen; then
    echo "Eigen was installed but find_package(Eigen3) still does not see it."
    find "${PREFIX}" -name 'Eigen3Config.cmake' 2>/dev/null || true
    exit 1
fi
echo "EIGEN-VERSION: ${EIGEN_VERSION}"
