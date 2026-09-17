#!/usr/bin/env bash
#
# One build-and-test run of PSOPT, identical on every distribution.
#
# Everything that differs between distributions lives in the Dockerfile that
# calls this script.  Nothing distribution-specific belongs here: if this file
# ever needs an "if Fedora" in it, the difference belongs in the Dockerfile
# instead, or it is a portability defect in PSOPT that should be fixed rather
# than worked around.
#
# The script is deliberately loud.  A distribution job that fails is read by
# somebody who does not have that distribution in front of them, so every step
# announces what it found before it uses it.
#
set -euo pipefail

SRC=${SRC:-/src}
BUILD=${BUILD:-/build}
PREFIX=${PREFIX:-/usr}
JOBS=${JOBS:-$(nproc)}

# Where the QP backends go.  Their own prefix rather than ${PREFIX}, so that they stay
# separable from the distribution's files and so that GNUInstallDirs does not send them
# to the multiarch directory, which it does for a CMake install whose prefix is /usr on
# a Debian derivative.
QP_PREFIX=${QP_PREFIX:-/opt/qp}
export QP_PREFIX

# The build tree's own library directory, both library directories under the
# prefix, and both under the QP prefix.  The examples run after the install, so
# they would find the library anyway on a distribution whose loader searches the
# prefix; naming them all costs nothing and removes one way for a job to fail for
# a reason that has nothing to do with PSOPT.  A source-built CppAD or IPOPT lands
# in ${PREFIX}/lib even where the distribution's own libraries are in lib64, and
# GALAHAD's plugin needs ${QP_PREFIX}/lib at load time because the backend it links
# is a shared library.
export LD_LIBRARY_PATH="${BUILD}/lib:${PREFIX}/lib:${PREFIX}/lib64:${QP_PREFIX}/lib:${QP_PREFIX}/lib64:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig:${PREFIX}/lib64/pkgconfig:${PKG_CONFIG_PATH:-}"

say() { printf '\n\033[1m== %s\033[0m\n' "$*"; }

# ---------------------------------------------------------------- environment
say "Environment"
( . /etc/os-release && echo "distribution : $PRETTY_NAME" ) || true
echo "kernel       : $(uname -r)"
echo "cmake        : $(cmake --version | head -1)"
echo "c++          : $(${CXX:-c++} --version | head -1)"
echo "fortran      : $(${FC:-gfortran} --version 2>/dev/null | head -1 || echo 'not present')"
echo "libc         : $(ldd --version 2>&1 | head -1)"
echo "python       : $(python3 --version 2>&1)"

# IPOPT is the dependency most likely to be missing or to differ in version
# between distributions, and its version is the one most likely to move an
# example's answer.  Report it before anything is built.
say "IPOPT"
# Which IPOPT this is, and where it came from, decides more than anything else
# here: a different interior-point path lands on a different mesh.  ensure_ipopt.sh
# has already run in the Dockerfile and has printed IPOPT-SOURCE and
# IPOPT-VERSION; this repeats the version so it sits beside the compiler in one
# block when somebody reads a failing job.
if pkg-config --exists ipopt; then
    echo "ipopt        : $(pkg-config --modversion ipopt)"
    echo "libs         : $(pkg-config --libs ipopt)"
else
    echo "ipopt        : NOT FOUND by pkg-config, although the image build should"
    echo "               have guaranteed it.  Check PKG_CONFIG_PATH in the Dockerfile."
    exit 1
fi

# Eigen is the other dependency whose provenance is worth recording: PSOPT's
# installed header includes <Eigen/Dense>, so a consumer of the installed
# package needs Eigen too, and its version has moved PSOPT's build before now.
say "Eigen"
if pkg-config --exists eigen3; then
    echo "eigen3       : $(pkg-config --modversion eigen3)"
else
    echo "eigen3       : no eigen3.pc; find_package(Eigen3 NO_MODULE) may still"
    echo "               succeed from Eigen3Config.cmake.  Configure will say."
fi

# ------------------------------------------------------------------ CppAD
# Built from source on every distribution, deliberately.  PSOPT links both the
# CppAD headers and the compiled libcppad_lib, and a distribution package that
# ships only the headers satisfies the compiler and then fails at link time.
# Building it the same way everywhere also keeps this matrix a test of PSOPT's
# portability rather than a test of ten distributions' CppAD packaging.
say "CppAD from source"
git clone --depth 1 https://github.com/coin-or/CppAD.git /tmp/cppad
cmake -S /tmp/cppad -B /tmp/cppad/build -D cppad_prefix="${PREFIX}"
cmake --build /tmp/cppad/build -j"${JOBS}"
cmake --install /tmp/cppad/build
rm -rf /tmp/cppad

# ------------------------------------------------------- QP backends for the SQP
# PSOPT's own SQP solver has no QP of its own: every subproblem goes to a backend, so
# WITH_SQP without one is not a build. Four are installed here through PSOPT's own
# installer scripts, which the matrix therefore tests as well.
#
# This is the expensive part of an image by a wide margin, GALAHAD above all, and it is
# what buys the SQP any coverage at all across distributions. Without it ctest runs two
# tests where it could run three, and PSOPT's own solver -- as opposed to IPOPT -- is
# tested on no distribution whatever.
"${SRC}/containers/common/ensure_qp_backends.sh"

# Register the QP prefix with the dynamic loader, as well as putting it on
# LD_LIBRARY_PATH above. The two are not the same thing and the difference is the point:
# LD_LIBRARY_PATH is inherited and any process in between can drop or replace it, while
# /etc/ld.so.conf.d is a property of the image and applies to every lookup, including the
# transitive ones that a RUNPATH does NOT cover -- DT_RUNPATH resolves only an object's
# own direct dependencies, not its dependencies' dependencies.
#
# Only GALAHAD's plugin needs this. OSQP and Clarabel are static archives and PIQP and
# ProxQP are header-only, so those plugins carry no external dependency at all. GALAHAD
# is linked as a shared library, its plugin was the one that failed to load on the first
# run that built it, and belt and braces costs two lines here.
#
# if/fi and not [ ] && ..., for legibility only. The patch that added this said the
# AND-list form would trip set -e on a prefix with no lib64; that was asserted without
# being tested and it is wrong. Bash exempts the left operand of && from set -e.
mkdir -p /etc/ld.so.conf.d
: > /etc/ld.so.conf.d/psopt-qp.conf
if [ -d "${QP_PREFIX}/lib" ];   then echo "${QP_PREFIX}/lib"   >> /etc/ld.so.conf.d/psopt-qp.conf; fi
if [ -d "${QP_PREFIX}/lib64" ]; then echo "${QP_PREFIX}/lib64" >> /etc/ld.so.conf.d/psopt-qp.conf; fi
ldconfig || true
echo "registered with the loader:"
cat /etc/ld.so.conf.d/psopt-qp.conf

# What was actually built, as cmake arguments. ensure_qp_backends.sh works this out from
# the files that exist and writes it down, so that a backend that failed to build leaves
# the others working instead of failing the configure.
QP_ARGS="$(cat "${QP_PREFIX}/psopt-qp-backends.args")"
say "PSOPT will be configured with: ${QP_ARGS}"

# ------------------------------------------------------------------ configure
say "Configure"
cmake -S "${SRC}" -B "${BUILD}" \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
      -DBUILD_EXAMPLES=ON \
      -DBUILD_TESTS=ON \
      -DHEADLESS=ON \
      -DPSOPT_AD_BACKEND=CPPAD \
      -DCMAKE_PREFIX_PATH="${PREFIX};${QP_PREFIX}" \
      ${QP_ARGS}

say "Build"
cmake --build "${BUILD}" -j"${JOBS}"

# --------------------------------------------------- the QP plugins that were built
# Loading a backend plugin sits between a successful build and a working SQP, and it is
# invisible in a build log until it fails. Three of the four plugins link their backend
# statically or are header-only, so they carry no dependency and cannot fail this way.
# GALAHAD's is built against a SHARED libgalahad_double and is the only one that can be
# present and unloadable, which is exactly what happened on the first run that built it.
# Listing them with their dependencies costs a second and is the difference between
# reading such a failure and guessing at it.
say "QP backend plugins"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
ls -l "${BUILD}/qp_plugins/" 2>/dev/null || echo "no qp_plugins directory: no backend was enabled"
for so in "${BUILD}"/qp_plugins/*.so; do
    [ -e "${so}" ] || continue
    echo "--- ${so}"
    ldd "${so}" 2>&1 | sed 's/^/    /' || true
done

# ------------------------------------------------------------------ unit tests
say "Unit tests"
ctest --test-dir "${BUILD}" --output-on-failure

# ------------------------------------------------------------------ install
say "Install"
cmake --install "${BUILD}"
ldconfig || true

# A consumer build is the only thing that tests the install rules, and those
# are what a distribution packager and a first-time user meet first.  Nothing
# else in the test suite reaches them: every other target builds against the
# source tree.
say "Consumer build against the installed package"
mkdir -p /tmp/consumer
cat > /tmp/consumer/CMakeLists.txt <<'EOF'
cmake_minimum_required(VERSION 3.16)
project(psopt_consumer CXX)
find_package(PSOPT REQUIRED)
add_executable(consumer consumer.cxx)
target_link_libraries(consumer PRIVATE PSOPT)
EOF
cat > /tmp/consumer/consumer.cxx <<'EOF'
// The smallest program that proves the installed headers and library are
// usable together: it constructs the library's own types and links.
#include "psopt.h"
int main()
{
    Alg  algorithm;
    Sol  solution;
    Prob problem;
    problem.name        = "installed-package smoke test";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);
    return 0;
}
EOF
cmake -S /tmp/consumer -B /tmp/consumer/build \
      -DCMAKE_PREFIX_PATH="${PREFIX}" -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/consumer/build
/tmp/consumer/build/consumer
echo "consumer built and ran against the installed package"

# ------------------------------------------------------------------ examples
# Only examples whose reference value is fixed by arithmetic or by an external
# algorithm, never by a previous PSOPT run.  Distributions ship different IPOPT
# and MUMPS versions, a different interior-point path lands on a different
# mesh, and a reference re-centred on one machine's output would make this
# matrix report the distribution's linear algebra as a defect in PSOPT.  The
# full thirteen-example regression stays in ci.yml, on one reference platform.
#
#   mineng_di     minimum-energy double integrator.  J* = 6 exactly, and the
#                 transcription is exact on it: cubic states, linear control,
#                 quadratic integrand, all inside what a Legendre interpolant
#                 on 40 nodes represents and Gauss-Lobatto quadrature
#                 integrates without error.  Arithmetic, not a measurement.
#   lts_costates  linear tangent steering.  tf = 419.58841 s follows from the
#                 analytical solution of Bryson and Ho, confirmed to nine
#                 figures by an independent quadrature.
#   geodesic      the exact geodesic on this spheroid is 5549.6819 km by
#                 Karney's algorithm, computed outside PSOPT entirely.
say "Examples with externally fixed answers"
export EXAMPLES="mineng_di,lts_costates,geodesic"
export REF_COSTS="6.000000e+00,4.195884e+02,5.549682e+03"
python3 "${SRC}/.github/scripts/run_examples.py" \
        --exe-dir "${BUILD}/examples" \
        --summary /tmp/distro_summary.json

say "PASS"
