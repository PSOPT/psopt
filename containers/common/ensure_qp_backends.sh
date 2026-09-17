#!/usr/bin/env bash
#
# Install the QP backends PSOPT's own SQP solver sends its subproblems to.
#
# This does not reimplement the installation. It calls scripts/build_qp_backends.sh and
# scripts/build_galahad.sh, which are the scripts the manual tells a user to run, so the
# matrix tests those too: if they stop working on a distribution, a job here says so
# before a user does. That is also why this file has to run after the source tree has
# been copied in, and so lives in build_and_test.sh's half of the work rather than in the
# Dockerfile's.
#
# Both are given --skip-deps, deliberately. Left to themselves they would apt-get or dnf
# their own prerequisites through sudo, which a container has no reason to allow and
# frequently has no sudo for; and neither knows zypper, so openSUSE would fall through
# their package-manager detection in any case. The prerequisites are the Dockerfile's
# job, which is where every other distribution difference lives.
#
# What each costs, since between them they dominate the running time of an image:
#
#   OSQP      C with CMake, about a minute.
#   PIQP      header-only C++ over Eigen, nothing to compile.
#   Clarabel  Rust, through the C interface of Clarabel.cpp. Needs a toolchain, and
#             installs rustup here when the distribution has not supplied one.
#   GALAHAD   a large Fortran library built with meson, and the slow one by a wide
#             margin. It is also the one PSOPT's SQP was measured against, which is why
#             it is worth the wait.
#
set -euo pipefail

SRC=${SRC:-/src}
QP_PREFIX=${QP_PREFIX:-/opt/qp}
JOBS=${JOBS:-$(nproc)}

say() { printf '\n\033[1m== %s\033[0m\n' "$*"; }

mkdir -p "${QP_PREFIX}"

# --------------------------------------------------------------------- Rust, for Clarabel
# Only Clarabel needs it. A distribution package is preferred when the Dockerfile has
# installed one, because it is quicker and it is what a user of that distribution would
# have; rustup is the fallback, and saying which happened matters, since a distribution
# Rust can be old enough for Clarabel.cpp to refuse it.
say "Rust toolchain (Clarabel only)"
if command -v cargo >/dev/null 2>&1; then
    echo "RUST-SOURCE: distribution package"
    echo "RUST-VERSION: $(cargo --version 2>&1 | head -1)"
else
    echo "RUST-SOURCE: rustup (no distribution package installed)"
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
        | sh -s -- -y --default-toolchain stable --profile minimal
    . "${HOME}/.cargo/env"
    echo "RUST-VERSION: $(cargo --version 2>&1 | head -1)"
fi
export PATH="${HOME}/.cargo/bin:${PATH}"

# ------------------------------------------------------------ OSQP, PIQP and Clarabel
say "OSQP, PIQP and Clarabel"
"${SRC}/scripts/build_qp_backends.sh" \
    --osqp --piqp --clarabel \
    --prefix "${QP_PREFIX}" \
    --jobs "${JOBS}" \
    --skip-deps --yes

# ------------------------------------------------------------------------- GALAHAD
say "GALAHAD"
"${SRC}/scripts/build_galahad.sh" \
    --prefix "${QP_PREFIX}" \
    --jobs "${JOBS}" \
    --skip-deps --yes

# ---------------------------------------------------------------------------- report
# One line per backend saying whether it is actually there, because the whole point of
# building four is to find out which of them a given distribution can carry. A job whose
# SQP tests pass on one backend and silently skip the other three has not tested what it
# appears to have tested.
say "What this image has"

osqp=$(find "${QP_PREFIX}" -name 'libosqpstatic.a'      2>/dev/null | head -1)
piqp=$(find "${QP_PREFIX}" -name 'piqp.hpp'             2>/dev/null | head -1)
clar=$(find "${QP_PREFIX}" -name 'libclarabel_c.a'      2>/dev/null | head -1)
gala=$(find "${QP_PREFIX}" -name 'libgalahad_double.so' 2>/dev/null | head -1)

printf 'QP-BACKEND osqp     : %s\n' "${osqp:-NOT BUILT}"
printf 'QP-BACKEND piqp     : %s\n' "${piqp:-NOT BUILT}"
printf 'QP-BACKEND clarabel : %s\n' "${clar:-NOT BUILT}"
printf 'QP-BACKEND galahad  : %s\n' "${gala:-NOT BUILT}"

# GALAHAD is the only backend installed as a SHARED library, so it is the only one whose
# own dependencies can be incomplete, and a plugin linking it then fails to load for a
# reason that has nothing to do with the plugin. Its dependencies are listed here, where
# they can be read, rather than left to be inferred from a load failure later.
if [ -n "${gala}" ]; then
    echo "--- ldd ${gala}"
    ldd "${gala}" 2>&1 | sed 's/^/    /' || true
fi

# if/fi rather than [ ] && ..., for legibility only. An earlier comment here claimed the
# AND-list form would trip set -e when the test failed, and that is not true: bash
# exempts every command in a && list except the one after the final &&, so the left
# operand of && never triggers it, even when the list is the last command of a
# redirected group. The claim was made without being tested and is corrected rather than
# quietly deleted, because it was also written into the message of the patch that
# introduced it.
n=0
for f in "${osqp}" "${piqp}" "${clar}" "${gala}"; do
    if [ -n "${f}" ]; then n=$((n+1)); fi
done
echo "QP-BACKENDS-BUILT: ${n} of 4"

# The cmake arguments that follow from what is actually here, written once so that
# build_and_test.sh does not repeat the detection and cannot disagree with it.
ARGS="${QP_PREFIX}/psopt-qp-backends.args"
{
    printf '%s' "-DWITH_SQP=ON"
    if [ -n "${osqp}" ]; then printf '%s' " -DWITH_OSQP=ON"; fi
    if [ -n "${piqp}" ]; then printf '%s' " -DWITH_PIQP=ON"; fi
    if [ -n "${clar}" ]; then printf '%s' " -DWITH_CLARABEL=ON -DCLARABEL_DIR=${QP_PREFIX}"; fi
    if [ -n "${gala}" ]; then printf '%s' " -DWITH_GALAHAD=ON -DGALAHAD_DIR=${QP_PREFIX}"; fi
    printf '\n'
} > "${ARGS}"
echo "cmake arguments written to ${ARGS}:"
cat "${ARGS}"

if [ "${n}" -eq 0 ]; then
    echo "No QP backend was built, so WITH_SQP cannot be configured and the SQP would go"
    echo "untested on this distribution. That is a failure rather than a result: the"
    echo "scripts above are PSOPT's own and are supposed to work here."
    exit 1
fi
