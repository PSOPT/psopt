#!/usr/bin/env bash
#
# build_qpoases.sh -- fetch, build and install qpOASES for PSOPT's SQP solver,
#                     on macOS or Linux.
#
# PSOPT's SQP currently requires qpOASES even when the quadratic programming
# subproblems are being sent to GALAHAD, because parts of the driver are written in
# qpOASES's vocabulary. That requirement is being removed; until it is, WITH_SQP=ON
# needs this library present. It has no bearing on which QP solver actually runs.
#
#   ./build_qpoases.sh                 # install under ~/qpoases-install
#   ./build_qpoases.sh --prefix /usr/local --sudo
#   ./build_qpoases.sh --real-blas     # link the platform BLAS instead (see below)
#   ./build_qpoases.sh --help
#
# The one thing that matters in this build
# ----------------------------------------
# qpOASES ships BLASReplacement.cpp and LAPACKReplacement.cpp -- small stand-ins for the
# handful of BLAS and LAPACK routines it needs -- and its CMakeLists globs src/*.cpp, so
# they are compiled in whether you want them or not. By default they define dgemm_ and
# dpotrf_ under exactly those names. Linked into a program that also has a real BLAS,
# they can capture those calls for everything in the image, including MUMPS inside IPOPT,
# which then crashes in its linear solver on problems that solved perfectly well before.
# The failure appears nowhere near qpOASES and is thoroughly unpleasant to find. PSOPT's
# CMake refuses to configure against a library in that state.
#
# This script builds with QPOASES_AVOID_LA_NAMING_CONFLICTS=ON, an option qpOASES
# provides for precisely this: the replacements are still compiled, but under the names
# qpOASES_gemm, qpOASES_dpotrf and so on, so they serve qpOASES and are invisible to
# everything else. Nothing has to be edited or stripped afterwards, and the result
# satisfies PSOPT's check by construction. The script verifies it before finishing.
#
# --real-blas instead deletes the two objects from the archive so qpOASES calls the
# platform's real BLAS. That is faster inside qpOASES and is the other supported answer.
# It makes no difference to PSOPT unless you intend to use qpOASES as the QP solver,
# which on any problem of interesting size you do not: it is a dense active-set method
# whose work per subproblem is cubic in the number of variables.
#
set -euo pipefail

PREFIX="${HOME}/qpoases-install"
SRCDIR="${HOME}/src/qpOASES"
QPOASES_REPO="https://github.com/coin-or/qpOASES.git"
QPOASES_REF="master"
JOBS=""
REAL_BLAS=0
SKIP_DEPS=0
USE_SUDO=0
ASSUME_YES=0

say()  { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
info() { printf '    %s\n' "$*"; }
die()  { printf '\n\033[1;31mError: %s\033[0m\n\n' "$*" >&2; exit 1; }

usage() {
    awk 'NR>1 { if ($0 !~ /^#/) exit; sub(/^# ?/, ""); print }' "$0"
    cat <<'EOF'

Options:
  --prefix DIR     where to install qpOASES      (default: ~/qpoases-install)
  --src DIR        where to clone the source     (default: ~/src/qpOASES)
  --ref REF        branch or tag to build        (default: master)
  --jobs N         parallel build jobs           (default: all cores)
  --real-blas      strip the bundled BLAS/LAPACK replacements and use the platform's
  --skip-deps      do not install any packages; assume cmake and a C++ compiler exist
  --sudo           use sudo for the install step (needed for a system prefix)
  --yes            do not prompt before installing packages
  --help           this message
EOF
    exit 0
}

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix)     PREFIX="$2";       shift 2 ;;
        --src)        SRCDIR="$2";       shift 2 ;;
        --ref)        QPOASES_REF="$2";  shift 2 ;;
        --jobs)       JOBS="$2";         shift 2 ;;
        --real-blas)  REAL_BLAS=1;       shift ;;
        --skip-deps)  SKIP_DEPS=1;       shift ;;
        --sudo)       USE_SUDO=1;        shift ;;
        --yes|-y)     ASSUME_YES=1;      shift ;;
        --help|-h)    usage ;;
        *)            die "unknown option '$1' (try --help)" ;;
    esac
done

OS="$(uname -s)"
PKGMGR=""
case "$OS" in
    Darwin)
        if   command -v port >/dev/null 2>&1; then PKGMGR="macports"
        elif command -v brew >/dev/null 2>&1; then PKGMGR="homebrew"
        fi ;;
    Linux)
        if   command -v apt-get >/dev/null 2>&1; then PKGMGR="apt"
        elif command -v dnf     >/dev/null 2>&1; then PKGMGR="dnf"
        elif command -v pacman  >/dev/null 2>&1; then PKGMGR="pacman"
        fi ;;
    *)  die "unsupported platform '$OS' (this script handles macOS and Linux)" ;;
esac

if [ -z "$JOBS" ]; then
    if   [ "$OS" = "Darwin" ] && command -v sysctl >/dev/null 2>&1; then JOBS="$(sysctl -n hw.ncpu)"
    elif command -v nproc >/dev/null 2>&1; then JOBS="$(nproc)"
    else JOBS=4
    fi
fi

say "qpOASES for PSOPT"
info "platform        $OS${PKGMGR:+ ($PKGMGR)}"
info "source          $SRCDIR  (ref: $QPOASES_REF)"
info "install prefix  $PREFIX"
info "parallel jobs   $JOBS"

# ---------------------------------------------------------------------------------
# Tools.  Unlike GALAHAD this needs nothing exotic: a C++ compiler, cmake and git.
# Apple clang is perfectly adequate here -- there is no Fortran and no OpenMP.
# ---------------------------------------------------------------------------------
confirm() {
    [ "$ASSUME_YES" = "1" ] && return 0
    printf '\n    %s\n    Proceed? [y/N] ' "$1"
    read -r reply </dev/tty || reply=""
    case "$reply" in [yY]*) return 0 ;; *) die "stopped at your request" ;; esac
}

install_deps() {
    case "$PKGMGR" in
        macports)
            if ! xcode-select -p >/dev/null 2>&1; then
                say "Xcode command line tools"
                info "These provide the C++ compiler and are needed before anything else."
                xcode-select --install || true
                die "rerun this script once the command line tools have finished installing"
            fi
            confirm "About to run: sudo port install cmake git"
            sudo port install cmake git ;;
        homebrew)
            confirm "About to run: brew install cmake git"
            brew install cmake git ;;
        apt)
            confirm "About to run: sudo apt-get install cmake g++ git"
            sudo apt-get update && sudo apt-get install -y cmake g++ git ;;
        dnf)
            confirm "About to run: sudo dnf install cmake gcc-c++ git"
            sudo dnf install -y cmake gcc-c++ git ;;
        pacman)
            confirm "About to run: sudo pacman -S cmake gcc git"
            sudo pacman -S --needed cmake gcc git ;;
        *)  info "no supported package manager found; skipping"
            info "you need: a C++ compiler, cmake, git" ;;
    esac
}

if [ "$SKIP_DEPS" = "1" ]; then
    say "Skipping package installation (--skip-deps)"
else
    say "Installing build tools"
    install_deps
fi

for t in cmake git; do
    command -v "$t" >/dev/null 2>&1 || die "'$t' not found. Install it, or rerun without --skip-deps."
done
say "Toolchain"
info "cmake     $(cmake --version | head -1)"

# ---------------------------------------------------------------------------------
# Source
# ---------------------------------------------------------------------------------
say "Fetching qpOASES"
if [ -d "$SRCDIR/.git" ]; then
    info "updating the existing clone at $SRCDIR"
    git -C "$SRCDIR" fetch --tags origin
    git -C "$SRCDIR" checkout "$QPOASES_REF"
    git -C "$SRCDIR" pull --ff-only origin "$QPOASES_REF" || true
else
    mkdir -p "$(dirname "$SRCDIR")"
    git clone "$QPOASES_REPO" "$SRCDIR"
    git -C "$SRCDIR" checkout "$QPOASES_REF"
fi
info "at $(git -C "$SRCDIR" log -1 --format='%h  %ci')"

# ---------------------------------------------------------------------------------
# Configure and build
#
# QPOASES_BUILD_EXAMPLES is off because PSOPT does not need them and they are the part
# of this project most likely to fail to build. A static library keeps the result a
# single file with no run-time loader path to arrange.
# ---------------------------------------------------------------------------------
BUILDDIR="$SRCDIR/build-psopt"

CMAKE_OPTS=(
    "-DCMAKE_BUILD_TYPE=Release"
    "-DCMAKE_INSTALL_PREFIX=$PREFIX"
    "-DCMAKE_POSITION_INDEPENDENT_CODE=ON"
    "-DQPOASES_BUILD_EXAMPLES=OFF"
    "-DBUILD_SHARED_LIBS=OFF"
)
if [ "$REAL_BLAS" = "1" ]; then
    info "BLAS: the platform's, with qpOASES's replacements stripped after the build"
else
    CMAKE_OPTS+=("-DQPOASES_AVOID_LA_NAMING_CONFLICTS=ON")
    info "BLAS: qpOASES's own, renamed so it cannot capture anyone else's"
fi

say "Configuring"
[ -d "$BUILDDIR" ] && { info "removing the previous build directory"; rm -rf "$BUILDDIR"; }
cmake -S "$SRCDIR" -B "$BUILDDIR" "${CMAKE_OPTS[@]}"

say "Building"
cmake --build "$BUILDDIR" -j "$JOBS"

say "Installing into $PREFIX"
if [ "$USE_SUDO" = "1" ]; then sudo cmake --install "$BUILDDIR"; else cmake --install "$BUILDDIR"; fi

# ---------------------------------------------------------------------------------
# Locate what was installed
# ---------------------------------------------------------------------------------
LIB=""
for d in "$PREFIX/lib" "$PREFIX/lib64"; do
    for e in a dylib so; do
        [ -f "$d/libqpOASES.$e" ] && LIB="$d/libqpOASES.$e" && break 2
    done
done
[ -n "$LIB" ] || die "libqpOASES was not installed under $PREFIX/lib. Read the build log above."
[ -f "$PREFIX/include/qpOASES.hpp" ] || die "qpOASES.hpp is missing from $PREFIX/include."

# ---------------------------------------------------------------------------------
# --real-blas: remove the replacement objects from the archive.
#
# This is what PSOPT's own error message recommends. macOS needs the archive's symbol
# table rebuilt afterwards, which is what ranlib is for; on Linux ar maintains it.
# ---------------------------------------------------------------------------------
if [ "$REAL_BLAS" = "1" ]; then
    say "Stripping the bundled BLAS/LAPACK replacements"
    case "$LIB" in
        *.a)
            for o in BLASReplacement.cpp.o LAPACKReplacement.cpp.o \
                     BLASReplacement.o LAPACKReplacement.o; do
                if ar t "$LIB" 2>/dev/null | grep -qx "$o"; then
                    info "removing $o"
                    ${USE_SUDO:+sudo} ar d "$LIB" "$o"
                fi
            done
            command -v ranlib >/dev/null 2>&1 && ${USE_SUDO:+sudo} ranlib "$LIB"
            ;;
        *)  info "not a static archive; nothing to strip. Rebuild without --real-blas if"
            info "PSOPT's configure step complains about dgemm_." ;;
    esac
fi

# ---------------------------------------------------------------------------------
# The check that matters: PSOPT refuses a qpOASES that defines dgemm_ or dpotrf_ as
# strong symbols. This is the same test its CMakeLists performs, run here so that a
# problem is found now rather than at PSOPT's configure step.
# ---------------------------------------------------------------------------------
say "Checking the install"
info "library   $LIB"
info "headers   $PREFIX/include/qpOASES.hpp"
if command -v nm >/dev/null 2>&1; then
    if nm "$LIB" 2>/dev/null | grep -qE " T _?(dgemm_|dpotrf_)$"; then
        die "the installed library still defines dgemm_ or dpotrf_ as strong symbols.
    PSOPT will refuse to configure against it, because those definitions can capture the
    real BLAS for the whole program and crash MUMPS inside IPOPT. Rerun this script
    without --real-blas to build with the renamed replacements instead."
    fi
    info "symbols   defines no dgemm_/dpotrf_ -- PSOPT's guard will accept it"
else
    info "symbols   nm not available; skipped the dgemm_/dpotrf_ check"
fi

ENVFILE="$PREFIX/qpoases-env.sh"
cat > "$ENVFILE" <<EOF
# qpOASES environment for PSOPT.  Source this, or copy it into your shell profile:
#     source $ENVFILE
export QPOASES_DIR="$PREFIX"
EOF

say "Done"
cat <<EOF

    qpOASES is installed in $PREFIX

    Configure PSOPT with both it and GALAHAD:

        source $ENVFILE
        source \${GALAHAD_PREFIX:-~/galahad-install}/galahad-env.sh

        cd /path/to/psopt
        cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=ON \\
              -DWITH_SQP=ON -DWITH_GALAHAD=ON \\
              -DQPOASES_DIR="$PREFIX" \\
              -DGALAHAD_DIR="\${GALAHAD_DIR:-\$HOME/galahad-install}"
        cmake --build build -j $JOBS

    qpOASES is a build-time requirement only. Leave algorithm.qp_solver = "GALAHAD";
    it is the sparse backend and the one the solver has been tuned against.

EOF
