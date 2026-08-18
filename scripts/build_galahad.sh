#!/usr/bin/env bash
#
# build_galahad.sh -- fetch, build and install GALAHAD for use as PSOPT's sparse QP
#                     backend, on macOS or Linux.
#
# PSOPT's SQP solver sends its quadratic programming subproblem to GALAHAD's QPA. This
# script installs the tools GALAHAD needs, clones it, configures it with the options
# PSOPT requires, builds it, and writes out the environment PSOPT and QPA expect. It
# does not touch PSOPT itself; it prints the cmake line to use when you are done.
#
# The meson options below are not a guess. They are the configuration of the GALAHAD
# build that every PSOPT SQP measurement in doc/SQP_ALL_EXAMPLES.md was made against.
#
#   ./build_galahad.sh                 # install under ~/galahad-install
#   ./build_galahad.sh --prefix /usr/local --sudo
#   ./build_galahad.sh --openblas      # link an optimised BLAS instead of GALAHAD's own
#   ./build_galahad.sh --skip-deps     # you have gfortran, meson, ninja and pkg-config
#   ./build_galahad.sh --help
#
# Two things about this build are worth knowing before you start.
#
# OpenMP is not optional. QPA's linear solver uses OpenMP *cancellation*, so GALAHAD
# must be compiled with -Dopenmp=true, and OMP_CANCELLATION=TRUE must be in the
# environment before a PSOPT program starts -- the OpenMP runtime reads that variable
# once, when it initialises, and nothing inside the process can set it afterwards.
# Without it the QP subproblems fail and the SQP looks simply broken. The env file this
# script writes sets it.
#
# The C interface is not optional either: PSOPT includes <galahad_qpa.h>, which only
# exists when GALAHAD is built with -Dciface=true.
#
set -euo pipefail

# ---------------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------------
PREFIX="${HOME}/galahad-install"
SRCDIR="${HOME}/src/GALAHAD"
GALAHAD_REPO="https://github.com/ralna/GALAHAD.git"
GALAHAD_REF="master"
JOBS=""
USE_OPENBLAS=0
SKIP_DEPS=0
USE_SUDO=0
WITH_TESTS=0
ASSUME_YES=0

say()  { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
info() { printf '    %s\n' "$*"; }
die()  { printf '\n\033[1;31mError: %s\033[0m\n\n' "$*" >&2; exit 1; }

usage() {
    # The header comment, up to the first line that is not a comment. Done by content
    # rather than by line number so that editing the header cannot silently spill the
    # script's own code into its help text.
    awk 'NR>1 { if ($0 !~ /^#/) exit; sub(/^# ?/, ""); print }' "$0"
    cat <<'EOF'

Options:
  --prefix DIR     where to install GALAHAD      (default: ~/galahad-install)
  --src DIR        where to clone the source     (default: ~/src/GALAHAD)
  --ref REF        branch or tag to build        (default: master)
  --jobs N         parallel build jobs           (default: all cores)
  --openblas       link OpenBLAS rather than building GALAHAD's own reference BLAS
  --with-tests     build and run GALAHAD's own test suite (slow, but thorough)
  --skip-deps      do not install any packages; assume the tools are present
  --sudo           use sudo for `meson install` (needed for a system prefix)
  --yes            do not prompt before installing packages
  --help           this message
EOF
    exit 0
}

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix)     PREFIX="$2";       shift 2 ;;
        --src)        SRCDIR="$2";       shift 2 ;;
        --ref)        GALAHAD_REF="$2";  shift 2 ;;
        --jobs)       JOBS="$2";         shift 2 ;;
        --openblas)   USE_OPENBLAS=1;    shift ;;
        --with-tests) WITH_TESTS=1;      shift ;;
        --skip-deps)  SKIP_DEPS=1;       shift ;;
        --sudo)       USE_SUDO=1;        shift ;;
        --yes|-y)     ASSUME_YES=1;      shift ;;
        --help|-h)    usage ;;
        *)            die "unknown option '$1' (try --help)" ;;
    esac
done

# ---------------------------------------------------------------------------------
# Platform and package manager
# ---------------------------------------------------------------------------------
OS="$(uname -s)"
PKGMGR=""

case "$OS" in
    Darwin)
        if   command -v port  >/dev/null 2>&1; then PKGMGR="macports"
        elif command -v brew  >/dev/null 2>&1; then PKGMGR="homebrew"
        fi
        ;;
    Linux)
        if   command -v apt-get >/dev/null 2>&1; then PKGMGR="apt"
        elif command -v dnf     >/dev/null 2>&1; then PKGMGR="dnf"
        elif command -v pacman  >/dev/null 2>&1; then PKGMGR="pacman"
        fi
        ;;
    *)  die "unsupported platform '$OS' (this script handles macOS and Linux)" ;;
esac

if [ -z "$JOBS" ]; then
    if   command -v sysctl  >/dev/null 2>&1 && [ "$OS" = "Darwin" ]; then JOBS="$(sysctl -n hw.ncpu)"
    elif command -v nproc   >/dev/null 2>&1; then JOBS="$(nproc)"
    else JOBS=4
    fi
fi

say "GALAHAD for PSOPT"
info "platform        $OS${PKGMGR:+ ($PKGMGR)}"
info "source          $SRCDIR  (ref: $GALAHAD_REF)"
info "install prefix  $PREFIX"
info "parallel jobs   $JOBS"

# ---------------------------------------------------------------------------------
# Tools
#
# GALAHAD is Fortran built with meson, so it needs a Fortran compiler, meson, ninja and
# pkg-config, plus a C compiler for the interface PSOPT uses. On macOS the C compiler
# comes from the Xcode command line tools; the Fortran compiler does not exist on the
# system at all and arrives with gcc from MacPorts or Homebrew.
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
                info "These provide the C compiler and are needed before anything else."
                xcode-select --install || true
                die "rerun this script once the command line tools have finished installing"
            fi
            confirm "About to run: sudo port install gcc14 meson ninja pkgconfig"
            sudo port install gcc14 meson ninja pkgconfig
            # MacPorts installs the compiler as gfortran-mp-14 and only creates a plain
            # 'gfortran' if you ask for it. Doing so keeps every later build simple.
            if ! command -v gfortran >/dev/null 2>&1; then
                info "selecting mp-gcc14 so that a plain 'gfortran' exists"
                sudo port select --set gcc mp-gcc14 || true
            fi
            ;;
        homebrew)
            confirm "About to run: brew install gcc meson ninja pkg-config"
            brew install gcc meson ninja pkg-config
            ;;
        apt)
            confirm "About to run: sudo apt-get install gfortran meson ninja-build pkg-config git"
            sudo apt-get update
            sudo apt-get install -y gfortran g++ meson ninja-build pkg-config git
            ;;
        dnf)
            confirm "About to run: sudo dnf install gcc-gfortran meson ninja-build pkgconfig git"
            sudo dnf install -y gcc-gfortran gcc-c++ meson ninja-build pkgconfig git
            ;;
        pacman)
            confirm "About to run: sudo pacman -S gcc-fortran meson ninja pkgconf git"
            sudo pacman -S --needed gcc-fortran meson ninja pkgconf git
            ;;
        *)
            info "no supported package manager found; skipping the install step"
            info "you need: a Fortran compiler, a C compiler, meson, ninja, pkg-config, git"
            ;;
    esac
}

if [ "$SKIP_DEPS" = "1" ]; then
    say "Skipping package installation (--skip-deps)"
else
    say "Installing build tools"
    install_deps
fi

# ---------------------------------------------------------------------------------
# Locate the compilers.
#
# GALAHAD is a mixed Fortran and C project and meson chooses each compiler
# independently: $FC or 'gfortran' for the Fortran, $CC or 'cc' for the C. Setting only
# FC is not enough, and on macOS it fails in a way that looks unrelated to the omission:
# 'cc' is Apple clang, Apple clang has no -fopenmp, and the build dies part way through
# the C interface with
#
#     clang: error: unsupported option '-fopenmp'
#
# -fopenmp cannot be dropped, because QPA's linear solver needs OpenMP cancellation.
# Apple clang can be made to do OpenMP with a separately installed libomp and a
# -Xpreprocessor incantation, but that is the LLVM runtime, and gfortran's half of
# GALAHAD would still be using GCC's libgomp: two OpenMP runtimes in one process, which
# is a well-known way to get hangs and crashes rather than an answer.
#
# So the C compiler is chosen to match the Fortran one -- the same GCC installation
# provides both -- and it is chosen by *testing* that it compiles an OpenMP program
# rather than by trusting its name. On macOS 'gcc' is usually Apple clang wearing a
# different hat, so the name proves nothing.
# ---------------------------------------------------------------------------------
find_fortran() {
    if [ -n "${FC:-}" ] && command -v "$FC" >/dev/null 2>&1; then echo "$FC"; return; fi
    if command -v gfortran >/dev/null 2>&1; then echo "gfortran"; return; fi
    for v in 15 14 13 12 11 10; do
        for n in "gfortran-mp-$v" "gfortran-$v"; do
            if command -v "$n" >/dev/null 2>&1; then echo "$n"; return; fi
        done
    done
    echo ""
}

# Does this compiler actually build and link an OpenMP program with plain -fopenmp?
supports_openmp() {
    local cc="$1" tmp rc
    command -v "$cc" >/dev/null 2>&1 || return 1
    tmp="$(mktemp -d)"
    printf '#include <omp.h>\nint main(void){ return omp_get_max_threads() > 0 ? 0 : 1; }\n' > "$tmp/omptest.c"
    "$cc" -fopenmp "$tmp/omptest.c" -o "$tmp/omptest" >/dev/null 2>&1
    rc=$?
    rm -rf "$tmp"
    return $rc
}

# Candidate C compilers, best first: an explicit $CC, then the gcc belonging to the same
# installation as $FC (gfortran-mp-14 -> gcc-mp-14, gfortran-14 -> gcc-14), then the gcc
# of the same version as $FC -- which is what catches a plain 'gfortran' put there by
# 'port select', where the name carries no version but the compiler still reports one --
# then any versioned gcc, then the plain names.
c_candidates() {
    local fcbase="${1##*/}" fcver="${2:-}" v
    [ -n "${CC:-}" ] && echo "$CC"
    case "$fcbase" in
        gfortran-mp-*) echo "gcc-mp-${fcbase##*-}" ;;
        gfortran-*)    echo "gcc-${fcbase##*-}" ;;
    esac
    if [ -n "$fcver" ]; then echo "gcc-mp-$fcver"; echo "gcc-$fcver"; fi
    for v in 15 14 13 12 11 10; do echo "gcc-mp-$v"; echo "gcc-$v"; done
    echo gcc
    echo cc
}

FC_FOUND="$(find_fortran)"
[ -n "$FC_FOUND" ] || die "no Fortran compiler found. Install one (MacPorts: sudo port install gcc14; Homebrew: brew install gcc) and rerun, or set FC."
export FC="$FC_FOUND"

say "Toolchain"
info "Fortran   $FC  ($("$FC" --version 2>/dev/null | head -1))"

FC_VER="$("$FC" -dumpversion 2>/dev/null || true)"
FC_VER="${FC_VER%%.*}"

CC_FOUND=""
while read -r cand; do
    [ -n "$cand" ] || continue
    if supports_openmp "$cand"; then CC_FOUND="$cand"; break; fi
done <<EOF
$(c_candidates "$FC" "$FC_VER")
EOF

if [ -z "$CC_FOUND" ]; then
    if [ "$OS" = "Darwin" ]; then
        die "no C compiler that supports -fopenmp was found.

    This is the usual macOS problem: 'cc' and often 'gcc' are Apple clang, which does not
    accept -fopenmp, and GALAHAD's QPA needs OpenMP. Install a real GCC and rerun:

        MacPorts:  sudo port install gcc14
        Homebrew:  brew install gcc

    If you have one under a name this script did not try, set CC and FC yourself, e.g.

        CC=gcc-mp-14 FC=gfortran-mp-14 $0 --skip-deps"
    fi
    die "no C compiler that supports -fopenmp was found. Install GCC, or set CC to one that does."
fi
export CC="$CC_FOUND"
info "C         $CC  ($("$CC" --version 2>/dev/null | head -1))"
info "          (verified: compiles and links -fopenmp)"

# A mismatched pair -- Apple clang for C, MacPorts gfortran for Fortran -- links two
# different OpenMP runtimes into one library. Worth saying out loud if it happens.
case "${CC##*/}" in
    gcc*) : ;;
    *)    info "note: $CC is not a gcc; if GALAHAD misbehaves at run time, try a matching GCC pair" ;;
esac
CC_VER="$("$CC" -dumpversion 2>/dev/null || true)"
if [ -n "$FC_VER" ] && [ -n "${CC_VER%%.*}" ] && [ "$FC_VER" != "${CC_VER%%.*}" ]; then
    info "note: GCC $CC_VER for C and GCC $FC_VER for Fortran. Both use libgomp so this"
    info "      normally works; set CC and FC yourself if you would rather they matched."
fi

for t in meson ninja pkg-config git; do
    command -v "$t" >/dev/null 2>&1 || die "'$t' not found. Install it, or rerun without --skip-deps."
done
info "meson     $(meson --version)"
info "ninja     $(ninja --version)"

# ---------------------------------------------------------------------------------
# Source
# ---------------------------------------------------------------------------------
say "Fetching GALAHAD"
if [ -d "$SRCDIR/.git" ]; then
    info "updating the existing clone at $SRCDIR"
    git -C "$SRCDIR" fetch --tags origin
    git -C "$SRCDIR" checkout "$GALAHAD_REF"
    git -C "$SRCDIR" pull --ff-only origin "$GALAHAD_REF" || true
else
    mkdir -p "$(dirname "$SRCDIR")"
    git clone "$GALAHAD_REPO" "$SRCDIR"
    git -C "$SRCDIR" checkout "$GALAHAD_REF"
fi
info "at $(git -C "$SRCDIR" log -1 --format='%h  %ci')"

# ---------------------------------------------------------------------------------
# Configure
#
# Every option here earns its place:
#
#   openmp=true    QPA's linear solver needs OpenMP cancellation. Without this the
#                  runtime OMP_CANCELLATION setting has nothing to act on.
#   ciface=true    PSOPT includes <galahad_qpa.h>, which is generated only for this.
#   double=true    PSOPT links libgalahad_double.
#   single, quadruple, pythoniface, examples, binaries = false
#                  PSOPT uses none of them and each one costs build time.
#
# BLAS and LAPACK default to GALAHAD's own reference implementation. That is the slower
# choice and the reliable one: it is what the PSOPT SQP benchmarks were run against,
# and it avoids the question of whether the platform's optimised BLAS was built with
# matching integer widths. --openblas overrides it.
# ---------------------------------------------------------------------------------
BUILDDIR="$SRCDIR/builddir/psopt"

MESON_OPTS=(
    "--prefix=$PREFIX"
    "--buildtype=release"
    "-Ddefault_library=shared"
    "-Dopenmp=true"
    "-Dciface=true"
    "-Dmodules=true"
    "-Ddouble=true"
    "-Dsingle=false"
    "-Dquadruple=false"
    "-Dpythoniface=false"
    "-Dexamples=false"
    "-Dbinaries=false"
)
if [ "$WITH_TESTS" = "1" ]; then MESON_OPTS+=("-Dtests=true"); else MESON_OPTS+=("-Dtests=false"); fi

if [ "$USE_OPENBLAS" = "1" ]; then
    MESON_OPTS+=("-Dlibblas=openblas" "-Dliblapack=openblas")
    # MacPorts and Homebrew are not on meson's default library search path.
    for d in /opt/local/lib /opt/homebrew/lib /usr/local/lib; do
        [ -d "$d" ] && MESON_OPTS+=("-Dlibblas_path=$d" "-Dliblapack_path=$d") && break
    done
    info "BLAS/LAPACK: OpenBLAS (meson falls back to GALAHAD's own if it is not found)"
else
    MESON_OPTS+=("-Dlibblas=" "-Dliblapack=")
    info "BLAS/LAPACK: GALAHAD's own reference build"
fi

say "Configuring"
cd "$SRCDIR"
# An existing build directory is removed outright rather than reconfigured. meson caches
# its choice of compiler, and 'meson setup --wipe' keeps the options it was first given,
# so a directory left behind by a run that picked the wrong compiler will quietly go on
# using it -- which is precisely the situation anyone rerunning this script after the
# clang '-fopenmp' failure is in. --wipe rebuilds everything from scratch in any case,
# so nothing is lost by being unambiguous about it.
if [ -d "$BUILDDIR" ]; then
    info "removing the previous build directory at $BUILDDIR"
    rm -rf "$BUILDDIR"
fi
meson setup "$BUILDDIR" "${MESON_OPTS[@]}"

say "Building  (this takes a while -- GALAHAD is large)"
meson compile -C "$BUILDDIR" -j "$JOBS"

if [ "$WITH_TESTS" = "1" ]; then
    say "Running GALAHAD's test suite"
    meson test -C "$BUILDDIR" || info "some GALAHAD tests failed; this is not always fatal for QPA -- read the log above"
fi

say "Installing into $PREFIX"
if [ "$USE_SUDO" = "1" ]; then
    sudo meson install -C "$BUILDDIR"
else
    meson install -C "$BUILDDIR"
fi

# ---------------------------------------------------------------------------------
# Verify that what PSOPT actually looks for is there
# ---------------------------------------------------------------------------------
say "Checking the install"
LIBFOUND=""
for e in so dylib a; do
    for d in "$PREFIX/lib" "$PREFIX/lib64"; do
        [ -f "$d/libgalahad_double.$e" ] && LIBFOUND="$d/libgalahad_double.$e" && break 2
    done
done
[ -n "$LIBFOUND" ] || die "libgalahad_double was not installed under $PREFIX/lib. Read the build log above."
[ -f "$PREFIX/include/galahad_qpa.h" ] || die "galahad_qpa.h is missing: GALAHAD was built without the C interface. Rerun; -Dciface=true is set by this script."
info "library   $LIBFOUND"
info "header    $PREFIX/include/galahad_qpa.h"

# ---------------------------------------------------------------------------------
# The environment
# ---------------------------------------------------------------------------------
ENVFILE="$PREFIX/galahad-env.sh"
if [ "$OS" = "Darwin" ]; then LIBVAR="DYLD_LIBRARY_PATH"; else LIBVAR="LD_LIBRARY_PATH"; fi

cat > "$ENVFILE" <<EOF
# GALAHAD environment for PSOPT.  Source this, or copy it into your shell profile:
#     source $ENVFILE

# Where PSOPT's cmake should look for GALAHAD.
export GALAHAD_DIR="$PREFIX"

# So the loader finds libgalahad_double at run time.
export $LIBVAR="$PREFIX/lib\${$LIBVAR:+:\$$LIBVAR}"

# GALAHAD's QPA uses OpenMP cancellation. The OpenMP runtime reads these once, when it
# initialises, so they cannot be set from inside a PSOPT program -- they must already be
# in the environment when it starts. Without them the QP subproblems fail and the SQP
# makes no progress at all.
export OMP_CANCELLATION=TRUE
export OMP_PROC_BIND=TRUE
EOF

say "Done"
cat <<EOF

    GALAHAD $(grep -m1 "version:" "$SRCDIR/meson.build" | sed "s/.*'\(.*\)'.*/\1/") is installed in $PREFIX

    1.  Load the environment (add this line to your ~/.zshrc or ~/.bashrc):

            source $ENVFILE

    2.  Configure PSOPT against it:

            source $ENVFILE
            cd /path/to/psopt
            cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=ON \\
                  -DWITH_SQP=ON -DWITH_GALAHAD=ON \\
                  -DGALAHAD_DIR="$PREFIX"
            cmake --build build -j $JOBS

    3.  In your problem:

            algorithm.nlp_method   = "SQP";
            algorithm.hessian      = "exact";
            algorithm.qp_solver    = "GALAHAD";
            algorithm.sqp_strategy = "FM";

    If the SQP reports that every QP subproblem failed, the first thing to check is that
    OMP_CANCELLATION=TRUE was exported *before* the program started.

EOF
