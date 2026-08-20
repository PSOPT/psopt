#!/usr/bin/env bash
#
# build_piqp_clarabel.sh -- fetch, build and install PIQP and Clarabel for use as PSOPT
#                           QP backends, on macOS or Linux.
#
# PSOPT's SQP solver sends its quadratic programming subproblem to a backend. GALAHAD's
# QPA is the default and scripts/build_galahad.sh installs it; this script installs the
# two interior-point backends, either or both, into one prefix and writes out the
# environment PSOPT's cmake wants. It does not touch PSOPT itself; it prints the cmake
# line to use when it is done.
#
#   ./build_piqp_clarabel.sh                    # both, under ~/qp-backends
#   ./build_piqp_clarabel.sh --piqp-only        # PIQP alone -- no Rust toolchain needed
#   ./build_piqp_clarabel.sh --prefix /usr/local --sudo
#   ./build_piqp_clarabel.sh --skip-deps        # you have cmake, a compiler and cargo
#   ./build_piqp_clarabel.sh --help
#
# What each of them costs you:
#
# PIQP is header-only C++14 over Eigen, which PSOPT already requires, so there is nothing
# to compile and no library to load at run time. It is installed here with its template
# instantiation turned off, which keeps it header-only: the plugin then carries no
# dependency of its own, which is the whole point of the plugin boundary.
#
# Clarabel is Rust, reached through the C interface of Clarabel.cpp, so it needs a Rust
# toolchain -- and if you have not got one this script offers to install rustup. Its
# build produces a static library that the plugin links, so again nothing has to be found
# at run time. Clarabel.cpp has no install target of its own, so this script does that
# part itself: the headers and the static library are copied into the prefix.
#
set -euo pipefail

# ---------------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------------
PREFIX="${HOME}/qp-backends"
SRCDIR="${HOME}/src"
PIQP_REPO="https://github.com/PREDICT-EPFL/piqp.git"
PIQP_REF="main"
CLARABEL_REPO="https://github.com/oxfordcontrol/Clarabel.cpp.git"
CLARABEL_REF="main"
JOBS=""
DO_PIQP=1
DO_CLARABEL=1
SKIP_DEPS=0
USE_SUDO=0
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
  --prefix DIR      where to install both       (default: ~/qp-backends)
  --src DIR         where to clone the sources  (default: ~/src)
  --piqp-only       build PIQP and skip Clarabel
  --clarabel-only   build Clarabel and skip PIQP
  --piqp-ref REF    branch or tag for PIQP      (default: main)
  --clarabel-ref REF  branch or tag for Clarabel  (default: main)
  --jobs N          parallel build jobs         (default: all cores)
  --skip-deps       do not install anything; assume cmake and cargo are present
  --sudo            use sudo when installing (needed for a system prefix)
  --yes             do not prompt before installing anything
  --help            this message
EOF
    exit 0
}

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix)        PREFIX="$2";        shift 2 ;;
        --src)           SRCDIR="$2";        shift 2 ;;
        --piqp-ref)      PIQP_REF="$2";      shift 2 ;;
        --clarabel-ref)  CLARABEL_REF="$2";  shift 2 ;;
        --jobs)          JOBS="$2";          shift 2 ;;
        --piqp-only)     DO_CLARABEL=0;      shift ;;
        --clarabel-only) DO_PIQP=0;          shift ;;
        --skip-deps)     SKIP_DEPS=1;        shift ;;
        --sudo)          USE_SUDO=1;         shift ;;
        --yes|-y)        ASSUME_YES=1;       shift ;;
        --help|-h)       usage ;;
        *)               die "unknown option '$1' (try --help)" ;;
    esac
done

[ "$DO_PIQP" = "1" ] || [ "$DO_CLARABEL" = "1" ] || die "--piqp-only and --clarabel-only together leave nothing to build"

# ---------------------------------------------------------------------------------
# Platform
# ---------------------------------------------------------------------------------
OS="$(uname -s)"
PKGMGR=""

case "$OS" in
    Darwin)
        if   command -v port >/dev/null 2>&1; then PKGMGR="macports"
        elif command -v brew >/dev/null 2>&1; then PKGMGR="homebrew"
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
    if   [ "$OS" = "Darwin" ] && command -v sysctl >/dev/null 2>&1; then JOBS="$(sysctl -n hw.ncpu)"
    elif command -v nproc >/dev/null 2>&1; then JOBS="$(nproc)"
    else JOBS=4
    fi
fi

RUN_INSTALL=(cmake --install)
[ "$USE_SUDO" = "1" ] && RUN_INSTALL=(sudo cmake --install)

say "PIQP and Clarabel for PSOPT"
info "platform        $OS${PKGMGR:+ ($PKGMGR)}"
info "sources         $SRCDIR"
info "install prefix  $PREFIX"
info "parallel jobs   $JOBS"
info "building        $([ "$DO_PIQP" = 1 ] && printf 'PIQP ')$([ "$DO_CLARABEL" = 1 ] && printf 'Clarabel')"

confirm() {
    [ "$ASSUME_YES" = "1" ] && return 0
    printf '\n    %s\n    Proceed? [y/N] ' "$1"
    read -r reply </dev/tty || reply=""
    case "$reply" in [yY]*) return 0 ;; *) die "stopped at your request" ;; esac
}

# ---------------------------------------------------------------------------------
# Tools
#
# Both want cmake, git and a C++ compiler; PSOPT needs all three already, so this is a
# check rather than an install in nearly every case. Clarabel additionally wants cargo,
# which is the one thing a PSOPT machine has no reason to have.
# ---------------------------------------------------------------------------------
say "Checking tools"

for t in git cmake; do
    command -v "$t" >/dev/null 2>&1 || die "$t is not on the PATH, and both builds need it"
done
info "git             $(git --version | awk '{print $3}')"
info "cmake           $(cmake --version | head -1 | awk '{print $3}')"

if [ "$DO_CLARABEL" = "1" ]; then
    if ! command -v cargo >/dev/null 2>&1; then
        # rustup puts cargo here and does not touch the current shell, so a previous
        # install in another session looks like no install at all until it is sourced.
        if [ -f "$HOME/.cargo/env" ]; then
            # shellcheck disable=SC1091
            . "$HOME/.cargo/env"
        fi
    fi
    if ! command -v cargo >/dev/null 2>&1; then
        if [ "$SKIP_DEPS" = "1" ]; then
            die "cargo is not on the PATH and --skip-deps was given. Install Rust, or use --piqp-only."
        fi
        confirm "Clarabel is written in Rust and there is no Rust toolchain here.
    About to run the official installer:  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    (--piqp-only skips Clarabel altogether if you would rather not.)"
        curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
        # shellcheck disable=SC1091
        . "$HOME/.cargo/env"
    fi
    command -v cargo >/dev/null 2>&1 || die "cargo still is not on the PATH; open a new shell and rerun, or use --piqp-only"
    info "cargo           $(cargo --version | awk '{print $2}')"
fi

mkdir -p "$SRCDIR"

# A clone that is already there is updated rather than replaced, so that a rerun after a
# failure does not spend the download again.
clone_or_update() {
    local repo="$1" dir="$2" ref="$3" recurse="$4"
    if [ -d "$dir/.git" ]; then
        info "updating $dir"
        git -C "$dir" fetch --quiet origin
        git -C "$dir" checkout --quiet "$ref"
        git -C "$dir" pull --quiet --ff-only origin "$ref" || true
        [ "$recurse" = "1" ] && git -C "$dir" submodule update --quiet --init --recursive
    else
        info "cloning $repo"
        if [ "$recurse" = "1" ]; then
            git clone --quiet --recurse-submodules --branch "$ref" "$repo" "$dir"
        else
            git clone --quiet --branch "$ref" "$repo" "$dir"
        fi
    fi
}

# ---------------------------------------------------------------------------------
# PIQP
# ---------------------------------------------------------------------------------
if [ "$DO_PIQP" = "1" ]; then
    say "PIQP"
    PIQP_SRC="$SRCDIR/piqp"
    clone_or_update "$PIQP_REPO" "$PIQP_SRC" "$PIQP_REF" 0

    # BUILD_WITH_TEMPLATE_INSTANTIATION=OFF is the one option that matters here. With it
    # on, piqp::piqp is a shared library of pre-instantiated templates and the PSOPT
    # plugin would acquire a run-time dependency on it; with it off the target is a pure
    # interface and the plugin stays self-contained.
    cmake -S "$PIQP_SRC" -B "$PIQP_SRC/build" \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX="$PREFIX" \
          -DBUILD_WITH_TEMPLATE_INSTANTIATION=OFF \
          -DBUILD_C_INTERFACE=OFF \
          -DBUILD_TESTS=OFF \
          -DBUILD_EXAMPLES=OFF \
          -DBUILD_BENCHMARKS=OFF

    info "installing into $PREFIX"
    "${RUN_INSTALL[@]}" "$PIQP_SRC/build" >/dev/null

    [ -f "$PREFIX/include/piqp/piqp.hpp" ] || die "piqp/piqp.hpp was not installed under $PREFIX/include. Read the log above."
    PIQP_CONFIG="$(find "$PREFIX" -name 'piqpConfig.cmake' -print -quit 2>/dev/null || true)"
    [ -n "$PIQP_CONFIG" ] || die "piqpConfig.cmake was not installed; PSOPT's find_package(piqp) will not see it."
    info "header          $PREFIX/include/piqp/piqp.hpp"
    info "cmake config    $PIQP_CONFIG"
fi

# ---------------------------------------------------------------------------------
# Clarabel
# ---------------------------------------------------------------------------------
if [ "$DO_CLARABEL" = "1" ]; then
    say "Clarabel"
    CLARABEL_SRC="$SRCDIR/Clarabel.cpp"
    clone_or_update "$CLARABEL_REPO" "$CLARABEL_SRC" "$CLARABEL_REF" 1

    cmake -S "$CLARABEL_SRC" -B "$CLARABEL_SRC/build" -DCMAKE_BUILD_TYPE=Release

    say "Building Clarabel (cargo will fetch its crates the first time)"
    cmake --build "$CLARABEL_SRC/build" -j "$JOBS"

    # Clarabel.cpp has no install target, so the two things PSOPT looks for are placed
    # by hand: the C headers, and the static library cargo leaves in its target tree.
    CLARABEL_LIB=""
    for d in "$CLARABEL_SRC/rust_wrapper/target/release" "$CLARABEL_SRC/build/rust_wrapper/target/release"; do
        [ -f "$d/libclarabel_c.a" ] && CLARABEL_LIB="$d/libclarabel_c.a" && break
    done
    [ -n "$CLARABEL_LIB" ] || die "libclarabel_c.a was not produced. Read the build log above."

    info "installing into $PREFIX"
    mkdir -p "$PREFIX/include" "$PREFIX/lib"
    if [ "$USE_SUDO" = "1" ]; then
        sudo cp -R "$CLARABEL_SRC/include/." "$PREFIX/include/"
        sudo cp "$CLARABEL_LIB" "$PREFIX/lib/"
    else
        cp -R "$CLARABEL_SRC/include/." "$PREFIX/include/"
        cp "$CLARABEL_LIB" "$PREFIX/lib/"
    fi

    [ -f "$PREFIX/include/clarabel.h" ]      || die "clarabel.h was not installed under $PREFIX/include"
    [ -f "$PREFIX/lib/libclarabel_c.a" ]     || die "libclarabel_c.a was not installed under $PREFIX/lib"
    info "header          $PREFIX/include/clarabel.h"
    info "library         $PREFIX/lib/libclarabel_c.a"
fi

# ---------------------------------------------------------------------------------
# The environment
# ---------------------------------------------------------------------------------
ENVFILE="$PREFIX/qp-backends-env.sh"
{
    echo "# PIQP and Clarabel environment for PSOPT.  Source this, or copy it into your"
    echo "# shell profile:"
    echo "#     source $ENVFILE"
    echo
    if [ "$DO_PIQP" = "1" ]; then
        echo "# Where PSOPT's find_package(piqp) should look."
        echo "export piqp_DIR=\"$(dirname "$PIQP_CONFIG")\""
        echo
    fi
    if [ "$DO_CLARABEL" = "1" ]; then
        echo "# Where PSOPT's cmake should look for Clarabel's header and static library."
        echo "export CLARABEL_DIR=\"$PREFIX\""
        echo
    fi
    echo "# Neither backend needs anything on the loader path: PIQP is header-only and the"
    echo "# Clarabel plugin links the static library."
} > "$ENVFILE"

say "Done"
cat <<EOF

    Installed in $PREFIX

    1.  Load the environment (add this line to your ~/.zshrc or ~/.bashrc):

            source $ENVFILE

    2.  Configure PSOPT against it. Keep GALAHAD in the same build -- it is the default
        backend and the one the solver is measured against:

            source $ENVFILE
            cd /path/to/psopt
            cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=ON -DBUILD_TESTS=ON \\
                  -DWITH_SQP=ON \\
                  -DWITH_GALAHAD=ON -DGALAHAD_DIR="\$GALAHAD_DIR" \\
$([ "$DO_PIQP" = 1 ] && printf '                  -DWITH_PIQP=ON \\\\\n')\
$([ "$DO_CLARABEL" = 1 ] && printf '                  -DWITH_CLARABEL=ON -DCLARABEL_DIR="$CLARABEL_DIR" \\\\\n')\
                  -DPSOPT_ALLOW_ENV_OVERRIDES=ON
            cmake --build build -j $JOBS

        The configure output should name each backend it enabled. Then run ctest once:
        with several backends built, qp_plugins_export_nothing_else is doing real work,
        and the sign-convention tests pin each new backend's multipliers against the QP
        every backend is pinned against.

    3.  In your problem, or from the environment if the build has the override option:

            algorithm.nlp_method = "SQP";
            algorithm.hessian    = "exact";
            algorithm.qp_solver  = "PIQP";        // or "Clarabel", "GALAHAD", "OSQP"

    Both are interior-point methods and need the subproblem's Hessian to be positive
    semidefinite, which an optimal control problem's rarely is; PSOPT raises its shift
    until the backend accepts the model, so this is handled, but it is why GALAHAD's QPA
    -- which takes an indefinite Hessian as posed -- remains the default.

EOF
