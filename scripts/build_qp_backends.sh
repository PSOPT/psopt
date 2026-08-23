#!/usr/bin/env bash
#
# build_qp_backends.sh -- fetch, build and install PSOPT's alternative QP backends --
#                         OSQP, PIQP and Clarabel -- on macOS or Linux.
#
# PSOPT's SQP solver sends its quadratic programming subproblem to a backend. GALAHAD's
# QPA is the default and scripts/build_galahad.sh installs it; this script installs the
# other three, any or all of them, into one prefix and writes out the environment PSOPT's
# cmake wants. It does not touch PSOPT itself; it prints the cmake line to use when it is
# done.
#
#   ./build_qp_backends.sh                    # all three, under ~/qp-backends
#   ./build_qp_backends.sh --osqp             # just OSQP
#   ./build_qp_backends.sh --osqp --piqp      # the two that need no Rust toolchain
#   ./build_qp_backends.sh --prefix /usr/local --sudo
#   ./build_qp_backends.sh --skip-deps        # you have cmake, a compiler and cargo
#   ./build_qp_backends.sh --help
#
# What each of them costs you:
#
# OSQP is C with CMake and no dependencies, and builds in under a minute. Only the static
# library is built here, because that is the one PSOPT's plugin links (osqp::osqpstatic)
# and a shared one would be an installed file nothing uses. The version matters: PSOPT's
# plugin is written against the 1.x API -- OSQPInt, OSQPCscMatrix_new, OSQPSettings_new --
# and will not compile against 0.6.x, so this script checks the installed headers and
# says so plainly rather than leaving you with a page of compiler errors.
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
OSQP_REPO="https://github.com/osqp/osqp.git"
OSQP_REF="v1.0.0"
PIQP_REPO="https://github.com/PREDICT-EPFL/piqp.git"
PIQP_REF="main"
CLARABEL_REPO="https://github.com/oxfordcontrol/Clarabel.cpp.git"
CLARABEL_REF="main"
JOBS=""
# With no backend named, all three are built. Naming any switches to naming them all.
DO_OSQP=1
DO_PIQP=1
DO_CLARABEL=1
CHOSE_ANY=0
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
  --prefix DIR      where to install them       (default: ~/qp-backends)
  --src DIR         where to clone the sources  (default: ~/src)
  --osqp            build OSQP     (naming any backend builds only the ones named)
  --piqp            build PIQP
  --clarabel        build Clarabel
  --osqp-ref REF    branch or tag for OSQP      (default: v1.0.0)
  --piqp-ref REF    branch or tag for PIQP      (default: main)
  --clarabel-ref REF  branch or tag for Clarabel  (default: main)
  --jobs N          parallel build jobs         (default: all cores)
  --skip-deps       do not install anything; assume cmake and cargo are present
  --sudo            use sudo when installing (needed for a system prefix)
  --yes             do not prompt before installing anything
  --help            this message

  --piqp-only and --clarabel-only are accepted as they were before this script grew
  OSQP, and mean the same as --piqp and --clarabel.
EOF
    exit 0
}

# The first named backend clears the "all of them" default, so that --osqp means OSQP
# and nothing else, and --osqp --piqp means those two.
pick() {
    if [ "$CHOSE_ANY" = "0" ]; then
        CHOSE_ANY=1; DO_OSQP=0; DO_PIQP=0; DO_CLARABEL=0
    fi
}

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix)        PREFIX="$2";        shift 2 ;;
        --src)           SRCDIR="$2";        shift 2 ;;
        --osqp-ref)      OSQP_REF="$2";      shift 2 ;;
        --piqp-ref)      PIQP_REF="$2";      shift 2 ;;
        --clarabel-ref)  CLARABEL_REF="$2";  shift 2 ;;
        --jobs)          JOBS="$2";          shift 2 ;;
        --osqp)          pick; DO_OSQP=1;     shift ;;
        --piqp)          pick; DO_PIQP=1;     shift ;;
        --clarabel)      pick; DO_CLARABEL=1; shift ;;
        --piqp-only)     pick; DO_PIQP=1;     shift ;;
        --clarabel-only) pick; DO_CLARABEL=1; shift ;;
        --skip-deps)     SKIP_DEPS=1;        shift ;;
        --sudo)          USE_SUDO=1;         shift ;;
        --yes|-y)        ASSUME_YES=1;       shift ;;
        --help|-h)       usage ;;
        *)               die "unknown option '$1' (try --help)" ;;
    esac
done

[ "$DO_OSQP" = "1" ] || [ "$DO_PIQP" = "1" ] || [ "$DO_CLARABEL" = "1" ] \
    || die "no backend selected"

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

say "QP backends for PSOPT"
info "platform        $OS${PKGMGR:+ ($PKGMGR)}"
info "sources         $SRCDIR"
info "install prefix  $PREFIX"
info "parallel jobs   $JOBS"
info "building        $([ "$DO_OSQP" = 1 ] && printf 'OSQP ')$([ "$DO_PIQP" = 1 ] && printf 'PIQP ')$([ "$DO_CLARABEL" = 1 ] && printf 'Clarabel')"

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
    command -v cargo >/dev/null 2>&1 || die "cargo still is not on the PATH; open a new shell and rerun, or build without --clarabel"

    # Being on the PATH is not the same as working, and the first cargo on the PATH is
    # not necessarily the one that works. rustup installs shims named cargo and rustc
    # that refuse to do anything until a default toolchain has been chosen:
    #
    #     error: rustup could not choose a version of cargo to run, because one wasn't
    #     specified explicitly, and no default is configured.
    #
    # A machine can carry those shims and a working cargo from a package manager at the
    # same time, in either order. Clarabel's CMake runs a bare `cargo build`, so what
    # matters is which one the PATH resolves -- and the only lever is the PATH. So every
    # cargo on it is tried in turn and the first that answers --version is the one used,
    # with its directory moved to the front so that the build resolves the same one.
    CARGO_WORKING=""
    CARGO_VERSION=""
    CARGO_REFUSED=""
    saved_ifs="$IFS"; IFS=":"
    for d in $PATH; do
        [ -n "$d" ] || continue
        [ -x "$d/cargo" ] || continue
        if v="$("$d/cargo" --version 2>&1)"; then
            CARGO_WORKING="$d/cargo"; CARGO_VERSION="$v"; break
        fi
        [ -n "$CARGO_REFUSED" ] || CARGO_REFUSED="$d/cargo
$v"
    done
    IFS="$saved_ifs"

    if [ -z "$CARGO_WORKING" ]; then
        die "every cargo on the PATH refuses to run. The first was:

    ${CARGO_REFUSED%%
*}
    ${CARGO_REFUSED#*
}

    If that is rustup saying no default toolchain is configured, either give it one:

        rustup default stable

    or install cargo from your package manager -- MacPorts: sudo port install cargo;
    Homebrew: brew install rust -- and rerun."
    fi

    # The build resolves cargo from the PATH itself, so the working one has to be in
    # front. On a machine with only rustup this changes nothing.
    CARGO_WORKING_DIR="$(dirname "$CARGO_WORKING")"
    case ":$PATH:" in
        "$CARGO_WORKING_DIR:"*) ;;
        *) PATH="$CARGO_WORKING_DIR:$PATH"; export PATH ;;
    esac
    info "cargo           ${CARGO_VERSION#cargo } ($CARGO_WORKING)"

    # Clarabel.cpp's build generates its C header by running cbindgen, which it installs
    # with `cargo install` and then calls by name. cargo install puts its binaries in
    # $CARGO_HOME/bin, and rustup adds that directory to the shell profile -- so on a
    # machine whose cargo came from rustup this is invisible and everything works. Where
    # cargo came from a package manager instead, MacPorts or Homebrew, nothing has put
    # that directory on the PATH, and the build gets most of the way through and then
    # fails with
    #
    #     Ignored package `cbindgen v0.29.4` is already installed
    #     /bin/sh: cbindgen: command not found
    #
    # which says both that the tool is installed and that it cannot be found, and does
    # not mention the PATH at all. So the directory goes on the PATH here, and cbindgen
    # is checked for before anything is built rather than in the middle of it.
    # Appended, not prepended, and the distinction matters. All that is wanted from this
    # directory is cbindgen; the cargo that has already been found is the one to keep
    # using. A machine can carry rustup's shims here *and* a working cargo from a package
    # manager, and putting this directory first replaces the cargo that works with a shim
    # that refuses to run until a default toolchain is chosen -- which is a worse failure
    # than the one being fixed, and one this script caused on a Mac where cargo came from
    # MacPorts.
    CARGO_BIN="${CARGO_HOME:-$HOME/.cargo}/bin"
    case ":$PATH:" in
        *":$CARGO_BIN:"*) ;;
        *) if [ -d "$CARGO_BIN" ]; then
               PATH="$PATH:$CARGO_BIN"; export PATH
               info "PATH            $CARGO_BIN appended (cargo install puts its binaries there)"
           fi ;;
    esac

    if ! command -v cbindgen >/dev/null 2>&1; then
        info "cbindgen        not found; installing it (Clarabel's build generates its C header with it)"
        cargo install cbindgen
        hash -r 2>/dev/null || true
    fi
    command -v cbindgen >/dev/null 2>&1 || die "cbindgen is not on the PATH even after installing it.
    cargo install puts it in $CARGO_BIN; put that directory on your PATH and rerun:

        export PATH=\"$CARGO_BIN:\$PATH\""
    info "cbindgen        $(cbindgen --version 2>/dev/null | awk '{print $2}') ($(command -v cbindgen))"
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
# OSQP
# ---------------------------------------------------------------------------------
if [ "$DO_OSQP" = "1" ]; then
    say "OSQP"
    OSQP_SRC="$SRCDIR/osqp"
    clone_or_update "$OSQP_REPO" "$OSQP_SRC" "$OSQP_REF" 1

    # Only the static library: PSOPT's plugin links osqp::osqpstatic, and a shared one
    # would be an installed file that nothing loads. Everything else is left at OSQP's
    # own defaults, which is the configuration every OSQP measurement in this repository
    # was made against.
    cmake -S "$OSQP_SRC" -B "$OSQP_SRC/build" \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX="$PREFIX" \
          -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
          -DOSQP_BUILD_STATIC_LIB=ON \
          -DOSQP_BUILD_SHARED_LIB=OFF

    say "Building OSQP"
    cmake --build "$OSQP_SRC/build" -j "$JOBS"

    info "installing into $PREFIX"
    "${RUN_INSTALL[@]}" "$OSQP_SRC/build" >/dev/null

    # OSQP 1.x puts its headers in <prefix>/include/osqp, which is what its exported
    # cmake target points at; 0.6.x put them a directory up. The plugin tries both
    # spellings, so the layout is not the thing to check -- the API is.
    OSQP_HEADER=""
    for h in "$PREFIX/include/osqp/osqp.h" "$PREFIX/include/osqp.h"; do
        [ -f "$h" ] && OSQP_HEADER="$h" && break
    done
    [ -n "$OSQP_HEADER" ] || die "osqp.h was not installed under $PREFIX/include. Read the log above."

    OSQP_API="$(dirname "$OSQP_HEADER")/osqp_api_functions.h"
    if ! grep -q "OSQPCscMatrix_new" "$OSQP_API" 2>/dev/null; then
        die "the OSQP installed in $PREFIX is not a 1.x release: its headers do not
    declare OSQPCscMatrix_new, which PSOPT's plugin calls. PSOPT needs OSQP 1.0 or
    later. Rerun with --osqp-ref v1.0.0, or point --prefix somewhere that does not
    already hold an older OSQP."
    fi

    OSQP_LIB=""
    for d in "$PREFIX/lib" "$PREFIX/lib64"; do
        [ -f "$d/libosqpstatic.a" ] && OSQP_LIB="$d/libosqpstatic.a" && break
    done
    [ -n "$OSQP_LIB" ] || die "libosqpstatic.a was not installed under $PREFIX/lib. PSOPT's plugin links osqp::osqpstatic."

    OSQP_CONFIG="$(find "$PREFIX" -name 'osqp-config.cmake' -print -quit 2>/dev/null || true)"
    [ -n "$OSQP_CONFIG" ] || die "osqp-config.cmake was not installed; PSOPT's find_package(osqp) will not see it."
    info "header          $OSQP_HEADER"
    info "library         $OSQP_LIB"
    info "cmake config    $OSQP_CONFIG"
fi

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
    echo "# QP backend environment for PSOPT.  Source this, or copy it into your"
    echo "# shell profile:"
    echo "#     source $ENVFILE"
    echo
    if [ "$DO_OSQP" = "1" ]; then
        echo "# Where PSOPT's find_package(osqp) should look."
        echo "export osqp_DIR=\"$(dirname "$OSQP_CONFIG")\""
        echo
    fi
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
    echo "# None of them needs anything on the loader path: PIQP is header-only, and the OSQP"
    echo "# and Clarabel plugins link static libraries."
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
$([ "$DO_OSQP" = 1 ] && printf '                  -DWITH_OSQP=ON \\\\\n')\
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
            algorithm.qp_solver  = "OSQP";        // or "PIQP", "Clarabel", "GALAHAD"

    All three need the subproblem's Hessian to be positive semidefinite, which an optimal
    control problem's rarely is: OSQP is an operator-splitting method and the other two
    are interior-point methods, and none of them can take indefinite curvature. PSOPT
    raises its shift until the backend accepts the model, so this is handled, but it is
    why GALAHAD's QPA -- which takes an indefinite Hessian as posed -- remains the
    default.

EOF
