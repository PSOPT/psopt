# Installing PSOPT

Pick your platform:

| platform | page | tested weekly by |
|---|---|---|
| Ubuntu 26.04 LTS, Ubuntu 24.04 LTS | [ubuntu.md](ubuntu.md) | `containers/ubuntu-26.04.Dockerfile`, `containers/ubuntu-24.04.Dockerfile` |
| Debian 13 | [debian.md](debian.md) | `containers/debian-13.Dockerfile` |
| Fedora 44 | [fedora.md](fedora.md) | `containers/fedora-44.Dockerfile` |
| openSUSE Leap 16.0, Tumbleweed | [opensuse.md](opensuse.md) | `containers/opensuse-leap-16.Dockerfile`, `containers/opensuse-tumbleweed.Dockerfile` |
| Arch Linux, Manjaro | [arch.md](arch.md) | `containers/arch.Dockerfile`, `containers/manjaro.Dockerfile` |
| macOS (Apple Silicon and Intel) | [macos.md](macos.md) | not automated; see the page |

Distributions that share a package manager share a page, because they share the
commands: Ubuntu 24.04 and 26.04 differ only in which IPOPT they ship, Arch and
Manjaro not at all.

**Every Linux page has an executable counterpart.** The Dockerfile named beside
it installs the same packages and then builds PSOPT, runs the tests, installs
it, and builds a separate program against the installed package — every week, on
GitHub's runners, and on every push. If a page and its Dockerfile ever disagree,
the Dockerfile is the one that has been run. `containers/README.md` describes
what those jobs do and what they deliberately do not check.

There is no page for a distribution that is not in that matrix. The nearest one
by package manager is usually right, and adding an image is a small job that
`containers/README.md` describes.

## What PSOPT needs

Three libraries, and only the first is likely to cost you any effort:

**IPOPT** is the nonlinear programming solver PSOPT uses by default. Where a
distribution packages it, use the package; where it does not, it is built from
source with [coinbrew](https://github.com/coin-or/coinbrew), which takes a few
minutes and is what the platform pages show. IPOPT needs a Fortran compiler for
MUMPS, its default linear solver. PSOPT has been built and tested against IPOPT
releases from 3.11.9 to 3.14.19.

**Eigen 3 or 5** is a header-only linear algebra library. Every distribution in
the matrix packages it, and PSOPT builds and passes its tests against both the
3.4 series and Eigen 5.0.1. Nothing has to be done beyond installing the
package.

**CppAD** supplies the automatic differentiation. It is *not* header-only:
PSOPT links `libcppad_lib` as well as including the headers, so a package that
ships only headers will satisfy the compiler and then fail at link time. On
Debian and Ubuntu the `libcppad-dev` package carries both and is the easy route;
elsewhere the platform pages build it from source, which takes about a minute.

Optionally, **GNUplot** for the plotting helpers, and **BLAS** and **LAPACK**,
which IPOPT wants and which every distribution packages.

Boost and ColPack were dependencies of older PSOPT releases and are not needed
by any current build. If you are following instructions that mention them, those
instructions are out of date.

## Building and installing PSOPT

Identical on every platform once the dependencies are in place:

```
git clone https://github.com/PSOPT/psopt.git
cd psopt
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=ON
cmake --build build -j
sudo cmake --install build
```

For a debug build, `-DCMAKE_BUILD_TYPE=Debug`. On a machine with no display, or
in a container, add `-DHEADLESS=ON` so that the plotting helpers write files
instead of trying to open a window.

Then run an example, which is the only check that matters:

```
cd build/examples/launch && ./launch
```

To run the unit tests as well, configure with `-DBUILD_TESTS=ON` and run
`ctest --test-dir build --output-on-failure`. They check costates against
closed-form adjoints, stationarity residuals, the constancy of the Hamiltonian
and much else, and they are the strongest evidence that a build is sound.

## If the configure step cannot find something

**`pkg-config` cannot find IPOPT.** PSOPT locates IPOPT through `pkg-config`, so
`ipopt.pc` must be on `PKG_CONFIG_PATH`. A source build of IPOPT puts it under
its own prefix, which is why the pages that build IPOPT set that variable. Check
with `pkg-config --modversion ipopt` before configuring PSOPT; if that prints a
version, PSOPT will find it.

A source build of anything installs into `<prefix>/lib` even on distributions
whose own libraries live in `<prefix>/lib64`, and neither `pkg-config` nor the
dynamic loader necessarily looks in both. Naming both costs nothing:

```
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:/usr/local/lib64/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH
```

If your IPOPT genuinely has no `.pc` file, one can be written by hand and placed
in a directory on `PKG_CONFIG_PATH`. The paths in it depend on where IPOPT was
installed:

```
prefix=/usr/local
exec_prefix=${prefix}
libdir=${exec_prefix}/lib
includedir=${prefix}/include/coin-or
Name: IPOPT
Description: Interior Point Optimizer
URL: https://github.com/coin-or/Ipopt
Version: 3.14.19
Cflags: -I${includedir}
Libs: -L${libdir} -lipopt
```

**CppAD not found.** The message names both halves, because both are needed:

```
PSOPT_AD_BACKEND=CPPAD but CppAD not found. Set -DCPPAD_INCLUDE_DIR=<dir with
cppad/cppad.hpp> and -DCPPAD_LIBRARY=<libcppad_lib>, or CPPAD_DIR.
```

If you installed CppAD to a prefix of your own, `export CPPAD_DIR=/your/prefix`
before configuring, or name the two paths directly.

**Eigen not found.** `find_package(Eigen3 NO_MODULE)` reads
`Eigen3Config.cmake`, which the distribution packages install. A source build of
Eigen installs it too; if it is somewhere unusual, add that prefix to
`CMAKE_PREFIX_PATH`.

**A linker error mentioning `dmumps_c`** belongs to PSOPT's own SQP solver, not
to an ordinary build. The main README's section on the SQP solver explains it.

## Other ways to install

A **Docker container** avoids installing anything on the host; the main README
describes the shipped `Dockerfile`. The images under `containers/` are built for
testing rather than for use, but any of them is also a working recipe, and each
one's first stanza is the package list for its distribution.

A **Python interface** is available for building models without writing C++; see
the user manual.
