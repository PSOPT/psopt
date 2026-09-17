# Installing PSOPT on macOS

Covers **Apple Silicon** (M1 to M4) and **Intel** Macs, using
[MacPorts](https://www.macports.org/install.php) for most dependencies. The
procedure below has been used successfully on an M2 Max, an M4 Pro and on Intel
hardware, most recently on macOS Tahoe 26.4.1.

macOS is not in the weekly container matrix — there is no macOS container — so
this page is maintained by hand and is the one place in the installation
documentation with no automated counterpart. If something here has gone stale,
that is why.

**Do not install IPOPT from MacPorts.** On Apple Silicon the MacPorts `ipopt`
port is built against a *parallel* (MPICH) build of MUMPS which calls
`MPI_Init` at library-load time and crashes as soon as an example is run.
Build IPOPT and MUMPS yourself, as a sequential solver, as step 3 describes.

## 1. Install MacPorts

From https://www.macports.org/install.php.

## 2. Install the dependencies

```
sudo port install cmake eigen3 git gnuplot pkgconfig
sudo port install gcc15        # for gfortran, at /opt/local/bin/gfortran-mp-15
```

The MacPorts `ipopt` port is deliberately omitted; it is built in step 3.

`gcc15` is wanted only for its Fortran compiler, which MUMPS needs. It installs
as `/opt/local/bin/gfortran-mp-15`; if you install a different GCC version,
adjust the `-mp-NN` suffix in step 3 to match.

## 3. Build IPOPT and MUMPS with coinbrew

```
git clone https://github.com/coin-or/coinbrew ~/coinbrew
cd ~/coinbrew
./coinbrew fetch Ipopt --no-prompt

export CC=/usr/bin/clang
export CXX=/usr/bin/clang++
export FC=/opt/local/bin/gfortran-mp-15

./coinbrew build Ipopt --prefix=$HOME/coin/dist --no-prompt \
      ADD_FFLAGS=-fallow-argument-mismatch
```

Why each of those settings matters:

- **`CC` and `CXX` set to Apple clang** make IPOPT use the `libc++` C++ standard
  library, matching PSOPT and the MacPorts libraries. Building IPOPT with the
  MacPorts `g++` instead links `libstdc++`, whose `std::string` is binary
  incompatible with `libc++`, and PSOPT then segfaults as soon as it passes
  options to IPOPT.
- **`FC` set to gfortran** compiles MUMPS, and
  `ADD_FFLAGS=-fallow-argument-mismatch` lets recent gfortran accept MUMPS's
  legacy Fortran.
- coinbrew builds MUMPS with its **sequential MPI stub**, so there is no MPICH
  and no load-time `MPI_Init` — the root cause of the MacPorts crash.
- Apple's **Accelerate** framework is detected automatically and used as a fast
  BLAS and LAPACK, which is excellent on Apple Silicon. No flag is needed.

If the build stops at IPOPT's **Java** unit test, that is harmless: it only
fails when the system `java` is an Intel JVM that cannot load an arm64 library.
Finish the install by hand:

```
cd ~/coinbrew/build/Ipopt/*/ && make install
```

## 4. Verify that IPOPT is sequential and on the right C++ library

```
otool -L ~/coin/dist/lib/libipopt.3.dylib | grep -iE 'mpi|c\+\+|stdc'
```

You should see `/usr/lib/libc++.1.dylib` and **no** `libmpi`, `libpmpi` or
`libstdc++`. Both halves matter, and both correspond to a failure described in
step 3. It is worth running: each of those failures shows up much later, as a
crash rather than as a build error.

## 5. Build CppAD

```
git clone --depth 1 https://github.com/coin-or/CppAD.git cppad.git
cd cppad.git
cmake -B build -D cppad_prefix=/usr/local
cmake --build build -j
sudo cmake --install build
cd ..
```

That puts the headers in `/usr/local/include/cppad/` and the library in
`/usr/local/lib/`, both of which PSOPT's CMake finds on the default search path.
If you install it to a prefix of your own, `export CPPAD_DIR=/your/prefix`
before configuring PSOPT.

## 6. Build PSOPT against your IPOPT

```
export PKG_CONFIG_PATH=$HOME/coin/dist/lib/pkgconfig:$PKG_CONFIG_PATH

git clone https://github.com/PSOPT/psopt.git
cd psopt
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=ON \
      -DCMAKE_PREFIX_PATH=$HOME/coin/dist \
      -DCMAKE_BUILD_RPATH=$HOME/coin/dist/lib \
      -DCMAKE_INSTALL_RPATH=$HOME/coin/dist/lib
cmake --build build -j
```

Add the `PKG_CONFIG_PATH` line to your `~/.zshrc`, placing
`~/coin/dist/lib/pkgconfig` **before** `/opt/local/lib/pkgconfig`, so that later
reconfigures keep finding this IPOPT rather than the MacPorts one.

## 7. Run an example

```
cd build/examples/launch && ./launch
```

Optionally, configure with `-DBUILD_TESTS=ON` and run the unit tests:

```
ctest --test-dir build --output-on-failure
```

## If you also want PSOPT's own SQP solver

The SQP solver is off by default and needs a QP backend; the main README has the
whole story, including whether you want it at all. Two things are specific to
macOS.

`scripts/build_galahad.sh` works with either MacPorts or Homebrew. On MacPorts
it also runs `port select` so that a plain `gfortran` exists, since MacPorts
installs the compiler as `gfortran-mp-NN` and meson looks for the plain name.

The SQP reads the inertia of the KKT matrix from MUMPS, and macOS will not
resolve a symbol through an indirect dylib. A coinbrew IPOPT records MUMPS
inside `libipopt`, so if a link fails with an undefined `dmumps_c`, point
`MUMPS_LIBRARY` at the library that holds it — `libcoinmumps` for a coinbrew
build. Since you built IPOPT and MUMPS into `~/coin/dist` above, the
`CMAKE_PREFIX_PATH` in step 6 is normally enough to find both the header and the
library, and there is nothing further to do.
