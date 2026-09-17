# Installing PSOPT on openSUSE

Covers **openSUSE Leap 16.0** and **openSUSE Tumbleweed**. Both are in the
weekly distribution matrix, as `containers/opensuse-leap-16.Dockerfile` and
`containers/opensuse-tumbleweed.Dockerfile`, which install exactly the packages
below and then build PSOPT, run the tests, install it and build a separate
program against the installed package.

openSUSE is the longest of the Linux routes, because **it has no usable IPOPT
package**: both Leap 16.0 and Tumbleweed build IPOPT from source in CI, and so
must you. That also makes this page the one to read on any distribution that
does not package IPOPT, whatever its package manager — the coinbrew steps below
are not specific to openSUSE.

## 1. Install the dependencies

```
sudo zypper install -y gcc gcc-c++ gcc-fortran make cmake git pkg-config \
    eigen3-devel blas-devel lapack-devel gnuplot \
    findutils diffutils patch wget unzip file which tar gzip
```

`findutils` and `diffutils` are in that list deliberately. Neither openSUSE base
image carries them, and their absence is not loud: IPOPT's configure reports
`cmp: command not found` and carries on to produce a working IPOPT that has
taken a different path, and later steps that call `find` read "no such program"
as "no such file". They cost nothing to install and they closed a failure that
took three rounds to diagnose.

## 2. Build IPOPT from source

[coinbrew](https://github.com/coin-or/coinbrew) fetches IPOPT with MUMPS and
builds both:

```
git clone https://github.com/coin-or/coinbrew ~/coinbrew
cd ~/coinbrew
./coinbrew fetch Ipopt --no-prompt
./coinbrew build Ipopt --prefix=$HOME/coin/dist --no-prompt \
      ADD_FFLAGS=-fallow-argument-mismatch
```

`ADD_FFLAGS=-fallow-argument-mismatch` lets recent gfortran accept MUMPS's
legacy Fortran. That build gives IPOPT 3.14.19, which is what CI measures on
both openSUSE images.

Then put its `pkg-config` file where PSOPT will look, and add the line to your
shell profile so that later reconfigures keep finding it:

```
export PKG_CONFIG_PATH=$HOME/coin/dist/lib/pkgconfig:$PKG_CONFIG_PATH
pkg-config --modversion ipopt        # should print 3.14.19
```

Installing into `/usr/local` instead of `$HOME/coin/dist` works equally well;
name both `lib` and `lib64` on `PKG_CONFIG_PATH` if you do, since a source build
lands in `lib` while openSUSE's own libraries are in `lib64`.

## 3. Build CppAD from source

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

## 4. Build and install PSOPT

```
git clone https://github.com/PSOPT/psopt.git
cd psopt
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=ON \
      -DCMAKE_PREFIX_PATH=$HOME/coin/dist
cmake --build build -j
sudo cmake --install build
```

Drop `-DCMAKE_PREFIX_PATH` if you installed IPOPT into `/usr/local`.

## 5. Check it works

```
cd build/examples/launch && ./launch
```

Optionally, configure with `-DBUILD_TESTS=ON` and run the unit tests:

```
ctest --test-dir build --output-on-failure
```

If an example fails to start because it cannot find `libipopt`, the library is
in a prefix the loader does not search:

```
export LD_LIBRARY_PATH=$HOME/coin/dist/lib:$LD_LIBRARY_PATH
```

## If you also want PSOPT's own SQP solver

The SQP solver is off by default and needs a QP backend; the main README has the
whole story, including whether you want it at all. What openSUSE needs beyond
the list above is only the tooling its backends are built with:

```
sudo zypper install -y meson ninja curl
```

Then `scripts/build_galahad.sh` and `scripts/build_qp_backends.sh` install the
backends themselves. Both take `--skip-deps`, which is worth knowing here: left
to themselves they try to install their own prerequisites through a package
manager they detect, and neither knows `zypper`, so on openSUSE you install the
two packages above yourself and pass `--skip-deps`. `curl` is there for
`rustup`, which only the Clarabel backend needs.
