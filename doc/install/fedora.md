# Installing PSOPT on Fedora

Covers **Fedora 44**. It is in the weekly distribution matrix as
`containers/fedora-44.Dockerfile`, which installs exactly the packages below and
then builds PSOPT, runs the tests, installs it and builds a separate program
against the installed package.

Fedora packages IPOPT and Eigen but not CppAD, so CppAD is built from source
here — about a minute's work. Fedora is also where a new compiler and a new
Eigen reach users first: it is usually a release of GCC ahead of the LTS
distributions, and it ships Eigen 5.0.1 where most others are still on 3.4.x.
PSOPT builds and passes its tests against both.

## 1. Install the dependencies

```
sudo dnf install -y gcc gcc-c++ gcc-gfortran make cmake git \
    pkgconf-pkg-config coin-or-Ipopt-devel MUMPS-devel \
    eigen3-devel blas-devel lapack-devel gnuplot \
    diffutils which file patch
```

Two of those are easy to miss and both have cost this project time:

- **`MUMPS-devel`** is separate from `coin-or-Ipopt-devel` on Fedora, which does
  not carry the MUMPS headers. Without it the PSOPT configure step cannot find
  `dmumps_c.h`.
- **`diffutils`** is absent from the Fedora base image. Nothing in an ordinary
  Fedora build needs `cmp` or `diff`, since IPOPT comes from a package, but a
  source build of IPOPT does, and its configure scripts carry on silently
  without them.

Fedora 44 supplies IPOPT 3.14.16 and Eigen 5.0.1, as last measured in CI.

## 2. Build CppAD from source

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

## 3. Build and install PSOPT

```
git clone https://github.com/PSOPT/psopt.git
cd psopt
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=ON
cmake --build build -j
sudo cmake --install build
```

If the configure step cannot find CppAD or IPOPT, note that a source build
installs into `/usr/local/lib` while Fedora's own libraries are in
`/usr/lib64`, and that neither `pkg-config` nor the loader necessarily looks in
both:

```
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:/usr/local/lib64/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH
```

## 4. Check it works

```
cd build/examples/launch && ./launch
```

Optionally, configure with `-DBUILD_TESTS=ON` and run the unit tests:

```
ctest --test-dir build --output-on-failure
```

## If you also want PSOPT's own SQP solver

The SQP solver is off by default and needs a QP backend; the main README has the
whole story, including whether you want it at all. What Fedora needs beyond the
list above is only the tooling its backends are built with:

```
sudo dnf install -y meson ninja-build curl wget unzip
```

Then `scripts/build_galahad.sh` and `scripts/build_qp_backends.sh` install the
backends themselves. `curl` is there for `rustup`, which only the Clarabel
backend needs.
