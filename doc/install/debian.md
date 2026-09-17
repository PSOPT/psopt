# Installing PSOPT on Debian

Covers **Debian 13**. It is in the weekly distribution matrix as
`containers/debian-13.Dockerfile`, which installs exactly the packages below and
then builds PSOPT, runs the tests, installs it and builds a separate program
against the installed package.

Debian packages everything PSOPT needs, IPOPT and CppAD included, so there is
nothing to build from source. The commands are the Ubuntu ones; the reason
Debian has a page of its own is that its IPOPT is usually a different version,
which is also why both are in the matrix.

## 1. Install the dependencies

As root, or with `sudo`:

```
apt-get update
apt-get install -y build-essential gfortran cmake git pkg-config \
    coinor-libipopt-dev libeigen3-dev libcppad-dev \
    libblas-dev liblapack-dev gnuplot
```

Debian 13 supplies IPOPT 3.14.17 and Eigen 3.4.0, as last measured in CI.

`build-essential` supplies gcc, g++ and make; `gfortran` is wanted by IPOPT's
MUMPS. `gnuplot` is optional and only affects PSOPT's plotting helpers.

## 2. Build and install PSOPT

```
git clone https://github.com/PSOPT/psopt.git
cd psopt
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=ON
cmake --build build -j
sudo cmake --install build
```

## 3. Check it works

```
cd build/examples/launch && ./launch
```

Optionally, configure with `-DBUILD_TESTS=ON` and run the unit tests:

```
ctest --test-dir build --output-on-failure
```

## If you also want PSOPT's own SQP solver

The SQP solver is off by default and needs a QP backend; the main README has the
whole story, including whether you want it at all. What Debian needs beyond the
list above is only the tooling its backends are built with:

```
apt-get install -y meson ninja-build curl patch wget unzip file
```

Then `scripts/build_galahad.sh` and `scripts/build_qp_backends.sh` install the
backends themselves. `curl` is there for `rustup`, which only the Clarabel
backend needs.

One Debian-specific note, if a link ever fails with an undefined `dmumps_c`:
`pkg-config --libs ipopt` lists `-ldmumps_seq` here, so the MUMPS library the
SQP solver reads inertia from resolves without help. That is not true of a
coinbrew-built IPOPT, which records the dependency inside `libipopt` instead.
