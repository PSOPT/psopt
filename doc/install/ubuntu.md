# Installing PSOPT on Ubuntu

Covers **Ubuntu 26.04 LTS** and **Ubuntu 24.04 LTS**. Both are in the weekly
distribution matrix, as `containers/ubuntu-26.04.Dockerfile` and
`containers/ubuntu-24.04.Dockerfile`, which install exactly the packages below
and then build PSOPT, run the tests, install it and build a separate program
against the installed package.

Ubuntu is the easiest of the platforms: everything PSOPT needs is packaged,
including IPOPT and CppAD, so there is nothing to build from source.

## 1. Install the dependencies

```
sudo apt-get update
sudo apt-get install -y build-essential gfortran cmake git pkg-config \
    coinor-libipopt-dev libeigen3-dev libcppad-dev \
    libblas-dev liblapack-dev gnuplot
```

`build-essential` supplies gcc, g++ and make; `gfortran` is wanted by IPOPT's
MUMPS. `gnuplot` is optional and only affects PSOPT's plotting helpers.

The IPOPT you get differs between releases, and it is worth knowing which,
because a different interior-point path reaches a different mesh and can move an
example's cost in the last few digits. As last measured in CI:

| release | IPOPT | Eigen |
|---|---|---|
| Ubuntu 26.04 LTS | 3.14.19 | 3.4.0 |
| Ubuntu 24.04 LTS | 3.11.9 | 3.4.0 |

Ubuntu 24.04's IPOPT is from 2013 and is the oldest in the matrix. PSOPT builds
and passes against it, which is the point of keeping 24.04 in the matrix, but if
you have a choice, a newer IPOPT is a better solver. Building a current one from
source is the procedure on the [openSUSE page](opensuse.md), which applies
unchanged here.

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
whole story, including whether you want it at all. What Ubuntu needs beyond the
list above is only the tooling its backends are built with:

```
sudo apt-get install -y meson ninja-build curl patch wget unzip file
```

Then `scripts/build_galahad.sh` and `scripts/build_qp_backends.sh` install the
backends themselves. `curl` is there for `rustup`, which only the Clarabel
backend needs.
