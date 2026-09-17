# Installing PSOPT on Arch Linux and Manjaro

Covers **Arch Linux** and **Manjaro**. Both are in the weekly distribution
matrix, as `containers/arch.Dockerfile` and `containers/manjaro.Dockerfile`,
which are byte-for-byte identical below their `FROM` line — deliberately, so
that one passing where the other fails is itself a diagnosis. Manjaro is Arch
with updates held back a few weeks, so a failure on Arch alone says the cause
arrived in the last fortnight.

Arch packages Eigen but neither IPOPT nor CppAD in the official repositories.
There are two ways round that and both are legitimate:

- **From the AUR**, which is quicker and is what the published `psopt-ci`
  Docker image uses.
- **From source**, which is what the weekly CI images do, and is therefore the
  route that is actually verified every week.

## 1. Install what is packaged

```
sudo pacman -Syu
sudo pacman -S --needed base-devel gcc-fortran cmake git pkgconf \
    eigen blas lapack gnuplot wget unzip python
```

`base-devel` supplies make, patch, file, which, findutils and diffutils along
with the compiler, which is why this list is shorter than the others.

Arch and Manjaro ship Eigen 5.0.1 where most distributions are still on 3.4.x.
PSOPT builds and passes its tests against both.

## 2a. IPOPT and CppAD from the AUR

```
sudo pacman -S --needed yay
yay -S coin-or-ipopt cppad
```

`yay` needs AUR support enabled. On ARM64 the AUR IPOPT build has been known to
fail; [Anaconda](https://www.anaconda.com/download) provides a working IPOPT
there, or use the source route below.

## 2b. IPOPT and CppAD from source

This is what CI does, and it gives IPOPT 3.14.19.

```
git clone https://github.com/coin-or/coinbrew ~/coinbrew
cd ~/coinbrew
./coinbrew fetch Ipopt --no-prompt
./coinbrew build Ipopt --prefix=$HOME/coin/dist --no-prompt \
      ADD_FFLAGS=-fallow-argument-mismatch
export PKG_CONFIG_PATH=$HOME/coin/dist/lib/pkgconfig:$PKG_CONFIG_PATH
```

```
git clone --depth 1 https://github.com/coin-or/CppAD.git ~/cppad.git
cd ~/cppad.git
cmake -B build -D cppad_prefix=/usr/local
cmake --build build -j
sudo cmake --install build
```

Add the `PKG_CONFIG_PATH` line to your shell profile so that later reconfigures
keep finding this IPOPT.

## 3. Build and install PSOPT

```
git clone https://github.com/PSOPT/psopt.git
cd psopt
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=ON
cmake --build build -j
sudo cmake --install build
```

If you built IPOPT into `$HOME/coin/dist`, add
`-DCMAKE_PREFIX_PATH=$HOME/coin/dist` to the configure line.

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
whole story, including whether you want it at all. What Arch needs beyond the
list above is only the tooling its backends are built with:

```
sudo pacman -S --needed meson ninja curl
```

Then `scripts/build_galahad.sh` and `scripts/build_qp_backends.sh` install the
backends themselves. `curl` is there for `rustup`, which only the Clarabel
backend needs.

## Building the container images yourself

If you build `containers/arch.Dockerfile` or `containers/manjaro.Dockerfile`
locally, note that both set `DisableSandbox` in `pacman.conf`. pacman 7 fetches
packages as an unprivileged user behind a seccomp and landlock filter, which a
container cannot install, and the failure — `error restricting syscalls via
seccomp: 22` — happens at the first `pacman -Syu`. That is a property of pacman
inside a container and has nothing to do with installing PSOPT on a real Arch
machine.
