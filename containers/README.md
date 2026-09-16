# Distribution tests

Does PSOPT build, pass its tests, install, and get consumed by somebody else's
project, on the Linux distributions people actually run? That is the only
question these images answer, and it is not the question
`.github/workflows/ci.yml` answers. The continuous integration workflow runs
thirteen examples against reference costs on one platform, and it is right to
call a change in any of them a regression. Here a change in an example's cost
is usually not a regression at all: it is a different IPOPT version taking a
different interior-point path to a different mesh. Asking the same question
twice, once per platform, would produce a wall of red that means nothing.

## What one job does

Every job runs `common/build_and_test.sh`, which is the same file everywhere:

1. report the distribution, the compilers, glibc, CMake, Python, IPOPT and
   Eigen, before anything is built;
2. build CppAD from source;
3. configure and build PSOPT with `-DBUILD_EXAMPLES=ON -DBUILD_TESTS=ON
   -DHEADLESS=ON -DPSOPT_AD_BACKEND=CPPAD`;
4. run `ctest`;
5. install PSOPT;
6. configure, build and run a **separate** project that does
   `find_package(PSOPT REQUIRED)` and links against the installed package;
7. run three examples whose answers are fixed outside PSOPT.

Step 6 is there because nothing else reaches the install rules. Every example
and every test in this repository is built inside the source tree and inherits
its include paths from it, so the installed package can be broken for years
without a single job going red. It was: `find_package(PSOPT)` had never worked
until September 2026, in two separate ways, and both were found by writing this
step and not by any test.

## What is checked, and what deliberately is not

Only three examples run, and the reason for each is that its answer comes from
somewhere other than a previous PSOPT run:

| example | reference | where the number comes from |
|---|---|---|
| `mineng_di` | `6.000000e+00` | the minimum-energy double integrator has the closed-form solution `u* = 6 - 12t`, `J* = 6` exactly, and the transcription is exact on it: cubic states, linear control, quadratic integrand, all inside what a Legendre interpolant on 40 nodes represents and Gauss-Lobatto quadrature integrates without error |
| `lts_costates` | `4.195884e+02` | linear tangent steering, whose minimum time follows from the analytical solution of Bryson and Ho and is confirmed to nine figures by an independent quadrature |
| `geodesic` | `5.549682e+03` | the exact geodesic on this spheroid, by Karney's algorithm, computed outside PSOPT entirely |

A departure from one of these is unambiguous: arithmetic does not vary by
distribution. The other ten examples in `ci.yml` have reference values that
some earlier build produced, which is exactly right for detecting a regression
on a fixed platform and exactly wrong here.

`ctest` is the stronger check in this matrix and it runs in full. Its tests
check costates against closed-form adjoints, stationarity residuals, the
constancy of the Hamiltonian, constraint coverage and the rest, all against
values that do not depend on which IPOPT solved the problem.

## Where dependencies come from

Two of them are not the same everywhere, and the difference is a result rather
than a nuisance, so it is reported instead of being smoothed away.
`common/ensure_ipopt.sh` and `common/ensure_eigen.sh` each use the
distribution's package when there is one and build from source when there is
not, printing one of

```
IPOPT-SOURCE: distribution package
IPOPT-SOURCE: built from source (no usable distribution package found)
```

and the version with it. A job that merely died where a distribution has no
IPOPT package would report "PSOPT fails on openSUSE", when the truth is "not
packaged there, and everything works once it is built". Those two outcomes have
to look different in the log, so they do.

CppAD is built from source on every distribution, without exception. PSOPT
links both the CppAD headers and the compiled `libcppad_lib`, and a package
that ships only the headers satisfies the compiler and then fails at link time.
Building it identically everywhere also keeps this a test of PSOPT's
portability and not a test of nine distributions' CppAD packaging.

Google Test is downloaded by `tests/CMakeLists.txt` when the distribution does
not provide it, so no image installs it and `ctest` runs the same tests
everywhere.

## The images

| file | base image | notes |
|---|---|---|
| `ubuntu-26.04.Dockerfile` | `ubuntu:26.04` | current LTS; IPOPT packaged |
| `ubuntu-24.04.Dockerfile` | `ubuntu:24.04` | previous LTS, still supported; the older toolchain is where a use of something later than C++17 shows up |
| `debian-13.Dockerfile` | `debian:13` | Debian's IPOPT is usually a different version from Ubuntu's, which is the point of having both |
| `fedora-44.Dockerfile` | `fedora:44` | newest GCC in the Red Hat family; IPOPT packaged as `coin-or-Ipopt-devel` |
| `rhel-10.Dockerfile` | `redhat/ubi10` | see the caveat below |
| `opensuse-leap-16.Dockerfile` | `opensuse/leap:16.0` | same sources as SUSE Linux Enterprise; conservative package set |
| `opensuse-tumbleweed.Dockerfile` | `opensuse/tumbleweed` | rolling; an early warning about dependencies that will reach the LTS distributions later |
| `arch.Dockerfile` | `archlinux:latest` | the distribution the published `psopt-ci` image is built on |
| `manjaro.Dockerfile` | `manjarolinux/base:latest` | Arch with updates held back a few weeks; when Arch fails and Manjaro passes, the cause changed in the last fortnight |

**The RHEL caveat.** `redhat/ubi10` is RHEL 10: the same glibc, compilers and
runtime. What it is not is the same set of *repositories*. The UBI repository
set is deliberately narrow, so CodeReady Builder and EPEL have to be enabled and
may not be reachable. That Dockerfile therefore allows those steps to fail and
lets the `ensure_` scripts build what is missing. A green job here says nothing
in PSOPT depends on what the wider repository set adds; it does not prove that a
subscribed RHEL 10 installation behaves identically.

**IPOPT on Arch and Manjaro** is built from source here, although `coin-or-ipopt`
is in the AUR and the published CI image uses it. An AUR build needs an
unprivileged user, an AUR helper and a chain of `makepkg` steps whose outcome
depends on what the AUR looked like that morning, none of which is a property of
PSOPT. The AUR route keeps being exercised by the image published from the
`Dockerfile` at the top of this repository.

## Running one locally

From the top of the repository, not from this directory: the build context is
the whole repository, because the image copies the source tree it is about to
build.

```
docker build --progress=plain -f containers/ubuntu-26.04.Dockerfile .
```

The build succeeds only if everything above succeeded, so there is nothing to
run afterwards; the log is the result. To get a shell in a working image with
PSOPT already installed, add `-t psopt-ubuntu` and then `docker run --rm -it
psopt-ubuntu`.

Podman works in place of Docker throughout.

## Adding a distribution

Copy the Dockerfile that is closest in package manager, change the base image
and the package names, and add one entry to the matrix in
`.github/workflows/distros.yml`. Nothing else should need to change.

If a new distribution seems to need a change to `common/build_and_test.sh`,
stop. That file has one rule and it is stated at the top of it: everything that
differs between distributions lives in the Dockerfile. A condition on the
distribution inside the shared script means either that the difference belongs
in the Dockerfile, or that PSOPT has a portability defect that should be fixed
in PSOPT.

## What the first local run found

The nine images were tried on a Mac Studio before any of this reached GitHub,
in the dependency-only mode described at the end of this file. Two things came
out of it that are worth keeping.

**pacman 7 will not download inside a container build unless its sandbox is
turned off.** Arch failed at the first `pacman -Syu` with `error restricting
syscalls via seccomp: 22`, followed by `switching to sandbox user 'alpm'
failed`. pacman 7 fetches packages as an unprivileged user behind a seccomp and
landlock filter, 22 is EINVAL, and the container could not install the filter at
all. Both Arch-family Dockerfiles now set `DisableSandbox` in the `[options]`
section of `pacman.conf` before anything else. Manjaro did not hit this, because
its base image starts with an older pacman and upgrades to 7 part-way through
the transaction, which is a good illustration of why both images are in the
matrix.

**Arch and Manjaro ship Eigen 5.** The Manjaro image reported `eigen3 5.0.1`,
with `Eigen3Config.cmake` present, where every other distribution here is on
3.4.x and PSOPT has only ever been built against 3.4. `ensure_eigen.sh` will not
intervene, since `find_package(Eigen3)` succeeds. So the first full build on
either of those two images is also the first time PSOPT has met Eigen 5, and a
compile failure there is a real finding about PSOPT and not about the
distribution.

## How far the package names have been checked

None of these images has been built yet, because no Docker daemon was available
where they were written. The package names are at different levels of
confidence and it is worth being exact about which is which.

Every name in every image has now resolved at least once. The seven images that
build natively on an arm64 Mac were checked there, and the two Arch-family
images were checked under x86_64 emulation, where Manjaro installed the whole
list and Arch reached the same list once the pacman sandbox was turned off.

That is a real result but it is not the same as the workflow passing. It was
measured on one machine, on one day, with two of the nine emulated, and a rolling
distribution can rename a package next week. The Red Hat UBI repository
identifiers remain the least certain thing in the set, since `rhel-10` is written
so that CodeReady Builder and EPEL may fail without failing the build.

So if a job in the first GitHub run dies inside `dnf install` or `pacman -S`,
that is still a statement about package names and not about PSOPT. Fix the name
and run it again. Only once a job has reached the configure step does what it
reports become a statement about the software.

One result is already worth noting in advance. Ubuntu 24.04 packages IPOPT
3.11.9, released in 2014, while Fedora 44 packages 3.14.16 and the source build
here uses 3.14.19. That eight-year spread across distributions is not a defect
in any of them and it is the main reason this matrix reports which IPOPT each
job used before it reports anything else.
