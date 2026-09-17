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
3. install the QP backends the SQP solver needs, through PSOPT's own
   `scripts/build_qp_backends.sh` and `scripts/build_galahad.sh`;
4. configure and build PSOPT with `-DBUILD_EXAMPLES=ON -DBUILD_TESTS=ON
   -DHEADLESS=ON -DPSOPT_AD_BACKEND=CPPAD`, plus `-DWITH_SQP=ON` and whichever
   backends step 3 actually produced;
5. run `ctest`;
6. install PSOPT;
7. configure, build and run a **separate** project that does
   `find_package(PSOPT REQUIRED)` and links against the installed package;
8. run three examples whose answers are fixed outside PSOPT.

Step 7 is there because nothing else reaches the install rules. Every example
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
portability and not a test of eight distributions' CppAD packaging.

Google Test is downloaded by `tests/CMakeLists.txt` when the distribution does
not provide it, so no image installs it and `ctest` runs the same tests
everywhere.

**The QP backends are built through PSOPT's own installer scripts**, not through
anything written for the containers. `scripts/build_qp_backends.sh` does OSQP,
PIQP and Clarabel; `scripts/build_galahad.sh` does GALAHAD. Those are the scripts
the manual tells a user to run, so the matrix tests them as well, and a job here
says they have stopped working on some distribution before a user does. Both are
given `--skip-deps`: left alone they would `apt-get` or `dnf` their own
prerequisites through `sudo`, which a container has no reason to allow, and
neither knows `zypper`, so openSUSE would fall through their detection anyway.
The prerequisites are the Dockerfile's job, like every other distribution
difference.

This matters more than it sounds. Until it was added, `WITH_SQP` was off in every
image, `ctest` ran two tests rather than three, and PSOPT's own solver, as
opposed to IPOPT, was tested on no distribution at all. The thirteen `SQPSolver`
tests and the second `ctest` pass that GALAHAD's OpenMP requirement needs now run
everywhere.

Rust, which only Clarabel needs, comes from `rustup` rather than a distribution
package, for the same reason CppAD is built from source: a distribution's Rust
can be old enough for Clarabel.cpp to refuse it, and one current toolchain
everywhere keeps this a test of PSOPT instead of a test of eight distributions'
Rust packaging. `ensure_qp_backends.sh` prefers a packaged `cargo` if a
Dockerfile ever installs one, and says which it used.

Not every backend need succeed. `ensure_qp_backends.sh` reports each as
`QP-BACKEND <name>: <path or NOT BUILT>`, writes the cmake arguments that follow
from what is actually present, and fails only when none of the four was built,
because at that point `WITH_SQP` cannot be configured and the SQP would go
untested while the job stayed green.

## The images

| file | base image | notes |
|---|---|---|
| `ubuntu-26.04.Dockerfile` | `ubuntu:26.04` | current LTS; IPOPT packaged |
| `ubuntu-24.04.Dockerfile` | `ubuntu:24.04` | previous LTS, still supported; the older toolchain is where a use of something later than C++17 shows up |
| `debian-13.Dockerfile` | `debian:13` | Debian's IPOPT is usually a different version from Ubuntu's, which is the point of having both |
| `fedora-44.Dockerfile` | `fedora:44` | newest GCC in the Red Hat family; IPOPT packaged as `coin-or-Ipopt-devel` |
| `opensuse-leap-16.Dockerfile` | `opensuse/leap:16.0` | same sources as SUSE Linux Enterprise; conservative package set |
| `opensuse-tumbleweed.Dockerfile` | `opensuse/tumbleweed` | rolling; an early warning about dependencies that will reach the LTS distributions later |
| `arch.Dockerfile` | `archlinux:latest` | the distribution the published `psopt-ci` image is built on |
| `manjaro.Dockerfile` | `manjarolinux/base:latest` | Arch with updates held back a few weeks; when Arch fails and Manjaro passes, the cause changed in the last fortnight |

**RHEL 10 is not in the matrix**, and `containers/rhel-10.Dockerfile` is kept
all the same. It was in the first run and it was the only job that did not pass;
it is out because the effort of chasing it further was not matched by what it
would establish, and not because PSOPT is known to have a problem on RHEL.

The difficulty is not RHEL, it is the image. `redhat/ubi10` really is RHEL 10,
with the same glibc, compilers and runtime, but it does not have the same set of
*repositories*: the UBI set is deliberately narrow, so CodeReady Builder and EPEL
have to be enabled from inside the container and that is the part that does not
reliably work. The Dockerfile is written to tolerate their absence and let the
`ensure_` scripts build what is missing, which is a fair test of PSOPT against a
bare RHEL, and a useful one for anyone inside an air-gapped or
subscription-limited estate. It can be built by hand, and the file says how.

Worth keeping straight if it is ever revived: a green job there would say that
nothing in PSOPT depends on what the wider repository set adds. It would not say
that a subscribed RHEL 10 installation behaves identically, which this image
cannot show either way.

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

## What the first runs established

The images were tried on a Mac Studio first, natively on arm64 except the two
Arch-family ones which were run under x86_64 emulation, and then on GitHub's
x86_64 runners. **Eight of the nine passed on GitHub; RHEL 10 was the exception
and is no longer in the matrix.** Four also passed a full local build on arm64,
which is the only evidence there is that PSOPT works on aarch64 Linux at all.

The versions those four reported are the argument for everything in the section
on what is checked:

| distribution | IPOPT | Eigen |
|---|---|---|
| Ubuntu 26.04 | 3.14.19, packaged | 3.4.0, packaged |
| Ubuntu 24.04 | 3.11.9, packaged | 3.4.0, packaged |
| Debian 13 | 3.14.17, packaged | 3.4.0, packaged |
| Fedora 44 | 3.14.16, packaged | 5.0.1, packaged |

Eleven years of IPOPT releases, and all three externally-fixed examples matched
their references in every one. A matrix comparing costs against values an
earlier PSOPT run produced would have been reporting those version differences
as regressions.

**PSOPT builds and passes against Eigen 5.** Fedora 44 ships 5.0.1, as do Arch
and Manjaro, where every other distribution here is on 3.4.x. This was the open
question when the Arch probe first reported it, and the answer is that it makes
no difference.

Two things the runs found that needed fixing or recording:

**pacman 7 will not download inside a container build unless its sandbox is
turned off.** Arch failed at the first `pacman -Syu` with `error restricting
syscalls via seccomp: 22`, followed by `switching to sandbox user 'alpm'
failed`. pacman 7 fetches packages as an unprivileged user behind a seccomp and
landlock filter, 22 is EINVAL, and the container could not install the filter at
all. Both Arch-family Dockerfiles now set `DisableSandbox` in the `[options]`
section of `pacman.conf` before anything else. Manjaro did not hit it, because
its base image starts on an older pacman and upgrades to 7 part-way through the
same transaction. Below the `FROM` line those two Dockerfiles are byte for byte
identical, so one passing and the other failing was itself the diagnosis: the
package names had to be fine and the difference had to be the base image or
pacman. That is the argument for carrying both images, arriving sooner than
expected.

**A full local run wants more disk than a Mac usually has spare.** Seven full
builds in one pass exhausted the Docker Desktop virtual disk, and an out-of-space
write left the containerd content store inconsistent, which then looked like an
unrelated I/O error. Build one image at a time with
`docker image prune -af && docker builder prune -af` between them. Note also
that deleting Docker's disk image while Docker Desktop is still running frees
nothing, because the process holds the descriptor open.

## How far the package names have been checked

Every name in every image has resolved, on two architectures and two machines.
That question is closed for now.

It is worth remembering what it is not. These are eight distributions on one
day, three of them rolling releases that can rename a package next week, which
is what the weekly schedule is for. If a job dies inside `apt-get`, `dnf`,
`zypper` or `pacman`, that remains a statement about package names and not about
PSOPT: fix the name and run it again. Only once a job has reached the configure
step does what it reports become a statement about the software.
