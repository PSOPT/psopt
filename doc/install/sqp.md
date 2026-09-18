# Building and configuring PSOPT with its own sparse SQP solver

*What you need beyond a working PSOPT build*

| dependency | why | where CMake looks |
|---|---|---|
| MUMPS | the SQP reads the inertia of the KKT matrix from MUMPS, which IPOPT already links as its default linear solver -- so this is almost always a matter of pointing at what you have, not installing anything | `MUMPS_DIR`, `CMAKE_PREFIX_PATH`, pkg-config's IPOPT dirs; `MUMPS_LIBRARY` to name the library directly |
| GALAHAD | the sparse QP backend | `GALAHAD_DIR` |

If you built IPOPT and MUMPS yourself with coinbrew, as the [macOS](doc/install/macos.md),
[openSUSE](doc/install/opensuse.md) and [Arch](doc/install/arch.md) pages describe, MUMPS is
already in your `~/coin/dist` prefix and the `CMAKE_PREFIX_PATH` those pages set is enough to
find both the header and the library. Nothing further to do.

Both halves are needed and they are found separately. On Debian and Ubuntu
`pkg-config --libs ipopt` lists `-ldmumps_seq` itself, so the library resolves whether or
not CMake looks for it; a coinbrew IPOPT records the dependency inside `libipopt` instead,
and macOS will not resolve a symbol through an indirect dylib. If a link fails with an
undefined `dmumps_c`, point `MUMPS_LIBRARY` at the library holding it -- `libcoinmumps`
for a coinbrew build.

*GALAHAD*

`scripts/build_galahad.sh` does the whole of this: it installs the build tools, clones
GALAHAD, configures it with the options PSOPT needs, builds and installs it, and writes
an environment file to source.

```
./scripts/build_galahad.sh                        # installs under ~/galahad-install
./scripts/build_galahad.sh --prefix /opt/galahad --sudo
./scripts/build_galahad.sh --help
```

It works on macOS with either MacPorts or Homebrew, and on Debian/Ubuntu, Fedora and
Arch. On MacPorts it also runs `port select` so that a plain `gfortran` exists, since
MacPorts installs the compiler as `gfortran-mp-14` and meson looks for the plain name.

If you would rather do it by hand, follow the instructions at
https://github.com/ralna/GALAHAD; the options that matter are `-Dopenmp=true` (QPA's
linear solver needs OpenMP cancellation) and `-Dciface=true` (PSOPT includes
`galahad_qpa.h`). GALAHAD is a Fortran package, so it needs `gfortran` (MacPorts:
`sudo port install gcc14`). PSOPT links the Fortran runtime by asking CMake's Fortran
compiler for its own implicit link line, so a MacPorts or Homebrew gcc in a versioned
directory is found without help.

GALAHAD's QPA uses OpenMP cancellation, which the OpenMP runtime reads **once**, when it
initialises. It cannot be set from inside the process, so it has to be in the environment
before the program starts:

```
export OMP_CANCELLATION=TRUE
export OMP_PROC_BIND=TRUE
```

Without these the QP subproblems fail and the solver makes no progress. Put them in your
shell profile.

*Configuring*

```
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=ON \
      -DWITH_SQP=ON -DWITH_GALAHAD=ON \
      -DGALAHAD_DIR=/path/to/galahad/prefix
cmake --build build -j
```

Each backend is built as a separate loadable module under `build/qp_plugins` and opened at
run time with `RTLD_LOCAL`. That is not tidiness: every one of these libraries carries its
own AMD/COLAMD ordering code under the same C symbol names, and linked into one image they
bind to each other's and corrupt the result. `include/psopt_qp_plugin.h` has the details.
A CTest case, `qp_plugins_export_nothing_else`, checks that each module exports only the
four ABI entry points and nothing more; run it with `ctest -R qp_plugins` after building
with `-DBUILD_TESTS=ON`.

*Running an example under a different solver without editing it*

Comparing solvers across many examples means running one binary many ways, which
otherwise means editing each example's source. Configuring with
`-DPSOPT_ALLOW_ENV_OVERRIDES=ON` -- off by default -- lets the environment override
`algorithm` settings instead:

```
PSOPT_NLP_METHOD=SQP PSOPT_HESSIAN=exact PSOPT_QP_SOLVER=GALAHAD ./brac1
```

`PSOPT_SQP_STRATEGY`, `PSOPT_QP_RESTORATION`, `PSOPT_ELASTIC_PENALTY` and
`PSOPT_QP_ITER_MAX` work the same way. Every override is announced on stdout, naming the
setting the source asked for and the one being used instead: a program that quietly
disregards its own source is an unpleasant thing to debug, and worse than the convenience
is worth. In a build without the option the variables are ignored entirely.

*Using it*

```cpp
algorithm.nlp_method  = "SQP";
algorithm.hessian     = "exact";       // sparse exact Hessian of the Lagrangian
algorithm.derivatives = "automatic";   // required by "exact"
```

`qp_solver` defaults to `"GALAHAD"` and `sqp_strategy` to `"FM"`, so neither needs setting
unless you want something else. The other options -- `qp_restoration`, `elastic_penalty`
and `qp_iter_max` -- have defaults that are the measured best across the example set;
`include/psopt.h` documents each of them and says what is known about when to change it.
