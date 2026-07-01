# PSOPT Python binding (CasADi front end)

Define an optimal-control problem entirely in Python — symbolic dynamics, costs,
events, linkages and observations written in CasADi — and solve it through the
native PSOPT/IPOPT core. No change is made to the C++ engine.

## How it works

Two layers:

* **Layer A — driver** (`src/_psopt_driver.cpp`): a pybind11 extension that builds
  the PSOPT `Prob`/`Alg` objects in the required order, registers the user
  functions, solves, and returns the trajectory as NumPy arrays.
* **Layer B — codegen** (`psopt/codegen.py`, `psopt/_emitter.py`): each problem's
  CasADi maths is traced to an `adouble`-templated C++ source, JIT-compiled to a
  small math `.so`, and cached by a hash of its source. The driver `dlopen`s it
  and registers its `extern "C"` functions as PSOPT function pointers.

The math `.so` performs `adouble` arithmetic only and links **CppAD**. PSOPT's own
symbols (e.g. `auto_link`) are resolved at load time: the `_psopt` extension
whole-archives `libPSOPT` and is imported with `RTLD_GLOBAL`, so those symbols are
visible to each dlopened math `.so`. This requires a position-independent
`libPSOPT`, the extension to be linked with `--whole-archive` + `--export-dynamic`
(`-rdynamic`), and default symbol visibility on the extension.

## Building

In-tree (recommended) — from the top-level PSOPT build:

```
cmake -S . -B build -DPSOPT_PYTHON=ON     # makes libPSOPT PIC and builds the binding
cmake --build build
```

Standalone, against a prebuilt PIC `libPSOPT.a`:

```
cmake -S python -B build/python \
      -DPSOPT_ROOT=<psopt source root> \
      -DPSOPT_LIB=<.../libPSOPT.a> \
      -Dpybind11_DIR=$(python -m pybind11 --cmakedir)
cmake --build build/python
```

Either way CMake builds `psopt/_psopt*.so` and **generates** `psopt/_toolchain.py`
(compiler, flags, includes and CppAD lib) from the same configuration used to
build PSOPT — there is no hand-maintained, machine-specific toolchain file. The
JIT cache defaults to `$XDG_CACHE_HOME/psopt` (override with `PSOPT_CACHE_DIR`).

`pip install ./python` drives the same CMake build (scikit-build-core); for a
standalone install pass the PSOPT location via
`--config-settings=cmake.define.PSOPT_ROOT=...` /
`--config-settings=cmake.define.PSOPT_LIB=...`.

## Examples and validation

Each example reproduces its native C++ baseline:

| example            | what it exercises                         | result            | native baseline   |
|--------------------|-------------------------------------------|-------------------|-------------------|
| `bryson_denham.py` | single phase                              | 3.9995386676223115| 3.999539e+00 (bit-identical) |
| `launch.py`        | 4 phases, 24 linkages (Delta-III ascent)  | −7529.6612513     | −7.529661e+03 (within NLP tol) |
| `cracking.py`      | parameters + observation (estimation)     | 4.319519e-03      | 4.319519e-03 (bit-identical) |
| `bryson_ir.py`     | integrated-residual transcription         | 5.049038e-06      | 5.049038e-06      |
| `bryson_mesh.py`   | hp mesh refinement (10→25 nodes)          | 3.9999969178      | 3.999997e+00      |
| `rv2oe_casadi.py`  | CasADi orbital-element dynamics helper     | (used by launch)  | —                 |

## Provenance

Co-developed with AI assistance (Claude). All results validated against native
PSOPT baselines.
