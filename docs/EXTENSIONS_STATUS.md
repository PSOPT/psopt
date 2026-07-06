# PSOPT extensions — status & known issues (branch `ecbrown`, 2026-07)

Closeout record for the extensions added on this branch. Each item lists what is
done, what is deferred, and any known issues.

## Done (built, validated, committed)

### CasADi NLP backend — exact derivatives + exact Hessian
- `algorithm.nlp_method="CASADI"`, `casadi_solver` plugin (ipopt / sqpmethod / blocksqp).
- Exact ADOL-C first derivatives via `Callback::get_jacobian` (kills the FD penalty).
- Exact Lagrangian Hessian via reverse-mode adjoint overrides (`hessian="exact"`);
  ~½ the IPOPT iterations (chain 236→18, catmix 46→16).
- sqpmethod indefinite-Hessian crash fixed (convex QP + `error_on_fail=false`).
- **Known issue**: exact-Hessian is not yet robust across the full suite — objectives
  match 28/32 solved-by-both, with segfaults on a few degenerate problems (e.g. `mpec`).
  The constraint-Hessian re-tape path (`ConRevJac`) needs hardening. Default remains
  limited-memory BFGS on exact gradients, which is robust.
- Files: `src/CASADI_interface.cxx`, `src/NLP_interface.cxx`, `src/validate.cxx`.

### Radau (LGR) and Gauss (LG) pseudospectral collocation
- `collocation_method="Radau"/"Gauss"`; nodes/weights/diff-matrix in
  `src/pseudospectral.cxx` (validated in `docs/ps_nodes_test.cpp`); endpoint
  interpolation (`Lt0/Ltf`) in `src/NLP_constraints.cxx`.
- bryson_denham = 3.999539 for Legendre/Radau/Gauss (Legendre unregressed).
- **Known issue**: Gauss is less robust than Radau on stiff/boundary-layer problems
  (hypersensitive) — as in the literature. Radau is the recommended default.

### Ross-Fahroo / Bellman `ps_method`
- `algorithm.ps_method="Ross-Fahroo"` (default) / `"Bellman"`; Bellman requires
  automatic mesh refinement. Verified: identical optima on the tested examples.
- Files: `src/util.cxx`, `src/psopt.cxx` (bellman bootstrap), `src/mesh.cxx`, `src/setup.cxx`.

### Linear solvers
- **SPRAL** (`ssids`, free, OpenMP-threaded) built into IPOPT (`PSOPT_WITH_SPRAL`);
  `ipopt_linear_solver="spral"` verified. Requires an OpenMP environment
  (`OMP_CANCELLATION/OMP_NESTED/OMP_PROC_BIND/OMP_STACKSIZE`) or SSIDS aborts (flag -53).
- **PARDISO** build plumbing (`PSOPT_PARDISO_LFLAGS` → `--with-pardiso`).
- Selector is already pass-through (any solver IPOPT is built with is selectable).

### Threading / parallelism (how to actually get parallel runs)
PSOPT's own code is serial — it contains **no OpenMP pragmas and no MPI calls**, and
calls **no BLAS directly** (it uses Eigen kernels). All parallelism therefore lives
*below* PSOPT and is enabled entirely by environment variables + the linear-solver
choice:
- **Threaded linear solver**: pick one built with OpenMP — `ipopt_linear_solver="spral"`
  (free) or `"ma97"` (HSL) — controlled by `OMP_NUM_THREADS` (MUMPS is sequential).
- **Dense BLAS/LAPACK** used inside IPOPT/MUMPS/SPRAL/HSL has its **own** thread knob,
  separate from `OMP_NUM_THREADS`:
  - OpenBLAS (the default here; `otool -L libipopt` → `libopenblas`) → `OPENBLAS_NUM_THREADS`
  - Intel MKL (used by `ipopt_linear_solver="pardisomkl"`) → `MKL_NUM_THREADS`

Standard recipe — set the count in ALL of them so it is honoured whichever BLAS/solver
is linked (otherwise OpenBLAS/MKL spin up their own default pool and "serial" isn't):
```sh
N=8
export OMP_NUM_THREADS=$N OPENBLAS_NUM_THREADS=$N MKL_NUM_THREADS=$N
export OMP_CANCELLATION=TRUE OMP_NESTED=TRUE OMP_PROC_BIND=TRUE OMP_STACKSIZE=64M   # SPRAL/SSIDS
# then: algorithm.ipopt_linear_solver = "spral";   (serial: N=1, all of the above = 1)
```
MPI is intentionally unsupported: PSOPT never invokes distributed MUMPS, and the MPI
build is documented to crash on macOS 26.

### Native constraint types (NEW)
- **Interior-point** constraints (`ninterior`, `interior_point_constraints`) — state/
  control interpolated to an interior time without a phase split; reuses the
  Radau/Gauss endpoint interpolation. Example `interior_ptc`.
- **Integral / isoperimetric** constraints (`nintegral`, `integral_constraints`) —
  `I_L ≤ ∫q dt ≤ I_U` via the collocation quadrature. Example `integral_ctc`.
- Both default to 0 count → zero regression on the 47 existing examples.
- Reformulation examples for the rest: `control_rate`, `path_window`,
  `chance_constraint`. See `docs/CONSTRAINTS.md`.

### Mixed-integer optimal control via SCIP (NEW)
- SCIP v11 + SoPlex built from source; standalone MILP MIOC demonstrator
  `scip-mioc/mioc_scip.c` (integer thrust, global bang-off-bang). In-core backend
  path documented in `docs/MIXED_INTEGER_SCIP.md`.

### AD-backend shootout
- Five backends (ADOL-C/CppAD/CppADCodeGen/Sacado/Enzyme) built + benchmarked in
  `ad-backends/` — the decision input for a future PSOPT AD provider.

## Deferred (identified, not implemented)
- **AD-backend integration into PSOPT** — `derivatives="enzyme"` / `"cppadcodegen"`
  providers (Enzyme first; CppADCodeGen needs a scalar-templated problem API). Highest-
  leverage next item; see `ad-backends/README.md`.
- **PARDISO / HSL ma97 validation** — need the licensed libraries (Panua PARDISO / MKL are
  unavailable on arm64 macOS; HSL is non-free). Plumbing is ready.
- **hp-adaptive p-order** refinement on top of the existing h-refinement.
- **Multiple shooting**, **indirect methods (PMP/BVP)**, **mixed-integer (SCIP/Bonmin)**,
  **robust/stochastic**, **homotopy/continuation** — each its own effort.

## Known environment issue
- The MacPorts **gcc-15** toolchain has an **intermittent assembler ICE** on macOS 26
  under load (`internal compiler error … signal terminated program as`). It is
  nondeterministic — retries usually succeed. The Apple-clang build (`build-ct`) is
  unaffected; use it for native-IPOPT work when the GNU superbuild is not required.

## Build
Self-contained GNU superbuild: `cmake -B build -DPSOPT_SUPERBUILD=ON [-DPSOPT_WITH_CASADI=ON]
[-DPSOPT_WITH_SPRAL=ON] [-DPSOPT_WITH_HSL=ON] [-DPSOPT_PARDISO_LFLAGS=...]`. See
`cmake/superbuild.cmake`.
</content>
