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
