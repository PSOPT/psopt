# PSOPT native IPOPT vs the CasADi backend — comparison and "net new"

*Branch `ecbrown`, 2026-07-05. Measured on macOS 26 / arm64 with the fully
self-contained GNU superbuild (`-DPSOPT_SUPERBUILD=ON -DPSOPT_WITH_CASADI=ON`).*

## TL;DR
The CasADi backend (`algorithm.nlp_method="CASADI"`) reuses PSOPT's **same**
transcription and hands the **same** NLP to IPOPT — just through
`casadi::nlpsol` instead of the native `IpoptApplication`. Empirically it
produces **identical objectives** on every test problem but is **5–25× slower**
(finite-difference + `MatrixXd↔DM` marshaling overhead). As currently built it
adds **no performance or capability benefit**. Its *potential* value —
alternative NLP solvers and exact-Hessian — is not yet realized (see "Net new").

## Architecture
Both paths share everything up to the NLP:

```
user problem (ADOL-C adouble dae/cost)
      │  PSOPT transcription (Legendre/Chebyshev/trapezoidal/Hermite-Simpson)
      ▼
   sparse NLP in Workspace  +  f/g callbacks  +  ADOL-C derivatives
      │
      ├── nlp_method="IPOPT"  → IPOPT_PSOPT (IpTNLP) → IpoptApplication      [native]
      └── nlp_method="CASADI" → casadi::nlpsol(<plugin>, PsoptNlp callbacks) [this]
```
`src/CASADI_interface.cxx` wraps PSOPT's `f`/`g` as a CasADi `Callback` and
solves with `nlpsol(algorithm.casadi_solver, ...)`. The linear solver
(`algorithm.ipopt_linear_solver`, e.g. mumps/ma97) still applies inside the
IPOPT plugin.

## Benchmark (native IPOPT vs CASADI-ipopt)
Same problems, same collocation, same tolerances; wall time is total run time.

| problem        | backend | status | secs | iters | objective     |
|----------------|---------|--------|------|-------|---------------|
| bryson_denham  | IPOPT   | SOLVED |  1   |  45   | 3.999539e+00  |
| bryson_denham  | CASADI  | SOLVED | 24   |  36   | 3.999539e+00  |
| hypersensitive | IPOPT   | SOLVED |  1   |  28   | 1.350105e+00  |
| hypersensitive | CASADI  | SOLVED |  6   |  28   | 1.350105e+00  |
| chain          | IPOPT   | SOLVED |  1   | 236   | 5.068480e+00  |
| chain          | CASADI  | SOLVED | 15   | 340   | 5.068480e+00  |
| catmix         | IPOPT   | SOLVED |  1   |  46   | -4.805320e-02 |
| catmix         | CASADI  | SOLVED |  7   |  37   | -4.805320e-02 |
| brac1          | IPOPT   | SOLVED |  2   |  52   | 8.247591e-01  |
| brac1          | CASADI  | SOLVED | 18   |  52   | 8.247591e-01  |
| obstacle       | IPOPT   | SOLVED |  6   |  34   | 4.571044e+00  |
| obstacle       | CASADI  | SOLVED | 55   |  37   | 4.571044e+00  |

Objectives match to displayed precision on all six; iteration counts are
comparable (same solver). The 5–25× slowdown is the finite-difference Jacobian
(CasADi perturbs the callback `n+1` times per evaluation) plus per-eval Eigen↔DM
copying.

## Capability matrix
| | native IPOPT | CASADI backend (as built) |
|---|---|---|
| Transcription | PSOPT | PSOPT (same) |
| Derivatives | ADOL-C exact (or numerical) | CasADi finite-difference over PSOPT callbacks |
| Hessian | limited-memory *or* exact (`sparse_hess`) | limited-memory only |
| Linear solver | mumps / ma57 / ma97 | same (inside IPOPT plugin) |
| Alternative NLP solvers | IPOPT, SNOPT | only `ipopt` (CasADi built WITH_IPOPT only) |
| Speed (these problems) | baseline | 5–25× slower |

## Net new — the honest verdict
- **As built: none.** Same transcription, same IPOPT, identical solutions, and
  slower. It is a compatibility/plumbing layer.
- **Potential, not yet realized:**
  1. *Alternative solvers* (`sqpmethod`, `blockSQP`) — the genuine capability gain
     over PSOPT's IPOPT/SNOPT-only set. **Now available**: CasADi is built with
     qpOASES + blockSQP (`WITH_QPOASES`/`WITH_BLOCKSQP`), so
     `algorithm.casadi_solver="sqpmethod"` runs (qpOASES solves the QP
     subproblems). But with the finite-difference backend + limited-memory
     Hessian it is impractically slow (bryson_denham did not converge in 400 s vs
     IPOPT's 1 s), so it is only useful once exact derivatives are supplied (2).
     (OSQP/fatrop omitted — their bundled builds fail on macOS 26.)
  2. *Exact derivatives* — the finite-difference penalty must be removed by
     feeding CasADi PSOPT's exact ADOL-C gradient/Jacobian (a documented
     follow-up; calling `IPOPT_PSOPT` standalone crashed and was reverted).
  3. *Exact Hessian / JIT / true symbolic transcription* — only possible by
     re-expressing the problem in CasADi MX/SX, which PSOPT's `adouble` callback
     API does not allow. That is a different problem front-end.

## How to use
```bash
cmake -B build -DPSOPT_SUPERBUILD=ON -DPSOPT_WITH_CASADI=ON   # builds CasADi + stack
```
```cpp
algorithm.nlp_method  = "CASADI";     // vs default "IPOPT"
algorithm.casadi_solver = "ipopt";    // plugin (only ipopt available as built)
algorithm.derivatives = "automatic";  // PSOPT's own f/g derivative mode (unchanged)
```

## Reproduce
`/opt/src/psopt_casadi_benchmark.sh` (rebuilds each example under both
backends from the superbuild tree; results in
`/opt/psopt/casadi-benchmark/RESULTS.tsv`).
</content>
