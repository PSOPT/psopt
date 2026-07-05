# AD-backend shootout for PSOPT

A head-to-head of **five automatic-differentiation backends** on a representative
PSOPT-style optimal-control kernel, to inform which AD tool(s) to integrate into
PSOPT (whose native AD is ADOL-C). Measured on macOS 26 / arm64, 2026-07-05.

## What is compared
The same scalar-**templated** problem (`ocp_problem.hpp`) — a nonlinear coupled
DAE residual `R^n→R^n` plus a scalar Lagrange cost `R^n→R`, using the
sin/cos/exp/products typical of PSOPT DAEs — is differentiated by each backend.
Every backend's Jacobian and gradient are validated against central finite
differences (`common.hpp`); all five agree to ~1e-10.

| backend | kind | how it differentiates the problem |
|---|---|---|
| **ADOL-C** | operator-overload, taped | PSOPT's native AD; tape once, `jacobian()`/`gradient()` |
| **CppAD** | operator-overload, taped | `AD<double>` tape, `ADFun::Jacobian()` |
| **CppADCodeGen** | source codegen | records with CppAD, **generates + compiles C** for the Jacobian, calls the compiled function |
| **Sacado** | operator-overload, forward | `Sacado::Fad::DFad` vector forward mode (one sweep, all directions) |
| **Enzyme** | LLVM-IR AD | differentiates the **compiled `double`** functions — no source retemplating; fwd for Jacobian columns, reverse for the gradient |

## Results (µs per full derivative eval: Jacobian n×n + gradient)

| n | ADOL-C | CppAD | CppADCodeGen | Sacado | Enzyme |
|---|--------|-------|--------------|--------|--------|
| 10 | 5.9 | 22.5 | **0.1** | 4.1 | 0.9 |
| 20 | 10.3 | 78.2 | **0.3** | 13.3 | 3.2 |
| 50 | 34.2 | 444.7 | **1.0** | 74.4 | 20.5 |
| 100 | 109.6 | 1734.8 | **3.5** | 304.3 | 92.6 |
| 200 | 414.6 | 6847.5 | **11.0** | 1213.1 | 395.2 |

All backends returned identical derivatives (`maxErr(J)` equal across the row).

## Takeaways for PSOPT integration
- **CppADCodeGen** is the runtime winner by 1–2 orders of magnitude — it emits
  compiled derivative code. Cost: a one-time codegen+compile per problem, and it
  needs the problem recorded via CppAD (operator-overload), i.e. a scalar-templated
  problem API. Best fit when a problem is solved many times / very large.
- **Enzyme** is the only backend that needs **no source retemplating** — it
  differentiates PSOPT's compiled callbacks directly — and it beats ADOL-C
  (~1.1–1.5×). It requires building with the LLVM/Enzyme toolchain (clang plugin),
  not GNU/gcc. Best fit to bolt onto PSOPT's existing `double`-level callbacks.
- **ADOL-C** (baseline) is solid and already integrated; its optimized sparse
  drivers scale better than plain CppAD/Sacado dense forward modes.
- **Sacado** forward vector mode is competitive at small n but O(n²) for the
  Jacobian; **CppAD**'s per-call reverse `Jacobian()` is slowest here.
- Two integration routes emerge: a **scalar-templated problem API** (unlocks
  CppAD/CppADCodeGen/Sacado) or an **Enzyme/LLVM path** (reuses the compiled
  callbacks). CppADCodeGen and Enzyme are the two worth pursuing for speed.

## Reproduce
`bash build_all.sh` (edit the tool paths at the top). Backends are fetched/built
from source: coin-or/CppAD, joaoleal/CppADCodeGen, Trilinos/Sacado (C++20, Kokkos
OFF), EnzymeAD/Enzyme (vs llvm-20); ADOL-C 2.7.2 from the PSOPT superbuild.
Operator-overload drivers compile with `g++-mp-15`; Enzyme with `llvm-20 clang`.
Full run captured in `RESULTS.txt`.
</content>
