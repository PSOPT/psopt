# Mixed-integer optimal control via SCIP — status & integration path

## Status
- **SCIP v11 + SoPlex built from source** (`/opt/claude/src/build_scip.sh` →
  `/opt/claude/scip`), arm64 macOS, GMP from MacPorts.
- **Working MIOC demonstrator**: `scip-mioc/mioc_scip.c` — a double integrator with
  integer thrust `u_k∈{−1,0,1}` transcribed to a MILP and solved to global optimality
  (bang-off-bang, total actuation 8.0). Standalone (links SCIP, not PSOPT) so the
  PSOPT build never depends on SCIP.

## Why standalone (not yet an in-core backend)
PSOPT transcribes to a **continuous** NLP and hands it to IPOPT/CasADi via callbacks
(ADOL-C derivatives). A full in-core SCIP backend (`algorithm.nlp_method="SCIP"`)
requires more than dropping in a solver:
1. **Integer-variable declaration** — per-phase integer/binary control & parameter
   index lists (new `phases_str` fields + bounds), so PSOPT knows which transcribed
   NLP variables SCIP must branch on.
2. **NLPI expression handoff** — SCIP does not consume an opaque IPOPT-style TNLP; it
   wants either its expression graph (cons_nonlinear) or an NLPI plugin. PSOPT's
   ADOL-C-taped objective/constraints must be exposed to SCIP's NLP interface (or the
   problem restricted to the MILP/MIQP subclass, which covers linear/quadratic
   transcriptions like the demonstrator).
3. **Solve loop** — replace the IPOPT `solve()` in `psopt()` with a SCIP path that
   builds the MINLP, sets integrality, solves, and unpacks the incumbent into the
   PSOPT `Sol` structure (states/controls per phase).

## Recommended path
- **Phase A (done)**: SCIP built + standalone MILP/MIQP MIOC demonstrator (`scip-mioc/`).
- **Phase B (done)**: integer-declaration plumbing — `problem.phase[i].integer_controls`
  and `integer_parameters` in `psopt.h`/`setup.cxx`; inert for existing problems.
- **Phase C (done)**: `src/SCIP_interface.cxx` — an in-core `algorithm.nlp_method="SCIP"`
  backend for the **linear-transcription subclass**. It linearises PSOPT's transcription
  about the feasible guess `x0` (exact when the free-variable model is linear; the raw
  transcription is bilinear in the free times `t0,tf`, so linearising at `x0`—not 0—is
  essential), marks the declared integer variables, and hands SCIP an explicit MILP;
  the incumbent is unpacked back into the PSOPT solution. A linearity self-check warns
  if the dynamics/cost are nonlinear. Example: `examples/scip_miop` — integer thrust
  `u∈{-1,0,1}`, global bang-off-bang, Σ|u|=8, solved via `nlp_method="SCIP"`.
  Build: `-DPSOPT_WITH_SCIP=ON` (compiles `SCIP_interface.cxx`, links `libscip`).
- **Phase D (done)**: NONLINEAR dynamics via a Sequential-Linear-Programming loop in
  `SCIP_interface.cxx` — re-linearise about the incumbent, let SCIP re-optimise the MILP
  within a trust region on the continuous variables, plus a second phase that pins the
  integer schedule and polishes the continuous states. A **second-order correction (SOC)**
  — evaluate the true constraints at the trial point, shift the model by that curvature
  residual, and re-solve — is the key heuristic: it converts the stalling first-order SLP
  into a **converging** solver. On the `scip_nlmiop` example (quadratic drag `0.03·vel²`)
  it drives the violation to **2.8e-8** (vs a first-order stall at ~7e-3). Reported in
  three tiers: converged / near-feasible / not-converged. Strong nonlinearity can still
  land in the near-feasible tier (drag `0.1·vel²` → ~1.3e-2, integer-schedule instability),
  better than first-order (~2.3e-2); keep the transcription linear for exact/global results.
- **Refinement (Phase D+)**: two second-order/OA pieces improve the nonlinear case:
  1. **Second-order correction (SOC)** in the SLP loop (curvature-aware re-solve) —
     converges mild–moderate nonlinearity (drag `0.03·vel²` → 2.8e-8).
  2. **Gauss-Newton feasibility polish** — after the SLP, pin the incumbent integer
     schedule and drive the continuous states to feasibility with a trust-region-free
     Gauss-Newton + SOC loop (no objective, pure restoration). Safe (skipped once
     converged); tightens cases where the schedule admits a feasible completion.
  3. **Outer-Approximation (OA) driver** — when residual infeasibility remains, alternate
     a SCIP **master** (integers + an objective epigraph + *accumulated* OA linearisation
     cuts `g(p)+J_p(x−p)`) with an **IPOPT NLP subproblem** that pins the master's integer
     schedule and solves the continuous problem to true nonlinear feasibility
     (`psopt_ipopt_resolve()`), adding each subproblem point as a new cut. This searches
     *different* integer schedules (not just polishing one), improving strongly-nonlinear
     cases (drag `0.1·vel²`: 1.3e-2 → 8.8e-3 through the IPOPT subproblem; the master then
     proves no better schedule under the cuts). For a feasible schedule the subproblem
     drives the violation to ~0.
  - The nested IPOPT call was unblocked by a real bug fix: the SCIP branch now allocates
    the Jacobian-sparsity workspace (`workspace.cxx`), which `IPOPT_PSOPT::get_nlp_info`
    writes — leaving it NULL (as the SCIP path did) segfaulted, which had looked like a
    reentrancy limit but was just an unallocated buffer.
- **Remaining frontier**: OA is exact/global only for *convex* MINLP; for nonconvex
  dynamics it is a strong heuristic (no-good/spatial-branching cuts, or SCIP's native
  nonlinear expression interface, would be needed for a global guarantee). A genuinely
  infeasible discretisation (no integer schedule reaches the bounds) is reported as the
  least-infeasible near-feasible point.

## In-core usage
```cpp
problem.phase[1].integer_controls.resize(1,1);
problem.phase[1].integer_controls(0) = 0;   // control 0 is integer at every node
algorithm.nlp_method = "SCIP";              // + linear dynamics & a linear objective
```
Keep the transcription linear in the free variables: linear dynamics, and a linear
objective (use an auxiliary control `a>=|u|` with `INT a dt` instead of `INT|u|dt`).
The SCIP backend does a one-shot MILP solve — use `mesh_refinement="manual"`.

## Build reference
Superbuild (SoPlex + SCIP as ExternalProjects; GMP + readline from
`PSOPT_SCIP_DEP_PREFIX`, MacPorts `/opt/local` by default):
```
cmake -S . -B build -DPSOPT_SUPERBUILD=ON -DPSOPT_WITH_SCIP=ON \
      -DPSOPT_SCIP_DEP_PREFIX=/opt/local
cmake --build build -j4            # builds ep_soplex, ep_scip (uses ep_ipopt for MINLP NLP)
```
Standalone (the demonstrator here):
```
bash /opt/claude/src/build_scip.sh      # SoPlex + SCIP -> /opt/claude/scip
cd scip-mioc && bash build.sh && \
  DYLD_LIBRARY_PATH=/opt/claude/scip/lib:/opt/local/lib ./mioc_scip
```
</content>
