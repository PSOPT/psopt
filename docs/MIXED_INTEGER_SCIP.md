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
- **Phase D (done, heuristic)**: NONLINEAR dynamics via a Sequential-Linear-Programming
  loop in `SCIP_interface.cxx` — re-linearise about the incumbent, let SCIP re-optimise
  the MILP within a trust region on the continuous variables, then a second phase pins
  the integer schedule and polishes the continuous states. First-order, so it returns a
  **near-feasible** integer solution, not guaranteed tight/global feasibility (violation
  scales with the nonlinearity; e.g. quadratic drag `0.03·vel²` → ~7e-3). Example:
  `examples/scip_nlmiop`. For exact/global results keep the transcription linear.
- **Beyond Phase D**: tight nonlinear MINLP would need SCIP's native nonlinear
  expression interface (not reachable from PSOPT's opaque ADOL-C callbacks) or a full
  SQP/outer-approximation with an NLP subproblem solver — a separate, larger effort.

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
