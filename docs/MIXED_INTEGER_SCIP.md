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
- **Phase A (done)**: SCIP built + standalone MILP/MIQP MIOC demonstrator.
- **Phase B**: add the integer-declaration plumbing to `psopt.h`/`setup.cxx`
  (`integer_controls`, `integer_parameters` per phase) — inert for existing problems.
- **Phase C**: a `SCIP_interface.cxx` backend for the **MILP/MIQP transcription
  subclass** (linear/quadratic dynamics & cost), routing the transcribed variables +
  integrality to SCIP and unpacking the solution. This covers switched-linear and
  many practical MIOC problems without needing SCIP's full nonlinear NLPI.
- **Phase D**: general MINLP via SCIP's expression interface (largest effort;
  optionally via CasADi's SCIP/bonmin plugins as an intermediate).

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
