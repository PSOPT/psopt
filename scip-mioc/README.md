# Mixed-integer optimal control with SCIP

PSOPT's NLP backends (IPOPT, CasADi) solve **continuous** optimal control. Problems
with **integer/binary controls** — on/off actuators, gear selection, mode switching,
minimum-actuation bang-off-bang — are mixed-integer optimal control (MIOC) and need
a MINLP/MILP solver. This directory demonstrates MIOC with **SCIP** (a free
constraint-integer-programming solver), independent of the PSOPT NLP core.

## Example: `mioc_scip.c`
A double integrator moved rest-to-rest (pos 0→1, vel 0→0) with an **integer thrust**
`u_k ∈ {−1,0,1}`, minimizing total actuation `Σ|u_k|`. Transcribed with explicit
Euler over N=30 stages (T=3) into a **MILP** and solved to **global optimality**:

```
optimal total actuation = 8.0   (SCIP: optimal solution found)
schedule: +1 +1 +1 +1 (accelerate) ... 0 0 0 (coast) ... −1 −1 −1 −1 (decelerate)
reaches pos=1.000, vel=0.000 at t=3.0
```

The result is the expected **bang-off-bang** integer control — global, not a rounded
continuous relaxation.

## Build & run
```
bash build.sh        # links the SCIP built into /opt/claude/scip (SoPlex LP backend)
DYLD_LIBRARY_PATH=/opt/claude/scip/lib:/opt/local/lib ./mioc_scip
```
SCIP (v11) + SoPlex are built from source (github.com/scipopt) via
`/opt/claude/src/build_scip.sh`.

## Relationship to PSOPT
This is a **standalone** demonstrator (links SCIP, not PSOPT) — deliberately not wired
into the PSOPT `examples/` build so the standard build never depends on SCIP. It shows
the transcription→MILP/MINLP→SCIP pattern that a future in-core PSOPT backend
(`algorithm.nlp_method="SCIP"` with per-phase integer control/parameter declarations,
routed through SCIP's NLPI) would use. See `docs/MIXED_INTEGER_SCIP.md` for that
integration path. For **continuous** optimal control, keep using PSOPT/IPOPT/CasADi.
</content>
