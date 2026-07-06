# PSOPT constraint types — catalogue & how to express them

PSOPT is a direct-collocation solver; every constraint below becomes an NLP
constraint. This catalogues what PSOPT supports and which example demonstrates
each pattern. Several "advanced" constraints are not *native declarations* but
are expressed exactly via the core mechanisms (path / events / linkages /
multi-phase / auxiliary states), as shown by the examples.

## Native mechanisms
| mechanism | declares | form | example |
|---|---|---|---|
| variable bounds | `bounds.{lower,upper}.{states,controls,parameters}` | `x_L ≤ x(t) ≤ x_U` at every node | most examples |
| path constraints | `npath`, `path[]` in `dae()`, `bounds.*.path` | `p_L ≤ g(x,u,p,t) ≤ p_U` at every node (equality if `p_L=p_U`) | `obstacle` |
| event constraints | `nevents`, `events()`, `bounds.*.events` | boundary/endpoint conditions | `bryson_denham` |
| phase linkages | `nlinkages`, `linkages()` | multi-phase coupling/continuity | `twophase_schwartz` |
| time bounds | `bounds.*.{StartTime,EndTime}` | fixed/free/bounded t₀,t_f | `bryson_denham` |

## Advanced constraints (expressed via the mechanisms above)
| constraint | technique | example |
|---|---|---|
| **integral / isoperimetric** `∫ q dt ⋛ C` | **NATIVE**: `problem.phase[i].nintegral`, `integral_constraints()`, `bounds.*.integral`; *or* the aux-state pattern | **`integral_ctc`** (native), `isoperimetric` (aux-state) |
| **interior-point** `g(x(t₁))⋛c` at an interior time | **NATIVE**: `problem.phase[i].ninterior`, `interior_time`, `interior_point_constraints()` (single phase); *or* the phase-split pattern | **`interior_ptc`** (native), `interior_point` (phase-split) |
| **control-rate** `|u̇| ≤ R` | promote `u` to a state, new control `v=u̇`, bound `|v|≤R` | **`control_rate`** (new) |
| **time-windowed path** `g(x)≤0` only for `t∈[t_a,t_b]` | smooth gate: `path = g − BIG·(1−gate(t))`, `gate≈1` in-window | **`path_window`** (new) |
| **chance / robust** `P{g≤0} ≥ 1−ε` | deterministic tightening `g + k·σ ≤ 0`, `k=Φ⁻¹(1−ε)` | **`chance_constraint`** (new) |
| **complementarity (MPEC)** `0 ≤ a ⊥ b ≥ 0` | smoothed NCP / relaxation `a·b ≤ μ`, `μ→0` | `mpec` |
| **pure state** `g(x)≤0` (no control) | path constraint on states only (watch constraint order/index) | `obstacle` (mixed) |

## New illustrative examples (added on branch `ecbrown`)
Built + verified against native IPOPT (Apple-clang build):
- **`control_rate`** — rest-to-rest double integrator, `|u̇| ≤ 4`; cost 1.846.
  Shows the control-as-state reformulation for control-rate limits.
- **`path_window`** — velocity limited to 1.30 only for `t∈[0.35,0.65]`; the
  constraint binds around the mid-move peak, raising cost 12.0 → 12.98.
- **`chance_constraint`** — `P{vel ≤ 1.30} ≥ 0.95` via tightening
  `vel + 1.6449·σ ≤ 1.30`; velocity caps at **1.251** (= the tightened limit),
  cost 12.0 → 13.87. A full version propagates `σ(t)` via a Lyapunov ODE as
  extra states.

## Native interior-point constraints (NEW)
`problem.phase[i].ninterior`, `interior_time` (normalized [0,1]), and the
`interior_point_constraints(g, xinterp, uinterp, params, t, index, iphase, ws)`
hook enforce a constraint on the state/control interpolated to an interior time
**without splitting the phase**. It reuses the Radau/Gauss endpoint-interpolation
(`lagrange_endpoint_weights`) and slots into the per-phase constraint block just
before the t0≤tf slot; `ninterior==0` (the default) leaves the layout unchanged,
so existing examples are byte-for-byte unregressed. Example: **`interior_ptc`**
(waypoint `pos(0.3·T)=0.5`, which raises the min-energy cost 12→38.1).

## Native integral constraints (NEW)
`problem.phase[i].nintegral` + the `integral_constraints(q, states, controls,
params, t, iphase, ws)` hook + `bounds.*.integral` declare `I_L ≤ ∫q dt ≤ I_U`
directly. PSOPT integrates each integrand with the **collocation quadrature**
(the same `workspace->w` weights as the running cost) and bounds the result; it
sits in the per-phase tail `[…|interior|integral|t0≤tf]`, and `nintegral==0`
(default) leaves the layout unchanged (zero regression). Example: **`integral_ctc`**
(area `∫pos dt=0.40` vs unconstrained 0.5 → cost 12→19.2).

## Genuinely missing (need new machinery / backends — deferred)
- **Uncertainty propagation** for true chance/robust control (covariance states).
- **Mixed-integer / switching** constraints — MINLP backend (SCIP). See
  `docs/MIXED_INTEGER_SCIP.md` for the integration status.
- **Conic (SOC/SDP)** structure (only smooth-nonlinear reformulations today).

See `docs/EXTENSIONS_STATUS.md` for the overall branch status.
</content>
