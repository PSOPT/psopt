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
| **integral / isoperimetric** `∫ q dt ⋛ C` | auxiliary state `I' = q`, bound `I(t_f)` via an event | `isoperimetric` |
| **interior-point** `g(x(t₁))⋛c` at an interior time | split the phase at `t₁`, constrain via an event + linkage | `interior_point` |
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

## Genuinely missing (need native support / new backends — deferred)
- **Native declarations** for integral and interior-point constraints (today they
  require the reformulations above). A native `n_interior_points` could reuse the
  Radau/Gauss endpoint-interpolation code (`lagrange_endpoint_weights`) to
  interpolate the state to any interior time; a native integral constraint could
  reuse the cost-quadrature. Both are invasive NLP-offset changes.
- **Uncertainty propagation** for true chance/robust control (covariance states).
- **Mixed-integer / switching** constraints (need a MINLP backend, e.g. Bonmin).
- **Conic (SOC/SDP)** structure (only smooth-nonlinear reformulations today).

See `docs/EXTENSIONS_STATUS.md` for the overall branch status.
</content>
