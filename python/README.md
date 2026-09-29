# PSOPT Python binding (CasADi front end)

Define an optimal-control problem entirely in Python — symbolic dynamics, costs,
events, linkages and observations written in CasADi — and solve it through the
native PSOPT/IPOPT core. No change is made to the C++ engine.

## How it works

Two layers:

* **Layer A — driver** (`src/_psopt_driver.cpp`): a pybind11 extension that builds
  the PSOPT `Prob`/`Alg` objects in the required order, registers the user
  functions, solves, and returns the trajectory as NumPy arrays.
* **Layer B — codegen** (`psopt/codegen.py`, `psopt/_emitter.py`): each problem's
  CasADi maths is traced to an `adouble`-templated C++ source, JIT-compiled to a
  small math `.so`, and cached by a hash of its source. The driver `dlopen`s it
  and registers its `extern "C"` functions as PSOPT function pointers.

The math `.so` performs `adouble` arithmetic only and links **CppAD**. PSOPT's own
symbols (e.g. `auto_link`) are resolved at load time: the `_psopt` extension
whole-archives `libPSOPT` and is imported with `RTLD_GLOBAL`, so those symbols are
visible to each dlopened math `.so`. This requires a position-independent
`libPSOPT`, the extension to be linked with `--whole-archive` + `--export-dynamic`
(`-rdynamic`), and default symbol visibility on the extension.

## Building

In-tree (recommended) — from the top-level PSOPT build:

```
cmake -S . -B build -DPSOPT_PYTHON=ON     # makes libPSOPT PIC and builds the binding
cmake --build build
```

Standalone, against a prebuilt PIC `libPSOPT.a`:

```
cmake -S python -B build/python \
      -DPSOPT_ROOT=<psopt source root> \
      -DPSOPT_LIB=<.../libPSOPT.a> \
      -Dpybind11_DIR=$(python -m pybind11 --cmakedir)
cmake --build build/python
```

Either way CMake builds `psopt/_psopt*.so` and **generates** `psopt/_toolchain.py`
(compiler, flags, includes and CppAD lib) from the same configuration used to
build PSOPT — there is no hand-maintained, machine-specific toolchain file. The
JIT cache defaults to `$XDG_CACHE_HOME/psopt` (override with `PSOPT_CACHE_DIR`).

`pip install ./python` drives the same CMake build (scikit-build-core); for a
standalone install pass the PSOPT location via
`--config-settings=cmake.define.PSOPT_ROOT=...` /
`--config-settings=cmake.define.PSOPT_LIB=...`.

## Discrete-valued controls and parameters

A control or a static parameter may be restricted to a finite admissible set. The
dynamics are written once, with the quantity treated as an ordinary control or
parameter; one declaration per quantity then makes it discrete:

```python
ph.declare_integer_control(0, [0.0, 1.0])        # control index, admissible values
ph.declare_integer_parameter(0, [0.0, 1.0, 2.0, 3.0])
```

Several integer controls may be declared in the same phase; PSOPT convexifies over
the Cartesian product of their admissible sets, so the weight controls carry an
SOS1 constraint over all admissible combinations, and each declared control is
recovered separately by sum-up rounding. Integer parameters take a different route
— a static parameter cannot chatter, so the relaxed optimum is generally not
realisable — and are solved exactly by enumeration; declaring one switches the
driver to `psopt_solve_integer` automatically.

When any integer control is declared, `sol.controls` holds the convexified weight
controls. The rounded trajectories come back separately:

```python
ic = sol.integer_controls[k]     # k-th declared integer control
ic.control, ic.time, ic.interval_widths, ic.n_switches
sol.integer_parameters[j]        # {'index': ..., 'value': ...}
```

Multiphase problems expose the same fields per phase on `MultiSolution`.

## What the solver gives back

`Problem.solve` returns a `Solution` (single phase) or a `MultiSolution` (several).
Both carry the trajectory and, since the interface was brought level with the C++
one, everything else the solve produced.

**Did it work.** Two different questions, and the interface answers them
differently.

A problem PSOPT *refuses* — an invalid option, a bad dimension, a guess of the
wrong shape — raises `psopt.PSOPTError`, carrying PSOPT's own diagnostic and the
Solution it came from:

```python
try:
    sol = prob.solve(alg)
except psopt.PSOPTError as e:
    print(e)                  # e.g. algorithm.ms_integrator must be "RK4", ...
    print(e.solution.status)  # the return codes; e.details has the full banner
```

This is the Python interface departing from PSOPT's own default on purpose.
`algorithm.on_error` defaults to `"fail-soft"` here, where the C++ default is
`"fail-fast"` and calls `exit()`. In C++ that is defensible; inside a Python
process it terminates the interpreter — no traceback, nothing to catch, a dead
kernel in a notebook — and with `print_level=0` it does so in silence. Passing
`on_error="fail-fast"` restores PSOPT's own behaviour, `exit()` included.

A solve that *runs* but does not converge is a result rather than an error: it
returns normally and says so, so that a non-converged trajectory can still be
looked at. `sol.objective` comes back whatever happened, so it cannot answer this
on its own:

```python
sol = prob.solve(alg)
if not sol.status.success:
    raise RuntimeError(sol.status.error_msg or
                       "NLP return code %d" % sol.status.nlp_return_code)
print(sol.status)     # success, nlp_return_code, error_flag, cpu_time
```

`nlp_return_code` means different things in the two NLP solvers and `success`
accounts for it: from IPOPT, 0 is "solved" and 1 is "solved to acceptable level",
both successes; from PSOPT's own SQP, 1 is the iteration limit reached and is a
failure.

**Multipliers and diagnostics.** `sol.costates` is the discrete adjoint — the thing
a solution is checked against the maximum principle with — and the rest sit on
`sol.duals`:

```python
sol.costates                     # (nstates, nnodes); list of these on a MultiSolution
sol.duals.hamiltonian            # constant on an autonomous free-final-time problem
sol.duals.dual_path              # one multiplier per path constraint, if any
sol.duals.dual_events            # one per event constraint, if any
sol.duals.relative_local_error   # what mesh refinement drives down
sol.duals.terminal_state         # Gauss collocation only; empty otherwise
sol.duals.terminal_costate       # the transversality condition, likewise
sol.dual_linkages                # MultiSolution only
```

**Work done.** `sol.status.mesh_stats` is one dict per mesh-refinement iteration —
method, nodes, variables, constraints, evaluation counts, the discretisation error
reached and the CPU time — which is the table PSOPT prints at the end of a run.

**Estimation statistics.** For a problem with an observation function, and with
`parameter_statistics="yes"` (the default):

```python
ps = sol.parameter_statistics        # None if PSOPT could not form them
ps.covariance, ps.standard_errors, ps.confidence_low, ps.confidence_high
ps.sigma_hat, ps.residuals
```

Residuals can be weighted and the parameter vector regularised, per phase:

```python
ph.residual_weights = 1.0 / sigma          # (nobserved, nsamples)
ph.regularization_factor = 1.0e-4          # adds this multiple of ||p||^2
```

## Algorithm options

`psopt.Algorithm` accepts every field of the C++ `Alg` structure. Beyond the core
four (`collocation_method`, `nlp_method`, `derivatives`, `scaling`) that means:

* **mesh refinement** — `mesh_refinement`, `mr_max_iterations`, `ode_tolerance`,
  `mr_max_growth_factor`, `mr_min_order`, `mr_max_order`, `mr_kappa`, `mr_M1`,
  `mr_switch_detection`, `switch_order`
* **integrated residuals** — `transcription_method="integrated-residual"` plus the
  `ir_*` family, including `ir_flexible_mesh` and `ir_element_local_controls`
* **multiple shooting** — `transcription_method="multiple-shooting"` plus the `ms_*`
  family: integrator, steps per segment, control parameterisation, path sampling,
  flexible segments and the index-1 DAE settings
* **PSOPT's own SQP** — `nlp_method="SQP"` plus `qp_solver`, `sqp_strategy`,
  `qp_restoration`, `qp_iter_max`, `trust_region`, `trust_region_radius`,
  `elastic_penalty`
* **the rest** — `hessian`, `objective_form`, `defect_scaling`, `diff_matrix`,
  `ipopt_linear_solver`, `ipopt_max_cpu_time`, `constraint_scaling`,
  `jac_sparsity_ratio`, `hess_sparsity_ratio`, `save_sparsity_pattern`,
  `nsteps_error_integration`, `parameter_statistics`, `parameter_estimation_norm`,
  `hessian_verify`, `on_error`, `max_integer_combinations`, `print_level`,
  `diagnostic_level`

An option left at `None` is not sent, so PSOPT's own default applies.

## What the code generator will not emit

The maths is traced from CasADi and emitted as C++ that CppAD tapes, so an
operation has to be one CppAD can differentiate. Branches on a symbolic value --
`ca.if_else`, comparisons -- are refused with a message saying so, because a tape
records the branch taken when it was made and then uses it everywhere. Use a smooth
transition, or split the problem into phases at the switching instant. `floor`,
`ceil`, `fmod`, `sign` and `copysign` are refused for the same reason: they are not
differentiable. Anything else that is missing is a gap rather than a refusal, and
the error says which of the two you have met.

## Examples and validation

Twenty-two examples, in `examples/`. Run one directly, or all of them:

```
cd python/examples
python3 obstacle.py
python3 run_all.py            # every example, with a pass/fail table
python3 run_all.py bryson     # just the ones whose name matches
```

Each runs from a source checkout with nothing installed: `_common.py` puts the
package on `sys.path` if it is not already importable. The examples that check
themselves exit non-zero when the answer does not match their reference, which is
what `run_all.py` reports.

| example | what it exercises | checked against |
|---|---|---|
| `brachistochrone.py` | free final time, costates, the Hamiltonian | the cycloid, in closed form (agrees to 3e-10); C++ `brac1` |
| `breakwell.py` | a state bound, and the costate jump at a boundary arc | 4/(9l) exactly; C++ `breakwell` |
| `obstacle.py` | path constraints and their multipliers | C++ `obstacle`, 9.970637e-01 |
| `multiple_shooting.py` | multiple shooting against collocation | J* = 6 and u* = 6 - 12t, in closed form |
| `integrated_residual.py` | integrated residuals on a singular arc | J* = 10/3, tf* = 4, in closed form |
| `hypersensitive.py` | hp mesh refinement and `mesh_stats` | C++ `hypersensitive`, 1.330826e+00 |
| `weighted_estimation.py` | residual weights, regularisation, covariance | C++ `cracking`, 4.319519e-03 |
| `cracking.py` | parameters and an observation function | C++ `cracking`, 4.319519e-03 |
| `bryson_denham.py` | the simplest single phase | C++ `bryson_denham`, 3.999539e+00 |
| `launch.py` | 4 phases, 24 linkages (Delta-III ascent) | C++ `launch`, -7.529661e+03 |
| `bryson_mesh.py` | hp mesh refinement, 10 -> 29 nodes | 3.999997e+00 |
| `bryson_ir.py` | integrated-residual pass-through | see the note below |
| `lotka_integer.py` | a binary integer control, sum-up rounding | 1.348104 / 1.351850, 4 switches |
| `integer_parameter.py` | an integer static parameter, by enumeration | closed form p = 2, J = 0.09 |
| `robust_arm.py` | robust optimal control by scenario augmentation | the nominal design must fail out of sample and the generated one must not |
| `robust_driver_arm.py` | the same problem through `psopt.robust` | the driver must certify the design over the whole set |
| `robust_driver_vdp.py` | two uncertain parameters, a state constraint, an expected cost | the robust design must hold the barrier for every plant |
| `robust_driver_risk.py` | expectation, mean-variance and CVaR on one problem | each objective is recomputed from the per-scenario costs |
| `robust_driver_estimate.py` | estimate the plant, then design against its covariance | the design must hold at the true plant, which it never saw |
| `robust_driver_tube.py` | an ancillary feedback gain against open loop | the closed loop is checked against an independent implementation |
| `robust_driver_gain.py` | who chooses the ancillary gain: given, scheduled, co-designed | a measured negative result --- the co-designed gain must be cheaper and must fail to certify |
| `robust_driver_cvar.py` | a tail measure on a set dense enough to resolve its tail | CVaR must be WORSE in the tail than the mean-minimising design at sixteen scenarios and better at forty-eight |
| `rv2oe_casadi.py` | orbital-element helper used by `launch.py` | not an example |

The first fourteen were run together on 17 September 2026 and passed;
`robust_arm.py` and five driver examples were added on 27 September 2026, and
`robust_driver_gain.py` and `robust_driver_cvar.py` on 28 September 2026. All twenty-two
were run together on 28 September 2026 and passed. `run_all.py` takes a name filter if you want a
subset.

The robust examples are the slow ones, because each calls `prob.solve` once per
iteration of a scenario-generation loop and integrates thousands of trajectories to
verify the result. Most take between forty seconds and three and a half minutes;
`robust_driver_cvar.py` is about two and a half, being four solves of a
forty-eight-scenario problem. `robust_driver_gain.py` is the outlier at three to eight
minutes depending on the build: it co-designs a feedback gain, which is a harder NLP than
any other example here, and the point of the example is what that buys. The C++
counterpart of the first, `examples/robust_arm/`, is the fuller study and takes minutes.

## Robust optimal control

`psopt.robust` turns a problem whose dynamics, path constraints and events depend
on an uncertain parameter into one deterministic problem --- M copies of the
state against one copy of the control, which is what makes the design
non-anticipative --- and wraps an outer loop around it. Nothing in the library is
modified; the driver assembles, calls `Problem.solve`, verifies, and calls it
again.

```python
from psopt.robust import RobustProblem, Gaussian

rp = RobustProblem(name="arm")
ph = rp.add_phase(nstates=4, ncontrols=2, nevents=8)
ph.dynamics = lambda x, u, p, t, th: ...        # th is the uncertain parameter
ph.events   = lambda xi, xf, p, t0, tf, th: ...
rp.uncertainty   = Gaussian(mean=[0.5], cov=[[0.15 ** 2]], truncate=3.0)
rp.initial_state = [0.0, 0.0, 0.5, 0.0]

out = rp.solve(alg, slack=0.0, risk="expectation", generate=True)
out.certificate      # the worst violation found anywhere in the set, and where
out.out_of_sample    # scored on parameters that took no part in the design
out.wait_and_see     # E[min] <= min E[.], the value of knowing theta in advance
```

`slack` is how much violation **of the bounds you declared** a design may leave,
not the constraint itself. A terminal ball of radius 0.02 belongs in the event
bounds; passing 0.02 as `slack` as well would quietly ask for a ball of radius
0.04. Zero demands the declared bounds hold everywhere in the set, which the
inward `margin` on the design is what makes attainable.

Scenario sets come from an unscented rule (`scenarios="sigma-points"`), a
low-discrepancy sequence (`"qmc"`), or a list (`Explicit`). With `generate=True`
the driver then adds scenarios where the design is actually failing, which is a
cutting plane on the semi-infinite constraint and is what lets it finish with a
certificate rather than a statistic.

`risk=` selects what is minimised: `"nominal"`, `"expectation"`,
`"mean-variance"` or `"cvar"`. The first two are weighted sums of per-scenario
costs and go straight into the integrand. The other two need each scenario's cost
as a quantity in its own right — the variance of a Lagrange cost across scenarios
is not the integral of anything — so the driver carries an extra state per
scenario for it, and needs `.cost_bounds`. CVaR is written by the
Rockafellar–Uryasev device with a slack static parameter per scenario rather than
a smoothed hinge, so the constraints are exact.

**Feedback.** By default a design is open loop: one control history, committed
before the uncertainty is revealed, serving every plant in the set unaided. That
is honest and it is expensive — on the arm it costs a factor of about three in
final time. Setting `.feedback` to a gain `K` turns the design into a tube,

    u_k(t) = u_bar(t) + K ( x_k(t) - x_ref(t) )

where `u_bar` is still the only decision variable and the reference is scenario 0
of the rule, which runs open loop by construction. The design remains
here-and-now: both `u_bar` and `K` are fixed before the uncertainty is revealed.
On the arm the price of robustness falls from +185% to +13% and fewer than half as
many scenarios are needed to certify it.

Three things come with it. The realised control differs from scenario to
scenario, so its bounds become path constraints and the driver adds them — set
`ms_path_samples` so they hold *inside* the segments too, or the correction
breaches the actuator limits between the nodes, which `certificate["control_excess"]`
reports. The verification integrates the closed loop, carrying the reference as
an extra column of the same sweep, so it assumes the controller regenerates the
reference by integrating the nominal model rather than storing it at the design's
node spacing. And the gain needs the state to be measurable, which the open-loop
design does not.

**Where the gain comes from, and who should not choose it.** `.feedback` takes four
kinds of value: a matrix, for a given constant gain; a callable of `t`, for a given
schedule `K(t)`; `"co-design"`, for a constant gain whose entries are optimised as
static parameters alongside the trajectory; and `"co-design-schedule"`, for a
time-varying gain optimised as extra *controls*, which is the representation to use
because the transcription already gives a control a time profile at the resolution
the trajectory has. The two co-designed forms need `.feedback_bounds = (lo, hi)`,
where each may be a scalar or an `(ncontrols, nstates)` array, and take
`.feedback_guess` as a starting point.

Co-design is easy to set up and does not work, and `robust_driver_gain.py` is the
measurement. Every co-designed variant tried on the arm — constant and scheduled,
bounds from wide to a box around a gain that certifies, scenario sets of three, nine
and sixteen points — came out cheaper on the objective than the given LQR gain and
failed its certificate, by between one and five orders of magnitude against a slack
of 1e-3. The reason is structural. A scenario set enters the design as a constraint
set, so the gain is rewarded for making those M plants cheap and charged nothing for
what it does to the rest of the family; and a gain multiplies a deviation that is
itself a function of the uncertain parameter, so its leverage on an unsampled plant
is bounded by nothing the design can see. Density does not fix it and a tighter
bound does not fix it. A gain has to be chosen by a criterion that quantifies over
the whole family — which is what solving a Riccati equation does and what minimising
over a finite sample cannot. `riccati_gain(out, Q, R)` is provided as the diagnostic
rather than as a design tool: it sweeps the Riccati equation backwards along a design
and reports how much the ideal gain actually varies, together with the error of
fitting it as a polynomial in `t`.

The driver therefore offers co-design and then measures it honestly.
`certificate["gain_on_bound"]` counts the entries that finished on their bound, which
is where an unconstrained co-design goes, and the report says what that means;
`certificate["integrator_drift"]` comes out as large as the violation on a
co-designed design, because an aggressive gain makes the closed loop stiff and a step
count adequate for the nominal problem is not adequate for that.

**Where the uncertainty should come from.** A covariance somebody chose is the
weakest part of a robust design. PSOPT's own parameter estimation returns one —
`sol.parameter_statistics.covariance` — and that is the distribution the design
ought to run against. `robust_driver_estimate.py` closes that loop inside one
model: estimate two arm parameters from noisy observations of a prescribed
manoeuvre, hand the covariance to the driver, and check the result on the true
plant, which nothing in the estimate or the design ever saw. The estimated
posterior is strongly correlated, and the example shows that discarding the
off-diagonal term is not a conservative simplification — it points the
uncertainty set along the wrong axes and leaves the design exposed in the one
direction the data did not determine.

Two cautions about the risk measures, both measured rather than asserted in
`robust_driver_risk.py`. A **tail measure needs a scenario set that resolves the
tail**: CVaR over the five-point unscented rule is a worst-case measure over
those five plants, because the rule matches a mean and a covariance and says
nothing about a tail. And the scenario set is limited by the transcription — see
below.

Three properties are worth knowing before relying on it.

*The verification does not use the design's integrator.* It builds CasADi
functions from the same user equations and integrates them separately: a
vectorised fixed-step RK4 for the search, SciPy's adaptive DOP853 for the
reported numbers, and the disagreement between the two is measured and printed.

*The objective a run returns is an upper bound, not the value of the problem.* The
warm chain the loop needs while the scenario set is small and growing also
conditions where it finishes: the loop follows one homotopy through a nonconvex
problem and lands on a local minimum the final scenario set does not require. Both
drivers therefore re-solve the final set cold, from the guess you already supplied,
and keep that design when it certifies and is cheaper (`polish=True` here,
`spec.polish` in C++). Measured on the arm at 25 nodes: the C++ loop finishes at
t_f = 8.9713 and the cold re-solve of its own twelve scenarios reaches **7.6387 with
no violation anywhere in the set**, 15% faster with a cleaner certificate, verified
over 4001 payloads by an integrator independent of the transcription.

The rule never trades a certificate for a cheaper number. In `robust_driver_arm.py`
the cold solve comes out at 8.6063 with a worst violation of 7.7e-03, seven times
the slack, so polish refuses it and the warm design at 8.9670 stands. That is also
why the two drivers end 17% apart on the same problem: on one scenario set the cold
solve certifies and on the other it does not. Both answers are certified designs.
**Read the certificate, not the objective**, and where the final time genuinely
matters, run the loop at more than one `seed` and keep the cheapest design that
certifies.

*A certificate is a claim about a search.* For one uncertain parameter the seeding
is effectively exhaustive; in more than one it is not, and `out.certificate` says
"nothing worse was found" rather than "nothing worse exists". The inner problem is
not concave.

*Not every problem can be robustified open-loop.* If the sensitivity of the
constrained state to the uncertain parameter obeys a linear equation whose
homogeneous part the control cannot shape --- a linear plant with an uncertain
gain and a pinned terminal state, or a kinematic vehicle with an uncertain speed
--- then that sensitivity is fixed by the boundary conditions and no control
reduces it. `robust_driver_vdp.py` sets out the argument. Check for it before
reaching for the driver.

**How many scenarios will fit.**

    degrees of freedom  =  free initial states + ncontrols * nodes + free times
                           + free parameters - equality events

A scenario brings `nstates` values at t0 that the defect equations do not
determine, and `nstates` pinned initial conditions that determine them. **Those
cancel**, so a scenario is free of itself, and what costs one is an equality event
*beyond* its initial condition — a pinned terminal state above all, which one
open-loop control cannot meet for several different plants anyway. Relaxing such an
event to a tolerance costs nothing at all, an inequality event taking no degree of
freedom. The count is a lower bound, because an equality event can be linearly
dependent on the defects, so it is used to explain a failure and never to refuse in
advance; left to IPOPT the failure arrives as return code -10 after the whole
problem has been assembled and taped.

There is one exception, and it is real: a transcription that collocates **every**
stored node holds one defect condition per stored state value, and the initial
condition is one more. The global Lobatto schemes — `Legendre` and `Chebyshev` — are
in that position, so a replicated state costs them `nstates` per scenario whatever
the mesh. On the arm they are refused above twelve scenarios at any node count,
where multiple shooting and `trapezoidal` collocation carry eighty. Use one of those
for a large scenario set; the driver warns if you do not.

**What this said until 28 September 2026, and why it was wrong.** The free initial
states were missing from the sum, and the formula nevertheless predicted the
observed wall exactly — the arm at 25 nodes refused above twelve scenarios, a
two-state problem at 31 nodes above fifteen. Both were phantom. PSOPT's defect block
holds `nstates * nodes` rows and a scheme with `nodes - 1` intervals fills only
`nstates * (nodes - 1)` of them; the rest were written as zeros with bounds
`[0, 0]`, and IPOPT counts equality rows against variables. The arm had 51 degrees
of freedom at every scenario count and was refused for having -9. `Alg`'s
`free_padded_defect_rows` frees those rows; it is off by default in PSOPT, because
it moves the iterate path on problems that are already rank-deficient, and **the
driver turns it on**, the phantom wall being the whole of what limited it. Measured
on the arm at 25 nodes, one cold solve per cell:

| transcription | before | after |
|---|---|---|
| multiple shooting, RK4 x 12, linear | refused above 12 | 80 scenarios, 305 s |
| `trapezoidal` | refused above 12 | 80 scenarios, 125 s |
| `Hermite-Simpson` | refused above 20 | 80 scenarios, 97 s |
| `Legendre` | refused above 12 | refused above 12 — genuinely |

Hermite-Simpson reached twenty before because its midpoint controls pay for the
phantom rows. Trapezoidal is the cheap choice for a large set.

**What control the verifier integrates.** A transcription decides not only where the
control is a decision variable but what the control *is* between those points, and the
verification integrator has to use the same reading or it is measuring a different
controller. The driver now reads it the way the transcription means it:

| transcription | reading | exact? |
|---|---|---|
| multiple shooting, `"constant"` | held across the segment | by construction |
| multiple shooting, `"linear"` | the ramp the integrator used | by construction |
| multiple shooting, `"quadratic"` | the parabola through node, midpoint, node | by construction |
| `Hermite-Simpson` | the same parabola | by construction |
| `trapezoidal` | straight lines | by convention — the scheme's own quadrature |
| `Legendre`, `Chebyshev` | straight lines | no: a degree-N polynomial read as chords |

The midpoint values come from `solution.controls_full`, the complete control history,
which PSOPT fills for exactly the two discretizations that carry a midpoint control. The
driver checks that its even columns are the nodal table before trusting them, and falls
back to straight lines with a warning if they are not.

What it is worth, on the arm at 25 nodes with five sigma points — one design per
transcription, each verified twice, once with the matched reading and once with the
blanket straight line the driver used before:

| transcription | matched | as a line | ratio |
|---|---|---|---|
| multiple shooting, `"constant"` | 1.511e+00 | 3.475e+00 | 2.3 |
| multiple shooting, `"linear"` | 1.380e+00 | 1.380e+00 | **1.0** |
| multiple shooting, `"quadratic"` | 2.947e+00 | 1.843e+00 | **0.6** |
| `Hermite-Simpson` | **8.301e-02** | 2.631e+00 | **31.7** |
| `trapezoidal` | 2.827e+00 | 2.827e+00 | **1.0** |

Hermite-Simpson designs were being charged thirty-two times the violation they have.
The two rows at 1.0 are the check that the change is a no-op where the old reading was
already right, which is every shipped example. And the row at **0.6** is the reason this
matters more than accuracy: reading a parabola as a chord *flatters* a design as easily
as it damns one, the parabola overshooting outside the chord, so the old reading was not
conservatively pessimistic — it was optimistic for the quadratic parameterisation. A
verifier that can flatter a design is worse than one that penalises it.

**What the density buys, on the risk measures.** The wall mattered most where a
scenario set has to resolve a *tail*, and `robust_driver_cvar.py` is the worked example:
CVaR at α = 0.95 on the van der Pol problem — the tail is the worst one plant in twenty —
scored on 2889 plants no design saw.

| scenarios | E[J] design's out-of-sample CVaR | CVaR design's | |
|---|---|---|---|
| 16 | 4.983 | 5.238 | CVaR is **+5.1% worse** |
| 48 | 5.103 | **4.635** | CVaR is **−9.2% better**, paying +19.2% on the mean |

At sixteen points a 5% tail holds eight tenths of a point, so minimising CVaR over it is
minimising which single scenario happens to be worst, and the design that comes out is
worse in the real tail than the one that simply minimised the mean. **Under the old
counting this problem carried twenty-two scenarios and was refused at twenty-four**,
CVaR's cost state and Rockafellar–Uryasev slacks included. The density at which the
measure starts to work sat just past the density the arithmetic refused.

The reported value is optimistic, as a sample average over the scenarios that were
optimised should be: 4.410 in sample against 4.635 realised, −4.8%. **The reported value
of a risk measure is never the certificate.**

At α = 0.75 the effect is much weaker and does not clearly favour CVaR at any density
tried, which is consistent with `robust_driver_risk.py`'s own account: at a gentle level
the tail is not where the action is.

Not implemented: scenario-dependent static parameters, and multi-phase problems,
which are refused rather than quietly mishandled. A given schedule `K(t)` is
implemented but is not recommended: the front end emits the maths for CppAD to tape
and refuses a branch on a symbolic value, so every interpolation is out and the
schedule has to be a polynomial in `t`. On the arm the Riccati gain needs degree 9
to fit to 3%, degree 3 misses by a factor of three, and at degree 9 the first solve
does not converge. `"co-design-schedule"` exists because a gain carried as extra
controls needs no basis and no branch; it is subject to the same objection as any
other co-designed gain.

One figure has moved and is flagged rather than quietly updated. `bryson_ir.py`
reports the integrated residual, which is a feasibility measure rather than a
cost, and this README previously recorded 3.498289481767e-04 for it, matching a
native C++ driver with the same options. It now returns 5.116118039065e-05. The
integrated-residual transcription has changed repeatedly since that figure was
taken, so the likeliest explanation is that the old number is simply stale -- but
the pass-through comparison against a native driver has not been repeated, so
that is an inference and not a measurement.

## Provenance

Co-developed with AI assistance (Claude). All results validated against native
PSOPT baselines.
