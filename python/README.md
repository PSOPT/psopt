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

Twenty examples, in `examples/`. Run one directly, or all of them:

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
| `rv2oe_casadi.py` | orbital-element helper used by `launch.py` | not an example |

The first fourteen were run together on 17 September 2026 and passed;
`robust_arm.py` and the five driver examples were added on 27 September 2026 and
all twenty were run together then. The suite now takes about nine minutes, most
of it in the robust examples; `run_all.py` takes a name filter if you want a
subset.

The three robust examples are the slow ones, at ten to sixty seconds each,
because each calls `prob.solve` once per iteration of a scenario-generation loop
and integrates thousands of trajectories to verify the result. The C++
counterpart, `examples/robust_arm/`, is the fuller study and takes minutes.

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
On the arm the price of robustness falls from +192% to +13% and half as many
scenarios are needed to certify it.

Three things come with it. The realised control differs from scenario to
scenario, so its bounds become path constraints and the driver adds them — set
`ms_path_samples` so they hold *inside* the segments too, or the correction
breaches the actuator limits between the nodes, which `certificate["control_excess"]`
reports. The verification integrates the closed loop, carrying the reference as
an extra column of the same sweep, so it assumes the controller regenerates the
reference by integrating the nominal model rather than storing it at the design's
node spacing. And the gain needs the state to be measurable, which the open-loop
design does not. The gain is given rather than co-designed: co-designing it makes
the problem bilinear, and a fixed ancillary gain is the standard first step.

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

**How many scenarios will fit.** Under multiple shooting the state at every node
is pinned by the defect equations, so the only free parameters in the augmented
problem are the shared control at the nodes, the free time endpoints and any free
static parameters. Every scenario brings its own *equality* events — its pinned
initial condition above all — out of that same budget:

    degrees of freedom  =  ncontrols * nodes + free times + free parameters
                           - equality events

Each scenario therefore costs as many degrees of freedom as it has pinned events,
and the scenario count is limited by the node count rather than by memory or time.
Left to IPOPT this arrives as return code -10 after the whole problem has been
assembled; the driver recognises it and prints the arithmetic with the two
remedies, which differ — raise the node count, or relax pinned events to a
tolerance, since an inequality event costs nothing here. The count is a lower
bound, because an equality event can be linearly dependent on the defects, so it
is used to explain a failure and never to refuse in advance.

Not implemented: scenario-dependent static parameters, a time-varying ancillary
gain from a Riccati sweep along the nominal trajectory (the obvious next step
after the fixed gain), co-design of the gain with the trajectory, and multi-phase
problems, which are refused rather than quietly mishandled.

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
