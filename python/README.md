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

Fourteen examples, in `examples/`. Run one directly, or all of them:

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
| `rv2oe_casadi.py` | orbital-element helper used by `launch.py` | not an example |

All fourteen were run together on 17 September 2026 and passed.

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
