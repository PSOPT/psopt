# The SQP across the whole example set

Every benchmark in this repository until now has used ten examples, chosen early and
kept for comparability. Ten is enough to steer a change and not enough to say what the
solver can do. This is the whole shipped example set, sixty-six of them, run once under
the SQP and once under IPOPT.

## How it was run

`algorithm.nlp_method = "SQP"`, `hessian = "exact"`, `qp_solver = "GALAHAD"`,
`sqp_strategy = "FM"`, against the same examples under IPOPT, each with a limit of 120
seconds and each example's own mesh refinement. `microgrid` is excluded because it is
not public; `x1` and one other have no built binary. The full table, with variable
counts, mesh counts, iteration counts, times and objectives, is in
`SQP_ALL_EXAMPLES.csv` beside this file.

An example counts as solved only if every one of its mesh refinements reached optimality
and it returned an objective. Reaching four meshes of five counts as a failure here,
which is strict but is the only line that does not require judgement.

## What it does

| | SQP | IPOPT |
|---|---|---|
| solved | **43** | 59 |
| timed out at 120 s | 17 | 2 |
| failed | 6 | 5 |

**Forty-three of sixty-six, against IPOPT's fifty-nine.** Of the forty-two both solve,
thirty-eight agree with IPOPT to within 1e-4 relative. The four that do not are
coulomb (1.7e-04), twophase_schwartz (3.0e-04, on an objective that is zero to machine
precision in both, so the relative measure means nothing), mpec (3.8e-04) and wheat
(2.3e-03). None is a wrong answer; all four are the same answer to fewer figures.

One example goes the other way: **lqr_radau is solved by the SQP and not by IPOPT.**

## Where the seventeen timeouts are, and what they mean

The limit is 120 seconds, and it is doing real work in this table. launch needs about 126
seconds in the configuration measured two commits ago, so it appears here as a timeout
having been a solve; low_thrust and zpm both need many minutes. Six of the seventeen had
finished some of their meshes when the clock stopped:

| example | meshes finished | IPOPT's time |
|---|---|---|
| f8 | 2 of 3 | 1 s |
| lotka | 3 of 4 | 20 s |
| lotka2 | 2 of 3 | 27 s |
| reorientation | 1 of 2 | 17 s |
| twoburn | 2 of 3 | 9 s |
| zpm | 0 of 2 | timed out too |

and eleven had not finished their first: bioreactor, breakwell, chance_constraint, crane,
launch, low_thrust, path_window, singular5, alpine. alpine and zpm defeat IPOPT within the
same limit as well.

The pattern across the sixteen where IPOPT succeeds and the SQP does not is uniform and
unsurprising: the SQP is between one and two orders of magnitude slower per problem, and
where that is enough to cross the limit it appears as a timeout. It is not that these
problems are unsolvable by this solver; it is that this solver costs what it costs.
`SQP_BACKEND_BENCHMARK.md` measures where the cost goes.

## The six failures

Four the SQP fails and IPOPT solves: bryson_max_range, isoperimetric, moon, user.
Two neither solves: delay_history and int_static_linear.

The second pair is the more interesting: an example that defeats both solvers is more
likely to be a problem with the formulation or the initial guess than with either NLP,
and both are worth a look on those grounds alone.

## What this says

The SQP is now a solver one could reasonably reach for, on two thirds of the problems
this library ships, at IPOPT's answers. That is a different statement from the one the
ten-example benchmark supported a fortnight ago, which was that it worked on the small
dense ones and fell over on the sparse ones.

It is not a replacement for IPOPT and this table does not pretend otherwise. IPOPT solves
sixteen problems within the limit that the SQP does not, and does the forty-two they share
considerably faster. What the SQP now offers is a second opinion that agrees to four
figures, an implementation whose every part is in this repository rather than behind a
third-party interface, and one problem that IPOPT does not solve at all.
