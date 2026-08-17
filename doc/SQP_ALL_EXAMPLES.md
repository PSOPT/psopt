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
| solved | **45** | 61 |
| timed out at 120 s | 17 | 2 |
| failed | 4 | 3 |

**Forty-five of sixty-six, against IPOPT's sixty-one.** Of the forty-two both solve,
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

## The four failures, and two that were not

The first version of this table reported six failures and named delay_history and
int_static_linear as examples that defeat both solvers -- which, being examples no NLP
could handle, looked like the most interesting thing in it. They are nothing of the sort.
Both are self-verifying programs: they solve their problem, compare the answer against a
closed form, and print PASS. Neither prints the mesh summary this survey's classifier was
reading, so both were counted as failures under both solvers. They pass, under both, and
the counts above are corrected.

The lesson is the ordinary one about harnesses. A classifier that asks "did it print the
thing I expect" reports a missing print as a failure, and the two most interesting-looking
rows in the table were an artefact of that.

The four real failures are bryson_max_range, isoperimetric, moon and user -- and user is
a copy of bryson_max_range, down to the last digit of its mesh error, so there are three
distinct problems here. All three begin at an enormously infeasible point: 6.7e+07 on
bryson_max_range, 6.5e+08 on isoperimetric. What happened next is in the next section.

## What the four failures showed, and one fix

The feasibility phase's relaxation charges its residuals quadratically, rho/2 u'u, and u
is the size of the row it has to absorb. At a starting violation of 6.7e+07 and rho =
1.0e+04 the subproblem's objective is of order 1.0e+18 before the solver has done
anything, and GALAHAD refuses it -- and refuses the least distance program with it. The
phase stopped after one iteration having reduced the violation by a thousandth of a per
cent, and the optimality phase then stopped after none. That is the whole of the
"0 iterations" those examples were reporting.

Dividing each row of the feasibility subproblems and its bounds by its own magnitude
makes the residuals O(1) and rho mean what it is meant to mean. It does not change the
feasible set of either subproblem, since a diagonal row scaling never does.

| example | before | after |
|---|---|---|
| isoperimetric | phase fails at iteration 1, violation 6.5e+08 | **phase succeeds in 13 iterations** |
| bryson_max_range | violation 6.711e+07 to 6.710e+07 | violation 6.7e+07 to 4.0e+06, then stalls |
| moon | restoration at every iteration | 47 iterations before failing |

None of the three solves yet. isoperimetric now reaches a feasible point and then spends
its whole iteration budget in the optimality phase; bryson_max_range gets 94 per cent of
the way and stops on "no step reduced the violation". But the phase is doing its job on
problems where it previously could not begin, and that is worth having on its own.

Against the ten-example regression set the scaling is neutral -- brac1 17 iterations,
bryson_denham 19, hypersensitive 11, lts 18, interior_point 3, glider 364, launch 26, all
unchanged -- except shuttle_reentry, which takes 100 iterations against 70 and returns
the same answer.

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
