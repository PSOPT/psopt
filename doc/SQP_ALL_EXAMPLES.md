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

An example counts as solved if it exited cleanly and either returned an objective or, in the
case of the self-verifying examples, printed PASS. Reaching four meshes of five and running
out of clock counts as a timeout, not a solve.

## What it does

| | SQP | IPOPT |
|---|---|---|
| solved | **48** | 60 |
| timed out at 120 s | 17 | 3 |
| failed | 1 | 3 |

**Forty-eight of sixty-six, against IPOPT's sixty.** Of the forty-four both solve and for
which both return an objective, forty-two agree with IPOPT to within 1e-4 relative. The two
that do not are mpec (1.3e-03) and wheat (2.3e-03); neither is a wrong answer, both are the
same answer to fewer figures.

Three examples go the other way: **delay_history, lqr_radau and path_window are solved by
the SQP and not by IPOPT** within the same limit.

Only one example now fails outright under the SQP: twolinkarm, and the section on the
constraint-scaling repair below explains why that particular row should not be read as a
result.

## Where the seventeen timeouts are, and what they mean

The limit is 120 seconds and it is doing real work in this table. launch, low_thrust and zpm
all need minutes rather than seconds. Five of the seventeen had finished some of their meshes
when the clock stopped:

| example | meshes finished | IPOPT's time |
|---|---|---|
| f8 | 2 | 2 s |
| lotka | 3 | 26 s |
| lotka2 | 2 | 33 s |
| reorientation | 1 | 21 s |
| twoburn | 2 | 12 s |

and twelve had not finished their first: alpine, bioreactor, breakwell, chance_constraint,
crane, fuller, launch, low_thrust, mineng_free_vf, mineng_target_set, singular5, zpm. alpine
and zpm defeat IPOPT within the same limit as well, and fuller defeats it outright.

The pattern across the examples where IPOPT succeeds and the SQP does not is uniform and
unsurprising: the SQP is between one and two orders of magnitude slower per problem, and
where that is enough to cross the limit it appears as a timeout. It is not that these
problems are unsolvable by this solver; it is that this solver costs what it costs.
`SQP_BACKEND_BENCHMARK.md` measures where the cost goes.

## The failures that were not failures, and the number that was not a number

The first version of this table reported six failures and named delay_history and
int_static_linear as examples that defeat both solvers -- which, being examples no NLP could
handle, looked like the most interesting thing in it. They are nothing of the sort. Both are
self-verifying programs: they solve their problem, compare the answer against a closed form,
and print PASS. Neither prints the mesh summary this survey's classifier was reading, so both
were counted as failures under both solvers. The classifier now accepts a PASS.

The second version reported four real failures -- bryson_max_range, isoperimetric, moon and
user -- and explained them by saying that all of them "begin at an enormously infeasible
point: 6.7e+07 on bryson_max_range, 6.5e+08 on isoperimetric". That explanation was wrong,
and wrong in a way worth recording. Those numbers were manufactured by PSOPT's automatic
constraint scaling, which gave a Jacobian row that is exactly zero at the initial guess a
factor of 1/sqrt(eps) = 6.7108864e+07. bryson_max_range's true initial violation is 1.0.
Thirty-four of the sixty-six examples had at least one such row. `CONSTRAINT_SCALING.md`
has the whole of it.

With that repaired, and with bryson_max_range's guess moved off the stationary point of its
own path constraint, bryson_max_range, user, isoperimetric and conic_soc all solve, and
bryson_max_range returns -1.712313 rather than the -0.251920 that a degenerate start had
been leading it into. twolinkarm moved the other way, and moon moved in, on a perturbation
of one part in a hundred million; those two are the measurement noise of this table rather
than results, and `CONSTRAINT_SCALING.md` quantifies that too.

## What this says

The SQP is now a solver one could reasonably reach for, on nearly three quarters of the
problems this library ships, at IPOPT's answers. That is a different statement from the one
the ten-example benchmark supported a fortnight ago, which was that it worked on the small
dense ones and fell over on the sparse ones.

It is not a replacement for IPOPT and this table does not pretend otherwise. IPOPT solves
fifteen problems within the limit that the SQP does not, and does the forty-four they share
considerably faster. What the SQP now offers is a second opinion that agrees to four
figures, an implementation whose every part is in this repository rather than behind a
third-party interface, and three problems that IPOPT does not solve within the same limit.

One caution about the counts. A handful of these examples take hundreds of iterations and
sit near the point where the QP subproblem stops being solvable, and two of them were shown
to change status under a perturbation of the constraint scaling of one part in a hundred
million. Forty-eight is the right number to quote, but it is not a number with three
significant figures behind it.
