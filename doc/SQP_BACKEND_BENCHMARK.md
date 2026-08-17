# QP backends for PSOPT's SQP solver: a comparison

This records one sweep, run in one configuration and one sitting, of the five QP
backends the SQP solver can use, over ten of the shipped examples. It is here because
the numbers are the reason the backend choice exists, and because they were expensive
to obtain: assembling them from separate runs, as I first tried, produced a table with
silent fallbacks in it and had to be thrown away.

## How it was run

Every run used `algorithm.nlp_method = "SQP"` with `algorithm.hessian = "exact"`, the
same mesh refinement the example ships with, and a limit of 240 seconds. The IPOPT
column is the same examples under their own solver, as the reference for both the
answer and the time. Two runs shared two cores throughout, so the times compare with
each other and not with a machine at rest. GALAHAD needs `OMP_CANCELLATION=TRUE` and
`OMP_PROC_BIND=TRUE` in the environment; both were set for every run so that the
comparison is uniform.

Entries are iterations, summed over the mesh refinements, and seconds. "fail" means the
solver stopped without reaching optimality on every mesh; "timeout" means it was still
going at 240 seconds; "crash" means the process died. An (a) marks an answer that
differs from IPOPT's by more than 1e-4 relative. The sizes are those of the first mesh,
which is the only one every backend reaches.

| example | variables | Jacobian nnz | density | IPOPT | qpOASES | ProxQP | QPALM | OSQP | GALAHAD |
|---|---|---|---|---|---|---|---|---|---|
| brac1 | 162 | 5133 | 25.1% | 1 s | 14 / 5 s | 14 / 1 s | 15 / 1 s | fail | 20 / 2 s |
| bryson_denham | 202 | 7813 | 24.8% | 1 s | 91 / 93 s (a) | 38 / 8 s | 230 / 10 s | fail | 35 / 7 s |
| hypersensitive | 102 | 2654 | 49.1% | 1 s | 10 / 0 s | 10 / 0 s | 10 / 0 s | 10 / 0 s | 10 / 0 s |
| lts | 152 | 3857 | 19.8% | 0 s | 19 / 1 s | 22 / 0 s | fail | fail | 21 / 3 s |
| interior_point | 84 | 895 | 22.7% | 0 s | 2 / 0 s | 2 / 0 s | 14 / 0 s | 2 / 0 s | 2 / 0 s |
| manutec | 326 | 3623 | 4.5% | 2 s | 103 / 134 s | 12 / 1 s | timeout | 59 / 8 s | 12 / 1 s |
| glider | 202 | 1257 | 3.7% | 3 s | timeout | timeout | fail | timeout | 388 / 19 s |
| shuttle_reentry | 482 | 4259 | 2.4% | 4 s | timeout | 186 / 35 s | timeout | 158 / 8 s | timeout |
| launch | 658 | 10283 | 2.8% | 6 s | timeout | timeout | crash | fail | timeout |
| low_thrust | 803 | 10295 | 2.0% | 36 s | timeout | timeout | timeout | fail | timeout |

The same data, unrounded and with the objective each run returned, is in
`SQP_BACKEND_BENCHMARK.csv` beside this file.

## What it shows

**The sparsity is where the interest is.** The first five examples are Legendre
pseudospectral and their Jacobians are 20 to 49 per cent dense; the last five are local
collocation at 2 to 5 per cent. Those two groups behave quite differently, and a
comparison drawn only from the first five -- which is what every earlier round of this
work had -- says almost nothing about the second.

**Two of the large problems are solved.** glider under GALAHAD in 388 iterations and 19
seconds, shuttle_reentry under ProxQP in 186 and 35 and under OSQP in 158 and 8, all at
IPOPT's answer. Nothing solved shuttle_reentry before. What made the difference was not
a faster backend but a correction to the inertia test that decides how much the Hessian
is shifted: it had been demanding that the KKT matrix have no zero eigenvalues, which a
shift of the Hessian cannot arrange when the Jacobian is rank-deficient, so on a problem
with a dependent constraint row the shift was pinned at its ceiling at every iteration
and the model reduced to a multiple of the identity.

**ProxQP and GALAHAD are the two to build on.** Between them they solve eight of the
ten cells that anybody solves outright, and they are the fastest two on five of the six
small examples -- manutec in 12 iterations and a second against qpOASES's 103 and 134.

**qpOASES times out on all five large problems**, which is the quadratic memory and
cubic work of its dense null-space factorisation arriving on schedule. It is still the
default, and the case for that is now only that it is the one backend with no external
dependency.

**QPALM crashes on one and fails or times out on five.** The crash is inside its own
factorisation and no wrapper can contain it.

**OSQP fails on five**, which is what a first-order method does when asked for 1e-8.
Where it does converge it is very fast -- shuttle_reentry in 8 seconds, the quickest
cell in the table.

## Two settings the table does not vary

The sweep above is one configuration: `algorithm.qp_restoration = "elastic"` and
`algorithm.elastic_penalty = "weights"`, which are the defaults. Both control what
happens when the linearised constraints are inconsistent and the subproblem has to be
relaxed, which on the small examples is occasional and on launch is every iteration.
They are settings rather than fixes because the evidence does not point one way.

`qp_restoration` chooses the shape of the relaxation. `"elastic"` gives every
constraint a pair of non-negative slacks, n + 2m variables, and can relax one row
without relaxing another. `"relaxation"` is Betts's, section 2.7: one variable in
[0, 1] that relaxes every row by the same fraction of its own infeasibility, n + 1
variables. Flexibility against size. On the small pseudospectral examples flexibility
wins and it is not close -- with `"relaxation"`, bryson_denham does not converge and
lts goes from 19 iterations to 28 under qpOASES and 74 under the multiplier pricing.
On launch, where 282 of 560 constraints are violated at the start, size wins by as
wide a margin: the elastic subproblem is 1778 variables, ProxQP took between five and
twenty seconds over each one, and timing the parts put 99 per cent of the solver's
wall clock there against 0.03 seconds for the plain subproblem and 0.02 for the
inertia test. The relaxed subproblem takes GALAHAD 20 iterations where the elastic one
takes the 1002 at which it gives up. Neither setting solves launch.

`elastic_penalty` chooses what the relaxation costs. `"weights"` prices it above the
merit function's penalty weights; `"multipliers"` prices it above the multipliers as
well, which is the classical condition on elastic mode and the theoretically correct
one. The default is the other one, for a reason worth recording: the rule was written
for the l1 merit function, whose weights have to dominate the multipliers anyway, so
pricing from the weights got the bound for free. Betts's augmented Lagrangian carries
the multipliers itself and its weights are least-norm -- machine epsilon for most of a
run -- so the price collapsed to its floor of ten, and on launch the same linearised
constraints a unit Hessian satisfied to a slack of 1.3e-02 came back one solve later
at 1.1, the objective having bought the difference.

Repairing it is a wash on this benchmark: 26 cells solved either way. `"multipliers"`
fixes bryson_denham under qpOASES, which goes from 91 iterations at an answer 0.09 per
cent from IPOPT's to 31 at IPOPT's exactly, and brac1 and lts each gain a backend;
against that OSQP loses interior_point and manutec, and ProxQP's bryson_denham answer
drifts by 0.3 per cent. Two things have to move with it or it is worse than useless:
the weights are no longer raised to the price after a restoration, which under an
augmented Lagrangian makes the merit stiff rather than exact, and the slacks'
regularisation scales with the price instead of sitting at 1.0e-08. Each of the three
alone breaks bryson_denham; they are one mechanism.

## The FM strategy, which is a different table

Everything above is `algorithm.sqp_strategy = "M"`: minimise from the starting guess,
taking the constraints and the objective together at every iteration. Betts's own
default is FM -- find a point feasible with respect to the constraints first, ignoring
the objective, the multipliers and the Hessian entirely, then minimise from there -- and
on the three backends worth keeping it looks like this. Entries are feasibility-phase
iterations plus optimality-phase iterations, and seconds.

| example | qpOASES | ProxQP | GALAHAD |
|---|---|---|---|
| brac1 | 5 + 5 / 2 s | 5 + 5 / 0 s | 5 + 5 / 1 s |
| bryson_denham | 6 + 8 / 14 s | 6 + 16 / 116 s (a) | fail |
| hypersensitive | 7 + 6 / 0 s | 7 + 6 / 0 s | 7 + 6 / 0 s |
| lts | 8 + 14 / 1 s | 8 + 14 / 0 s | 8 + 15 / 1 s |
| interior_point | 1 + 2 / 0 s | 1 + 2 / 0 s | 1 + 2 / 0 s |
| manutec | timeout | 11 + 11 / 1 s | 11 + 11 / 2 s |
| glider | timeout | timeout | 19 + 364 / 16 s |
| shuttle_reentry | timeout | 41 + 136 / 8 s | 39 + 70 / 13 s |
| **launch** | timeout | timeout | **24 + 23 / 57 s** |
| low_thrust | timeout | timeout | timeout |

**launch is solved, for the first time by anything here, at IPOPT's answer to seven
figures.** Its starting guess violates 282 of 560 constraints; under strategy M the
linearisation is inconsistent at every iterate, restoration runs every iteration, and 99
per cent of the wall clock goes into a relaxed subproblem that is three times the size of
the real one. The feasibility phase drives the maximum violation from 7.1e-01 to 2.1e-11
in 23 iterations of a subproblem with no objective in it, switching to the relaxation at
iteration 3 when the least distance program fails and back at iteration 17 when it takes
a full step, exactly as section 2.8.2 describes. The optimality phase then takes 23
iterations from a point its model can be trusted at.

**It also settles the backend question.** With FM, GALAHAD solves eight of the ten,
against ProxQP's six and qpOASES's five, and it is the only backend that reaches glider
or launch at all.

The small examples improve throughout, which was not the point of the exercise but is
worth recording: bryson_denham under qpOASES goes from 91 iterations and 93 seconds to 14
and 14, brac1 from 14 to 10, interior_point from 2 to 3 but with the first two spent
proving feasibility. Splitting the work in two makes each half easier than the whole.

Two cells go the other way. manutec under qpOASES was solved under M and times out under
FM. bryson_denham under GALAHAD stops during the optimality phase -- the feasibility
phase finds its point in six iterations, and the run then hangs inside a GALAHAD
subproblem, which is a backend fault rather than a strategy one. Neither is a reason to
prefer M: qpOASES is the backend being retired, and the hang is the same one that made
these sweeps need SIGKILL rather than SIGTERM.

The default remains "M", because changing a default changes behaviour for everyone who
has a working script, and because the backend question should be settled first. On the
evidence here the configuration to recommend is `sqp_strategy = "FM"` with
`qp_solver = "GALAHAD"`.

## Two harder problems, outside the sweep

`low_thrust` is in the table above and fails there. `zpm` -- the zero propellant manoeuvre
of the International Space Station, Bedrossian's problem and one of Betts's own showcases
-- is not in the sweep at all, because IPOPT alone needs 153 seconds and 725 iterations on
it and the ten-example sweep is long enough already. Both were run under FM with GALAHAD
and a 25-minute budget.

| | IPOPT | SQP, FM + GALAHAD |
|---|---|---|
| zpm | 6.680107e+06, 153 s | **6.680106e+06**, 1048 s |
| low_thrust | -2.203403e-01, 28 s | -2.061760e-01, 1016 s |

**zpm agrees with IPOPT to seven significant figures**, and its discretisation error comes
out slightly the better of the two, 9.27e-04 against 1.01e-03. Under strategy M it never
leaves the first mesh.

Neither run converged on every mesh, and both failed in the same place, which is the
useful part:

| mesh | zpm | low_thrust |
|---|---|---|
| 1 | feasible in 38, then the QP failed after 135 | feasibility hit its 1000-iteration limit |
| 2 | feasible in 3, then the iteration limit at 1000 | feasible in 6, then the limit at 1000 |
| 3 | feasible in 3, solved in 64 | feasible in 7, solved in 56 |
| 4 | feasible in 2, solved in 22 | feasible in 5, solved in 14 |
| 5 | feasible in 2, solved in 3 | feasible in 3, solved in 17 |

The last three meshes are easy on both. All the difficulty is on the first two, and once
there is a solution to start the next mesh from, the solver walks the rest of the way. zpm
still reaches IPOPT's answer because mesh refinement recovers; low_thrust does not, and its
6 per cent error is the arithmetic of building meshes 3 to 5 on a mesh 2 that never
converged.

Two specific things are responsible, and both are now visible in the traces rather than
inferred.

**The feasibility phase converges linearly on low_thrust's first mesh.** The steplength is
1.00 throughout the tail -- the line search is not cutting anything -- and the violation
still falls by only about a tenth per iteration, from 1.26e+01 to 3.3e-03 over 47
iterations and not to the tolerance within a thousand. A least distance step that satisfies
the linearised constraints exactly should do better than that near a feasible point, so
either the linearisation is a poor model at that scale or the step is being limited by the
variable bounds it has to respect.

**The optimality phase plateaus on the early meshes.** On low_thrust's second mesh the
violation reaches 2.7e-13 and the dual error freezes at 3.18e-05 -- a factor of thirty from
the tolerance -- with the objective static to nine significant figures, a full steplength,
and the Hessian shifted at every iteration.

Chasing this produced one fix and two dead ends, and the dead ends are worth recording so
that nobody spends the afternoon twice.

The fix: on zpm a relaxed subproblem came back from GALAHAD with multipliers of 2.1e+50,
and at the next iteration with exactly zero -- an overflow and the NaN cascade after it.
Those multipliers enter the Hessian of the Lagrangian, whose Gerschgorin bound then reads
-2.8e+50, and the Levenberg shift scaled by that bound is carried forward into every
iteration that follows, decaying by a fifth each time and needing seventy iterations to
come back to a sane value. Seventy iterations of a model that is a multiple of the
identity. A QP answer with a non-finite entry in it is now rejected, which is what it
deserves: the restoration that follows is a step the solver knows how to take. With the
guard, the largest multiplier zpm produces is 1.0e+06 and the largest shift 1.7e+08,
both in scale with the problem.

The first dead end: capping the carried shift at the current Gerschgorin ceiling, which is
correct in principle -- a shift larger than |sigma| + 1 makes the model's own curvature
irrelevant -- and which breaks bryson_denham under strategy M. Measured and withdrawn.

The second: the trust-region multipliers. A variable whose step is pinned by the trust
region rather than by a real bound has its bound multiplier set to zero, and it looked as
though the dual error was then picking up an unremovable residual at exactly those
variables. Instrumenting it says otherwise -- the dual error computed over the unpinned
variables alone is the same number, to four figures, at most iterations.

What the instrumentation does show is that zpm's first mesh is not converging quietly at
all: the dual error passes through 9.3e+11, 1.2e+15 and 1.7e+13 at successive iterations,
with between 77 and 187 of its 283 variables pinned by a trust region that has shrunk to
1.6e-02.

That looked like the thread to pull, and pulling it produced two more dead ends and one
finding. All four attempts are listed here because each cost an hour and none of them
needs repeating.

**Not the unconverged multipliers.** Instrumenting the backend's return status showed that
forty-one of zpm's hundred and nine subproblems stop rather than solve, and that the
1.2e+15 multipliers came from one of them, where the current multipliers were 5.2e-03.
Refusing the duals of a subproblem that did not converge is the obvious response and is
what Betts prescribes for the first iteration, where there is no estimate to be had
either. It costs bryson_denham its solution outright. Refusing them only when they jump by
more than four orders of magnitude is neutral on all eight examples that solve -- and takes
zpm from three meshes of five to none. Suppressing those multipliers makes zpm worse, not
better, so they are a symptom.

**Not the iteration budget.** GALAHAD's stopping reasons on zpm are 258 solved, 21 at its
iteration limit and 9 declaring the subproblem primal infeasible; the last is the honest
signal that sends the solver to restoration and is working as intended. Raising the budget
does not rescue the 21. In two hundred seconds at a budget of 200, QPA solves 276
subproblems and stops at its limit on 25; at 3000 it solves 34 and stops on 38. The
subproblems that hit the limit are ones it cannot solve, not ones it needs longer for, and
letting them grind costs the ones behind them.

**The finding, since repaired.** The GALAHAD plugin never set `control.maxit`, so the
`max_iter` the plugin interface carries never reached it: every subproblem ran to
GALAHAD's own default while PSOPT believed it had set a budget. Three of the four backends
honoured it; GALAHAD alone did not. Passing the value through is one line and is not safe
by itself, because PSOPT's default of 200 is tighter than GALAHAD's own and bryson_denham
under GALAHAD fails at it. The budget is now `algorithm.qp_iter_max`, and the default was
chosen by sweeping it -- 200, 500, 1000 and 2000, over eight examples and two backends:

| budget | GALAHAD | ProxQP |
|---|---|---|
| 200 | bryson_denham fails; shuttle_reentry 113, launch 35 | lts 22, otherwise as below |
| 500 | bryson_denham 38; shuttle_reentry 88, launch 28 | lts 24 |
| **1000** | **bryson_denham 35; shuttle_reentry 70, launch 23** | **lts 28** |
| 2000 | bryson_denham 35; shuttle_reentry 136, launch fails | lts 24 |

1000 is best for GALAHAD on every cell that moves and costs ProxQP six iterations on lts.
It is also GALAHAD's own default, so honouring the budget changes that backend's behaviour
not at all while making the option work for the other three. Above it, launch fails and
shuttle_reentry doubles: a subproblem allowed to grind is a subproblem taken away from the
ones behind it, which is Betts's argument in section 2.7 arriving as a number.

So the open question is narrower than it was but still open: on zpm's first mesh, seven per
cent of subproblems are ones QPA cannot solve at any budget, and the multipliers they
return are noise that the solver nonetheless needs. What is not yet known is whether those
subproblems are hard because of the problem, because of the relaxation that produced them,
or because of the trust region that shaped them.

## The line search, replaced

Halving was the crudest thing that works. Betts fits a quadratic and a cubic to the merit
function and imposes the Wolfe condition to stop steplengths becoming too small (section
2.6.1); what is implemented here is the standard safeguarded interpolation -- a quadratic
through phi(0), phi'(0) and the first rejected trial, a cubic through those and the
second, each new trial confined to [0.1, 0.5] of the last. Betts's Wolfe condition needs
the slope at the trial point, which for this merit function means a gradient and a whole
Jacobian there, one full derivative evaluation per trial; the purpose it serves is to stop
the steplength collapsing, and a floor does that for nothing.

| example | halving | interpolation |
|---|---|---|
| brac1 | 14 / 4 s | 17 / 5 s |
| **bryson_denham** | **91 / 74 s, 3.995918 (a)** | **18 / 12 s, 3.999539** |
| hypersensitive | 10 / 0 s | 11 / 0 s |
| lts | 19 / 1 s | 20 / 1 s |
| interior_point | 2 / 0 s | 2 / 0 s |
| glider (FM) | 364 / 16 s | 364 / 22 s |
| shuttle_reentry (FM) | 70 / 12 s | 70 / 13 s |
| launch (FM) | 23 / 56 s | 42 / 179 s |

bryson_denham is why this is kept. It has been the weak cell in this table for four
commits -- ninety-one iterations at an answer 0.09 per cent from IPOPT's -- and it is now
eighteen iterations at IPOPT's answer exactly. That is not a faster route to the same
place; it is a different and better place.

launch is what it costs: 42 iterations against 23, and three times the wall clock. It
still solves, at the same answer. The safeguard band is the reason and it is not a free
parameter -- at [0.25, 0.5] launch takes 30 iterations instead of 42 and bryson_denham
stops converging altogether. The wider band is kept because what it buys on bryson_denham
is an answer rather than a count.

## The initial trust region, measured

The radius the exact-Hessian trust region starts at was 1.0, on the reasoning that the
variables are scaled and an O(1) region is the natural first guess. That was reasoning
rather than measurement, and the measurement disagrees.

Counting GALAHAD's return codes over zpm's first mesh, by initial radius:

| radius | solved | primal infeasible | iteration limit | failure rate | subproblems in 150 s |
|---|---|---|---|---|---|
| 0.1 | 386 | 19 | 6 | 6% | 411 |
| 1.0 | 44 | 38 | 22 | 58% | 104 |
| 10 | 161 | 23 | 21 | 21% | 205 |

A region the linearisation cannot be satisfied inside is a subproblem the backend cannot
solve, and a subproblem it cannot solve is a restoration the solver did not need. At a
radius of 1.0 more than half the subproblems on that mesh were being refused. This is
the answer to the question left open two commits ago -- whether zpm's unsolvable
subproblems were hard because of the problem, the relaxation, or the trust region. It is
the trust region.

0.3 is the default rather than 0.1 because the small dense examples want the larger of
the two: at 0.1, brac1 takes 23 iterations against 17 and lts 28 against 20. At 0.3 they
are 17 and 18, and what the reduction buys elsewhere is kept:

| example | radius 1.0 | radius 0.3 |
|---|---|---|
| brac1 | 17 / 5 s | 17 / 5 s |
| bryson_denham | 18 / 13 s | 19 / 13 s |
| hypersensitive | 11 / 0 s | 11 / 0 s |
| lts | 20 / 1 s | 18 / 1 s |
| interior_point | 2 / 0 s | 3 / 0 s |
| **manutec** (qpOASES) | **103 / 106 s** | **13 / 15 s** |
| glider (FM) | 364 / 16 s | 364 / 17 s |
| shuttle_reentry (FM) | 70 / 12 s | 70 / 13 s |
| **launch** (FM) | **42 / 195 s** | **26 / 126 s** |

manutec is the largest single change in this table since the benchmark was started, and
it came from a constant that had never been questioned.

## What it does not show

No backend is close to IPOPT, which solves all ten in a few seconds each. The SQP is
not yet a replacement for it, and this sweep does not claim otherwise.

bryson_denham under qpOASES is the one cell that is worse than the previous sweep: 91
iterations where it took 29, at an answer 0.09 per cent from IPOPT's. It is feasible to
2.4e-07 and inside the mesh tolerance, so it is a coarser solution of the same problem
rather than a wrong one, but it is the weak spot on this group and the same change that
repaired glider caused it.

launch and low_thrust remain out of reach for every backend, as they were.
