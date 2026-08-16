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

## What it does not show

No backend is close to IPOPT, which solves all ten in a few seconds each. The SQP is
not yet a replacement for it, and this sweep does not claim otherwise.

bryson_denham under qpOASES is the one cell that is worse than the previous sweep: 91
iterations where it took 29, at an answer 0.09 per cent from IPOPT's. It is feasible to
2.4e-07 and inside the mesh tolerance, so it is a coarser solution of the same problem
rather than a wrong one, but it is the weak spot on this group and the same change that
repaired glider caused it.

launch and low_thrust remain out of reach for every backend, as they were.
