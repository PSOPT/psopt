# QP backends for PSOPT's SQP solver: a comparison

This records one sweep, run in one configuration and one sitting, of the five QP
backends the SQP solver can use, over ten of the shipped examples. It is here because
the numbers are the reason the backend choice exists, and because they were expensive
to obtain: assembling them from separate runs, as I first tried, produced a table with
silent fallbacks in it and had to be thrown away.

It replaces an earlier sweep taken before three changes to the SQP itself: the exact
Hessian's modification came under the control of the inertia of the KKT matrix, the l1
merit function was replaced by Betts's augmented Lagrangian, and the Levenberg
parameter began to be carried between iterations rather than restarted at every one.
The earlier numbers are not comparable and are not kept.

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
| brac1 | 162 | 5133 | 25.1% | 1 s | 14 / 3 s | 14 / 0 s | 15 / 1 s | fail | 20 / 1 s |
| bryson_denham | 202 | 7813 | 24.8% | 1 s | 29 / 26 s | 19 / 1 s | crash | 65 / 8 s (a) | 35 / 5 s |
| hypersensitive | 102 | 2654 | 49.1% | 0 s | 10 / 0 s | 10 / 0 s | 10 / 0 s | 10 / 0 s | 10 / 0 s |
| lts | 152 | 3857 | 19.8% | 0 s | 19 / 1 s | 22 / 0 s | fail | fail | 21 / 2 s |
| interior_point | 84 | 895 | 22.7% | 0 s | 2 / 0 s | 2 / 0 s | 14 / 0 s | 2 / 0 s | 2 / 0 s |
| manutec | 326 | 3623 | 4.5% | 1 s | 103 / 89 s | 12 / 1 s | fail | 59 / 5 s | 12 / 1 s |
| glider | 152 | 937 | 4.8% | 2 s | timeout | timeout | fail | fail | fail |
| shuttle_reentry | 482 | 4259 | 2.4% | 3 s | timeout | timeout | timeout | timeout | timeout |
| launch | 658 | 10283 | 2.8% | 4 s | timeout | timeout | crash | fail | timeout |
| low_thrust | 803 | 10295 | 2.0% | 21 s | timeout | timeout | timeout | timeout | timeout |

The same data, unrounded and with the objective each run returned, is in
`SQP_BACKEND_BENCHMARK.csv` beside this file.

## What it shows

**The sparsity is where the interest is.** The first five examples are Legendre
pseudospectral and their Jacobians are 20 to 49 per cent dense; the last five are local
collocation at 1 to 5 per cent. Those two groups behave quite differently, and a
comparison drawn only from the first five -- which is what every earlier round of this
work had -- says almost nothing about the second.

**The small examples are settled.** Every one of the first six is now solved by at least
three backends, at the same answer IPOPT gives, in a number of iterations that is no
longer embarrassing: brac1 in 14, manutec in 12, interior_point in 2. brac1 in
particular used to converge under qpOASES to a stationary point 2.8 per cent above
IPOPT's answer, and no longer does. The iteration counts on this group fell by roughly
half against the previous sweep, and the change responsible was the carry-over of the
Levenberg parameter, which had been written but never read.

**ProxQP and GALAHAD are the two to build on.** They agree with IPOPT on all six of the
examples anybody solves, and they are the fastest two on five of them -- manutec in 12
iterations and a second against qpOASES's 103 and eighty-nine. ProxQP is the better of
the two on the dense pseudospectral problems and GALAHAD the more consistent; neither
gets near the large sparse ones.

**qpOASES times out on all five large problems**, which is the quadratic memory and
cubic work of its dense null-space factorisation arriving on schedule. On the small
dense ones it is no longer the most reliable either -- bryson_denham takes it 26
seconds against ProxQP's one -- and the case for it as the default is now only that it
is the backend with no external dependency.

**QPALM crashes on two and fails on three.** The crashes are inside its own
factorisation and no wrapper can contain them.

**OSQP fails on four and gets bryson_denham wrong**, returning 4.93 where the answer is
4.00. It is a first-order method being asked for 1e-8, which is the wrong question to
put to it; nothing here is a defect in OSQP.

## What it does not show

No backend is close to IPOPT, which solves all ten in a few seconds each. The SQP is
not yet a replacement for it, and this sweep does not claim otherwise.

The five large problems are worse than they were. The previous sweep had GALAHAD
solving glider in 117 iterations and 8 seconds and QPALM solving shuttle_reentry in
269 and 14; both now fail. Reverting only the change recorded here does not bring them
back, so the loss belongs to the inertia control or to the augmented Lagrangian merit
function rather than to the Levenberg carry-over, and it is the next thing to find.
What the large problems gained in exchange is reach rather than success: shuttle_reentry
now gets through three mesh refinements to 1240 variables before the clock stops, where
before it did not leave the first.
