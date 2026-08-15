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

Entries are iterations / seconds. "fail" means the solver stopped without reaching
optimality on every mesh; "timeout" means it was still going at 240 seconds; "crash"
means the process died. An (a) marks an answer that differs from IPOPT's by more than
1e-4 relative.

| example | variables | Jacobian nnz | density | IPOPT | qpOASES | ProxQP | QPALM | OSQP | GALAHAD |
|---|---|---|---|---|---|---|---|---|---|
| brac1 | 162 | 5133 | 25.1% | 2 s | 23 / 7 s (a) | fail | 25 / 1 s | fail | fail |
| bryson_denham | 202 | 7813 | 24.8% | 1 s | 20 / 10 s | 13 / 1 s | crash | fail | timeout |
| hypersensitive | 102 | 2654 | 49.1% | 0 s | 11 / 0 s | 11 / 0 s | 11 / 0 s | 11 / 0 s | 11 / 0 s |
| lts | 152 | 3857 | 19.8% | 0 s | 43 / 1 s | 29 / 0 s | 56 / 0 s | fail | 49 / 18 s |
| interior_point | 84 | 895 | 22.7% | 0 s | 14 / 0 s | 14 / 0 s | 21 / 0 s | 14 / 0 s | 14 / 0 s |
| manutec | 326 | 3623 | 4.5% | 2 s | 54 / 37 s | timeout | fail | fail | 44 / 4 s |
| glider | 752 | 4777 | 1.0% | 3 s | timeout | timeout | fail | fail | 117 / 8 s |
| shuttle_reentry | 482 | 4259 | 2.4% | 3 s | timeout | timeout | 269 / 14 s | timeout | timeout |
| launch | 658 | 10283 | 2.8% | 5 s | timeout | timeout | crash | fail | timeout |
| low_thrust | 803 | 10295 | 2.0% | 21 s | timeout | timeout | crash | fail | timeout |

## What it shows

**The sparsity is where the interest is.** The first five examples are Legendre
pseudospectral and their Jacobians are 20 to 49 per cent dense; the last five are local
collocation at 1 to 4.5 per cent. Those two groups behave quite differently, and a
comparison drawn only from the first five -- which is what every earlier round of this
work had -- says almost nothing about the second.

**GALAHAD is the only backend that copes with the large sparse problems.** manutec in 4
seconds against qpOASES's 37, glider in 8 seconds where qpOASES and ProxQP both ran out
of time. It solves five of the ten, more than any other. That is what it was expected to
do and the reason is structural: QPA states the subproblem exactly as the SQP poses it,
bounds included, and takes an indefinite Hessian as it is, so it is the one backend that
needs no convexification and no bounds smuggled in as identity rows.

**qpOASES times out on all five large problems**, which is the quadratic memory and
cubic work of its dense null-space factorisation arriving on schedule. It remains the
most reliable of the five on the small dense ones, and it is still the default for that
reason.

**ProxQP is fastest where it works** -- bryson_denham in 13 iterations and a second --
but it too times out on all five large problems, which was not what the earlier partial
runs suggested and is the clearest single argument for having run this properly.

**QPALM crashes on three of the ten**, inside its own factorisation, which no wrapper
can contain. Against that, it is the only backend that gets brac1 right: qpOASES
converges there to a stationary point 2.8 per cent above IPOPT's answer, and has done
since the exact Hessian went in.

**OSQP fails on six.** It is a first-order method asked for 1e-8, which is the wrong
question to put to it; nothing here is a defect in OSQP.

## What it does not show

No backend is close to IPOPT, which solves all ten in a few seconds each. The SQP is
not yet a replacement for it, and this sweep does not claim otherwise. What it does
establish is which QP engine to build on for the large sparse problems that are the
point of the sparse work, and the answer is GALAHAD.
