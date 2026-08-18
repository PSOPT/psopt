# A degenerate row is not a large gradient

This note records a defect in PSOPT's automatic constraint scaling, how it was found, what
it was doing to a third of the shipped examples, and what changed. It is written up
separately from the SQP work because it is not an SQP defect: it predates the SQP, it
affects every solver PSOPT can call, and the SQP's only role was to be the first solver
sensitive enough to complain about it.

## What was found

`determine_constraint_scaling_factors` in `src/scaling.cxx` scales each constraint so that
the corresponding row of the scaled Jacobian has unit Euclidean norm. That is the right
idea, and the loop that implemented it read:

```cpp
double sqeps = sqrt(PSOPT_extras::GetEPS());

if ( jac_row_norm(i) < 1.e7 ) {
        (*workspace->constraint_scaling)(i) = 1.0/(jac_row_norm(i)+sqeps);
}
else  {
        if ( jac_row_norm(i) > 1.e7 ) {
                (*workspace->constraint_scaling)(i) = 1.0/1.e7;
        }
}
```

The `+ sqeps` is there to stop a division by zero. It stops the division by zero, and in
doing so it hands a row whose Jacobian is *exactly* zero the largest factor the loop is
capable of producing, `1/sqrt(eps)` = 67108864. A row that carries no information about its
own magnitude is thereby declared to be the most delicately scaled row in the problem.

Two smaller things are wrong with the same loop. A row norm of exactly `1.e7` satisfies
neither branch, so its factor is left at whatever the workspace matrix happened to contain
-- the previous mesh's value, or uninitialised memory on the first mesh, since Eigen's
`resize` does not initialise. And there is no lower clamp to match the upper one, so a row
norm of, say, `2.e-8` still receives a factor of nearly `5.e7`.

## What it was doing

`bryson_max_range` is the clearest case, because there the arithmetic is exact. The problem
asks for `u1^2 + u2^2 = 1` at each of fifty nodes, and its initial guess is `u = 0` -- the
one point in the whole control space where that constraint's gradient, `(2*u1, 2*u2)`,
vanishes. All fifty path rows therefore had a zero Jacobian at the guess, all fifty were
multiplied by 67108864, and a constraint violated by exactly 1.0 was presented to the NLP
as a violation of 6.7108864e+07. That number, which had been read as evidence that the
problem simply starts a very long way from feasibility, is `1/sqrt(eps)` and nothing else.

It is not rare. Running every shipped example and counting rows whose Jacobian norm is
below `sqrt(eps)` at the initial guess:

| | examples |
|---|---|
| with at least one degenerate row | **34 of 66** |
| with ten or more | bryson_max_range (53), user (53), conic_soc (40), twoburn (26), heat (10) |

A third of the example set was being handed a scale factor of 6.7e+07 on at least one row.

The reason this went unnoticed for so long is that IPOPT applies its own gradient-based
scaling on top of PSOPT's (`nlp_scaling_method` is set to `gradient-based` in
`src/NLP_interface.cxx`), and an interior-point method with a filter line search and its own
restoration absorbs a badly scaled equality row without much complaint. The SQP has no
second scaling layer, and its feasibility phase reports the violation it is actually given.
It was right to complain.

## The change

A row whose Jacobian norm is at or below `sqrt(eps)` now receives a factor of 1.0, and every
other factor is clamped symmetrically to `[1.e-7, 1.e7]`, which also closes the gap at
exactly `1.e7`.

Choosing 1.0 is not an invention. The two scaling routines either side of this one already
do exactly that. `determine_objective_scaling` ends with
`if (nrm_g != 0.0 && ...) scale = 1/nrm_g; else scale = 1.0;`, and
`determine_scaling_factors_for_variables` initialises every factor to one and overwrites it
only where the bounds are informative. The constraint routine was the only one of the three
that responded to an absence of information with a large number instead of a neutral one.

## What the change costs and what it buys

Every example with a degenerate row was run under both solvers, before and after.

Under IPOPT, no example changed from solved to failed or from failed to solved. Six
objectives moved in their trailing figures -- bioreactor 1.3e-04, mpec 1.6e-03, lotka
3.7e-06, lotka2 2.5e-05, low_thrust 1.0e-05, conic_soc 1.0e-05 -- which is what a different
scaling of the same problem is expected to do. The other twenty-six were identical.

Under the SQP, two of the thirty-four changed from failing to solving, conic_soc and
isoperimetric, and none of the thirty-four went the other way. isoperimetric's whole
difficulty was this: one degenerate row, a starting violation reported as 6.5e+08, and a
feasibility phase that could not begin. Across the full sixty-six the SQP goes from
forty-five solved to forty-eight, which counts the guess repair below and two examples that
moved for reasons that are not this change at all -- see the section after next.

## The one that did not come free

Two examples did change their answer materially, and since they are the same problem there
is really one: `bryson_max_range` and its copy `user` went from -1.712313 to -0.251920.
Both points are feasible, both are KKT points, and the second is far worse -- the problem is
a maximisation, so -1.71 is the better objective.

The explanation is not that the new scaling is worse. It is that the guess `u = 0` sits at
the stationary point of `u1^2 + u2^2`, so the linearisation there supplies no direction at
all, and which local solution the run falls into is an accident of the arithmetic. The old
scaling was not choosing the better optimum on any principle; it was landing there.

The proper repair is to the guess rather than to the scaling. `u1 = 1, u2 = 0` satisfies the
path constraint exactly and costs nothing to write, and with it all four combinations of
{old, new} scaling and {IPOPT, SQP} return the same answer:

| | IPOPT | SQP |
|---|---|---|
| old scaling | -1.712314 | -1.712327 |
| new scaling | -1.712314 | -1.712327 |

That the SQP now solves this problem under the *old* scaling as well shows how much of its
earlier failure here was the guess and how much was the scaling: on this example, the guess
was the whole of it. isoperimetric and conic_soc, whose degenerate rows are intrinsic and
not an artefact of a poorly chosen starting point, needed the scaling change.

The value -1.7123 was checked outside PSOPT altogether, by a trapezoidal collocation of the
same problem in Python solved with SLSQP from ten random multistart guesses. It finds one
feasible local optimum, x_f = 1.712313, agreeing with PSOPT to seven figures. -0.2519 is
not a second optimum worth reporting; it is a stationary point the degenerate start led
into.

## Two examples that flipped, and why they do not count

Across the whole example set two problems changed status in the opposite directions:
`moon` went from failing at 47 iterations to solving in 17, and `twolinkarm` went from
solving in 270 iterations to failing at 210. Both were rebuilt and re-run against the old
scaling to confirm the attribution, and both reproduce.

Neither is evidence about this change. Instrumenting the loop to compare the old and new
factor for every row of every example gives:

| example | degenerate rows | largest relative change in any factor |
|---|---|---|
| bryson_max_range | 3 | 6.7e+07 |
| isoperimetric | 1 | 6.7e+07 |
| moon | 0 | **1.5e-08** |
| twolinkarm | 0 | **1.1e-08** |
| brac1 | 0 | 1.1e-08 |

Where there is no degenerate row, the only difference between the two rules is the dropped
`+ sqeps`, which moves each factor by about `sqeps` relative -- a part in a hundred million,
a couple of bits of a double. moon and twolinkarm both take hundreds of iterations and both
sit close to the point where the QP subproblem stops being solvable, and a perturbation in
the last two bits of the scaling is enough to send them either way. Changing the compiler's
optimisation level would do as much.

The useful conclusion is not that the change fixed moon and broke twolinkarm. It is that a
count of solved examples has a genuine uncertainty of a problem or two at the margin,
because a handful of the examples are that sensitive. brac1 is perturbed by the same 1.1e-08
and does not move at all; that is the normal case.

## What to take from it

The initial violation of 6.7e+07 on bryson_max_range and 6.5e+08 on isoperimetric was
recorded, in the first version of `SQP_ALL_EXAMPLES.md`, as a property of those problems:
they "begin at an enormously infeasible point". They do not. They begin at an ordinary
point, and the number was manufactured by the scaling on the way in. It is worth being
suspicious of a quantity that is large, unexplained, and the same across unrelated problems
-- 6.7108864e+07 is not a physical quantity, and dividing it by `1/sqrt(eps)` gives exactly
one.
