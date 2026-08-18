# PSOPT's SQP against Betts's, component by component

> **Note (superseded in part).** qpOASES appears throughout the tables below as one of
> the backends measured. It was removed from PSOPT after this work: it is a dense
> active-set method, and the timings here are the evidence for that decision rather than
> a record of a supported configuration. Its rows are left as measured; nothing else in
> this document depends on them.

The SQP solver in `src/SQP_interface.cxx` was not built from Betts's chapter 2 and then
completed; it was built from a dense SQP and then had pieces of chapter 2 grafted onto
it, one increment at a time, each kept or discarded on what it did to the ten-example
benchmark. That is a reasonable way to work but it leaves an obvious question
unanswered, which this file answers: how much of the algorithm in the book is actually
there.

The reference is *Practical Methods for Optimal Control Using Nonlinear Programming*,
3rd edition, chapter 2, sections 2.3 to 2.8. Equation numbers below are the book's.

## The audit

| # | Component | Book | Status |
|---|---|---|---|
| 1 | Sparse QP subproblem, standard form | 2.3, (2.17)–(2.19) | delegated to the backend |
| 2 | Schur-complement active set over one LDL' of K0 | 2.3, (2.21)–(2.26) | absent; the backends do as they please |
| 3 | Warm start of the QP from the previous active set | 2.3 | absent; measured as not worth building — see below |
| 4 | Augmented Lagrangian merit function | 2.4, (2.27) | present for the constraints |
| 5 | Constraint slacks s | 2.4, (2.28) | present, `merit_slacks` |
| 6 | Bound slacks t, multipliers nu, weights Xi | 2.4, (2.27), (2.29) | **absent** |
| 7 | Step moves (x, lambda, s) together | 2.4, (2.30)–(2.34) | present |
| 8 | Bound-multiplier and bound-slack steps | 2.4, (2.32), (2.35) | **absent** |
| 9 | Least-norm penalty weights Theta | 2.4, (2.36)–(2.37) | present |
| 10 | Threshold psi0 raised when the min-norm solution is zero | 2.4 | absent; fixed at machine epsilon |
| 11 | Hessian modification H = HL + tau(|sigma|+1)I | 2.5, (2.39) | diverged, deliberately |
| 12 | Gerschgorin bound sigma | 2.5, (2.40) | present, sets the ceiling |
| 13 | tau from the rho1/rho2/rho3 trust-region rule | 2.5, (2.42)–(2.44) | removed, deliberately |
| 14 | Inertia test In(K0) = (nf, m, 0) | 2.5, (2.41) | present, rank-tolerant form |
| 15 | KKT termination test | 2.6.1, step 1(b) | present |
| 16 | Exact Hessian of the Lagrangian | 2.5, (2.38) | present, sparse, by AD |
| 17 | First iteration solves two QPs, the first with H = I | 2.6.1 | present, `multiplier_pass` |
| 18 | Inertia control rung (i): raise tau and retry | 2.6.1, step 2(c)i | present |
| 19 | Inertia control rung (ii): tau = 1 fails, set H = I and retry | 2.6.1, step 2(c)ii | **absent** |
| 20 | Inertia control rung (iii): H = I fails, go find a feasible point | 2.6.1, step 2(c)iii | absent; the phase exists, the hook does not |
| 21 | Line search on a quadratic and cubic model of the merit | 2.6.1 | absent; plain halving |
| 22 | Wolfe condition, to stop steplengths collapsing | 2.6.1 | **absent** |
| 23 | tau and the weights start at zero, merit starts as the Lagrangian | 2.6.1 | present |
| 24 | Strategy M, minimize from x0 | 2.6.2 | present |
| 25 | Strategy FM, find a feasible point first — *his default* | 2.6.2 | present |
| 26 | Strategy F | 2.6.2 | present |
| 26b | Strategy FME | 2.6.2 | absent |
| 27 | Segregating constraint difficulties from objective ones | 2.7 | present, via FM |
| 28 | LDP feasibility step, min half p'p | 2.8.1, (2.46)–(2.47) | present |
| 29 | Relaxation step, min half p'p + rho/2 u'u | 2.8.1, (2.49)–(2.50) | present |
| 30 | Constraint-violation merit Mv | 2.8.1, (2.51) | present |
| 31 | LDP / relaxation switching on alpha = 1 | 2.8.2 | present |
| 32 | Abandoning the primary step on fill, conditioning, inertia, QP count | 2.8.2 | inertia only |

Twenty present, two partial, eleven absent, after the feasibility phase went in; it was thirteen, four and fifteen when this file was first written. That count flatters the gap and
understates it at the same time, so it is worth saying where the weight sits.

## What is actually complete

**Betts's minimization process — his strategy M — is essentially implemented.** Steps 1
through 5 of section 2.6.1 are all there: the KKT test, the exact sparse Hessian of the
Lagrangian, the two-QP start that gets first-order multipliers from a unit Hessian, the
inertia-controlled modification, the augmented Lagrangian merit function with its slacks
and its least-norm weights, and the step that moves the multipliers and slacks alongside
the variables. Four of the entries marked "diverged" or "removed" are places where the
book was followed, measured, and departed from on evidence, each recorded in the source
and in the commit that made the change: the Levenberg parameter is now an absolute shift
escalated from the last one that worked rather than a fraction of the Gerschgorin bound
restarted at every iteration, because the bound moves with the multipliers by six orders
of magnitude within three iterations; and his rho1/rho2/rho3 rule for that parameter is
not reinstated because the inertia of the KKT matrix answers the same question without a
step having to be taken and judged first.

Three things in the implementation are *not* in Betts and are there because they were
needed: a trust region on the step (Betts gives the Levenberg parameter that role), a
second-order correction for the Maratos effect, and a damped BFGS model for when the
exact Hessian is not wanted.

## What is missing, in order of what it is costing

*Written before the feasibility phase was built; section 2.8 and the FM strategy headed
this list and have since gone in. What follows is the rest of it, with the original
entry on the feasibility phase kept below as an addendum because the reasoning in it is
what the implementation followed.*

**Second, the last two rungs of the inertia-control ladder.** When the shift reaches its
ceiling we currently give up on the model and go to restoration. Betts sets H = I and
tries once more, and only when *that* fails concludes that the constraints are the
problem. The middle rung is nearly free to add and is the correct reading of what a
failed shift means: it separates "the Hessian is unusable here" from "the constraints are
inconsistent here", which is the same segregation principle as FM at the scale of a
single iteration. Rung (iii) is a jump into the feasibility phase, so it arrives with the
first item.

**Third, the line search.** Halving is the crudest thing that works. Betts fits a
quadratic and a cubic to the merit function, and imposes the Wolfe condition specifically
to stop steplengths collapsing. Our traces are full of collapsing steplengths — 6.25e-02,
1.95e-03, 1.19e-07 on the harder examples — and each one is a full evaluation of the
objective and all the constraints thrown away. This is cheap to fix and would show up
everywhere at once.

**Fourth, the QP warm start.** Betts is explicit that the efficiency of his QP comes from
factorising the KKT matrix once and updating a small dense Schur complement as the active
set changes, and that as the NLP converges the active set stops changing and the QP
iteration count collapses. We call `init()` on a fresh problem for every subproblem, at
every iteration, so every subproblem pays a cold start. qpOASES and GALAHAD both support
hot starts; nothing in the plugin ABI currently carries an active set across a call, so
this is an interface change as well as a code change.

**Fifth, the bound terms in the merit function.** Betts's merit function carries slacks,
multipliers and penalty weights for the simple bounds as well as the general constraints;
ours carries only the general constraints. The justification recorded in the source is
that the iterates are clipped into the box, so (x − t) is zero — that is true at a KKT
point and when a bound is inactive, but not in general along a step where a bound
multiplier is nonzero and the variable is interior. This is a real simplification rather
than a no-op, and it is the least of the five.

## How far, in effort

The minimization process is done. What remains is one substantial piece of work and three
small ones.

The feasibility phase is the substantial one: an LDP subproblem, the relaxation fallback,
the constraint-violation merit function, the switching logic between the two, and the FM
strategy that decides when to run the phase at all. It touches the driver's structure
rather than sitting inside it, and it needs its own tests. Call it a week of careful work
and expect it to be the increment that either makes launch and low_thrust solvable or
proves they need something else.

The inertia ladder's middle rung is an afternoon. The line search model is a day. The
merit function's bound terms are a day, and would only be worth doing after the others,
since their effect is likely to be small and hard to separate from noise at the factor-of-
two sensitivity these iteration counts show.

The QP warm start is a day of code and a change to the plugin ABI, but its value depends
on settling the backend question first: it is worth doing once, for the backend that is
going to be *the* backend, and not five times.

## The honest summary

Of Betts's algorithm we have the optimisation half and not the feasibility half. That is
not a coincidence of what was easy — it is the order the increments were attempted in,
and each one was chosen by what the benchmark rewarded on the six examples that already
solved. The four that do not solve are all failures of the half that is missing.


## Addendum: the feasibility phase, built

Section 2.8 and the FM strategy went in the increment after this file was written, and
the prediction it makes above turned out to be right: launch is solved, by GALAHAD, in 24
feasibility iterations and 23 optimality ones, at IPOPT's answer. The details are in
SQP_BACKEND_BENCHMARK.md. What remains of the original list is the inertia ladder's
middle rung and its hook into the phase, the line search model and the Wolfe condition,
the QP warm start, the merit function's bound terms, and the FME strategy -- which now
looks more interesting than it did, because the first optimality step on bryson_denham
takes the violation from 2.3e-11 back to 1.3 and FME exists to stop exactly that.


## Addendum: two more components built, measured, and not kept

Both were implemented in full from the book, both are faithful, and both regress the
benchmark. They are recorded here rather than in the commit history alone, because "we
tried it and it made things worse" is the sort of finding that otherwise gets
rediscovered every eighteen months.

**Rung (ii) of the inertia ladder, section 2.6.1 step 2(c)ii.** When the shift reaches
its ceiling and the inertia is still wrong, replace the Hessian by the identity rather
than sending an unusable model to the subproblem. Tried four times, in four
configurations, and it fails in all of them. The mechanism is that replacing H by the
identity makes the step *longer*: the ceiling-valued shift it replaces had been acting as
a strong regulariser, since H + delta I with delta = |sigma| + 1 gives a step of
-g/delta, and at |sigma| of 8.8e+02 that is very short indeed. On brac1 under strategy M
the first step then drives the scaled objective of a minimum-time problem to 2e-14, and
the merit function accepts it, because at the first iteration the penalty weights are at
machine epsilon by construction and the merit is therefore the objective alone.

The interpolating line search does not save it -- the step is accepted, not truncated, so
there is nothing for the line search to refuse. Nor do the merit function's bound terms,
which was the next hypothesis, on the grounds that brac1's collapse is the final time
going to its bound. Under FM, where the optimality phase starts from a feasible point,
brac1 and bryson_denham and lts are all fine and launch fails instead.

**The merit function's bound terms, (2.27), (2.29), (2.32) and (2.35).** Bound slacks t,
bound multipliers nu, penalty weights Xi, and the least-norm weight calculation of (2.37)
extended from m to m + n components as Betts writes it. This is a straightforward piece
of work and the result is faithful; the one thing worth flagging for anyone repeating it
is the sign. Betts writes his Lagrangian as g - G'lambda - nu and PSOPT writes
grad f + J'lambda - z, so his lambda is minus PSOPT's -- but his nu is PSOPT's z with the
*same* sign, both being subtracted. The constraint term flips and the bound term does
not. Getting that wrong costs brac1 and bryson_denham their solutions and looks exactly
like the terms being harmful.

With the signs right, the bound terms are neutral on all five small examples -- brac1 17
iterations, bryson_denham 18, hypersensitive 11, lts 20, interior_point 2, cell for cell
-- and they break launch, which goes from 42 iterations and a solution to 17 and none.
That is the whole measured effect: no benefit anywhere and one clear loss. The earlier
assessment in this file, that the bound terms are "a real simplification rather than a
no-op", turns out to be right about the theory and wrong about which direction the
practice runs.

Both are kept as files rather than commits, and the audit table above still marks them
absent, which is the honest state: the code exists, it does what the book says, and it is
not in the library.


## Addendum: the warm start, measured and declined

The QP warm start was the last substantial item on the list above, and the argument for
it is Betts's own: the efficiency of his QP comes from factorising the KKT matrix once
and updating a small dense Schur complement as the active set changes, and as the NLP
converges the active set stops changing and the QP iteration count collapses (section
2.3). We call the equivalent of `init()` on a fresh problem for every subproblem, so
every one pays a cold start, and carrying an active set across the plugin boundary would
be an ABI change as well as a code change.

Before building it, the thing to check is whether the count is collapsing already. It is
printed in every trace, so this cost nothing.

| example | backend | SQP iterations | QP iterations, mean | first five | last five |
|---|---|---|---|---|---|
| glider | GALAHAD | 364 | 14 | 4, 4, 4, 4, 4 | 20, 20, 20, 20, 20 |
| shuttle_reentry | GALAHAD | 70 | 6 | 6, 6, 6, 6, 6 | 4, 4, 4, 4, 4 |
| launch | GALAHAD | 42 | 89 | 123, 4, 7, 13, 10 | 6, 4, 4, 6, 6 |
| brac1 | qpOASES | 17 | 517 | 537, 530, 771, 392, 441 | 633, 577, 625, 581, 571 |
| bryson_denham | qpOASES | 18 | 897 | 527, 1052, 851, 862, 1607 | 601, 548, 594, 585, 754 |
| lts | qpOASES | 20 | 202 | 128, 226, 326, 336, 214 | 152, 148, 149, 337, 374 |

The two backends could not be further apart. **GALAHAD is already doing what warm starting
is supposed to achieve**: four to twenty iterations a subproblem, and on launch the count
falls from 123 on the first to four by the last. There is no warm start that improves on
four iterations. **qpOASES is the opposite** -- two hundred to nine hundred iterations,
and flat from the first subproblem to the last, which is precisely the signature of
paying a cold start every time.

So the warm start is worth a great deal for the backend being retired and essentially
nothing for the one being kept. It is not built, and this table is the reason. If the
backend decision is ever revisited, the measurement to repeat is the last two columns:
flat means there is something to win, falling means there is not.

One implementation note for whoever does revisit it. qpOASES's `hotstart()` is not the
mechanism it looks like -- it accepts a new gradient and new bounds but not a new Hessian
or Jacobian, and in an SQP both change at every iteration. The applicable call is `init()`
with a guessed working set. GALAHAD's QPA takes a warm start through its entry status and
the x_stat and c_stat arrays, which would need either an opaque per-call-site context in
the plugin ABI or state held inside the plugin and keyed on the problem dimensions.
