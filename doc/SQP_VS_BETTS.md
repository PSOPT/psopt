# PSOPT's SQP against Betts's, component by component

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
| 3 | Warm start of the QP from the previous active set | 2.3 | **absent** |
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
| 20 | Inertia control rung (iii): H = I fails, go find a feasible point | 2.6.1, step 2(c)iii | **absent** |
| 21 | Line search on a quadratic and cubic model of the merit | 2.6.1 | absent; plain halving |
| 22 | Wolfe condition, to stop steplengths collapsing | 2.6.1 | **absent** |
| 23 | tau and the weights start at zero, merit starts as the Lagrangian | 2.6.1 | present |
| 24 | Strategy M, minimize from x0 | 2.6.2 | present |
| 25 | Strategy FM, find a feasible point first — *his default* | 2.6.2 | **absent** |
| 26 | Strategies FME and F | 2.6.2 | absent |
| 27 | Segregating constraint difficulties from objective ones | 2.7 | **absent** |
| 28 | LDP feasibility step, min half p'p | 2.8.1, (2.46)–(2.47) | **absent** |
| 29 | Relaxation step, min half p'p + rho/2 u'u | 2.8.1, (2.49)–(2.50) | partial |
| 30 | Constraint-violation merit Mv | 2.8.1, (2.51) | **absent** |
| 31 | LDP / relaxation switching on alpha = 1 | 2.8.2 | absent |
| 32 | Abandoning the primary step on fill, conditioning, inertia, QP count | 2.8.2 | inertia only |

Thirteen present, four partial, fifteen absent. That count flatters the gap and
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

**First, the whole of section 2.8 and the FM strategy — finding a feasible point before
optimising.** This is the single largest gap, and it is not a refinement: FM is the
default strategy in Betts's own software, and section 2.7 explains why in terms that
describe our benchmark exactly. His first premise is to *segregate difficulties caused by
the constraints from difficulties caused by the objective*, and his mechanism is to find
a feasible point first while ignoring the objective, the multipliers and the Hessian
entirely. What we do instead is carry the objective into every relaxed subproblem and
judge the resulting step with a merit function that is still trying to optimise.

The two examples the solver cannot touch are the two where this bites. On launch, 282 of
560 constraints are violated at the starting guess, the linearisation is inconsistent at
every iterate, and restoration therefore runs at every iteration; timing the parts put 99
per cent of the wall clock inside the relaxed subproblem. That is precisely the situation
section 2.7 says the FM strategy exists to avoid. A feasibility phase would replace those
subproblems with the much smaller LDP of (2.46) — no objective, no Hessian, no
multipliers, minimum-norm step — and judge it against the constraint-violation merit
(2.51) rather than an augmented Lagrangian that has no business being consulted yet.

I attempted a version of this in one increment and withdrew it. It failed for a reason
worth recording: I zeroed the objective in the relaxed subproblem but left the line search
judging the step with the augmented Lagrangian, so the search cut a perfectly good
feasibility step to a thirty-second because the objective had gone the wrong way. The
phase needs its own merit function, which is (2.51), and its own step, which is (2.46).
Half of it is worse than none.

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
