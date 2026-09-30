"""A tail measure needs a scenario set that resolves its tail, and now it can have one.

`robust_driver_risk.py` compares three risk measures on one problem and finds, at
alpha = 0.75, that sixteen quasi-random scenarios is about where CVaR starts doing
what it says. This example asks the same question at a level worth calling a tail --
alpha = 0.95, the worst one plant in twenty -- and the answer is not sixteen.

THE ARITHMETIC OF A TAIL

CVaR_alpha averages the worst 1 - alpha of the cost distribution, and the
sample-average approximation can only see the distribution the scenario set describes.
At alpha = 0.95 a set of sixteen points has eight tenths of a point in its tail.
Minimising CVaR over it is not minimising a tail measure; it is minimising something
whose value is decided by which single scenario happens to be worst, and the design
that comes out is worse in the real tail than the design that simply minimised the
mean. At forty-eight points the tail holds two and a half, and the measure starts to
earn its name.

WHY THIS EXAMPLE DID NOT EXIST UNTIL NOW

Because forty-eight scenarios were not reachable. Every transcription in PSOPT reuses
the collocation layout, whose defect block holds nstates*(nodes) rows while a scheme
with nodes-1 intervals can fill only nstates*(nodes-1) of them; the rest were written
as zeros with bounds [0,0], and IPOPT counts equality rows against variables. On an
augmented problem that is nstates PER SCENARIO of counted equalities holding no
dynamics. Measured on this problem: with the rows counted it carries **twenty-two**
scenarios and is refused at twenty-four with Not_Enough_Degrees_Of_Freedom, whatever
freedom the design still has.

So the density at which CVaR begins to beat the expectation design sits just past the
density at which the arithmetic refused it. `algorithm.free_padded_defect_rows` frees
those rows and the driver turns it on; the scenario count then costs nothing in
degrees of freedom at all, a scenario's pinned initial condition exactly cancelling
the state values at t0 that the defects leave undetermined.

WHICH TRANSCRIPTION, FOR A LARGE SET

Trapezoidal collocation, and the reason is cost rather than accuracy. Measured here at
forty-eight scenarios: 40 s against 94 s for multiple shooting with RK4 and eight
steps a segment, and 42 s for Hermite-Simpson, with all three agreeing on the
out-of-sample figures to four significant figures. This problem is smooth and has no
path constraint, so there is nothing for the integrator accuracy to buy; on a problem
with either, multiple shooting earns its extra cost. The verifier reads a trapezoidal
control as straight lines between the nodes, which is that scheme's own quadrature
assumption and the best-defined reading available there.
"""
import time
import numpy as np
import casadi as ca
from _common import psopt
from psopt.robust import RobustProblem, Gaussian

TF, NODES = 5.0, 41
X0 = np.array([1.0, 0.0])
MEAN = np.array([1.00, 1.00])
SD = np.array([0.20, 0.35])
RHO = -0.35
COV = np.array([[SD[0] ** 2, RHO * SD[0] * SD[1]],
                [RHO * SD[0] * SD[1], SD[1] ** 2]])
ALPHA = 0.95
# Sixteen is below the density the measure needs and forty-eight is above it. Both are
# solved for each risk measure, because the comparison that matters is not "CVaR against
# the mean" but "CVaR against the mean AT THE SAME SCENARIO COUNT" -- a denser set
# changes both designs, and attributing that to the risk measure would be the usual
# error.
COARSE, DENSE = 16, 48


def build():
    """The van der Pol oscillator of robust_driver_risk.py, same uncertainty."""
    rp = RobustProblem(name="robust_driver_cvar")
    ph = rp.add_phase(nstates=2, ncontrols=1, nevents=2)
    ph.nodes = [NODES]
    ph.dynamics = lambda x, u, p, t, th: ca.vertcat(
        x[1], -th[0] * x[0] + th[1] * (1.0 - x[0] ** 2) * x[1] + u[0])
    ph.integrand = lambda x, u, p, t: 0.5 * (x[0] ** 2 + x[1] ** 2 + u[0] ** 2)
    ph.events = lambda xi, xf, p, t0, tf, th: ca.vertcat(xi[0], xi[1])
    ph.bounds.lower.states = [-5.0, -5.0]
    ph.bounds.upper.states = [5.0, 5.0]
    ph.bounds.lower.controls = [-10.0]
    ph.bounds.upper.controls = [10.0]
    ph.bounds.lower.events = list(X0)
    ph.bounds.upper.events = list(X0)
    ph.bounds.t0 = (0.0, 0.0)
    ph.bounds.tf = (TF, TF)
    ph.guess.states = np.vstack([np.linspace(1.0, 0.0, NODES), np.zeros(NODES)])
    ph.guess.controls = np.zeros((1, NODES))
    ph.guess.time = np.linspace(0.0, TF, NODES).reshape(1, NODES)
    rp.uncertainty = Gaussian(mean=MEAN, cov=COV, truncate=2.5)
    rp.initial_state = X0
    # CVaR carries each scenario's cost as a state, and a state needs bounds. Asked for
    # rather than guessed: a bound that turned out to be active would silently change
    # the risk measure into something else.
    rp.cost_bounds = (0.0, 200.0)
    return rp


alg = psopt.Algorithm(collocation_method="trapezoidal",
                      scaling="automatic", nlp_tolerance=1.0e-7,
                      nlp_iter_max=3000, print_level=0)

# One test set, drawn once, used to score every design. No plant in it takes any part in
# any design, and it is shared so that comparing two designs is not also comparing two
# samples.
TEST = Gaussian(mean=MEAN, cov=COV, truncate=2.5).sample(
    3000, np.random.default_rng(20260928))
TEST = TEST[[Gaussian(mean=MEAN, cov=COV, truncate=2.5).contains(t) for t in TEST]]


def tail_mean(J, alpha):
    """CVaR, out of sample: the mean of the worst 1 - alpha of a sample."""
    k = max(1, int(round((1.0 - alpha) * len(J))))
    return float(np.sort(J)[-k:].mean())


def design(risk, M):
    rp = build()
    t0 = time.time()
    out = rp.solve(alg, slack=0.0, risk=risk, cvar_alpha=ALPHA,
                   scenarios="qmc", n_scenarios=M, generate=False,
                   out_of_sample=0, wait_and_see=0, verbose=False)
    # The realised cost of THIS design on plants it never saw, by the driver's own cost
    # sweep, which shares the equations with the transcription and nothing else.
    J = rp._costs_many(TEST, (out.time, out.controls),
                       np.asarray(out.design.parameters).ravel(), nsub=8)
    return out, J, time.time() - t0


print("van der Pol, CVaR at alpha = %.2f: the tail is the worst one plant in %d"
      % (ALPHA, round(1.0 / (1.0 - ALPHA))))
print("  scored on %d plants inside the set, none of which any design saw" % len(TEST))
print("  %d nodes, trapezoidal collocation" % NODES)

print("\n  %-13s %3s %10s %10s %10s %7s"
      % ("risk", "M", "reported", "E[J] out", "CVaR out", "secs"))
res = {}
for M in (COARSE, DENSE):
    for risk in ("expectation", "cvar"):
        out, J, el = design(risk, M)
        res[(risk, M)] = (out, J)
        print("  %-13s %3d %10.5f %10.5f %10.5f %7.0f"
              % (risk, M, out.objective, J.mean(), tail_mean(J, ALPHA), el))

# ---- what the density bought ------------------------------------------------------
def tails(M):
    return (tail_mean(res[("expectation", M)][1], ALPHA),
            tail_mean(res[("cvar", M)][1], ALPHA))


e16, c16 = tails(COARSE)
e48, c48 = tails(DENSE)
print("\n  at %d scenarios the CVaR design is %+.1f%% in the tail against the design"
      % (COARSE, 100.0 * (c16 - e16) / e16))
print("  that minimised the mean -- WORSE, because at alpha = %.2f a set of %d points"
      % (ALPHA, COARSE))
print("  has %.1f of a point in its tail, and minimising over that is minimising"
      % ((1.0 - ALPHA) * COARSE))
print("  which single scenario happens to be worst.")
print("\n  at %d it is %+.1f%%, and the trade is now a real one: it pays %+.1f%% on the"
      % (DENSE, 100.0 * (c48 - e48) / e48,
         100.0 * (res[("cvar", DENSE)][1].mean()
                  - res[("expectation", DENSE)][1].mean())
         / res[("expectation", DENSE)][1].mean()))
print("  mean to take %.1f%% off the tail. Whether that is worth it is the modelling"
      % (-100.0 * (c48 - e48) / e48))
print("  decision the measure exists to express; that it is AVAILABLE is the point.")

# The sample-average approximation, against the thing it approximates.
rep48 = res[("cvar", DENSE)][0].objective
print("\n  and the in-sample estimate is optimistic, as a sample average over the very")
print("  scenarios that were optimised should be: %.5f reported against %.5f realised,"
      % (rep48, c48))
print("  %+.1f%%. The gap is what a larger set buys next, and it is why the reported"
      % (100.0 * (rep48 - c48) / c48))
print("  value of a risk measure is never the certificate.")

print("\n  Under the old counting this problem carried twenty-two scenarios and was")
print("  refused at twenty-four. The density at which the measure starts to work sat")
print("  just past the density the arithmetic refused, which is why this example could")
print("  not have been written before algorithm.free_padded_defect_rows.")

# A test of the finding, not of a value: the measure must FAIL at the coarse set and
# WORK at the dense one, and the trade at the dense set must be a real trade.
ok = (all(res[k][0].converged for k in res)
      and c16 > e16                     # coarse: CVaR is worse in the tail
      and c48 < 0.95 * e48              # dense: better by at least five per cent
      and res[("cvar", DENSE)][1].mean() > res[("expectation", DENSE)][1].mean()
      and rep48 < c48)                  # the in-sample estimate is optimistic
print("\nobjective      : %.12g" % rep48)
print("RESULT: %s" % ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
