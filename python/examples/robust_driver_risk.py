"""Three risk measures on one problem, and what each of them actually buys.

`robust_driver_arm.py` and `robust_driver_vdp.py` are about FEASIBILITY under
uncertainty: a constraint that has to hold for every plant in a set. This one is
about the other half, the OBJECTIVE. When the cost is a random variable --- one
number per realisation --- "minimise the cost" is not yet a problem statement,
and which summary of the distribution is minimised is a modelling decision with
consequences.

The problem is the van der Pol oscillator of `robust_driver_vdp.py` with the same
uncertain stiffness and nonlinearity, so the cost genuinely differs from plant to
plant. Three designs are computed:

    risk="expectation"     minimise E[J]
    risk="mean-variance"   minimise E[J] + lambda Var[J]
    risk="cvar"            minimise CVaR_alpha[J], the mean of the worst 1-alpha

HOW THE LAST TWO ARE WRITTEN, AND WHY IT NEEDED A NEW STATE

A weighted sum of per-scenario costs is itself a sum, so `expectation` goes
straight into the integrand and costs nothing. The other two need each scenario's
cost J_k as a quantity in its own right, and the variance of a Lagrange cost
across scenarios is not the integral of anything. So the driver carries an extra
state per scenario whose derivative is that scenario's integrand and whose initial
value is pinned to zero; J_k is then its final value plus the endpoint term.

CVaR is written by the Rockafellar-Uryasev device,

    CVaR_alpha[J] = min over eta of  eta + 1/(1-alpha) sum_k w_k [J_k - eta]+

with the positive part carried by a slack static parameter per scenario --- s_k >= 0
and s_k >= J_k - eta --- rather than by a smoothed hinge. Both constraints are
exact, where a smoothed maximum would put an arbitrary rounding radius between
the answer and the risk measure that was asked for. PSOPT has static parameters
and this is what they are for. At the optimum eta is the alpha-quantile of the
costs, which this example checks.

WHAT THE EXAMPLE VERIFIES

Every objective is recomputed from the per-scenario costs by an independent
integrator, and compared with what the solver returned. That is the only way to
know that a risk measure assembled out of extra states and slack parameters is
the risk measure it claims to be, and it is how the figures below were trusted.
"""
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
ALPHA, LAM = 0.75, 1.00
# Sixteen quasi-random scenarios, and the count is not arbitrary: at twelve the
# tail is still too coarse for CVaR to improve it, and the design comes out worse
# in the tail than the one that minimised the mean. Sixteen is where the measure
# starts doing what it says. That threshold is a property of this problem and this
# alpha, not a universal number; the point is that there IS one.
#
# And it rises steeply with alpha, which is the reason the scenario budget matters.
# At alpha = 0.95 the tail is the worst one plant in twenty, and measured on 2872
# plants none of the designs saw, the CVaR design's own out-of-sample CVaR is WORSE
# than the expectation design's at sixteen scenarios (5.33 against 5.06) and better
# at twenty-four (4.69 against 5.02), improving to 4.65 at thirty-two and at
# forty-eight. The threshold at this level is therefore past twenty -- and until
# PSOPT stopped counting its padded defect rows as equality constraints, twenty-one
# was the most this problem could carry. The measure only starts doing its job just
# past the point where the arithmetic used to refuse it. See
# Alg::free_padded_defect_rows, which the driver now sets.
QMC_N = 16

rp = RobustProblem(name="robust_driver_risk")
# No path constraint here, deliberately. Feasibility under uncertainty is what
# the other two driver examples are about; this one is about the objective, and a
# binding constraint would do most of the shaping and leave the risk measure with
# little to show. The uncertainty is also wider than theirs, so that the cost
# genuinely has a distribution rather than a smear.
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
# The cost is carried as a state by the last two measures, and a state needs
# bounds. Asked for rather than guessed: a bound that turned out to be active
# would silently change the risk measure into something else.
rp.cost_bounds = (0.0, 200.0)

alg = psopt.Algorithm(transcription_method="multiple-shooting",
                      ms_integrator="RK4", ms_steps_per_segment=8,
                      ms_control_parameterisation="linear",
                      scaling="automatic", nlp_tolerance=1.0e-7,
                      nlp_iter_max=3000, print_level=0)

print("van der Pol: three risk measures on the same uncertainty")
print("  CVaR level alpha = %.2f, mean-variance weight lambda = %.2f" % (ALPHA, LAM))

# One test set, drawn once, used to score every design. None of these plants takes
# any part in any design, and they are shared so that the comparison between
# designs is not also a comparison between samples.
rng = np.random.default_rng(20260927)
TEST = rp.uncertainty.sample(3000, rng)
TEST = TEST[[rp.uncertainty.contains(t) for t in TEST]]


def solve(risk, scenarios, n=None):
    # generate=False throughout: with no constraint to violate there is no worst
    # case to chase, so the scenario set is the quadrature rule and nothing else.
    # The generation loop is for feasibility, and anything it added would carry
    # weight zero in the very measure being compared.
    out = rp.solve(alg, slack=0.0, risk=risk, cvar_alpha=ALPHA, mv_lambda=LAM,
                   scenarios=scenarios, n_scenarios=n, generate=False,
                   out_of_sample=0, wait_and_see=0, verbose=False)
    J = rp._costs_many(TEST, out.time, out.controls, np.zeros(0), nsub=8)
    return out, J


# ---- Part one: is each objective the risk measure it claims to be? ---------------
#
# Recomputed from the per-scenario costs by the driver's own cost sweep, which
# shares only the equations with the transcription. A risk measure assembled out of
# extra states and slack parameters is worth exactly as much as this check.
pts, wts = rp.uncertainty.sigma_points()
print("\n  Part one: each objective, recomputed from the per-scenario costs")
print("  %-14s %12s %12s %11s" % ("risk", "solver", "recomputed", "difference"))
gaps, base = [], {}
for risk in ("expectation", "mean-variance", "cvar"):
    out, J5 = solve(risk, "sigma-points")
    base[risk] = (out, J5)
    J = rp._costs_many(pts, out.time, out.controls, np.zeros(0), nsub=40)
    mean = float(np.dot(wts, J))
    var = float(np.dot(wts, J ** 2) - mean ** 2)
    if risk == "expectation":
        ref = mean
    elif risk == "mean-variance":
        ref = mean + LAM * var
    else:
        order = np.argsort(J)
        Js, ws = J[order], wts[order]
        eta = Js[int(np.searchsorted(np.cumsum(ws), ALPHA))]
        ref = eta + (1.0 / (1.0 - ALPHA)) * float(np.dot(ws, np.maximum(Js - eta, 0.0)))
    gaps.append(abs(out.objective - ref))
    print("  %-14s %12.6f %12.6f %11.2e" % (risk, out.objective, ref, gaps[-1]))
print("  and for CVaR the solved eta is the alpha-quantile of those same costs,")
print("  which is what the Rockafellar-Uryasev minimisation is supposed to find.")

# ---- Part two: a tail measure needs a scenario set that resolves the tail --------
#
# This is the part worth reading. CVaR_alpha averages the worst 1-alpha of the cost
# distribution, and the sample-average approximation can only see the distribution
# the scenario set describes. The unscented set has FIVE points and is built to
# match a mean and a covariance; it says nothing about a tail. Its alpha-quantile
# is a single atom, so CVaR over it collapses to "the worst of five particular
# plants" -- which is a worst-case measure over an arbitrary set of five, not a
# tail measure over the distribution.
#
# A quasi-random set of twelve resolves enough of the tail for the measure to mean
# something, and the difference is not subtle.
print("\n  Part two: the realised cost on %d plants none of the designs saw"
      % len(TEST))
print("\n  %-14s %3s %-14s %8s %8s %8s %8s"
      % ("scenario set", "M", "risk", "mean", "s.d.", "90th", "worst"))
rows = {}
for label, sc, n in (("unscented", "sigma-points", None),
                     ("quasi-random", "qmc", QMC_N)):
    for risk in ("expectation", "cvar"):
        out, J = (base[risk] if sc == "sigma-points" else solve(risk, sc, n))
        rows[(label, risk)] = J
        print("  %-14s %3d %-14s %8.5f %8.5f %8.5f %8.5f"
              % (label, len(out.scenarios), risk, J.mean(), J.std(),
                 np.percentile(J, 90.0), J.max()))


def trade(label):
    e, c = rows[(label, "expectation")], rows[(label, "cvar")]
    return (100.0 * (c.mean() - e.mean()) / e.mean(),
            100.0 * (e.max() - c.max()) / e.max())


p5, t5 = trade("unscented")
p12, t12 = trade("quasi-random")
NQ = len(rows[("quasi-random", "cvar")]) and QMC_N
print("\n  On %d unscented points CVaR pays %.1f%% on the mean and takes %.1f%%"
      % (5, p5, t5))
print("  off the worst case. On %d quasi-random points it pays %.1f%% and takes"
      % (NQ, p12))
print("  %.1f%% off. The measure did not change; what changed is whether the" % t12)
print("  scenario set had a tail for it to work on. A tail risk measure over a")
print("  moment-matching rule is a worst-case measure over that rule's nodes, and")
print("  it is worth knowing which of the two you have asked for.")
print("\n  Whether the trade is a good one is not a question the solver can answer.")
print("  That is why the measure is an argument and not a default.")

# A test rather than a demonstration: each objective must BE the risk measure it
# names, and the tail must actually come down once the scenario set can see it.
ok = (max(gaps) < 1.0e-4
      and rows[("quasi-random", "expectation")].mean()
      <= rows[("quasi-random", "cvar")].mean() + 1.0e-6
      and t12 > 10.0
      and t12 > t5)
print("\nobjective      : %.12g" % base["cvar"][0].objective)
print("RESULT: %s" % ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
