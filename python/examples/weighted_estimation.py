"""Parameter estimation with weighted residuals, and what the estimates are worth.

A catalytic cracking reaction, fitted to 21 noisy samples of two states:

    y1' = -(theta1 + theta3) y1^2,   y2' = theta1 y1^2 - theta2 y2,

with theta1, theta2, theta3 estimated from the data. This is the Python port of
examples/cracking in the C++ distribution; the unweighted least-squares residual
is 4.319519e-03, which the first solve below reproduces.

The point of the example is everything after that. An estimate without an
uncertainty is not an estimate, and PSOPT computes the covariance of the parameter
vector, its confidence intervals and the residual standard deviation. It also lets
each residual be weighted, which is what you do when the observations are not
equally trustworthy: the usual choice is the reciprocal of each measurement's
standard deviation, so that a noisy channel does not pull the fit around.

Here the second state is measured about five times less precisely than the first,
so the second solve weights it down and the two fits are compared. The unweighted
fit is not wrong -- it is the answer to a different question, in which both
channels are asserted to be equally reliable.
"""
import numpy as np
import casadi as ca
from _common import psopt, report

y1 = [1.0, 0.8105, 0.6208, 0.5258, 0.4345, 0.3903, 0.3342, 0.3034, 0.2735, 0.2405,
      0.2283, 0.2071, 0.1669, 0.153, 0.1339, 0.1265, 0.12, 0.099, 0.087, 0.077, 0.069]
y2 = [0.0, 0.2, 0.2886, 0.301, 0.3215, 0.3123, 0.2716, 0.2551, 0.2258, 0.1959,
      0.1789, 0.1457, 0.1198, 0.0909, 0.0719, 0.0561, 0.046, 0.028, 0.019, 0.014, 0.01]
tm = [0.0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25,
      0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.65, 0.75, 0.85, 0.95]
NS = len(tm)


def build(weights=None, regularization=None, observations=None):
    prob = psopt.Problem(name="cracking_weighted")
    ph = prob.add_phase(nstates=2, ncontrols=0, nparameters=3, nobserved=2, nsamples=NS)
    ph.nodes = [80]
    ph.dynamics = lambda x, u, p, t: ca.vertcat(-(p[0] + p[2]) * x[0] ** 2,
                                                p[0] * x[0] ** 2 - p[1] * x[1])
    ph.observation = lambda x, u, p, t: ca.vertcat(x[0], x[1])
    ph.bounds.lower.states = [0.0, 0.0]
    ph.bounds.upper.states = [2.0, 2.0]
    ph.bounds.lower.parameters = [0.0, 0.0, 0.0]
    ph.bounds.upper.parameters = [20.0, 20.0, 20.0]
    ph.bounds.t0 = (0.0, 0.0)
    ph.bounds.tf = (0.95, 0.95)
    ph.observation_nodes = np.array(tm).reshape(1, NS)
    ph.observations = (np.vstack([y1, y2]) if observations is None else observations)
    if weights is not None:
        ph.residual_weights = weights
    if regularization is not None:
        ph.regularization_factor = regularization
    sg = np.vstack([np.linspace(1.0, 0.069, 40), np.linspace(0.30, 0.01, 40)])
    ph.guess.states = sg
    ph.guess.time = np.linspace(0.0, 0.95, 40).reshape(1, 40)
    ph.guess.parameters = np.zeros((3, 1))
    return prob


def show(sol, label):
    print("\n---- %s" % label)
    print("objective      : %.9e" % sol.objective)
    print("status         : %s" % sol.status)
    p = sol.parameters
    ps = sol.parameter_statistics
    if ps is None:
        print("no parameter statistics were computed")
        return p, None
    se = ps.standard_errors
    print("sigma_hat      : %.6e   (estimated residual standard deviation)" % ps.sigma_hat)
    for k in range(len(p)):
        print("theta%d         : %8.5f +/- %7.5f    95%% interval [%8.5f, %8.5f]"
              % (k + 1, p[k], se[k], ps.confidence_low[k], ps.confidence_high[k]))
    return p, ps


alg = psopt.Algorithm(collocation_method="Hermite-Simpson", nlp_tolerance=1.0e-6,
                      nlp_iter_max=1000, print_level=0,
                      parameter_statistics="yes")

sol_u = build().solve(alg)
ok = report(sol_u, "Unweighted fit (C++ example cracking)",
            reference=4.319519e-03, tol=1.0e-5)
p_u, ps_u = show(sol_u, "Unweighted: both channels asserted equally reliable")

# What weighting is actually for. Three samples of the second state are corrupted
# -- trebled, as a stuck sensor might -- and the fit is run twice: once asserting
# that every sample is equally reliable, and once weighting those three down by a
# factor of a hundred. The first estimate moves away from the clean one; the second
# comes back to it. Nothing here needs to know the true parameters, only which
# samples are suspect.
BAD = [8, 9, 10]
y2_bad = list(y2)
for i in BAD:
    y2_bad[i] *= 3.0

sol_c = build(observations=np.vstack([y1, y2_bad])).solve(alg)
p_c, _ = show(sol_c, "Corrupted data, every sample weighted equally")

W = np.ones((2, NS))
for i in BAD:
    W[1, i] = 0.01
sol_w = build(observations=np.vstack([y1, y2_bad]), weights=W).solve(alg)
p_w, ps_w = show(sol_w, "Corrupted data, the three bad samples weighted down 100x")

print("\nestimates, against the fit to the clean data:")
print("  %-10s %10s %10s %10s" % ("", "clean", "corrupted", "downweighted"))
for k in range(3):
    print("  theta%-5d %10.5f %10.5f %10.5f" % (k + 1, p_u[k], p_c[k], p_w[k]))
d_c = np.linalg.norm(np.asarray(p_c) - np.asarray(p_u))
d_w = np.linalg.norm(np.asarray(p_w) - np.asarray(p_u))
print("  distance from the clean fit: %.4f corrupted, %.4f downweighted" % (d_c, d_w))
ok = ok and d_w < 0.5 * d_c

# Tikhonov regularisation, for completeness: it pulls the parameter vector towards
# zero and is what you reach for when the data do not identify all of it.
sol_r = build(regularization=1.0e-2).solve(alg)
print("\nwith regularization_factor = 1e-2: theta = %s" %
      np.array2string(sol_r.parameters, precision=5))
print("  ||theta||: %.4f unregularised -> %.4f regularised"
      % (np.linalg.norm(p_u), np.linalg.norm(sol_r.parameters)))

ok = ok and ps_u is not None and ps_w is not None
ok = ok and np.linalg.norm(sol_r.parameters) < np.linalg.norm(p_u)
print("\nRESULT: %s" % ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
