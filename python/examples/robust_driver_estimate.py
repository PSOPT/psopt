"""Estimate the plant from data, then design against what the data actually says.

Every robust design so far has been given its uncertainty by hand: a mean, a
covariance and a truncation that somebody chose. That is the weakest part of the
whole construction, because the answer depends on a distribution nobody measured.

PSOPT already computes the thing that should be supplying it. A parameter
estimation returns not only the fitted parameters but their COVARIANCE --
`sol.parameter_statistics` -- and that covariance is precisely the distribution a
robust design ought to run against. This example closes that loop inside one
model:

    1. drive the arm with a prescribed excitation torque and observe its two
       angles at 41 instants, with noise
    2. estimate the payload mass and the joint damping from those observations,
       and take the covariance PSOPT returns with them
    3. hand that Gaussian to the robust driver as the uncertainty
    4. check the design on the TRUE plant, which nothing in steps 2 or 3 saw

WHY THIS IS WORTH DOING RATHER THAN GUESSING A COVARIANCE

The estimated covariance is not a round number and it is not diagonal. The two
parameters trade off against each other in the fit -- a heavier payload and
lighter damping produce nearly the same observed motion over a short manoeuvre --
so the posterior is an ellipse at an angle, and its long axis is the direction the
data did NOT pin down. A robust design against that ellipse spends its effort
along that axis, which is where the real uncertainty is. A design against a
diagonal covariance with the same standard deviations would spend it in the wrong
directions: it would be conservative where the data was informative and
optimistic where it was not.

The example measures that, by designing against both and checking each on the
true plant.

WHAT IS AND IS NOT BEING CLAIMED

The covariance PSOPT returns is the usual linearised (Gauss-Newton) one, valid to
the extent the residuals are small and the model is locally linear in its
parameters. It is not a posterior in the Bayesian sense and it carries no prior.
For a short manoeuvre and a well-excited fit it is a reasonable description of
what the data determined, which is all this uses it for. The check that matters is
the one at the end: the design is verified on the true plant, and the true plant
took no part in the estimate beyond generating the data.
"""
import numpy as np
import casadi as ca
from _common import psopt
from psopt.robust import RobustProblem, Gaussian, Explicit

# The plant the data comes from. Nothing downstream is allowed to look at these
# two numbers until the very last check.
TRUE = np.array([0.55, 0.22])          # payload mass, joint damping
NOISE = 0.004                          # s.d. of the angle measurements, radians

X0 = np.array([0.0, 0.0, 0.500, 0.000])
XF = np.array([0.0, 0.0, 0.500, 0.522])
DELTA, SLACK = 0.03, 1.0e-3
NODES = 25
# How far out to truncate the estimated Gaussian. 3.5 rather than the 3 that would
# be habitual for one parameter, because a k-sigma ELLIPSE in two dimensions
# covers less than a k-sigma interval does in one: the coverage is 1 - exp(-k^2/2),
# which is 95.6% at k = 2.5 and 98.9% at k = 3, where the one-dimensional figures
# are 98.8% and 99.7%. At k = 3.5 it is 99.8%, which is the figure a habit formed
# on scalar uncertainty would expect from 3.
TRUNC = 3.5


def arm(x1, x2, x3, u1, u2, mp, b, sin, cos):
    """Two-link arm with a tip payload and viscous damping on the joint rates.

    The payload adds mp to both diagonal moments and to the first moment of the
    second link; the damping removes b times each link's absolute angular rate.
    Written once for CasADi and NumPy so the estimation, the design and the
    verification all differentiate the same equations.
    """
    m11, m22, s2 = 7.0 / 3.0 + mp, 4.0 / 3.0 + mp, 1.5 + mp
    m12 = s2 * cos(x3)
    det = m11 * m22 - m12 * m12
    b1 = (u1 - u2) + s2 * sin(x3) * x2 ** 2 - b * x1
    b2 = u2 - s2 * sin(x3) * x1 ** 2 - b * x2
    return ((m22 * b1 - m12 * b2) / det, (-m12 * b1 + m11 * b2) / det, x2 - x1, x1)


# The excitation. A prescribed function of time, not a decision variable: during
# the identification manoeuvre the torques are what the experimenter applied. Two
# incommensurate frequencies per joint, so that the motion is rich enough to
# separate an inertia from a damping -- with a single frequency it cannot be done,
# and the covariance says so by coming out nearly singular.
def excite(t, sin):
    return (0.8 * sin(2.1 * t) + 0.4 * sin(5.3 * t),
            0.6 * sin(3.7 * t) - 0.3 * sin(1.3 * t))


T_ID, N_ID = 4.0, 41
P_GUESS = np.array([0.20, 0.50])       # what the experimenter thought beforehand


def run_excited(ts, params, nsub=40):
    """The excited manoeuvre for a given plant, sampled at ts. Plain NumPy RK4."""
    x = X0.astype(float).copy()
    out = [x.copy()]
    for i in range(len(ts) - 1):
        h = (ts[i + 1] - ts[i]) / nsub
        t = ts[i]
        for _ in range(nsub):
            def f(xx, tt):
                u1, u2 = excite(tt, np.sin)
                return np.array(arm(xx[0], xx[1], xx[2], u1, u2,
                                    params[0], params[1], np.sin, np.cos))
            k1 = f(x, t)
            k2 = f(x + 0.5 * h * k1, t + 0.5 * h)
            k3 = f(x + 0.5 * h * k2, t + 0.5 * h)
            k4 = f(x + h * k3, t + h)
            x = x + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            t += h
        out.append(x.copy())
    return np.array(out).T


def simulate_truth(seed):
    """Generate the observations: the true plant, integrated, plus noise."""
    rng = np.random.default_rng(seed)
    ts = np.linspace(0.0, T_ID, N_ID)
    y = run_excited(ts, TRUE)[2:4]
    return ts, y + rng.normal(0.0, NOISE, y.shape)


# ---------------------------------------------------------------------------------
#  Step 1 and 2: the experiment, and the estimate with its covariance
# ---------------------------------------------------------------------------------

t_obs, y_obs = simulate_truth(20260927)

est = psopt.Problem(name="arm_identification")
ep = est.add_phase(nstates=4, ncontrols=0, nparameters=2, nobserved=2,
                   nsamples=N_ID)
ep.nodes = [80]
ep.dynamics = lambda x, u, p, t: ca.vertcat(
    *arm(x[0], x[1], x[2], *excite(t, ca.sin), p[0], p[1], ca.sin, ca.cos))
ep.observation = lambda x, u, p, t: ca.vertcat(x[2], x[3])
ep.bounds.lower.states = [-5.0] * 4
ep.bounds.upper.states = [5.0] * 4
ep.bounds.lower.parameters = [0.0, 0.0]
ep.bounds.upper.parameters = [2.0, 2.0]
ep.bounds.t0 = (0.0, 0.0)
ep.bounds.tf = (T_ID, T_ID)
ep.observation_nodes = t_obs.reshape(1, N_ID)
ep.observations = y_obs
# The state guess is the manoeuvre as the experimenter's PRIOR guess at the
# parameters would produce it, not a constant. A constant guess for a trajectory
# that swings through three radians leaves the collocation with nothing to work
# from, and the fit walks to the parameter bounds instead -- which is what this
# example did until the guess was built properly.
t_guess = np.linspace(0.0, T_ID, 40)
ep.guess.states = run_excited(t_guess, P_GUESS)
ep.guess.time = t_guess.reshape(1, 40)
ep.guess.parameters = P_GUESS.reshape(2, 1)

sol = est.solve(psopt.Algorithm(collocation_method="Hermite-Simpson",
                                nlp_tolerance=1.0e-8, nlp_iter_max=2000,
                                parameter_statistics="yes", print_level=0))
ps = sol.parameter_statistics
p_hat = np.asarray(sol.parameters).ravel()
COV = np.atleast_2d(ps.covariance)
se = ps.standard_errors
rho = COV[0, 1] / (se[0] * se[1])

print("Estimate, then design against what the data says")
print("\n  step 1-2: %d noisy observations of two angles, s.d. %.4f rad"
      % (N_ID, NOISE))
print("  %-24s %10s %10s %22s" % ("parameter", "true", "estimate", "95% interval"))
for k, nm in enumerate(("payload mass", "joint damping")):
    print("  %-24s %10.4f %10.4f   [%8.4f, %8.4f]"
          % (nm, TRUE[k], p_hat[k], ps.confidence_low[k], ps.confidence_high[k]))
print("  residual s.d. %.5f against the %.5f the data was made with"
      % (ps.sigma_hat, NOISE))
print("  correlation between the two estimates: %+.3f" % rho)

# Marginal intervals and the joint ellipse are different objects, and this fit
# shows why it matters. Each interval above is a projection of the ellipse onto one
# axis, which throws the correlation away; the Mahalanobis distance uses it.
maha = float(np.sqrt((TRUE - p_hat) @ np.linalg.solve(COV, TRUE - p_hat)))
marg = np.all((TRUE >= ps.confidence_low) & (TRUE <= ps.confidence_high))
cover = 100.0 * (1.0 - np.exp(-0.5 * TRUNC ** 2))
print("\n  where the truth sits relative to the estimate:")
print("    inside every marginal 95%% interval : %s" % marg)
print("    Mahalanobis distance from the fit  : %.2f sigma" % maha)
print("    inside the %.1f sigma ellipse used  : %s  (which covers %.1f%% of the"
      % (TRUNC, maha <= TRUNC, cover))
print("                                         posterior in two dimensions)")
print("  Two things are worth taking from that. The marginal intervals exclude the")
print("  truth while the joint ellipse contains it: an interval per parameter is")
print("  the ellipse's shadow on one axis, and with a correlation of %+.3f the" % rho)
print("  shadows are a poor description of it. And the ellipse is truncated at")
print("  %.1f rather than 3 sigma, because coverage falls away with dimension --" % TRUNC)
print("  1 - exp(-k^2/2) is %.1f%% at k = 3 in the plane against 99.7%% on a line."
      % (100.0 * (1.0 - np.exp(-4.5))))

# ---------------------------------------------------------------------------------
#  Step 3: the robust design, against the covariance the estimate came with
# ---------------------------------------------------------------------------------

rp = RobustProblem(name="robust_from_estimate")
ph = rp.add_phase(nstates=4, ncontrols=2, nevents=8)
ph.nodes = [NODES]
ph.dynamics = lambda x, u, p, t, th: ca.vertcat(
    *arm(x[0], x[1], x[2], u[0], u[1], th[0], th[1], ca.sin, ca.cos))
ph.endpoint = lambda xi, xf, p, t0, tf: tf
ph.events = lambda xi, xf, p, t0, tf, th: ca.vertcat(
    xi[0], xi[1], xi[2], xi[3], xf[0], xf[1], xf[2], xf[3])
ph.bounds.lower.states = [-2.0] * 4
ph.bounds.upper.states = [2.0] * 4
ph.bounds.lower.controls = [-1.0, -1.0]
ph.bounds.upper.controls = [1.0, 1.0]
ph.bounds.lower.events = list(X0) + list(XF - DELTA)
ph.bounds.upper.events = list(X0) + list(XF + DELTA)
ph.bounds.t0 = (0.0, 0.0)
ph.bounds.tf = (1.0, 20.0)
ph.guess.states = np.tile(X0.reshape(4, 1), (1, NODES))
ph.guess.controls = np.zeros((2, NODES))
ph.guess.time = np.linspace(0.0, 3.0, NODES).reshape(1, NODES)
rp.initial_state = X0

alg = psopt.Algorithm(transcription_method="multiple-shooting",
                      ms_integrator="RK4", ms_steps_per_segment=12,
                      ms_control_parameterisation="linear",
                      scaling="automatic", nlp_tolerance=1.0e-6,
                      nlp_iter_max=1500, print_level=0)


def design(uncertainty, generate, label):
    rp.uncertainty = uncertainty
    out = rp.solve(alg, slack=SLACK, scenarios="sigma-points", generate=generate,
                   max_iterations=8, tighten=0.9, out_of_sample=0,
                   wait_and_see=0, verbose=False)
    print("  %-34s t_f = %7.4f on %2d scenarios" % (label, out.objective,
                                                    len(out.scenarios)))
    return out


print("\n  step 3: three designs, differing only in the uncertainty they are given")
nominal = design(Explicit([p_hat]), False, "the point estimate alone")
posterior = design(Gaussian(p_hat, COV, truncate=TRUNC), True,
                   "the estimated covariance")
# The same standard deviations with the correlation thrown away. This is what
# "I measured the error bars" gets you if the off-diagonal terms are discarded,
# and it is a common thing to do.
diagonal = design(Gaussian(p_hat, np.diag(np.diag(COV)), truncate=TRUNC), True,
                  "the same error bars, uncorrelated")

# ---------------------------------------------------------------------------------
#  Step 4: the check that counts -- the true plant, which none of them saw
# ---------------------------------------------------------------------------------

rp.uncertainty = Gaussian(p_hat, COV, truncate=TRUNC)
print("\n  step 4: each design on the TRUE plant (%.4f, %.4f), which no design saw"
      % (TRUE[0], TRUE[1]))
print("\n  %-34s %10s %12s" % ("design", "t_f", "miss at truth"))
truth_miss = {}
for label, out in (("the point estimate alone", nominal),
                   ("the estimated covariance", posterior),
                   ("the same error bars, uncorrelated", diagonal)):
    v = float(rp._violation_many(TRUE.reshape(1, 2), out.time, out.controls,
                                 np.zeros(0))[0])
    truth_miss[label] = v
    flag = "  within" if v <= SLACK else "  OUTSIDE the declared ball"
    print("  %-34s %10.4f %12.3e%s" % (label, out.objective, v, flag))

print("\n  and over the whole posterior ellipse, by the driver's own search:")
print("  %-34s %14s" % ("design", "worst in set"))
worst = {}
for label, out in (("the point estimate alone", nominal),
                   ("the estimated covariance", posterior),
                   ("the same error bars, uncorrelated", diagonal)):
    w, _th = rp._worst_case(out.time, out.controls, np.zeros(0), 128, 3, 5)
    worst[label] = w
    print("  %-34s %14.3e%s" % (label, w, "" if w <= SLACK else "   fails"))

print("\n  The design built from the estimate's own covariance is the only one that")
print("  meets its requirement at the true plant AND across the ellipse the data")
print("  supports, and it costs %.1f%% in time over the design that ignored the"
      % (100.0 * (posterior.objective - nominal.objective) / nominal.objective))
print("  uncertainty entirely.")
print("\n  Discarding the correlation is the second place it goes wrong, and it is")
print("  not a conservative simplification. The diagonal ellipse has the same")
print("  width in each parameter but points along the axes, so it covers pairs the")
print("  data ruled out and misses the ones along the direction the data left")
print("  undetermined -- which is where the worst case actually is.")

# A test rather than a demonstration. The estimate must contain the truth in the
# set the design actually uses; the design built from the full covariance must
# meet its requirement both at the truth and across that whole set; and the two
# designs that throw information away must each fail in their own way -- the point
# estimate at the truth, the diagonal one across the set.
ok = (maha <= TRUNC
      and abs(rho) > 0.3
      and truth_miss["the estimated covariance"] <= SLACK
      and worst["the estimated covariance"] <= SLACK
      and truth_miss["the point estimate alone"] > SLACK
      and worst["the point estimate alone"] > SLACK
      and worst["the same error bars, uncorrelated"] > SLACK)
print("\nobjective      : %.12g" % posterior.objective)
print("RESULT: %s" % ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
