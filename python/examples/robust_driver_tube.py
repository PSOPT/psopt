"""What an ancillary feedback gain buys, and what it costs.

Every robust design in the other driver examples is OPEN LOOP: one control
history, committed before the uncertainty is revealed, and required to serve
every plant in the set unaided. That is an honest formulation and it is an
expensive one. On the two-link arm it costs a factor of about three in final
time, and the assessment this work follows predicted exactly that --- an
open-loop robust trajectory over a long horizon is usually so conservative as to
be useless.

The standard remedy is a TUBE: a nominal trajectory, plus an ancillary feedback
law that pulls each realisation back towards it.

    u_k(t) = u_bar(t) + K ( x_k(t) - x_ref(t) )

Both u_bar and K are fixed before the uncertainty is revealed, so the design is
still here-and-now and nothing adapts to theta. What changes is that one history
no longer has to serve every plant by itself. This example solves the same
problem both ways and measures the difference.

WHERE THE GAIN COMES FROM, AND WHY IT IS GIVEN RATHER THAN DESIGNED

Co-designing K with the trajectory makes the problem bilinear in the decision
variables, so the standard first step is a fixed ancillary gain computed
separately. Here it is an LQR gain for the arm linearised about its target with
the nominal payload, which is a regulator for the end of the manoeuvre rather
than for the whole of it. That is crude and it is enough: the gain does not have
to be optimal, only stabilising, because the optimiser is still choosing the
whole nominal trajectory around it.

A time-varying gain from a Riccati sweep along the nominal trajectory would be
better and is the obvious next step. It is not here.

THE COST, WHICH IS NOT IN THE OBJECTIVE

A tube is not free, and the price shows up in a place the final time does not
record. The open-loop minimum-time solution is bang-bang: both torques sit
against their bounds almost everywhere. A correction added to a control already
at its bound saturates, and a saturated correction is not the control that was
designed. So the nominal trajectory has to back away from its bounds to leave
the ancillary law room to work, and the driver enforces that by imposing the
control bounds on the REALISED control of every scenario.

The verification checks it independently: the realised control is sampled far
more finely than the design constrains it, and the excess is reported. A design
whose corrections quietly ask for more actuator than exists is one no plant can
execute, whatever its certificate says about the states.
"""
import numpy as np
import casadi as ca
from scipy.linalg import solve_continuous_are
from _common import psopt
from psopt.robust import RobustProblem, Gaussian

MU, SIGMA, DELTA, SLACK = 0.50, 0.15, 0.03, 1.0e-3
X0 = np.array([0.0, 0.0, 0.500, 0.000])
XF = np.array([0.0, 0.0, 0.500, 0.522])
NODES = 25


def arm_rhs(x, u, mp):
    x1, x2, x3 = x[0], x[1], x[2]
    m11, m22, s2 = 7.0 / 3.0 + mp, 4.0 / 3.0 + mp, 1.5 + mp
    m12 = s2 * ca.cos(x3)
    det = m11 * m22 - m12 * m12
    b1 = (u[0] - u[1]) + s2 * ca.sin(x3) * x2 ** 2
    b2 = u[1] - s2 * ca.sin(x3) * x1 ** 2
    return ca.vertcat((m22 * b1 - m12 * b2) / det,
                      (-m12 * b1 + m11 * b2) / det, x2 - x1, x1)


# ---- the ancillary gain ----------------------------------------------------------
# LQR about the target, on the nominal plant. Q weights the two angles ten times the
# two rates, because it is the angles the terminal ball is drawn around; R is unit,
# which keeps the corrections comparable in size to the torques themselves.
xs = ca.SX.sym("x", 4)
us = ca.SX.sym("u", 2)
fs = arm_rhs(xs, us, MU)
A = np.array(ca.Function("A", [xs, us], [ca.jacobian(fs, xs)])(XF, np.zeros(2)))
B = np.array(ca.Function("B", [xs, us], [ca.jacobian(fs, us)])(XF, np.zeros(2)))
Q, R = np.diag([1.0, 1.0, 10.0, 10.0]), np.eye(2)
GAIN = -np.linalg.solve(R, B.T @ solve_continuous_are(A, B, Q, R))

print("Open loop against a tube, on the arm with an uncertain payload")
print("\n  ancillary gain, LQR about the target on the nominal plant:")
for r in range(2):
    print("    [%s ]" % " ".join("%8.4f" % GAIN[r, c] for c in range(4)))
print("  closed-loop eigenvalues %s"
      % np.array2string(np.linalg.eigvals(A + B @ GAIN), precision=3))

# ---- the problem, built once ------------------------------------------------------
rp = RobustProblem(name="robust_driver_tube")
ph = rp.add_phase(nstates=4, ncontrols=2, nevents=8)
ph.nodes = [NODES]
ph.dynamics = lambda x, u, p, t, th: arm_rhs(x, u, th[0])
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
ph.bounds.tf = (1.0, 15.0)
ph.guess.states = np.tile(X0.reshape(4, 1), (1, NODES))
ph.guess.controls = np.zeros((2, NODES))
ph.guess.time = np.linspace(0.0, 3.0, NODES).reshape(1, NODES)
rp.uncertainty = Gaussian(mean=[MU], cov=[[SIGMA ** 2]], truncate=3.0)
rp.initial_state = X0

# ms_path_samples is the setting that makes the tube implementable, and it is
# worth understanding why. Multiple shooting imposes path constraints at the
# segment boundaries; between them nothing holds, and the bound that matters here
# -- the realised control staying inside the actuator's range -- is a path
# constraint. Left at the default the correction breaches it by 1.83e-02 between
# the nodes; at four samples per segment it breaches by 8.13e-04, a factor of 22,
# and the final time does not move. The facility exists for exactly this, and it
# samples only INEQUALITY components of the path vector, which is what these are.
alg = psopt.Algorithm(transcription_method="multiple-shooting",
                      ms_integrator="RK4", ms_steps_per_segment=12,
                      ms_control_parameterisation="linear",
                      ms_path_samples=4,
                      scaling="automatic", nlp_tolerance=1.0e-6,
                      nlp_iter_max=1500, print_level=0)


def run(gain, label):
    ph.feedback = gain
    out = rp.solve(alg, slack=SLACK, scenarios="sigma-points", generate=True,
                   max_iterations=10, tighten=0.9, out_of_sample=250,
                   wait_and_see=0, verbose=False)
    c, o = out.certificate, out.out_of_sample
    print("  %-12s t_f %7.4f  M %2d  worst %8.2e  within %5.1f%%  u beyond %7.2e"
          % (label, out.objective, len(out.scenarios), c["violation"],
             o["within"], c["control_excess"]))
    return out


# The nominal design, for the baseline that makes the comparison mean anything.
ph.feedback = None
rp.uncertainty = Gaussian(mean=[MU], cov=[[1.0e-12]], truncate=3.0)
nom = rp.solve(alg, slack=SLACK, scenarios="sigma-points", generate=False,
               out_of_sample=0, wait_and_see=0, verbose=False)
rp.uncertainty = Gaussian(mean=[MU], cov=[[SIGMA ** 2]], truncate=3.0)

print("\n  %-12s %7s     %2s  %8s  %6s  %9s"
      % ("design", "t_f", "M", "worst", "within", "u beyond"))
print("  %-12s t_f %7.4f  M %2d  (the deterministic design, for comparison)"
      % ("nominal", nom.objective, 1))
openloop = run(None, "open loop")
tube = run(GAIN, "tube")

# ---- the closed loop, against an independent implementation -----------------------
#
# The driver integrates the closed loop inside its own vectorised sweep. This is a
# separate NumPy implementation of the same law -- reference and scenario marched
# together, the correction formed from their difference -- sharing nothing with the
# driver but the equations of motion. Without it, a verification of a closed-loop
# design is only a claim that two parts of the same file agree.
def closed_loop_miss(t, u, mp_grid, gain, nsub=32):
    k = len(mp_grid)
    xr = np.tile(X0.reshape(4, 1), (1, k)).astype(float)     # the reference
    xk = xr.copy()                                           # the scenarios

    def f(xx, uu, mp):
        x1, x2, x3 = xx[0], xx[1], xx[2]
        m11, m22, s2 = 7.0 / 3.0 + mp, 4.0 / 3.0 + mp, 1.5 + mp
        m12 = s2 * np.cos(x3)
        det = m11 * m22 - m12 * m12
        b1 = (uu[0] - uu[1]) + s2 * np.sin(x3) * x2 ** 2
        b2 = uu[1] - s2 * np.sin(x3) * x1 ** 2
        return np.vstack([(m22 * b1 - m12 * b2) / det,
                          (-m12 * b1 + m11 * b2) / det, x2 - x1, x1])

    mu_col = np.full(k, MU)
    for i in range(len(t) - 1):
        h = (t[i + 1] - t[i]) / nsub
        for sp in range(nsub):
            def stage(xr_, xk_, w):
                ub = (u[:, i] + w * (u[:, i + 1] - u[:, i])).reshape(2, 1)
                uu = ub if gain is None else ub + gain @ (xk_ - xr_)
                return (f(xr_, np.tile(ub, (1, k)), mu_col),
                        f(xk_, np.tile(uu, (1, k))[:2] if uu.shape[1] == 1 else uu,
                          mp_grid))
            a1, b1 = stage(xr, xk, sp / nsub)
            a2, b2 = stage(xr + 0.5 * h * a1, xk + 0.5 * h * b1, (sp + 0.5) / nsub)
            a3, b3 = stage(xr + 0.5 * h * a2, xk + 0.5 * h * b2, (sp + 0.5) / nsub)
            a4, b4 = stage(xr + h * a3, xk + h * b3, (sp + 1.0) / nsub)
            xr = xr + (h / 6.0) * (a1 + 2 * a2 + 2 * a3 + a4)
            xk = xk + (h / 6.0) * (b1 + 2 * b2 + 2 * b3 + b4)
    return np.max(np.abs(xk - XF.reshape(4, 1)), axis=0)


grid = np.linspace(MU - 3 * SIGMA, MU + 3 * SIGMA, 601)
mine = closed_loop_miss(tube.time, tube.controls, grid, GAIN)
theirs = rp._violation_many(grid.reshape(-1, 1), tube.time, tube.controls,
                            np.zeros(0))
gap = float(np.max(np.abs(np.maximum(mine - DELTA, 0.0) - theirs)))

print("\n  the closed loop against an independent NumPy implementation, %d payloads:"
      % len(grid))
print("    largest disagreement in max(0, miss - delta) : %.2e" % gap)
print("    worst terminal miss anywhere in the set      : %.4f  (the ball is %.2f)"
      % (mine.max(), DELTA))

# ---- what it bought ---------------------------------------------------------------
op = 100.0 * (openloop.objective - nom.objective) / nom.objective
tb = 100.0 * (tube.objective - nom.objective) / nom.objective
print("\n  price of robustness, against the deterministic design:")
print("    open loop : %+6.1f%%   (%d scenarios to certify)"
      % (op, len(openloop.scenarios)))
print("    tube      : %+6.1f%%   (%d scenarios to certify)"
      % (tb, len(tube.scenarios)))
print("\n  The tube gives back %.0f%% of what the open-loop design paid, and needs"
      % (100.0 * (op - tb) / op))
print("  half as many scenarios to certify. It is not free: the nominal torques")
print("  have to leave room for the correction, so the bang-bang structure of the")
print("  deterministic solution is gone; the realised control has to be held")
print("  inside the segments and not only at their ends; and the gain needs the")
print("  state to be measurable, which the open-loop design does not.")

ok = (tube.converged and openloop.converged
      and tube.objective < openloop.objective
      and tube.certificate["violation"] <= SLACK
      and tube.certificate["control_excess"] < 2.0e-3
      and gap < 1.0e-3
      and mine.max() <= DELTA + 10.0 * SLACK)
print("\nobjective      : %.12g" % tube.objective)
print("RESULT: %s" % ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
