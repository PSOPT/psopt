"""Minimum-time path past two circular obstacles -- path constraints from Python.

A vehicle moves at a fixed speed V = 2.138 with the heading as control, from
(0,0) to (1.2, 1.6), and must stay clear of two circles of radius sqrt(0.1)
centred at (0.4, 0.5) and (0.8, 1.5). Each is one path constraint,

    (x - xc)^2 + (y - yc)^2 >= 0.1,

imposed by giving the path vector a lower bound. This is the Python port of
examples/obstacle in the C++ distribution, whose answer is tf = 9.970637e-01 s.

It is the example to read for path constraints and for their multipliers: an
obstacle the trajectory does not touch has a multiplier of zero, and one it
grazes does not, which is complementary slackness and is visible in
sol.duals.dual_path.
"""
import numpy as np
import casadi as ca
from _common import psopt, report

V = 2.138
N = 30
prob = psopt.Problem(name="obstacle")
ph = prob.add_phase(nstates=2, ncontrols=1, nevents=4, npath=2)
ph.nodes = [N]

ph.dynamics = lambda x, u, p, t: ca.vertcat(V * ca.cos(u[0]), V * ca.sin(u[0]))
# The path vector: squared distance from each obstacle centre, bounded below.
ph.path = lambda x, u, p, t: ca.vertcat((x[0] - 0.4) ** 2 + (x[1] - 0.5) ** 2,
                                        (x[0] - 0.8) ** 2 + (x[1] - 1.5) ** 2)
ph.endpoint = lambda xi, xf, p, t0, tf: tf
ph.events = lambda xi, xf, p, t0, tf: ca.vertcat(xi[0], xi[1], xf[0], xf[1])

ph.bounds.lower.states = [0.0, 0.0]
ph.bounds.upper.states = [2.0, 2.0]
ph.bounds.lower.controls = [-np.pi / 2.0]
ph.bounds.upper.controls = [np.pi]
ph.bounds.lower.events = [0.0, 0.0, 1.2, 1.6]
ph.bounds.upper.events = [0.0, 0.0, 1.2, 1.6]
ph.bounds.lower.path = [0.1, 0.1]
ph.bounds.upper.path = [100.0, 100.0]
ph.bounds.t0 = (0.0, 0.0)
ph.bounds.tf = (0.5, 3.0)

# The guess matters here, and not because the solve is delicate. The problem has
# more than one taut path: going up and to the right of both circles gives
# tf = 0.997064, and the next homotopy class gives 0.999929, which is also a local
# minimum and which a heading guess of zero finds. So the guess sweeps the heading
# from 0.5 to 1.4 radians and integrates the kinematics -- it goes up and to the
# right, which is the only sensible way round, and nothing about it is fitted to
# the answer. This is the same guess the C++ example uses, where it is reported to
# be insensitive: over eighteen of twenty perturbations of its two angles it
# returns the same minimum to eight figures.
NG = 30
th = np.linspace(0.5, 1.4, NG)
dt = 1.0 / (NG - 1)
px = np.concatenate([[0.0], np.cumsum(V * np.cos(th) * dt)[:-1]])
py = np.concatenate([[0.0], np.cumsum(V * np.sin(th) * dt)[:-1]])
ph.guess.states = np.vstack([px, py])
ph.guess.controls = th.reshape(1, NG)
ph.guess.time = np.linspace(0.0, 1.0, NG).reshape(1, NG)

alg = psopt.Algorithm(collocation_method="Hermite-Simpson", nlp_tolerance=1.0e-6,
                      nlp_iter_max=1000, print_level=0,
                      mesh_refinement="automatic", ode_tolerance=1.0e-8)
sol = prob.solve(alg)

ok = report(sol, "Obstacle avoidance (C++ example obstacle)",
            reference=9.970637e-01, tol=1.0e-5)

x, y = sol.states[0], sol.states[1]
d1 = np.sqrt((x - 0.4) ** 2 + (y - 0.5) ** 2)
d2 = np.sqrt((x - 0.8) ** 2 + (y - 1.5) ** 2)
r = np.sqrt(0.1)
print("closest to #1  : %.6f   (radius %.6f, clearance %+.2e)" % (d1.min(), r, d1.min() - r))
print("closest to #2  : %.6f   (radius %.6f, clearance %+.2e)" % (d2.min(), r, d2.min() - r))
ok = ok and d1.min() >= r - 1.0e-6 and d2.min() >= r - 1.0e-6

# Complementary slackness: the multiplier is non-zero only where the constraint is
# active. Summing |mu| over each obstacle says which one actually shapes the path.
mu = sol.duals.dual_path
print("path multipliers: obstacle #1 sum|mu| = %.4f, obstacle #2 sum|mu| = %.4f"
      % (np.abs(mu[0]).sum(), np.abs(mu[1]).sum()))
print("                : a constraint that is never active has multiplier zero")

# A straight line would take this long; the detour is the price of the obstacles.
straight = np.hypot(1.2, 1.6) / V
print("straight-line   : %.6f s (infeasible); detour costs %+.2f%%"
      % (straight, 100.0 * (sol.objective - straight) / straight))

print("\nRESULT: %s" % ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
