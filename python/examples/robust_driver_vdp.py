"""Two uncertain parameters, a state constraint, and an expected cost.

The companion to ``robust_driver_arm.py``, and deliberately unlike it in every
structural respect the driver has to cope with:

    arm                                 this
    ---                                 ----
    one uncertain parameter             two, with correlated uncertainty
    free final time                     fixed final time
    Mayer objective (t_f)               Lagrange objective, and an EXPECTATION
    terminal constraints only           a state PATH constraint, for every plant
    worst case found by scanning a line found by searching a plane

The plant is the van der Pol oscillator with neither of its two coefficients
known exactly,

    x1' = x2,    x2' = -w x1 + mu (1 - x1^2) x2 + u,

controlled from (1, 0) over five seconds at least expected cost, and required to
keep x1 above a barrier for EVERY plant in the uncertainty set rather than for
the nominal one. The stiffness w and the nonlinearity mu are negatively
correlated, which is what a covariance estimated from data usually looks like and
is why the uncertainty set is a tilted ellipse rather than a box.

A TRAP WORTH KNOWING ABOUT BEFORE REACHING FOR A DRIVER

The obvious second example is a vehicle with an uncertain speed, or a double
integrator with an uncertain control gain, and neither works. The reason is
structural.

Take x2' = -c x2 + b u with b uncertain and the terminal state pinned. The
sensitivity of the final state to b obeys d(dx2)/dt = -c dx2 + u, whose
homogeneous part is -c whatever the control does, so
dx2(t_f)/db = integral exp(-c(t_f - s)) u(s) ds, which the terminal condition
itself fixes at (x2(t_f) - exp(-c t_f) x2(0)) / b. No control can reduce it. The
same argument applies to a kinematic vehicle, where the derivative of the final
position with respect to a speed error IS the total displacement, fixed by the
two endpoints.

So for a LINEAR plant with an uncertain gain and a pinned terminal state, robust
open-loop design buys exactly nothing, and the augmented problem is infeasible for
any tolerance smaller than the unavoidable spread. What the arm and this
oscillator have in common, and what those examples lack, is that the sensitivity's
own homogeneous part depends on the trajectory --- so the control can shape it,
and there is something for a robust design to do. Worth checking for in a new
problem before running anything.

WHAT THE SCENARIO SET IS FOR, WHICH IS TWO THINGS

This example uses ``risk="expectation"``, so the scenario set is at once a
quadrature rule for the expected cost and a constraint set for feasibility. The
driver keeps those apart: the five unscented points carry the rule's weights, and
every scenario the generation loop adds afterwards enters the constraints with
weight zero. They are not quadrature nodes --- they are the places the design was
failing, which is a biased sample of the uncertainty by construction, and letting
them into the expectation would quietly replace the risk measure with a
worst-case-weighted one.
"""
import numpy as np
import casadi as ca
from _common import psopt
from psopt.robust import RobustProblem, Gaussian, Explicit

BARRIER = -0.40         # the state constraint, x1 >= BARRIER, for every plant
SLACK = 0.0             # violation of the declared barrier a design may leave.
                        # Not the barrier itself: BARRIER is declared in the path
                        # bounds, and slack is how far BELOW it the design may
                        # still dip somewhere in the uncertainty set. Zero here,
                        # because the margin below makes zero attainable.
TF, NODES = 5.0, 31
X0 = np.array([1.0, 0.0])

# w: stiffness. mu: nonlinearity. Negatively correlated, as an estimated pair
# usually is, which tilts the uncertainty ellipse off the coordinate axes and
# makes the worst case a genuinely two-dimensional search.
MEAN = np.array([1.00, 1.00])
SD = np.array([0.12, 0.20])
RHO = -0.35
COV = np.array([[SD[0] ** 2, RHO * SD[0] * SD[1]],
                [RHO * SD[0] * SD[1], SD[1] ** 2]])

rp = RobustProblem(name="robust_driver_vdp")
ph = rp.add_phase(nstates=2, ncontrols=1, nevents=2, npath=1)
ph.nodes = [NODES]

ph.dynamics = lambda x, u, p, t, th: ca.vertcat(
    x[1], -th[0] * x[0] + th[1] * (1.0 - x[0] ** 2) * x[1] + u[0])
ph.integrand = lambda x, u, p, t: 0.5 * (x[0] ** 2 + x[1] ** 2 + u[0] ** 2)
ph.path = lambda x, u, p, t, th: ca.vertcat(x[0])
# Only the initial state is an event, and it is shared by every scenario. The
# terminal state is free, which is what lets this problem escape the trap above:
# nothing forces a fixed terminal sensitivity.
ph.events = lambda xi, xf, p, t0, tf, th: ca.vertcat(xi[0], xi[1])

ph.bounds.lower.states = [-5.0, -5.0]
ph.bounds.upper.states = [5.0, 5.0]
ph.bounds.lower.controls = [-10.0]
ph.bounds.upper.controls = [10.0]
ph.bounds.lower.events = list(X0)
ph.bounds.upper.events = list(X0)
ph.bounds.lower.path = [BARRIER]
ph.bounds.upper.path = [5.0]
ph.bounds.t0 = (0.0, 0.0)
ph.bounds.tf = (TF, TF)

ph.guess.states = np.vstack([np.linspace(1.0, 0.0, NODES), np.zeros(NODES)])
ph.guess.controls = np.zeros((1, NODES))
ph.guess.time = np.linspace(0.0, TF, NODES).reshape(1, NODES)

rp.uncertainty = Gaussian(mean=MEAN, cov=COV, truncate=2.5)
rp.initial_state = X0

alg = psopt.Algorithm(transcription_method="multiple-shooting",
                      ms_integrator="RK4", ms_steps_per_segment=10,
                      ms_control_parameterisation="linear",
                      scaling="automatic", nlp_tolerance=1.0e-7,
                      nlp_iter_max=3000, print_level=0)

print("van der Pol with an uncertain plant, through psopt.robust")
print("  stiffness and nonlinearity, correlation %.2f; x1 >= %.2f for every plant"
      % (RHO, BARRIER))

out = rp.solve(alg, slack=SLACK, risk="expectation", scenarios="sigma-points",
               generate=True, max_iterations=10,
               # The barrier is ONE-SIDED, so there is no two-sided half-width for
               # `tighten` to take a fraction of and the driver would give it no
               # margin at all. A one-sided constraint left exactly at its bound is
               # where the miss between two scenarios appears, so it is set here.
               margin=dict(events=None, path=np.array([0.04])),
               n_seed=256, out_of_sample=500, wait_and_see=3)

# ---- the nominal design, for the comparison that makes the numbers mean anything -
rp.uncertainty = Explicit([MEAN])
nom = rp.solve(alg, slack=SLACK, scenarios="explicit", generate=False,
               out_of_sample=0, verbose=False)
rp.uncertainty = Gaussian(mean=MEAN, cov=COV, truncate=2.5)
nom_worst, nom_th = rp._worst_case(nom.time, nom.controls, np.zeros(0), 256, 3, 1)
nom_score = rp._score(nom.time, nom.controls, np.zeros(0), 500, SLACK, 20260927)

print("\n  %-26s %10s %13s %20s" % ("design", "J", "worst in set", "at (w, mu)"))
print("  %-26s %10.5f %13.3e   (%.4f, %.4f)"
      % ("nominal (mean plant only)", nom.objective, nom_worst, nom_th[0], nom_th[1]))
print("  %-26s %10.5f %13.3e   (%.4f, %.4f)"
      % ("robust (%d scenarios)" % len(out.scenarios), out.objective,
         out.certificate["violation"], out.certificate["theta"][0],
         out.certificate["theta"][1]))
# "Violation zero" is the right answer and an uninformative number. The quantity
# with meaning is how far x1 actually dips, over the whole set, under each design.
# Asking the violation measure for it needs no new machinery: a lower bound placed
# above anything reachable turns the measure into max(bound - x1), so the worst
# excursion falls straight out of it.
PROBE = 5.0
ph.bounds.lower.path = [PROBE]
probe_set = rp.uncertainty.quasi_random(256, seed=3)
dip_robust = PROBE - float(rp._violation_many(probe_set, out.time, out.controls,
                                              np.zeros(0)).max())
dip_nominal = PROBE - float(rp._violation_many(probe_set, nom.time, nom.controls,
                                               np.zeros(0)).max())
ph.bounds.lower.path = [BARRIER]
print("\n  lowest x1 reached anywhere in the uncertainty set")
print("    nominal design : %+.4f   (the barrier is %.2f, so it breaks it)"
      % (dip_nominal, BARRIER))
print("    robust design  : %+.4f   (clear of the barrier, and of the %.2f the"
      % (dip_robust, BARRIER + 0.04))
print("                              design was actually given)")

print("\n  within tolerance out of sample : nominal %.1f%%, robust %.1f%%"
      % (nom_score["within"], out.out_of_sample["within"]))
print("  price of robustness            : J %.5f -> %.5f, a factor of %.2f"
      % (nom.objective, out.objective, out.objective / nom.objective))

# Does the worst case really need a plane to be searched, or would a lattice of
# the same budget do? The comparison is made on the NOMINAL design, because that
# is where there is a worst case to find: the robust design has none anywhere in
# the set, so measuring two searches against it would compare two zeros.
g = int(np.sqrt(256))
grid = np.array([[a, b]
                 for a in np.linspace(MEAN[0] - 2.5 * SD[0], MEAN[0] + 2.5 * SD[0], g)
                 for b in np.linspace(MEAN[1] - 2.5 * SD[1], MEAN[1] + 2.5 * SD[1], g)])
grid = np.array([t for t in grid if rp.uncertainty.contains(t)])
grid_worst = float(rp._violation_many(grid, nom.time, nom.controls,
                                      np.zeros(0)).max())
print("\n  searching the NOMINAL design's worst case, on a budget of 256 points:")
print("    a %dx%d lattice, %d of whose points are in the set : %.4f"
      % (g, g, len(grid), grid_worst))
print("    the driver's oracle                               : %.4f" % nom_worst)
print("    the lattice misses %.1f%% of it."
      % (100.0 * (nom_worst - grid_worst) / nom_worst))

ok = (out.converged
      and out.certificate["violation"] <= SLACK
      and out.out_of_sample["within"] >= 95.0
      and dip_robust > BARRIER
      and nom_worst > 0.1
      and out.certificate["integrator_drift"] < 1.0e-4)
print("\nobjective      : %.12g" % out.objective)
print("RESULT: %s" % ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
