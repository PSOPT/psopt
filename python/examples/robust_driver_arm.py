"""The robust two-link arm again, through the driver instead of by hand.

``robust_arm.py`` writes out the whole robust design loop: the augmented problem,
the verification integrator, the worst-case scan, the tightening, the warm start.
This example is the same problem posed through ``psopt.robust``, and it exists to
check the driver against something already measured rather than to demonstrate
anything new about the problem.

What the user writes is the NOMINAL problem, once, with the uncertain parameter
as one extra argument to the dynamics and the events. Everything else --- the M
copies of the state, the one copy of the control that makes the design
non-anticipative, the replicated bounds, the independent verifier, the search for
the worst payload, the scenario that gets added and the warm start it is solved
from --- is the driver's.

Two things in the printout are worth reading against ``robust_arm.py``:

  * the starting scenario set. The unscented set for a scalar Gaussian is
    {mu, mu +- sqrt(3) sigma} with weights {2/3, 1/6, 1/6}, which is exactly the
    three-point Gauss-Hermite rule that example uses. The two agree by
    construction and not by coincidence, and it is worth knowing that the "sigma
    points" of the estimation literature and the "Gauss-Hermite nodes" of the
    quadrature literature are the same three numbers here.

  * the certificate. The driver reports the worst violation it FOUND and how many
    integrations went into finding it. For one uncertain parameter that search is
    effectively exhaustive; the point of saying it that way is that in more than
    one parameter it is not, and the wording should not have to change.
"""
import numpy as np
import casadi as ca
from _common import psopt
from psopt.robust import RobustProblem, Gaussian

MU, SIGMA, DELTA = 0.50, 0.15, 0.03
X0 = np.array([0.0, 0.0, 0.500, 0.000])
XF = np.array([0.0, 0.0, 0.500, 0.522])
# 25 rather than the 21 of robust_arm.py. The reason WAS the scenario budget: the
# driver used to run out of degrees of freedom at ten scenarios on a 21-node mesh,
# and this problem needs more. That budget turned out to be phantom -- PSOPT's
# defect block holds nstates*(nodes) rows and multiple shooting fills only
# nstates*(nodes-1) of them, the rest being equality rows of zeros that IPOPT
# counts against the variables, nstates of them PER SCENARIO. The driver now sets
# algorithm.free_padded_defect_rows, the problem has 51 degrees of freedom at every
# scenario count, and the arm carries eighty scenarios where it used to be refused
# above twelve. 25 nodes is left as it stands because every number printed below
# was measured on it; 21 would now do.
NODES = 25


def arm_rhs(x, u, mp):
    """Planar two-link arm in absolute coordinates, carrying a payload of mass mp.

    The payload enters only through the inertia: a point mass at the tip of the
    second link adds mp to both diagonal moments and to the first moment of the
    second link about the second joint. Setting mp = 0 recovers the coefficients
    of ``examples/twolinkarm``.
    """
    x1, x2, x3 = x[0], x[1], x[2]
    m11, m22, s2 = 7.0 / 3.0 + mp, 4.0 / 3.0 + mp, 1.5 + mp
    m12 = s2 * ca.cos(x3)
    det = m11 * m22 - m12 * m12
    b1 = (u[0] - u[1]) + s2 * ca.sin(x3) * x2 ** 2
    b2 = u[1] - s2 * ca.sin(x3) * x1 ** 2
    return ca.vertcat((m22 * b1 - m12 * b2) / det,
                      (-m12 * b1 + m11 * b2) / det,
                      x2 - x1,
                      x1)


rp = RobustProblem(name="robust_driver_arm")
ph = rp.add_phase(nstates=4, ncontrols=2, nevents=8)
ph.nodes = [NODES]

# theta is the uncertain parameter vector; here it has one component, the payload.
ph.dynamics = lambda x, u, p, t, th: arm_rhs(x, u, th[0])
ph.endpoint = lambda xi, xf, p, t0, tf: tf
ph.events = lambda xi, xf, p, t0, tf, th: ca.vertcat(xi[0], xi[1], xi[2], xi[3],
                                                     xf[0], xf[1], xf[2], xf[3])

ph.bounds.lower.states = [-2.0] * 4
ph.bounds.upper.states = [2.0] * 4
ph.bounds.lower.controls = [-1.0, -1.0]
ph.bounds.upper.controls = [1.0, 1.0]
# The initial states are pinned, which every scenario shares and so can be. The
# terminal states are held in a ball of radius DELTA, which is what makes the
# problem well posed for more than one scenario: one open-loop torque history
# cannot steer several different plants to the same point exactly.
ph.bounds.lower.events = list(X0) + list(XF - DELTA)
ph.bounds.upper.events = list(X0) + list(XF + DELTA)
ph.bounds.t0 = (0.0, 0.0)
ph.bounds.tf = (1.0, 15.0)

ph.guess.states = np.tile(X0.reshape(4, 1), (1, NODES))
ph.guess.controls = np.zeros((2, NODES))
ph.guess.time = np.linspace(0.0, 3.0, NODES).reshape(1, NODES)

rp.uncertainty = Gaussian(mean=[MU], cov=[[SIGMA ** 2]], truncate=3.0)
rp.initial_state = X0

alg = psopt.Algorithm(transcription_method="multiple-shooting",
                      ms_integrator="RK4", ms_steps_per_segment=12,
                      ms_control_parameterisation="linear",
                      scaling="automatic", nlp_tolerance=1.0e-6,
                      nlp_iter_max=1000, print_level=0)

print("Two-link arm with an uncertain payload, through psopt.robust")
# SLACK is not the terminal tolerance. DELTA is already declared in the event
# bounds above; slack is how far OUTSIDE that ball the design may still stray
# somewhere in the uncertainty set. Passing DELTA here would quietly ask for a
# ball of radius 2*DELTA -- which is what an earlier version of this example did,
# and what the cross-check at the foot of the file was written to catch.
SLACK = 1.0e-3

out = rp.solve(alg, slack=SLACK, scenarios="sigma-points",
               generate=True, max_iterations=12, tighten=0.9,
               out_of_sample=400, wait_and_see=2)

print("\n  starting scenarios were the unscented set %s"
      % np.array2string(rp.uncertainty.sigma_points()[0].ravel(), precision=4))
print("  three-point Gauss-Hermite gives  [%.4f %.4f %.4f]"
      % (MU - np.sqrt(3) * SIGMA, MU, MU + np.sqrt(3) * SIGMA))

# ---- the nominal design, for the comparison that makes the numbers mean anything -
# One scenario at the mean is the ordinary deterministic design. It is solved here
# through the same driver so that nothing but the scenario set differs.
from psopt.robust import Explicit                                    # noqa: E402

rp.uncertainty = Explicit([[MU]])
nom = rp.solve(alg, slack=SLACK, scenarios="explicit", generate=False,
               out_of_sample=0, verbose=False)
rp.uncertainty = Gaussian(mean=[MU], cov=[[SIGMA ** 2]], truncate=3.0)
nom_worst, nom_theta = rp._worst_case((nom.time, nom.controls), np.zeros(0),
                                      128, 3, 1)

print("\n  %-22s %9s %13s %11s" % ("design", "t_f", "worst in set", "at m_p"))
print("  %-22s %9.4f %13.3e %11.4f"
      % ("nominal", nom.objective, nom_worst, nom_theta[0]))
print("  %-22s %9.4f %13.3e %11.4f"
      % ("robust", out.objective, out.certificate["violation"],
         out.certificate["theta"][0]))

# ---- the driver's verifier, against an independent implementation ----------------
#
# The driver builds its verification integrator out of the user's CasADi
# equations. This checks it against the hand-written NumPy one in robust_arm.py,
# which shares nothing with it but the physics. The two measure different things
# by construction -- the driver reports how far the terminal state falls OUTSIDE
# the declared ball, the hand-written one reports the miss itself -- so the
# identity being tested is
#
#     driver violation  ==  max(0, miss - DELTA)
#
# and it is written out rather than eyeballed, because the first version of this
# example got exactly that relation wrong and read a design that missed by 0.053
# as one that had met a tolerance of 0.03.
def miss_numpy(t, u, mp_grid, nsub=32):
    """|x(t_f) - XF| at each payload, by an independent fixed-step RK4."""
    k = len(mp_grid)
    x = np.tile(X0.reshape(4, 1), (1, k)).astype(float)

    def f(xx, uu):
        x1, x2, x3 = xx[0], xx[1], xx[2]
        m11, m22, s2 = 7.0 / 3.0 + mp_grid, 4.0 / 3.0 + mp_grid, 1.5 + mp_grid
        m12 = s2 * np.cos(x3)
        det = m11 * m22 - m12 * m12
        b1 = (uu[0] - uu[1]) + s2 * np.sin(x3) * x2 ** 2
        b2 = uu[1] - s2 * np.sin(x3) * x1 ** 2
        return np.vstack([(m22 * b1 - m12 * b2) / det,
                          (-m12 * b1 + m11 * b2) / det, x2 - x1, x1])

    for i in range(len(t) - 1):
        h = (t[i + 1] - t[i]) / nsub
        for sstep in range(nsub):
            uA = u[:, i] + (sstep / nsub) * (u[:, i + 1] - u[:, i])
            uH = u[:, i] + ((sstep + 0.5) / nsub) * (u[:, i + 1] - u[:, i])
            uB = u[:, i] + ((sstep + 1.0) / nsub) * (u[:, i + 1] - u[:, i])
            k1 = f(x, uA)
            k2 = f(x + 0.5 * h * k1, uH)
            k3 = f(x + 0.5 * h * k2, uH)
            k4 = f(x + h * k3, uB)
            x = x + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return np.max(np.abs(x - XF.reshape(4, 1)), axis=0)


grid = np.linspace(MU - 3 * SIGMA, MU + 3 * SIGMA, 2001)
mine = miss_numpy(out.time, out.controls, grid)
theirs = rp._violation_many(grid.reshape(-1, 1), (out.time, out.controls), np.zeros(0))
gap = float(np.max(np.abs(np.maximum(mine - DELTA, 0.0) - theirs)))
print("\n  the driver's verifier against an independent NumPy RK4, over %d payloads:"
      % len(grid))
print("    largest disagreement in max(0, miss - delta) : %.2e" % gap)
print("    worst terminal miss anywhere in the set      : %.4f  (the ball is %.2f)"
      % (mine.max(), DELTA))

# A test rather than a demonstration: the driver must certify the robust design
# over the whole set, must bring essentially every sampled plant inside the
# declared ball, the two verifiers must agree, and the nominal design must fail.
ok = (out.converged
      and out.certificate["violation"] <= SLACK
      and out.out_of_sample["within"] >= 95.0
      and mine.max() <= DELTA + 10.0 * SLACK
      and gap < 1.0e-3
      and nom_worst > 3.0 * DELTA)
print("\nobjective      : %.12g" % out.objective)
print("RESULT: %s" % ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
