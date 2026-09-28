"""Who should choose the ancillary gain, and a measured answer: not the optimiser.

`robust_driver_tube.py` shows that a tube --- a nominal trajectory plus a feedback
law pulling each realisation back towards it --- recovers almost all of what an
open-loop robust design pays: on this arm the penalty over the deterministic design
falls from +192% to +13%. It uses the simplest possible gain, an LQR regulator for
the plant linearised about its TARGET, held constant through the manoeuvre. That
gain looks like a placeholder. The arm spends almost all of the manoeuvre nowhere
near the target, the gain is regulating about a point the plant is not at, and there
is +13% still on the table. Handing the gain to the optimiser --- CO-DESIGN --- is
the obvious next move.

It does not work, and the reason is structural rather than numerical. This example
is that measurement.

WHAT THE DRIVER OFFERS

    .feedback = K              a matrix: a given constant gain
    .feedback = lambda t: ...  a given SCHEDULE, K(t)
    .feedback = "co-design"    a constant gain, its entries optimised as static
                               parameters alongside the nominal trajectory
    .feedback = "co-design-schedule"
                               a TIME-VARYING gain, optimised as extra controls

with `.feedback_bounds = (lo, hi)` for the two co-designed forms, where lo and hi
are scalars or (ncontrols, nstates) arrays, and `.feedback_guess` for their starting
point.

WHY CO-DESIGN IS NOT HARD TO SET UP

The usual objection is that a gain multiplying a state deviation makes the problem
bilinear in the decision variables. True and irrelevant: the problem was already
nonconvex, the dynamics are already nonlinear in the state, and a product of two
decision variables is just another nonlinear term. PSOPT has static parameters, so a
constant gain goes in as eight more of them. A gain that varies in time goes in as
eight more CONTROLS, which is the representation to use: the transcription already
gives a control a time profile at exactly the resolution the trajectory has, with no
basis to choose and no interpolation to tape. (A polynomial in t was tried first.
It is worse in every way: the Riccati gain along this manoeuvre needs degree 9 to
fit to 3%, degree 3 misses by a factor of three, and at degree 9 the first solve
does not converge at all.)

So setting co-design up is easy. That is the trap.

WHY IT FAILS

A scenario set does two jobs at once, and for the gain they pull in opposite
directions. As a CONSTRAINT SET it is a list of plants the design must serve. As a
sample it is meant to stand for the whole family. The nominal control can be chosen
against it because one control history must serve every member as it is, and its
effect on an unsampled plant is close to its effect on a sampled one nearby.

A gain is not in that position. It multiplies a deviation which is itself a function
of theta, so its leverage on a plant the design never saw is bounded by nothing the
design can see. The optimiser is paid in final time for a large gain and charged
nothing for what that gain does off-sample.

The measurements, from this file and from the sweeps behind it. Slack is 1e-3
throughout, "within" is the fraction of 150 out-of-sample draws inside it, and the
scenario counts are the sets the designs were built on:

    given LQR gain, generated               t_f 3.558   worst 0.0e+00   100% within
    co-design, |K| <= 25,           3 pts   t_f 3.307   worst 2.1e+02     0% within
    co-design-schedule, |K| <= 10,  3 pts   t_f 3.300   worst 6.2e+01     0% within
    co-design-schedule, |K| <= 10,  9 pts   t_f 3.262   worst 3.3e-01     6% within
    co-design-schedule, |K| <= 10, 16 pts   t_f 3.396   worst 6.8e+01     0% within
    co-design, box around LQR,   generated  t_f 3.438   worst 2.4e-02    74% within
    co-design-schedule, box,     generated  t_f 3.418   worst 4.8e-02    46% within

Every co-designed design is cheaper on the objective than the design with the given
gain, and every one fails its certificate --- by between one and five orders of
magnitude. Two things the table rules out:

IT IS NOT THE SCENARIO SET BEING TOO THIN. Nine points are two orders of magnitude
better than three and sixteen are no better than three, which is not a trend but the
luck of which local minimum IPOPT found. Density is also what a scenario method pays
for: sixteen points cost 21 minutes.

IT IS NOT THE BOUND BEING TOO LOOSE. With a wide bound the gain does run to it: at
|K| <= 25 the gain comes out as

    [ -6.49   25.00   -3.49  -21.31 ]
    [ 21.91   25.00   -8.36  -10.19 ]

against an LQR gain whose largest entry is 4.38 --- two entries exactly on the bound,
two more within 15% of it --- and that design is 200 000 times over slack. But the
boxed SCHEDULE, confined to within half its own magnitude of a gain
that certifies, finishes with only 11 of its 200 values on the bound and a schedule
that rises smoothly from |K| 2.7 to 5.2 over the manoeuvre --- a thoroughly plausible
gain, not a bound artefact --- and it still misses by a factor of 48. The defect is
in what the objective measures, not in how much gain it is allowed.

The gain that works is chosen by a criterion that never looks at the scenario set at
all: stabilise the linearisation, penalise state and effort quadratically, solve the
Riccati equation. That criterion quantifies over the whole family implicitly, which
is the property the gain needs and an objective on a finite sample cannot supply.

The gain that works is chosen by a criterion that never looks at the scenario set at
all: stabilise the linearisation, penalise state and effort quadratically, solve the
Riccati equation. That criterion quantifies over the whole family implicitly, which
WHAT THE DRIVER DOES ABOUT IT

It offers co-design and then measures it honestly, which is the only defensible
posture for a feature that is easy to reach for and usually wrong:

  * the verification integrates the CLOSED loop with the gain that was DESIGNED, so
    a co-designed gain is checked as the controller it is. This example found a bug
    in which it was not: the solved gain was not read back, and every co-designed
    design was being verified against the LQR gain it had started from. Every one of
    them certified. The numbers above are the ones that appeared once it did not.
  * `certificate["gain_on_bound"]` counts the entries sitting on their bound, and
    the report says in as many words what that means.
  * `certificate["control_excess"]` measures what the gain asks of the actuator
    between the nodes, where the design constrains nothing.
  * `certificate["integrator_drift"]` measures the search integrator against the
    adaptive one, and on the co-designed design it comes out as large as the
    violation: an aggressive gain makes the closed loop stiff, and a step count
    adequate for the nominal problem is not adequate for that. Reported, not hidden.

WHAT THIS FILE RUNS

The two designs that make the point, at the same settings, checked the same way: the
given LQR gain, and a co-design confined to a box around it --- the best-behaved
co-design of the six, so the comparison is against co-design at its most favourable.
The test at the end asserts the finding rather than a value: the given gain must
certify, and the co-designed one must both undercut it on the objective and fail to
certify. If co-design ever starts certifying here, this file should fail and be
rewritten, because the conclusion above will have changed.
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


rp = RobustProblem(name="robust_driver_gain")
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

# The given gain: LQR for the plant linearised about its target, with the nominal
# payload. Nothing about it knows that the payload is uncertain, and that is the
# point -- it is chosen for a property of the whole family, not of a sample.
Q, R = np.diag([1.0, 1.0, 10.0, 10.0]), np.eye(2)
xsym, usym = ca.SX.sym("x", 4), ca.SX.sym("u", 2)
fsym = arm_rhs(xsym, usym, MU)
A = np.array(ca.Function("A", [xsym, usym], [ca.jacobian(fsym, xsym)])(XF, np.zeros(2)))
B = np.array(ca.Function("B", [xsym, usym], [ca.jacobian(fsym, usym)])(XF, np.zeros(2)))
FIXED = -np.linalg.solve(R, B.T @ solve_continuous_are(A, B, Q, R))

# The box the co-design gets: each entry within half its own magnitude, plus half a
# unit so that the entries near zero have somewhere to go. This is co-design at its
# most favourable -- it starts from a gain that certifies and cannot leave its
# neighbourhood -- and it is the variant that comes closest to working.
HALF = 0.5 * np.abs(FIXED) + 0.5
BOX = (FIXED - HALF, FIXED + HALF)

alg = psopt.Algorithm(transcription_method="multiple-shooting",
                      ms_integrator="RK4", ms_steps_per_segment=12,
                      ms_control_parameterisation="linear",
                      ms_path_samples=4,
                      scaling="automatic", nlp_tolerance=1.0e-6,
                      nlp_iter_max=3000, print_level=0)


def run(gain, label, bounds=None, guess=None, max_iterations=4):
    ph.feedback = gain
    ph.feedback_bounds = bounds
    ph.feedback_guess = guess
    out = rp.solve(alg, slack=SLACK, scenarios="sigma-points", generate=True,
                   max_iterations=max_iterations, tighten=0.9, out_of_sample=150,
                   wait_and_see=0, verbose=False)
    c, o = out.certificate, out.out_of_sample
    print("  %-22s t_f %7.4f  M %2d  worst %8.2e  within %5.1f%%  "
          "u beyond %8.2e  %s"
          % (label, out.objective, len(out.scenarios), c["violation"],
             o["within"], c["control_excess"],
             "certifies" if out.converged else "FAILS"))
    return out


print("Who should choose the ancillary gain, on the arm with an uncertain payload")

# The deterministic design, the floor the price of robustness is measured against. A
# one-point uncertainty, so nothing is being made robust.
ph.feedback = None
rp.uncertainty = Gaussian(mean=[MU], cov=[[1.0e-12]], truncate=3.0)
nom = rp.solve(alg, slack=SLACK, scenarios="sigma-points", generate=False,
               out_of_sample=0, wait_and_see=0, verbose=False)
rp.uncertainty = Gaussian(mean=[MU], cov=[[SIGMA ** 2]], truncate=3.0)
print("\n  deterministic design (no uncertainty): t_f = %.4f" % nom.objective)

print("\n  %-22s %11s %4s %10s %8s %9s"
      % ("gain", "t_f", "M", "worst", "within", "u beyond"))
given = run(FIXED, "given, LQR")
# Three iterations, not ten. The generation loop is there to close the gap between a
# design and the worst parameter it did not see, and it cannot close a gap the gain
# is reopening: this co-design reaches exactly the same five-scenario design at three
# iterations, at four and at ten (measured), and stopping at three saves eleven
# minutes of a quarter-hour example.
codes = run("co-design", "co-designed in a box", bounds=BOX, guess=FIXED,
            max_iterations=3)

G = np.asarray(codes.design.parameters).ravel()[:8].reshape(2, 4)
print("\n  the co-designed gain, against the gain it started from:")
for r in range(2):
    print("    co-designed [%s ]" % " ".join("%7.3f" % G[r, c] for c in range(4)))
    print("    given       [%s ]" % " ".join("%7.3f" % FIXED[r, c] for c in range(4)))
nb = codes.certificate.get("gain_on_bound", 0)
print("  %d of its %d entries finished on the edge of the box."
      % (nb, codes.certificate.get("gain_entries", 8)))

print("\n  price of robustness against the deterministic design:")
for label, out in (("given", given), ("co-designed", codes)):
    print("    %-12s %+6.1f%%   worst violation %8.2e against a slack of %.0e  %s"
          % (label, 100.0 * (out.objective - nom.objective) / nom.objective,
             out.certificate["violation"], SLACK,
             "certifies" if out.converged else "does not certify"))

print("\n  The co-designed design is not only uncertifiable, it is hard to verify.")
print("  The search integrator disagrees with the adaptive one by %.2e on it,"
      % codes.certificate["integrator_drift"])
print("  against %.2e on the given gain's design -- so on the co-designed one the"
      % given.certificate["integrator_drift"])
print("  disagreement is as large as the violation being reported. An aggressive")
print("  gain makes the closed loop stiff, and a step count adequate for the")
print("  nominal problem is not adequate for that. The driver reports the")
print("  disagreement rather than letting the certificate rest on it silently.")

print("\n  The co-designed gain is cheaper and does not work. It is cheaper BECAUSE")
print("  it does not work: the objective is evaluated on the scenarios, the gain's")
print("  damage is done off them, and nothing in the design charges it. A gain has")
print("  to be chosen for a property of the whole family -- which is what solving")
print("  a Riccati equation does and what minimising over a sample cannot.")

# A test of the finding, not of a value. The given gain must certify; the co-designed
# one must undercut it on the objective and fail. Were co-design ever to certify here
# this file would fail, which is the intended behaviour: the conclusion it documents
# would have changed and it should be rewritten.
ok = (given.converged
      and given.certificate["violation"] <= SLACK
      and given.certificate["control_excess"] < 5.0e-3
      and given.objective < 2.0 * nom.objective
      and codes.objective < given.objective          # cheaper
      and not codes.converged                        # and wrong
      and codes.certificate["violation"] > 5.0 * SLACK
      and nb > 0)                                    # by taking all the gain allowed
print("\nobjective      : %.12g" % given.objective)
print("RESULT: %s" % ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
