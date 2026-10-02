"""A robust optimal control problem, posed by scenario augmentation.

The arm is Weinreb and Bryson (1985), Section IV, whose equations carry the tip
mass symbolically as the ratio mu = M/m; Luus (2000), Section 12.4.2, is the mu = 1
case with the coefficients evaluated, which is what ``examples/twolinkarm`` solves.
The payload used here is mu = 1 + m_p. That identity is checked and not assumed:
substituting it into their equations reproduces this model exactly.

The two-link arm of ``examples/twolinkarm`` carries a payload whose mass is not
known exactly. One torque history and one final time must be committed BEFORE the
payload is revealed, and must bring every plant in the uncertainty set to the
target within a tolerance. That is a here-and-now decision, and it is what
separates robust optimal control from solving the problem once per sample:

    Solving one optimal control problem per sample and averaging the answers does
    NOT give a robust control. It gives the wait-and-see solution, in which every
    realisation is optimised with foreknowledge of its own uncertainty. Since
    E[min] <= min E[.], its cost is a lower bound no implementable control
    attains, and the average of the control histories solves nothing at all. What
    couples the samples into one problem is that they SHARE the control --
    non-anticipativity -- and here that is structural: one control vector and M
    copies of the state.

So the transcription is ONE deterministic problem with ``4*M`` states. Nothing in
PSOPT has to change, and in Python the scenario count is simply a variable that
the dynamics, the events and the bounds are built from -- which is the part of
this that the interface makes short.

Two ways of choosing the scenarios are compared, and they do not do equally well:

  * the Gauss-Hermite nodes of the payload distribution, which is quadrature and
    is the right tool for an EXPECTATION. Here the constraint has to hold
    everywhere, and quadrature pins the terminal miss at its own nodes while
    saying nothing about the gaps between them.

  * scenarios GENERATED where the current design is failing: solve, find the
    payload the resulting control serves worst over the whole uncertainty set,
    add it, re-solve. That is a cutting-plane method on a semi-infinite
    constraint, it needs one call to ``prob.solve`` per iteration, and it ends
    with a certificate rather than a statistic.

The verification integrator here is deliberately NOT PSOPT's: it is a few lines
of NumPy, vectorised over payloads, sharing only the equations of motion. A
design checked against a different implementation of the plant is checking the
implementation and not the design.

``examples/robust_driver`` is the fuller C++ study of the same problem, through
the library's own driver: a sweep in the quadrature order, out-of-sample
statistics for every design, and the measurement showing that the integrator
accuracy validated on the nominal slew is not enough for the robust one. This is
the short version, cut down to run in under a minute, and it is the one that
writes the loop out BY HAND: the augmented problem, the verification integrator
and the generation loop are all here, which is what to read if you want to see
what a driver does rather than use one.
"""
import numpy as np
import casadi as ca
from _common import psopt

MU, SIGMA, DELTA = 0.50, 0.15, 0.03     # payload mean, s.d., terminal tolerance
LO, HI = MU - 3.0 * SIGMA, MU + 3.0 * SIGMA
X0 = np.array([0.0, 0.0, 0.500, 0.000])
XF = np.array([0.0, 0.0, 0.500, 0.522])
NODES, STEPS = 21, 12

# ---------------------------------------------------------------------------------
# The arm. The payload enters only through the inertia: it adds m_p to the two
# diagonal moments and to the first moment of the second link, with unit lengths.
# ---------------------------------------------------------------------------------


def arm_rhs(x1, x2, x3, u1, u2, mp, sin, cos):
    """Equations of motion, written once for both CasADi and NumPy."""
    m11, m22, s2 = 7.0 / 3.0 + mp, 4.0 / 3.0 + mp, 1.5 + mp
    m12 = s2 * cos(x3)
    det = m11 * m22 - m12 * m12
    b1 = (u1 - u2) + s2 * sin(x3) * x2 ** 2
    b2 = u2 - s2 * sin(x3) * x1 ** 2
    return ((m22 * b1 - m12 * b2) / det,
            (-m12 * b1 + m11 * b2) / det,
            x2 - x1,
            x1)


def build(payloads, delta, guess=None):
    """The augmented problem: one copy of the arm per payload, one control."""
    m = len(payloads)
    prob = psopt.Problem(name="robust_arm")
    ph = prob.add_phase(nstates=4 * m, ncontrols=2, nevents=8 * m)
    ph.nodes = [NODES]

    def dynamics(x, u, p, t):
        return ca.vertcat(*[c for k, mp in enumerate(payloads)
                            for c in arm_rhs(x[4 * k], x[4 * k + 1], x[4 * k + 2],
                                             u[0], u[1], mp, ca.sin, ca.cos)])

    # The initial states, known exactly, then the terminal states, which are held
    # in a ball of radius delta rather than pinned: one open-loop torque history
    # cannot steer several different plants to the same point exactly, so pinning
    # them is infeasible for m > 1 and the tolerance is the modelling choice that
    # makes the problem well posed.
    ph.dynamics = dynamics
    ph.events = lambda xi, xf, p, t0, tf: ca.vertcat(
        *[xi[4 * k + j] for k in range(m) for j in range(4)],
        *[xf[4 * k + j] for k in range(m) for j in range(4)])
    ph.endpoint = lambda xi, xf, p, t0, tf: tf

    ph.bounds.lower.states = [-2.0] * (4 * m)
    ph.bounds.upper.states = [2.0] * (4 * m)
    ph.bounds.lower.controls = [-1.0, -1.0]
    ph.bounds.upper.controls = [1.0, 1.0]
    ph.bounds.lower.events = list(np.tile(X0, m)) + list(np.tile(XF - delta, m))
    ph.bounds.upper.events = list(np.tile(X0, m)) + list(np.tile(XF + delta, m))
    ph.bounds.t0 = (0.0, 0.0)
    # The robust slew is a slow one -- see the closing note -- so the bound of 10
    # that the deterministic problem uses is not enough. It must not be generous
    # either: a free final time with a wide bound lets the solver escape into a
    # distant local minimum on the harder scenario sets.
    ph.bounds.tf = (1.0, 15.0)

    if guess is None:
        ph.guess.states = np.tile(X0, m).reshape(4 * m, 1) * np.ones((1, NODES))
        ph.guess.controls = np.zeros((2, NODES))
        ph.guess.time = np.linspace(0.0, 3.0, NODES).reshape(1, NODES)
    else:
        # Dynamically consistent: each scenario's own plant, integrated through
        # the previous design's control. A constant state guess leaves the
        # augmented problem with m copies of an infeasible arc, and the solver
        # then wanders off to a distant local minimum.
        t, u = guess
        ph.guess.states = np.vstack([simulate(t, u, np.array([mp]), 4)[0][:, :, 0]
                                     for mp in payloads])
        ph.guess.controls = u
        ph.guess.time = t.reshape(1, NODES)
    return prob


def solve(payloads, delta, guess=None):
    prob = build(payloads, delta, guess)
    sol = prob.solve(psopt.Algorithm(
        transcription_method="multiple-shooting",
        ms_integrator="RK4",
        ms_steps_per_segment=STEPS,
        ms_control_parameterisation="linear",
        scaling="automatic", nlp_tolerance=1.0e-6, nlp_iter_max=1000,
        print_level=0))
    if not sol.status.success:
        return None
    return sol.time.reshape(NODES), sol.controls, sol.status


# ---------------------------------------------------------------------------------
# The verifier: fixed-step RK4 through every payload at once. The control returned
# by multiple shooting with a linear parameterisation IS piecewise linear through
# the node table, so interpolating linearly inside the interval being integrated
# reproduces the designed control exactly.
# ---------------------------------------------------------------------------------


def simulate(t, u, mp, nsub):
    """States at the nodes for each payload; returns (states, terminal miss)."""
    k = len(mp)
    x = np.tile(X0.reshape(4, 1), (1, k)).astype(float)
    hist = np.empty((4, len(t), k))
    hist[:, 0, :] = x

    def f(x, uu):
        d = arm_rhs(x[0], x[1], x[2], uu[0], uu[1], mp, np.sin, np.cos)
        return np.vstack(d)

    for i in range(len(t) - 1):
        h = (t[i + 1] - t[i]) / nsub
        for s in range(nsub):
            uA = u[:, i] + (s / nsub) * (u[:, i + 1] - u[:, i])
            uH = u[:, i] + ((s + 0.5) / nsub) * (u[:, i + 1] - u[:, i])
            uB = u[:, i] + ((s + 1.0) / nsub) * (u[:, i + 1] - u[:, i])
            k1 = f(x, uA)
            k2 = f(x + 0.5 * h * k1, uH)
            k3 = f(x + 0.5 * h * k2, uH)
            k4 = f(x + h * k3, uB)
            x = x + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        hist[:, i + 1, :] = x
    return hist, np.max(np.abs(x - XF.reshape(4, 1)), axis=0)


def worst_case(design, ngrid=241):
    """The worst payload in the set for a given control, by exhaustive scan.

    The uncertainty is one scalar here, so a scan is both exhaustive and cheap,
    and it is honest in a way a gradient search on a non-concave function would
    not be. For a higher-dimensional uncertainty this is the step that would have
    to become an optimisation, and it is what makes the scheme bilevel.
    """
    grid = np.linspace(LO, HI, ngrid)
    miss = simulate(design[0], design[1], grid, 8)[1]
    j = int(np.argmax(miss))
    return miss[j], grid[j]


def score(design, sample):
    miss = simulate(design[0], design[1], sample, 8)[1]
    return miss.mean(), miss.max(), 100.0 * np.mean(miss <= DELTA)


# ---------------------------------------------------------------------------------

print("Two-link arm with an uncertain payload")
print("  m_p ~ N(mu = %.3f, sigma = %.3f), terminal tolerance delta = %.3f" %
      (MU, SIGMA, DELTA))
print("  uncertainty set for verification: [%.3f, %.3f]" % (LO, HI))

rng = np.random.default_rng(20260927)
sample = rng.normal(MU, SIGMA, 2000)
sample = sample[(sample >= LO) & (sample <= HI)]     # the set the design was given

# ---- the nominal design, which robustness has to be measured against -------------
nominal = solve([MU], DELTA)
assert nominal is not None, "the nominal design failed"
print("\n  nominal design (the mean payload alone): t_f = %.6f" % nominal[0][-1])

# ---- quadrature scenarios --------------------------------------------------------
# Gauss-Hermite nodes for the standard normal, mapped onto the payload. Five of
# them: enough that the scenario set spans the interesting part of the set, and
# few enough that the failure below is not about the arithmetic.
z = np.polynomial.hermite_e.hermegauss(5)[0]
gh = solve(list(MU + SIGMA * z), DELTA, guess=nominal[:2])
assert gh is not None, "the Gauss-Hermite design failed"

# ---- generated scenarios ---------------------------------------------------------
# The design is given a tightened tolerance and verified against the full one. The
# scenarios are satisfied to delta exactly -- an active constraint is active -- so
# the miss BETWEEN two scenarios is necessarily a little larger than delta, and a
# loop that demands delta over the whole set could not otherwise terminate.
print("\n  generating scenarios where the design is failing")
print("  %2s %3s %10s %13s %9s %10s" % ("it", "M", "t_f", "worst in set",
                                        "at m_p", "within d"))
# Seeded with the two ends of the uncertainty set as well as its middle. The
# extremes are almost always among the active scenarios when the dependence on the
# uncertain parameter is monotone, as it is here, and seeding them saves the loop
# from discovering it one iteration at a time.
payloads = [LO, MU, HI]
design = solve(payloads, 0.9 * DELTA, guess=nominal[:2]) or nominal
for it in range(12):
    wc, at = worst_case(design)
    print("  %2d %3d %10.6f %13.3e %9.4f %9.1f%%" %
          (it, len(payloads), design[0][-1], wc, at, score(design, sample)[2]))
    if wc <= DELTA:
        print("\n  converged: no payload in the set misses by more than delta")
        break
    payloads.append(at)
    nxt = solve(payloads, 0.9 * DELTA, guess=design[:2])
    if nxt is None:
        print("  adding m_p = %.4f made the problem unsolvable" % at)
        break
    design = nxt

# ---- what it bought and what it cost ---------------------------------------------
print("\n  over %d sampled payloads inside the set\n" % len(sample))
print("  %-27s %3s %8s %11s %10s %9s" %
      ("design", "M", "t_f", "worst set", "mean miss", "within d"))
rows = [("nominal (mean payload only)", 1, nominal),
        ("Gauss-Hermite nodes", 5, gh),
        ("generated scenarios", len(payloads), design)]
stats = {}
for name, m, d in rows:
    mean, _, frac = score(d, sample)
    stats[name] = frac
    print("  %-27s %3d %8.4f %11.3e %10.3e %8.1f%%" %
          (name, m, d[0][-1], worst_case(d)[0], mean, frac))

print("\n  price of robustness: t_f %.4f -> %.4f, a factor of %.2f. The robust"
      % (nominal[0][-1], design[0][-1], design[0][-1] / nominal[0][-1]))
print("  slew is slow because the payload enters through the inertia: drive hard")
print("  and it matters, drive gently and the trajectory approaches a quasi-static")
print("  one on which it matters much less. That is where the insensitivity is")
print("  bought, and why the bound on t_f had to be raised above the nominal one.")

# A test and not a demonstration: the nominal design must fail out of sample and
# the generated one must not.
ok = (stats["nominal (mean payload only)"] < 25.0 and
      stats["generated scenarios"] >= 95.0)
print("\nobjective      : %.12g" % design[0][-1])
print("RESULT: %s" % ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
