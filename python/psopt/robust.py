"""Robust optimal control by scenario augmentation.

This module turns a problem whose dynamics, path constraints and events depend on
an uncertain parameter vector into ONE deterministic problem that PSOPT solves as
it stands, and wraps an outer loop around it. Nothing in the library is modified;
the driver assembles, calls ``Problem.solve``, verifies, and calls it again.

    import psopt
    from psopt.robust import RobustProblem, Gaussian

    rp = RobustProblem(name="arm")
    ph = rp.add_phase(nstates=4, ncontrols=2, nevents=8)
    ph.dynamics = lambda x, u, p, t, th: ...      # th is the uncertain parameter
    ph.events   = lambda xi, xf, p, t0, tf, th: ...
    ...
    rp.uncertainty   = Gaussian(mean=[0.5], cov=[[0.15 ** 2]], truncate=3.0)
    rp.initial_state = [0.0, 0.0, 0.5, 0.0]
    out = rp.solve(alg, slack=0.02)

    out.certificate     # worst violation found anywhere in the uncertainty set
    out.out_of_sample   # scored on payloads never used to design

WHAT IS BEING SOLVED, AND WHAT IS NOT

The control is chosen BEFORE the uncertainty is revealed and must serve every
realisation --- a here-and-now decision. Solving the problem once per sample and
averaging is a different computation: it gives the wait-and-see solution, whose
cost is a lower bound that no implementable control attains, because every
realisation was optimised with foreknowledge of its own uncertainty. Since
E[min] <= min E[.], that bound is genuinely useful as a measure of the value of
information, and ``wait_and_see=`` computes it, but it is a diagnostic and not a
design.

What couples the scenarios is that they share the control. In the augmented
problem that is structural: M copies of the state, one copy of the control.

THREE THINGS THIS DRIVER TAKES SERIOUSLY

*The verification is independent of the design.* The design is transcribed by
PSOPT and integrated by its segment integrator; the verification builds a CasADi
function from the same user equations and integrates it with SciPy's DOP853 at
tight tolerances, which is a different method, a different step control and a
different code path. Checking a design against the integrator that produced it
checks nothing.

*In-sample risk is not the answer.* The optimiser drives the constraint to its
bound at every scenario it is given, so the in-sample figure is the bound by
construction and carries no information. Every number this driver reports is
either measured over payloads that took no part in the design, or is a worst case
over the whole set.

*A certificate is a claim about a search.* For one uncertain parameter the worst
case can be found by exhaustive scan. In more than one it cannot, and the oracle
here is quasi-random seeding followed by local refinement. The result object
records how many evaluations went into it, and ``RobustSolution.certificate``
says "nothing worse was found" rather than "nothing worse exists". The inner
problem is not concave and no cheap method can say more.

WHAT IS NOT HERE YET

CVaR and other risk measures beyond expectation and mean-variance (the
Rockafellar-Uryasev device wants a static parameter per scenario, which is
expressible but not written); scenario-dependent static parameters; feedback of
any kind, so the designs are open-loop and therefore conservative over long
horizons. Multi-phase problems are refused rather than quietly mishandled.
"""
import warnings

import numpy as np
import casadi as ca

from . import Problem, Algorithm, _Bounds, _Guess

__all__ = ["RobustProblem", "RobustSolution", "Gaussian", "Uniform", "Explicit"]


# ----------------------------------------------------------------------------------
#  The uncertainty
# ----------------------------------------------------------------------------------
#
#  An uncertainty supplies three separate things, and they are separate on purpose:
#  a SET, over which the design must hold and over which the worst case is sought;
#  a DISTRIBUTION, from which out-of-sample payloads are drawn and against which a
#  quadrature rule is defined; and a scenario rule. Conflating the set with the
#  distribution is how a design ends up being scored on its own truncation.


class Uncertainty(object):
    """Base class. See Gaussian, Uniform and Explicit."""

    @property
    def dim(self):
        raise NotImplementedError

    def sigma_points(self):
        """Deterministic scenario set: (points, weights)."""
        raise NotImplementedError

    def quasi_random(self, m, seed=0):
        """m low-discrepancy points spread over the SET (not the distribution)."""
        raise NotImplementedError

    def sample(self, m, rng):
        """m draws from the DISTRIBUTION, untruncated."""
        raise NotImplementedError

    def contains(self, theta):
        raise NotImplementedError

    def box(self):
        """A bounding box of the set, as (lo, hi). Used to seed the oracle."""
        raise NotImplementedError

    def clip(self, theta):
        """The nearest point of the set, or theta itself if it is inside."""
        raise NotImplementedError

    def boundary_points(self):
        """Extreme points of the set, added to every worst-case search.

        A low-discrepancy sequence fills the interior well and reaches the
        boundary only by accident, and the worst parameter very often sits
        exactly on it -- a heavier payload than any considered, the strongest
        current. These are put in by hand so that the search cannot miss the
        boundary through bad luck in the seeding.
        """
        raise NotImplementedError

    def describe(self):
        raise NotImplementedError


class Gaussian(Uncertainty):
    """theta ~ N(mean, cov), with the SET the ellipsoid at `truncate` sigma.

    The set is {theta : (theta-mu)' inv(cov) (theta-mu) <= truncate^2}, which for
    one parameter is the familiar mu +- k sigma and in n dimensions covers rather
    less probability than the same k does in one. The coverage is the chi-square
    distribution with n degrees of freedom at k^2, which for k = 3 runs

        99.73% in one dimension, 98.89% in two, 97.07% in three,

    so a habit formed on scalar uncertainty truncates a three-parameter posterior
    at a set holding only 97% of it. The driver reports how much of the sample
    fell outside, because that is the part the design promises nothing about.
    """

    def __init__(self, mean, cov, truncate=3.0):
        self.mean = np.atleast_1d(np.asarray(mean, dtype=float))
        self.cov = np.atleast_2d(np.asarray(cov, dtype=float))
        self.truncate = float(truncate)
        if self.cov.shape != (self.dim, self.dim):
            raise ValueError("Gaussian: cov is %s, expected (%d, %d)"
                             % (self.cov.shape, self.dim, self.dim))
        # Cholesky of the covariance: the map z -> mean + L z sends the unit ball
        # to the set and the standard normal to the distribution, which is what
        # every method below is written in terms of.
        self._L = np.linalg.cholesky(self.cov)
        self._Linv = np.linalg.inv(self._L)

    @property
    def dim(self):
        return self.mean.size

    def sigma_points(self):
        """The unscented set: 2n+1 points matching the mean and the covariance.

        kappa is fixed at 3-n, the usual choice, which makes the fourth moment
        right for a scalar Gaussian. For n > 3 the central weight goes negative;
        that is admissible for a quadrature but not for a scenario set the design
        must satisfy, so the driver raises rather than pretending.
        """
        n = self.dim
        kappa = 3.0 - n
        if n + kappa <= 0:
            raise ValueError("Gaussian.sigma_points: n + kappa <= 0 for n = %d" % n)
        if kappa < 0:
            raise ValueError(
                "Gaussian.sigma_points: the central weight is negative for n = %d "
                "(kappa = 3 - n = %g). A negative weight is admissible for a "
                "quadrature but not for a scenario set every one of whose members "
                "is a constraint. Use scenarios='qmc' instead." % (n, kappa))
        c = np.sqrt(n + kappa)
        pts = [self.mean.copy()]
        wts = [kappa / (n + kappa)]
        for k in range(n):
            col = c * self._L[:, k]
            pts.append(self.mean + col)
            pts.append(self.mean - col)
            wts.extend([0.5 / (n + kappa), 0.5 / (n + kappa)])
        return np.array(pts), np.array(wts)

    def quasi_random(self, m, seed=0):
        from scipy.stats import qmc
        s = qmc.Sobol(d=self.dim, scramble=True, seed=seed)
        # Uniform in the unit BALL, mapped through L: inverse-CDF radius so the
        # points do not pile up at the centre in higher dimension. Sobol' only
        # balances on powers of two, so the count is rounded up to one.
        raw = s.random_base2(m=int(np.ceil(np.log2(max(2, m)))))
        z = _ball_from_unit_cube(raw)
        return self.mean + (self._L @ (self.truncate * z).T).T

    def sample(self, m, rng):
        return rng.multivariate_normal(self.mean, self.cov, size=m)

    def _radius(self, theta):
        z = self._Linv @ (np.asarray(theta, dtype=float) - self.mean)
        return float(np.linalg.norm(z))

    def contains(self, theta):
        return self._radius(theta) <= self.truncate + 1e-12

    def box(self):
        h = self.truncate * np.sqrt(np.diag(self.cov))
        return self.mean - h, self.mean + h

    def clip(self, theta):
        r = self._radius(theta)
        if r <= self.truncate:
            return np.asarray(theta, dtype=float)
        z = self._Linv @ (np.asarray(theta, dtype=float) - self.mean)
        return self.mean + self._L @ (z * (self.truncate / r))

    def boundary_points(self):
        pts = []
        for k in range(self.dim):
            col = self.truncate * self._L[:, k]
            pts.append(self.mean + col)
            pts.append(self.mean - col)
        return np.array(pts)

    def describe(self):
        sd = np.sqrt(np.diag(self.cov))
        return ("Gaussian, mean %s, s.d. %s, set = ellipsoid at %g sigma"
                % (_fmt(self.mean), _fmt(sd), self.truncate))


class Uniform(Uncertainty):
    """theta uniform on the box [lo, hi], which is also the set."""

    def __init__(self, lo, hi):
        self.lo = np.atleast_1d(np.asarray(lo, dtype=float))
        self.hi = np.atleast_1d(np.asarray(hi, dtype=float))
        if self.lo.shape != self.hi.shape:
            raise ValueError("Uniform: lo and hi have different shapes")
        if np.any(self.hi < self.lo):
            raise ValueError("Uniform: hi < lo in component %d"
                             % int(np.argmin(self.hi - self.lo)))

    @property
    def dim(self):
        return self.lo.size

    def sigma_points(self):
        """Centre plus the 2n face centres, at the variance-matching distance.

        For a uniform distribution the standard deviation of each component is
        (hi-lo)/sqrt(12), and the unscented spread sqrt(n+kappa) applied to it is
        what makes the set match the second moment rather than merely span the box.
        """
        n = self.dim
        kappa = 3.0 - n
        if kappa < 0:
            raise ValueError("Uniform.sigma_points: negative central weight for "
                             "n = %d; use scenarios='qmc'" % n)
        mean = 0.5 * (self.lo + self.hi)
        sd = (self.hi - self.lo) / np.sqrt(12.0)
        c = np.sqrt(n + kappa)
        pts = [mean.copy()]
        wts = [kappa / (n + kappa)]
        for k in range(n):
            step = np.zeros(n)
            step[k] = c * sd[k]
            pts.append(np.minimum(mean + step, self.hi))
            pts.append(np.maximum(mean - step, self.lo))
            wts.extend([0.5 / (n + kappa), 0.5 / (n + kappa)])
        return np.array(pts), np.array(wts)

    def quasi_random(self, m, seed=0):
        from scipy.stats import qmc
        s = qmc.Sobol(d=self.dim, scramble=True, seed=seed)
        raw = s.random_base2(m=int(np.ceil(np.log2(max(2, m)))))
        return self.lo + raw * (self.hi - self.lo)

    def sample(self, m, rng):
        return rng.uniform(self.lo, self.hi, size=(m, self.dim))

    def contains(self, theta):
        th = np.asarray(theta, dtype=float)
        return bool(np.all(th >= self.lo - 1e-12) and np.all(th <= self.hi + 1e-12))

    def box(self):
        return self.lo.copy(), self.hi.copy()

    def clip(self, theta):
        return np.clip(np.asarray(theta, dtype=float), self.lo, self.hi)

    def boundary_points(self):
        """Every corner of the box up to eight dimensions, face centres beyond."""
        n = self.dim
        if n <= 8:
            import itertools
            return np.array([[self.hi[k] if b else self.lo[k] for k, b in
                              enumerate(bits)]
                             for bits in itertools.product([0, 1], repeat=n)])
        mid = 0.5 * (self.lo + self.hi)
        pts = []
        for k in range(n):
            for v in (self.lo[k], self.hi[k]):
                q = mid.copy()
                q[k] = v
                pts.append(q)
        return np.array(pts)

    def describe(self):
        return "Uniform on the box %s to %s" % (_fmt(self.lo), _fmt(self.hi))


class Explicit(Uncertainty):
    """A finite list of parameter vectors, with weights. The set is the list.

    Use this when the uncertainty is genuinely discrete --- a handful of known
    operating modes, or a posterior already reduced to particles. The worst-case
    oracle degenerates to evaluating every member, which makes its certificate
    exhaustive rather than a search.
    """

    def __init__(self, points, weights=None):
        self.points = np.atleast_2d(np.asarray(points, dtype=float))
        if weights is None:
            self.weights = np.full(len(self.points), 1.0 / len(self.points))
        else:
            w = np.asarray(weights, dtype=float)
            self.weights = w / w.sum()

    @property
    def dim(self):
        return self.points.shape[1]

    def sigma_points(self):
        return self.points.copy(), self.weights.copy()

    def quasi_random(self, m, seed=0):
        return self.points.copy()

    def sample(self, m, rng):
        idx = rng.choice(len(self.points), size=m, p=self.weights)
        return self.points[idx]

    def contains(self, theta):
        return bool(np.any(np.all(np.isclose(self.points, theta), axis=1)))

    def box(self):
        return self.points.min(axis=0), self.points.max(axis=0)

    def clip(self, theta):
        d = np.linalg.norm(self.points - np.asarray(theta, dtype=float), axis=1)
        return self.points[int(np.argmin(d))]

    def boundary_points(self):
        return self.points.copy()

    def describe(self):
        return "Explicit, %d points" % len(self.points)


def _ball_from_unit_cube(raw):
    """Map points of the unit cube to the unit ball, uniformly by volume."""
    m, n = raw.shape
    # Direction from a normal deviate built out of the cube by the Box-Muller
    # transform on pairs, then normalised; radius by the inverse CDF u^(1/n).
    from scipy.special import ndtri
    g = ndtri(np.clip(raw, 1e-12, 1.0 - 1e-12))
    nrm = np.linalg.norm(g, axis=1, keepdims=True)
    nrm[nrm == 0.0] = 1.0
    direction = g / nrm
    # Radii stratified over (0, 1], with 1 attained: the worst parameter very
    # often sits ON the boundary of the set, and a seeding that never reaches it
    # would miss exactly the point the certificate is about.
    r = np.power(np.linspace(1.0 / m, 1.0, m), 1.0 / n)
    return direction * r[:, None]


def _fmt(v):
    return "[" + ", ".join("%.4g" % x for x in np.atleast_1d(v)) + "]"


# ----------------------------------------------------------------------------------
#  The problem statement
# ----------------------------------------------------------------------------------


class RobustPhase(object):
    """One phase of a robust problem.

    Identical to psopt.Phase except that dynamics, path, events and observation
    take the uncertain parameter vector as a final argument, and the bounds and
    the guess are stated at the NOMINAL size --- the driver replicates them.
    """

    def __init__(self, nstates, ncontrols, nevents=0, npath=0, nparameters=0):
        self.nstates = nstates
        self.ncontrols = ncontrols
        self.nevents = nevents
        self.npath = npath
        self.nparameters = nparameters
        self.nodes = [20]
        self.dynamics = None    # f(x, u, p, t, theta)
        self.path = None        # g(x, u, p, t, theta)
        self.integrand = None   # L(x, u, p, t)        -- shared, theta-free
        self.endpoint = None    # phi(xi, xf, p, t0, tf) -- shared, theta-free
        self.events = None      # e(xi, xf, p, t0, tf, theta)
        self.bounds = _Bounds()
        self.guess = _Guess()
        # Range the per-scenario cost is certainly inside. Required by the risk
        # measures that carry it as a state -- mean-variance and CVaR -- because a
        # state needs bounds, and a bound that turned out to be active would
        # silently change the risk measure into something else.
        self.cost_bounds = None
        # Ancillary feedback gain. One of:
        #   None                  open loop
        #   array (ncontrols, nstates)   a constant gain
        #   callable(t) -> rows   a gain scheduled on time, returned as nested
        #                         sequences of CasADi expressions
        #   "co-design"           the entries become static parameters and are
        #                         optimised with the trajectory
        # See the FEEDBACK section of RobustProblem._augment.
        self.feedback = None
        # (lo, hi) for a co-designed gain's entries. Each may be a scalar, or an
        # (ncontrols, nstates) array, which is how a co-design is confined to a
        # neighbourhood of a gain already known to stabilise the family.
        self.feedback_bounds = None
        self.feedback_guess = None     # a gain to start a co-design from


class RobustSolution(object):
    """What a robust solve produced, and how much of it can be believed."""

    def __init__(self):
        self.design = None          # the psopt Solution of the final augmented solve
        self.scenarios = None       # the scenario set it was built on
        self.time = None            # 1 x N node times of the design
        self.controls = None        # ncontrols x N control table
        self.objective = None
        self.certificate = None     # dict: worst violation found over the set
        self.out_of_sample = None   # dict: scored on samples that took no part
        self.wait_and_see = None    # lower bound on the objective, or None
        self.history = []           # one row per outer iteration
        self.n_solves = 0
        self.n_verifications = 0
        self.converged = False
        # True when the loop stopped because the transcription could not carry
        # another scenario, rather than because the design was good enough.
        self.budget_exhausted = False
        # A co-designed gain schedule, (ncontrols*nstates, N), or None.
        self.gain_schedule = None

    def report(self, printer=print):
        """Print what was found, in the order it should be read."""
        c = self.certificate
        printer("  scenarios in the final design : %d" % len(self.scenarios))
        printer("  objective                     : %.6f" % self.objective)
        printer("  worst violation over the set  : %.3e  at theta = %s"
                % (c["violation"], _fmt(c["theta"])))
        printer("  found over                    : %d evaluations of an independent"
                % c["evaluations"])
        printer("                                  integrator; nothing worse was")
        printer("                                  found, which is not a proof that")
        printer("                                  nothing worse exists")
        if c.get("feedback"):
            printer("  realised control beyond bounds: %.2e  (the ancillary gain's"
                    % c["control_excess"])
            printer("                                  demand on the actuator,")
            printer("                                  sampled between the nodes")
            printer("                                  where nothing constrains it)")
        if c.get("gain_on_bound"):
            printer("  CO-DESIGNED GAIN ON ITS BOUND : %d of %d entries, largest"
                    % (c["gain_on_bound"], c["gain_entries"]))
            printer("                                  |K| %.3g. The bound, not the"
                    % c["gain_norm"])
            printer("                                  problem, is choosing the")
            printer("                                  gain: the design is buying")
            printer("                                  cost with stability it is not")
            printer("                                  being charged for, and the")
            printer("                                  worst violation above is")
            printer("                                  where that is charged")
        printer("  search vs adaptive integrator : %.2e  (%d fixed steps per node"
                % (c["integrator_drift"], c["nsub"]))
        printer("                                  interval against DOP853; if this")
        printer("                                  approached the slack the")
        printer("                                  certificate would be about the")
        printer("                                  wrong trajectory)")
        o = self.out_of_sample
        if o is not None:
            printer("  out of sample, inside the set : mean %.3e, worst %.3e, "
                    "%.1f%% within slack" % (o["mean"], o["worst"], o["within"]))
            printer("  drawn outside the set         : %d of %d (%.2f%%), worst %.3e"
                    % (o["n_out"], o["n"], 100.0 * o["n_out"] / o["n"], o["worst_out"]))
            if o.get("cost"):
                c2 = o["cost"]
                printer("  realised cost out of sample   : mean %.5f, s.d. %.5f, "
                        "90th %.5f, worst %.5f"
                        % (c2["mean"], c2["sd"], c2["p90"], c2["worst"]))
        if self.wait_and_see is not None:
            printer("  wait-and-see lower bound      : %.6f  (the value of knowing"
                    % self.wait_and_see)
            printer("                                  theta in advance: %.1f%%)"
                    % (100.0 * (self.objective - self.wait_and_see)
                       / abs(self.wait_and_see)))
        if self.budget_exhausted:
            printer("  STOPPED EARLY                 : the transcription could not")
            printer("                                  carry another scenario; see")
            printer("                                  the message above for the two")
            printer("                                  remedies")
        printer("  calls to psopt()              : %d" % self.n_solves)


class RobustProblem(object):
    """A robust optimal control problem, solved by scenario augmentation."""

    def __init__(self, name="robust"):
        self.name = name
        self._phases = []
        self.uncertainty = None
        # The state the verification integrator starts from. Either an array, or a
        # callable of theta when the initial condition itself depends on it. It is
        # asked for rather than inferred: the driver cannot in general read x(t0)
        # off the event bounds, and guessing it would put a silent error in the one
        # place the whole verification rests on.
        self.initial_state = None
        # Per-constraint scale factors for the violation measure, so that a metre
        # and a radian are not added together. One entry per event and per path
        # constraint; left None they are all 1 and the measure is the raw infinity
        # norm, which is right only when the constraints share units.
        self.event_scale = None
        self.path_scale = None
        # Accepted here as well as on the phase, for the same reason cost_bounds is.
        self.feedback = None
        self.feedback_bounds = None
        self.feedback_guess = None
        # Accepted here as well as on the phase, because the phase is where it
        # belongs and the problem is where a reader reaches for it first. A plain
        # Python object accepts any attribute silently, so offering only one of the
        # two spellings would mean the other quietly did nothing.
        self.cost_bounds = None

    def add_phase(self, nstates, ncontrols, nevents=0, npath=0, nparameters=0):
        if self._phases:
            raise NotImplementedError(
                "RobustProblem: multi-phase problems are not supported. The "
                "augmentation would have to replicate the linkages as well, and "
                "the scenario copies of a phase boundary are not independent. "
                "Refusing rather than assembling something that looks right.")
        ph = RobustPhase(nstates, ncontrols, nevents, npath, nparameters)
        self._phases.append(ph)
        return ph

    # -- assembly ------------------------------------------------------------------

    def _augment(self, thetas, margins, guess=None, weights=None, risk="nominal",
                 cvar_alpha=0.9, mv_lambda=1.0, pguess=None):
        """Build the deterministic psopt.Problem for a given scenario set.

        The scenario set does two jobs and the driver keeps them apart. As a
        CONSTRAINT SET every member is a plant the design must serve, and every
        member counts equally --- there is no such thing as a constraint that holds
        with weight one sixth. As a QUADRATURE RULE it estimates the risk measure,
        and there the weights are the rule's and matter.

        So the scenarios added by the generation loop enter the constraints and
        carry weight zero in the objective. They are not quadrature nodes; they are
        the places the design was failing, which is a biased sample of the
        uncertainty by construction. Letting them into the risk measure would be
        quietly replacing it with a worst-case-weighted one.

        CARRYING THE COST OF EACH SCENARIO

        "nominal" and "expectation" are sums of per-scenario costs, so they can be
        written straight into the integrand and cost nothing. Mean-variance and
        CVaR cannot: both need each scenario's cost J_k as a quantity in its own
        right, and the variance of a Lagrange cost across scenarios is not the
        integral of anything. So for those two the augmentation carries an extra
        state per scenario whose derivative is that scenario's integrand and whose
        initial value is pinned to zero. J_k is then its final value plus the
        scenario's endpoint term, available to any function of the final state.

        CVaR is written by the Rockafellar-Uryasev device,

            CVaR_alpha = min over eta of  eta + 1/(1-alpha) * sum_k w_k [J_k - eta]+

        with the positive part carried by a slack static parameter per scenario
        rather than by a smoothed hinge: s_k >= 0 and s_k >= J_k - eta, both exact,
        where a smoothed max would put an arbitrary rounding radius between the
        answer and the risk measure that was asked for. PSOPT has static parameters
        and this is what they are for.
        """
        rp = self._phases[0]
        M = len(thetas)
        n, m = rp.nstates, rp.ncontrols
        ne, npth = rp.nevents, rp.npath
        npu = rp.nparameters                       # the USER's static parameters

        if risk not in ("nominal", "expectation", "mean-variance", "cvar"):
            raise ValueError(
                "RobustProblem: risk=%r. Implemented: 'nominal', 'expectation', "
                "'mean-variance', 'cvar'." % (risk,))
        use_cost = risk in ("mean-variance", "cvar")
        if rp.cost_bounds is None and self.cost_bounds is not None:
            rp.cost_bounds = self.cost_bounds
        if rp.feedback_bounds is None and self.feedback_bounds is not None:
            rp.feedback_bounds = self.feedback_bounds
        if rp.feedback_guess is None and self.feedback_guess is not None:
            rp.feedback_guess = self.feedback_guess
        if use_cost and rp.cost_bounds is None:
            raise ValueError(
                "RobustProblem: risk=%r carries each scenario's cost as a state, "
                "and a state needs bounds. Set .cost_bounds = (lo, hi) to a range "
                "the per-scenario cost is certainly inside. It is asked for rather "
                "than guessed because a bound that turns out to be active silently "
                "changes the risk measure into something else." % (risk,))
        if rp.feedback is None and self.feedback is not None:
            rp.feedback = self.feedback
        fb = _feedback_kind(rp.feedback, m, n) if M > 1 else ("none", None)
        fb_kind, fb_value = fb
        # A co-designed gain's entries ARE decision variables: m*n static
        # parameters, optimised alongside the nominal control. Nothing about that
        # is structurally hard -- the problem was already nonconvex, and a product
        # of two decision variables is just another nonlinear term -- but the gain
        # must be bounded, because nothing else stops it growing.
        n_gain = m*n if fb_kind == "co-design" else 0
        # A TIME-VARYING co-designed gain is carried as extra CONTROLS, not as
        # parameters and not as a polynomial in t. The transcription already gives
        # a control a time profile -- piecewise linear within a segment, under
        # multiple shooting -- so the gain inherits one for free, at exactly the
        # resolution the trajectory itself has, with no basis to choose and no
        # branch to refuse. A polynomial in t was tried first and is not practical:
        # see the note in examples/robust_driver_gain.py.
        n_gainu = m*n if fb_kind == "co-design-schedule" else 0
        mu_ = m + n_gainu                      # controls the augmented problem has
        if (n_gain or n_gainu) and rp.feedback_bounds is None \
                and self.feedback_bounds is None:
            raise ValueError(
                "RobustProblem: a co-designed gain needs .feedback_bounds = "
                "(lo, hi). The entries are decision variables and nothing else "
                "bounds them; an unbounded gain runs away into saturation, where "
                "the realised control is not the control that was designed.")

        nx = n * M + (M if use_cost else 0)
        # Parameter layout: the user's, then any co-designed gain, then CVaR's eta
        # and slacks. Fixed here once, because three places read it back.
        gain_off = npu
        cvar_off = npu + n_gain
        npar = cvar_off + (M + 1 if risk == "cvar" else 0)
        # Under feedback the REALISED control differs from scenario to scenario, so
        # its bounds are no longer the decision variable's bounds and have to be
        # imposed as path constraints. They are inequalities, so they cost no
        # degrees of freedom -- see _dof.
        nu_rows = m * (M - 1) if fb_kind != "none" else 0

        w = (np.ones(M) / M if weights is None else np.asarray(weights, dtype=float))
        if len(w) != M:
            raise ValueError("RobustProblem: %d weights for %d scenarios"
                             % (len(w), M))

        nev = ne * M + (M if use_cost else 0) + (M if risk == "cvar" else 0)
        prob = Problem(name=self.name)
        ph = prob.add_phase(nstates=nx, ncontrols=mu_, nevents=nev,
                            npath=npth * M + nu_rows, nparameters=npar)
        ph.nodes = list(rp.nodes)

        th_list = [np.asarray(t, dtype=float) for t in thetas]
        xs = lambda x, k: x[n * k:n * (k + 1)]          # noqa: E731  scenario k
        cs = lambda x, k: x[n * M + k]                  # noqa: E731  its cost state

        # FEEDBACK
        #
        # Open loop, every scenario is driven by the same control history and the
        # design has to find one history that serves them all. That is what makes
        # an open-loop robust design expensive: the arm pays a factor of nearly
        # three in final time for it.
        #
        # With an ancillary gain the scenarios are driven by
        #
        #     u_k(t) = u_bar(t) + K ( x_k(t) - x_ref(t) )
        #
        # where u_bar is still the only decision variable and x_ref is the
        # REFERENCE trajectory -- scenario 0, the centre of every scenario rule the
        # driver offers. Scenario 0 therefore runs open loop by construction, since
        # its deviation from itself is zero, and every other scenario is corrected
        # towards it. The design is still here-and-now: u_bar and K are both fixed
        # before the uncertainty is revealed, and nothing adapts to theta. What the
        # feedback does is stop one history from having to serve every plant
        # unaided.
        #
        # The gain, in whichever of its three forms was asked for. It is a function
        # of time and of the static parameters so that one expression covers all of
        # them: a constant ignores both, a scheduled gain reads the first, and a
        # co-designed gain reads the second.
        def gain_at(u, p, t):
            if fb_kind == "constant":
                return ca.DM(fb_value)
            if fb_kind == "scheduled":
                return ca.reshape(ca.vertcat(*[ca.vertcat(*r)
                                               for r in _rows_of(fb_value(t))]),
                                  m, n)
            if fb_kind == "co-design-schedule":
                return ca.reshape(u[m:m + n_gainu], m, n)
            return ca.reshape(p[gain_off:gain_off + n_gain], m, n)

        def uk(x, u, p, t, k):
            ub = u[0:m]
            if fb_kind == "none" or k == 0:
                return ub
            return ub + ca.mtimes(gain_at(u, p, t), xs(x, k) - xs(x, 0))

        def dynamics(x, u, p, t):
            rows = [rp.dynamics(xs(x, k), uk(x, u, p, t, k), p[0:npu], t, th_list[k])
                    for k in range(M)]
            if use_cost:
                rows += [rp.integrand(xs(x, k), uk(x, u, p, t, k), p[0:npu], t)
                         if rp.integrand is not None else ca.SX(0.0)
                         for k in range(M)]
            return ca.vertcat(*rows)

        ph.dynamics = dynamics
        if npth or nu_rows:
            def path(x, u, p, t):
                rows = [rp.path(xs(x, k), uk(x, u, p, t, k), p[0:npu], t, th_list[k])
                        for k in range(M)] if npth else []
                # The realised control of every corrected scenario, so that its own
                # bounds can be imposed on it. Without this the ancillary gain is
                # free to ask for torque the actuator has not got, and the design
                # would be one no plant could execute. It matters more for a
                # co-designed gain than for a given one, because there the
                # optimiser is actively pushing the gain around.
                if fb_kind != "none":
                    rows += [uk(x, u, p, t, k) for k in range(1, M)]
                return ca.vertcat(*rows)

            ph.path = path

        def J_of(xi, xf, p, t0, tf, k):
            """Scenario k's total cost, as a function of the final state."""
            term = (rp.endpoint(xs(xi, k), xs(xf, k), p[0:npu], t0, tf)
                    if rp.endpoint is not None else ca.SX(0.0))
            return term + (cs(xf, k) if use_cost else ca.SX(0.0))

        def events(xi, xf, p, t0, tf):
            rows = []
            if ne:
                rows += [rp.events(xs(xi, k), xs(xf, k), p[0:npu], t0, tf, th_list[k])
                         for k in range(M)]
            if use_cost:
                rows += [cs(xi, k) for k in range(M)]        # each cost starts at 0
            if risk == "cvar":
                eta = p[cvar_off]
                rows += [p[cvar_off + 1 + k] - (J_of(xi, xf, p, t0, tf, k) - eta)
                         for k in range(M)]                  # s_k >= J_k - eta
            return ca.vertcat(*rows)

        if nev:
            ph.events = events

        if risk == "nominal":
            # The cost of the first scenario, which is the quadrature rule's centre
            # for every rule the driver offers. Feasibility is robust; the objective
            # is the nominal one, and saying so is the whole of the honesty here.
            if rp.endpoint is not None:
                ph.endpoint = lambda xi, xf, p, t0, tf: rp.endpoint(
                    xs(xi, 0), xs(xf, 0), p[0:npu], t0, tf)
            if rp.integrand is not None:
                ph.integrand = lambda x, u, p, t: rp.integrand(
                    xs(x, 0), uk(x, u, p, t, 0), p[0:npu], t)
        elif risk == "expectation":
            # A weighted sum of per-scenario costs is itself a sum, so it goes
            # straight into the integrand and the endpoint and needs no cost state.
            if rp.endpoint is not None:
                ph.endpoint = lambda xi, xf, p, t0, tf: sum(
                    float(w[k]) * rp.endpoint(xs(xi, k), xs(xf, k), p[0:npu], t0, tf)
                    for k in range(M))
            if rp.integrand is not None:
                ph.integrand = lambda x, u, p, t: sum(
                    float(w[k]) * rp.integrand(xs(x, k), uk(x, u, p, t, k),
                                               p[0:npu], t)
                    for k in range(M))
        elif risk == "mean-variance":
            lam = float(mv_lambda)

            def mv(xi, xf, p, t0, tf):
                J = [J_of(xi, xf, p, t0, tf, k) for k in range(M)]
                mean = sum(float(w[k]) * J[k] for k in range(M))
                second = sum(float(w[k]) * J[k] * J[k] for k in range(M))
                return mean + lam * (second - mean * mean)

            ph.endpoint = mv
        else:                                                   # cvar
            a = float(cvar_alpha)
            if not 0.0 <= a < 1.0:
                raise ValueError("RobustProblem: cvar_alpha must be in [0, 1)")

            def cvar(xi, xf, p, t0, tf):
                return p[cvar_off] + (1.0 / (1.0 - a)) * sum(
                    float(w[k]) * p[cvar_off + 1 + k] for k in range(M))

            ph.endpoint = cvar

        ph.bounds.lower.states = _tile(rp.bounds.lower.states, M)
        ph.bounds.upper.states = _tile(rp.bounds.upper.states, M)
        if use_cost:
            clo, chi = rp.cost_bounds
            ph.bounds.lower.states = list(ph.bounds.lower.states) + [clo] * M
            ph.bounds.upper.states = list(ph.bounds.upper.states) + [chi] * M
        ph.bounds.lower.controls = list(rp.bounds.lower.controls)
        ph.bounds.upper.controls = list(rp.bounds.upper.controls)
        if n_gainu:
            glo, ghi = _gain_box(rp.feedback_bounds if rp.feedback_bounds is not None
                                 else self.feedback_bounds, m, n)
            ph.bounds.lower.controls = ph.bounds.lower.controls + list(glo)
            ph.bounds.upper.controls = ph.bounds.upper.controls + list(ghi)
        ph.bounds.lower.parameters = rp.bounds.lower.parameters
        ph.bounds.upper.parameters = rp.bounds.upper.parameters
        plo = list(rp.bounds.lower.parameters or [])
        phi = list(rp.bounds.upper.parameters or [])
        if n_gain:
            glo, ghi = _gain_box(rp.feedback_bounds if rp.feedback_bounds is not None
                                 else self.feedback_bounds, m, n)
            plo = plo + list(glo)
            phi = phi + list(ghi)
            ph.bounds.lower.parameters = plo
            ph.bounds.upper.parameters = phi
        if risk == "cvar":
            clo, chi = rp.cost_bounds
            # eta lives on the same scale as the cost; each slack is a positive
            # part of a difference of two costs, so it cannot exceed their range.
            ph.bounds.lower.parameters = plo + [clo] + [0.0] * M
            ph.bounds.upper.parameters = phi + [chi] + [chi - clo] * M
        ph.bounds.t0 = rp.bounds.t0
        ph.bounds.tf = rp.bounds.tf

        elo, ehi = [], []
        if ne:
            lo, hi = _shrink(rp.bounds.lower.events, rp.bounds.upper.events,
                             margins["events"])
            elo, ehi = _tile(lo, M), _tile(hi, M)
        if use_cost:
            elo, ehi = elo + [0.0] * M, ehi + [0.0] * M
        if risk == "cvar":
            clo, chi = rp.cost_bounds
            elo, ehi = elo + [0.0] * M, ehi + [chi - clo] * M
        if nev:
            ph.bounds.lower.events, ph.bounds.upper.events = elo, ehi
        plo, phi = [], []
        if npth:
            lo, hi = _shrink(rp.bounds.lower.path, rp.bounds.upper.path,
                             margins["path"])
            plo, phi = _tile(lo, M), _tile(hi, M)
        if nu_rows:
            # The user's own control bounds, applied to each corrected scenario's
            # realised control. No margin: these are the actuator's limits, not a
            # requirement the design is being asked to meet with room to spare.
            plo = plo + list(np.tile(np.asarray(rp.bounds.lower.controls,
                                                dtype=float), M - 1))
            phi = phi + list(np.tile(np.asarray(rp.bounds.upper.controls,
                                                dtype=float), M - 1))
        if plo:
            ph.bounds.lower.path, ph.bounds.upper.path = plo, phi

        N = ph.nodes[-1]
        if guess is not None:
            t_g, u_g, x_per = guess
            ph.guess.time = np.asarray(t_g).reshape(1, N)
            ph.guess.controls = _widen_controls(np.asarray(u_g), mu_, N,
                                                self._gain_guess(rp, m, n))
            rows = list(x_per)
        else:
            ph.guess.time = (np.asarray(rp.guess.time).reshape(1, N)
                             if rp.guess.time is not None
                             else np.linspace(0.0, 1.0, N).reshape(1, N))
            base_u = (np.asarray(rp.guess.controls).reshape(m, N)
                      if rp.guess.controls is not None else np.zeros((m, N)))
            ph.guess.controls = _widen_controls(base_u, mu_, N,
                                                self._gain_guess(rp, m, n))
            base = (np.asarray(rp.guess.states) if rp.guess.states is not None
                    else np.zeros((n, N)))
            rows = [base] * M
        if use_cost:
            rows = list(rows) + [self._cost_guess(rows[k],
                                                  np.asarray(ph.guess.controls),
                                                  np.asarray(ph.guess.time).ravel())
                                 for k in range(M)]
        ph.guess.states = np.vstack(rows)
        pg = (np.asarray(rp.guess.parameters, dtype=float).ravel()
              if rp.guess.parameters is not None else np.zeros(npu))
        if n_gain:
            # Start from a gain that already works. A co-designed gain is a
            # nonconvex problem in its own right, and starting it at zero starts it
            # at the open-loop design, which is the expensive local minimum this is
            # meant to escape. Once the generation loop has a gain of its own,
            # start from THAT instead: the scenario being added is one the previous
            # gain nearly served, so its value is a far better guess than the LQR
            # one, and re-starting each iteration from the LQR gain throws away
            # everything the loop has learned about the gain.
            pv = (None if pguess is None
                  else np.asarray(pguess, dtype=float).ravel())
            if pv is not None and len(pv) >= npu + n_gain:
                g0 = pv[npu:npu + n_gain].reshape(m, n)
            else:
                g0 = (rp.feedback_guess if rp.feedback_guess is not None
                      else self.feedback_guess)
                g0 = (np.zeros((m, n)) if g0 is None
                      else np.atleast_2d(np.asarray(g0, dtype=float)))
            pg = np.concatenate([pg, np.asarray(g0, dtype=float).reshape(-1)])
        if risk == "cvar":
            pg = np.concatenate([pg, np.zeros(M + 1)])
        ph.guess.parameters = pg.reshape(-1, 1) if npar else None
        return prob

    @staticmethod
    def _dof(ph):
        """Degrees of freedom the transcription has left, as a LOWER bound.

        Under multiple shooting the state at every node is pinned by the defect
        equations, so the state count cancels and what remains is

            dof = ncontrols * N + free times + free static parameters
                  - equality events

        Every scenario brings its own equality events -- its pinned initial
        condition above all -- and they come out of a budget the shared control
        fixes. Inequality events cost nothing, which is a second reason to relax a
        terminal condition rather than pin it.

        It is a LOWER bound because equality events can be linearly dependent on
        the defect equations, and a dependent constraint removes no freedom. So a
        non-positive count does not prove the problem is over-determined, and this
        is used to EXPLAIN a failure rather than to refuse in advance -- refusing
        on it would block problems that solve perfectly well.
        """
        N = ph.nodes[-1]
        nfree_t = ((1 if ph.bounds.t0[0] != ph.bounds.t0[1] else 0)
                   + (1 if ph.bounds.tf[0] != ph.bounds.tf[1] else 0))
        npar_free = 0
        if ph.bounds.lower.parameters is not None:
            plo = np.asarray(ph.bounds.lower.parameters, dtype=float)
            phi = np.asarray(ph.bounds.upper.parameters, dtype=float)
            npar_free = int(np.sum(phi > plo))
        neq = 0
        if ph.bounds.lower.events is not None:
            elo = np.asarray(ph.bounds.lower.events, dtype=float)
            ehi = np.asarray(ph.bounds.upper.events, dtype=float)
            neq = int(np.sum(ehi <= elo))
        return (ph.ncontrols * N + nfree_t + npar_free - neq, nfree_t, npar_free,
                neq, N)

    def _dof_message(self, ph, M, risk):
        """Why a solve with this many scenarios had nothing left to move."""
        dof, nfree_t, npar_free, neq, N = self._dof(ph)
        per = max(1, neq // max(1, M))
        room = ph.ncontrols * N + nfree_t + npar_free
        most = max(1, (room - 1) // per)
        extra = ("\n      risk=%r adds one pinned event per scenario of its own, "
                 "for the cost state's zero initial value." % risk
                 if risk in ("mean-variance", "cvar") else "")
        return (
            "%d scenarios leave the transcription about %d degrees of freedom: %d "
            "control(s) at %d nodes, %d free time endpoint(s) and %d free static "
            "parameter(s), against %d equality events -- some %d per scenario.\n"
            "      Multiple shooting pins the state at every node through its "
            "defect equations, so the shared control is the only thing left to "
            "satisfy them with, and each scenario spends its own pinned events out "
            "of that budget.\n"
            "      Two remedies, and they differ: raise the node count, which buys "
            "%d per node and would carry roughly %d scenarios at this mesh shape; "
            "or relax pinned events to a tolerance, which costs nothing per "
            "scenario because an inequality event takes no degree of freedom at "
            "all.%s"
            % (M, dof, ph.ncontrols, N, nfree_t, npar_free, neq, per,
               ph.ncontrols, most, extra))

    def _cost_guess(self, x_rows, u_rows, t_row):
        """A running-cost guess, by integrating the integrand along a guessed arc.

        Zero would do and would be worse: the cost state's terminal value IS the
        objective under mean-variance and CVaR, so a guess of zero starts the
        solver with an objective estimate that is wrong by the whole of the cost.
        """
        if self._vL is None:
            return np.zeros((1, len(t_row)))
        npu = self._phases[0].nparameters
        # The guess table handed in may carry a co-designed gain schedule in its
        # trailing rows; the user's integrand knows only the control.
        u_rows = np.asarray(u_rows)[0:self._phases[0].ncontrols, :]
        L = np.array([float(self._vL(x_rows[:, j], u_rows[:, j], np.zeros(npu),
                                     t_row[j]))
                      for j in range(len(t_row))])
        c = np.concatenate([[0.0], np.cumsum(0.5 * (L[1:] + L[:-1]) * np.diff(t_row))])
        return c.reshape(1, -1)

    # -- the verification, at two levels -------------------------------------------
    #
    #  The worst-case search asks for the violation at thousands of parameter
    #  vectors, and the reported numbers ask for a handful of them to be right. Those
    #  are different requirements and the driver meets them with different
    #  integrators.
    #
    #  The SEARCH uses a fixed-step RK4 evaluated at every candidate parameter
    #  simultaneously: CasADi's map() turns the user's equations into one function
    #  over K columns, so a thousand trajectories cost about what one costs in
    #  Python-call overhead, which is what the cost actually is. The REPORTED figures
    #  -- the certificate at the worst parameter, and the agreement check -- use
    #  SciPy's DOP853 with adaptive steps at tolerances far tighter than the design's:
    #  a different method, a different step control and a different code path.
    #
    #  Both share the user's equations and nothing else with the transcription, which
    #  is the property that makes either of them a verification at all. `check_steps`
    #  measures the fixed-step integrator against the adaptive one and reports the
    #  disagreement, so the cheap one is not trusted merely because it is cheap.

    def _build_verifier(self):
        """CasADi functions for the user's equations, scalar and mapped."""
        rp = self._phases[0]
        n, m, npar = rp.nstates, rp.ncontrols, rp.nparameters
        d = self.uncertainty.dim
        x = ca.SX.sym("x", n)
        u = ca.SX.sym("u", m)
        p = ca.SX.sym("p", npar)
        t = ca.SX.sym("t", 1)
        th = ca.SX.sym("th", d)
        xi = ca.SX.sym("xi", n)
        xf = ca.SX.sym("xf", n)
        t0 = ca.SX.sym("t0", 1)
        tf = ca.SX.sym("tf", 1)
        f = ca.Function("f", [x, u, p, t, th],
                        [ca.vertcat(rp.dynamics(x, u, p, t, th))])
        g = (ca.Function("g", [x, u, p, t, th],
                         [ca.vertcat(rp.path(x, u, p, t, th))])
             if rp.npath else None)
        e = (ca.Function("e", [xi, xf, p, t0, tf, th],
                         [ca.vertcat(rp.events(xi, xf, p, t0, tf, th))])
             if rp.nevents else None)
        L = (ca.Function("L", [x, u, p, t], [rp.integrand(x, u, p, t)])
             if rp.integrand is not None else None)
        phi = (ca.Function("phi", [xi, xf, p, t0, tf],
                           [rp.endpoint(xi, xf, p, t0, tf)])
               if rp.endpoint is not None else None)
        return f, g, e, L, phi

    def _maps(self, K):
        """Mapped versions of the verifier functions, cached by width."""
        if self._map_cache.get("K") != K:
            self._map_cache = dict(
                K=K,
                f=self._vf.map(K),
                g=self._vg.map(K) if self._vg is not None else None,
                e=self._ve.map(K) if self._ve is not None else None,
                L=self._vL.map(K) if self._vL is not None else None,
                phi=self._vphi.map(K) if self._vphi is not None else None)
        return self._map_cache

    def _costs_many(self, thetas, t_nodes, u_nodes, params, nsub=None):
        """The realised cost J at every parameter in `thetas`, in one sweep.

        The same vectorised RK4 as the violation sweep, carrying the integrand
        alongside the state, plus the endpoint term at the finish. This is what
        makes the out-of-sample COST distribution reportable rather than merely
        the out-of-sample feasibility, and the difference between a design chosen
        for its mean and one chosen for its tail shows up here and nowhere else.
        """
        rp = self._phases[0]
        thetas = np.atleast_2d(np.asarray(thetas, dtype=float))
        Kg = self._gain()
        nreal = len(thetas)
        if Kg is not None:
            thetas = np.vstack([np.atleast_2d(self._theta_ref), thetas])
        K = len(thetas)
        nsub = self._nsub if nsub is None else nsub
        mp = self._maps(K)
        TH = thetas.T
        P = np.tile(np.asarray(params, dtype=float).ravel()[:rp.nparameters]
                    .reshape(-1, 1), (1, K))
        X = self._x0_of(thetas)
        X0 = X.copy()
        J = np.zeros(K)
        N = len(t_nodes)

        mu = rp.ncontrols

        def realise(Ufull, Xc, tt):
            ub = Ufull[0:mu]
            if Kg is None:
                return ub
            with np.errstate(invalid="ignore", over="ignore"):
                return ub + Kg(tt, params, Ufull[:, 0]) @ (Xc - Xc[:, 0:1])

        def LL(Xc, Uc, tt):
            if mp["L"] is None:
                return np.zeros(K)
            return np.asarray(mp["L"](Xc, Uc, P, np.full((1, K), tt))).ravel()

        for i in range(N - 1):
            ta, tb = t_nodes[i], t_nodes[i + 1]
            h = (tb - ta) / nsub
            ua, ub = u_nodes[:, i:i + 1], u_nodes[:, i + 1:i + 2]
            for sstep in range(nsub):
                w0, wh, w1 = (sstep / float(nsub), (sstep + 0.5) / float(nsub),
                              (sstep + 1.0) / float(nsub))
                BA = np.tile(ua + w0 * (ub - ua), (1, K))
                BH = np.tile(ua + wh * (ub - ua), (1, K))
                BB = np.tile(ua + w1 * (ub - ua), (1, K))
                tA = ta + sstep * h
                tH = ta + (sstep + 0.5) * h
                tB = ta + (sstep + 1) * h
                UA = realise(BA, X, tA)
                k1 = np.asarray(mp["f"](X, UA, P, np.full((1, K), tA), TH))
                l1 = LL(X, UA, tA)
                X2 = X + 0.5 * h * k1
                UH2 = realise(BH, X2, tH)
                k2 = np.asarray(mp["f"](X2, UH2, P, np.full((1, K), tH), TH))
                l2 = LL(X2, UH2, tH)
                X3 = X + 0.5 * h * k2
                UH3 = realise(BH, X3, tH)
                k3 = np.asarray(mp["f"](X3, UH3, P, np.full((1, K), tH), TH))
                l3 = LL(X3, UH3, tH)
                X4 = X + h * k3
                UB4 = realise(BB, X4, tB)
                k4 = np.asarray(mp["f"](X4, UB4, P, np.full((1, K), tB), TH))
                l4 = LL(X4, UB4, tB)
                X = X + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
                J = J + (h / 6.0) * (l1 + 2 * l2 + 2 * l3 + l4)
        if mp["phi"] is not None:
            J = J + np.asarray(mp["phi"](X0, X, P, np.full((1, K), t_nodes[0]),
                                         np.full((1, K), t_nodes[-1]))).ravel()
        return J[1:] if Kg is not None else J[:nreal]

    def _x0_of(self, thetas):
        """Initial states, one column per parameter vector."""
        if callable(self.initial_state):
            return np.array([np.asarray(self.initial_state(th), dtype=float)
                             for th in thetas]).T
        x0 = np.asarray(self.initial_state, dtype=float)
        return np.tile(x0.reshape(-1, 1), (1, len(thetas)))

    def riccati_gain(self, out, Q, R, Qf=None, degree=9, nsub=6):
        """A time-varying ancillary gain, from a Riccati sweep along a design.

        The fixed gain of the first tube is a regulator about one point -- usually
        the target -- and the manoeuvre spends most of its time elsewhere. This
        linearises the plant along the REFERENCE trajectory a design actually
        produced,

            A(t) = df/dx,  B(t) = df/du   at (x_ref(t), u_bar(t), theta_ref),

        integrates the Riccati equation backwards from P(t_f) = Qf,

            -dP/dt = A'P + P A - P B inv(R) B' P + Q,

        and returns K(t) = -inv(R) B'(t) P(t). It is therefore an OUTER iteration:
        design, linearise, re-design. Two passes are usually enough, because the
        second design moves the trajectory much less than the first.

        The gain comes back as a callable of t whose entries are POLYNOMIALS of the
        given degree, fitted to the Riccati solution. A table interpolated between
        nodes would be the obvious alternative and cannot be used: the emitter
        refuses a branch on a symbolic value, and every interpolation is a branch.
        A polynomial is elementary arithmetic, so it survives the trip to C++ --
        and a gain is an approximation being approximated, which is a comfortable
        place to fit a polynomial.

        Qf DEFAULTS TO THE INFINITE-HORIZON VALUE FUNCTION, and that matters.
        With Qf = Q the gain has a terminal boundary layer: P relaxes towards Q
        over the last part of the horizon and the gain collapses with it --
        measured on the arm, |K| falls from 7.0 to 3.0 over the final tenth. The
        feedback is then weakest exactly where the terminal ball has to be met, and
        a global polynomial cannot represent the layer either. Setting Qf to the
        solution of the algebraic Riccati equation at the terminal linearisation
        removes it: the gain ends at the steady-state value it would hold on an
        infinite horizon, and |K| runs 7.08, 7.53, 7.19 instead.

        Two things this does not do. The reference between nodes is interpolated
        linearly rather than integrated: it is a gain, and the expense is not
        warranted. And the fit is over the interval the design produced, so if the
        final time moves a long way at the next solve the gain is being used
        slightly outside it -- which is why the outer iteration re-fits.

        `.fit_error` on the returned callable is the largest discrepancy between
        the fitted polynomial and the Riccati gain it approximates. Read it: what
        the design uses is the polynomial, not the sweep.
        """
        rp = self._phases[0]
        n, m, npu = rp.nstates, rp.ncontrols, rp.nparameters
        d = self.uncertainty.dim
        Q = np.atleast_2d(np.asarray(Q, dtype=float))
        R = np.atleast_2d(np.asarray(R, dtype=float))
        Qf = None if Qf is None else np.atleast_2d(np.asarray(Qf, dtype=float))

        xs_ = ca.SX.sym("x", n)
        us_ = ca.SX.sym("u", m)
        ps_ = ca.SX.sym("p", npu)
        ts_ = ca.SX.sym("t", 1)
        th_ = ca.SX.sym("th", d)
        f = ca.vertcat(rp.dynamics(xs_, us_, ps_, ts_, th_))
        Af = ca.Function("A", [xs_, us_, ps_, ts_, th_], [ca.jacobian(f, xs_)])
        Bf = ca.Function("B", [xs_, us_, ps_, ts_, th_], [ca.jacobian(f, us_)])

        t_nodes = np.asarray(out.time).ravel()
        u_nodes = np.asarray(out.controls).reshape(m, len(t_nodes))
        x_nodes = np.asarray(out.design.states)[0:n, :]        # scenario 0
        th = np.asarray(self._theta_ref, dtype=float)
        pv = np.zeros(npu)
        Rinv = np.linalg.inv(R)

        def at(tt):
            x = np.array([np.interp(tt, t_nodes, x_nodes[i]) for i in range(n)])
            u = np.array([np.interp(tt, t_nodes, u_nodes[i]) for i in range(m)])
            return (np.asarray(Af(x, u, pv, tt, th)),
                    np.asarray(Bf(x, u, pv, tt, th)))

        if Qf is None:
            # The infinite-horizon value function at the terminal linearisation.
            # See the note above on why Q itself is the wrong terminal weight.
            from scipy.linalg import solve_continuous_are
            AT = np.asarray(Af(x_nodes[:, -1], u_nodes[:, -1], pv, t_nodes[-1], th))
            BT = np.asarray(Bf(x_nodes[:, -1], u_nodes[:, -1], pv, t_nodes[-1], th))
            Qf = solve_continuous_are(AT, BT, Q, R)

        from scipy.integrate import solve_ivp

        # The Riccati equation runs BACKWARDS from P(t_f) = Qf:
        #
        #     -dP/dt = A'P + P A - P B inv(R) B' P + Q
        #
        # and solve_ivp marches forwards, so it is integrated in reversed time
        # s = t_f - t, where ds = -dt and the sign flips once:
        #
        #      dP/ds = A'P + P A - P B inv(R) B' P + Q,   P(s=0) = Qf.
        #
        # Written out because the sign was wrong the first time, and the Riccati
        # equation integrated in the unstable direction does not fail politely: it
        # blows up in finite time, the fit returns NaN coefficients, and what
        # surfaces is a compiler error about a non-finite constant in the emitted
        # C++, a long way from here.
        def rhs_reversed(ss, pflat):
            P = pflat.reshape(n, n)
            A, B = at(t_nodes[-1] - ss)
            return (A.T @ P + P @ A - P @ B @ Rinv @ B.T @ P + Q).ravel()

        grid = np.linspace(0.0, t_nodes[-1] - t_nodes[0], max(20, nsub*len(t_nodes)))
        r = solve_ivp(rhs_reversed, (0.0, grid[-1]), Qf.ravel(), t_eval=grid,
                      method="LSODA", rtol=1e-8, atol=1e-10)
        if not r.success or not np.all(np.isfinite(r.y)):
            raise RuntimeError(
                "riccati_gain: the Riccati sweep did not stay finite (%s). The "
                "equation is being integrated in reversed time from P(t_f) = Qf; "
                "if it diverges, the usual causes are a Q or R that is not "
                "positive definite, or a linearisation taken about a trajectory "
                "the design never actually follows." % r.message)

        times = t_nodes[-1] - grid
        Ks = np.empty((len(grid), m, n))
        for j in range(len(grid)):
            P = r.y[:, j].reshape(n, n)
            P = 0.5*(P + P.T)                       # symmetry, lost to round-off
            _A, B = at(times[j])
            Ks[j] = -(Rinv @ B.T @ P)

        order = np.argsort(times)
        tt, Ks = times[order], Ks[order]

        # Fit in NORMALISED time. A polynomial of this degree against absolute
        # seconds is badly conditioned, and the fit IS the gain: nothing
        # downstream ever sees the sweep.
        t0, tf = float(tt[0]), float(tt[-1])
        span = tf - t0 if tf > t0 else 1.0
        tau = (tt - t0)/span
        coeffs = [[np.polyfit(tau, Ks[:, i, j], degree) for j in range(n)]
                  for i in range(m)]
        if not np.all(np.isfinite(np.asarray(coeffs, dtype=float))):
            raise RuntimeError(
                "riccati_gain: the polynomial fit of degree %d returned "
                "non-finite coefficients. Lower the degree." % degree)
        err = 0.0
        for i in range(m):
            for j in range(n):
                err = max(err, float(np.max(np.abs(
                    np.polyval(coeffs[i][j], tau) - Ks[:, i, j]))))

        def K_of_t(t):
            z = (t - t0)/span
            return [[_horner(coeffs[i][j], z) for j in range(n)] for i in range(m)]

        K_of_t.coefficients = coeffs
        K_of_t.samples = (tt, Ks)
        K_of_t.fit_error = err
        K_of_t.interval = (t0, tf)
        return K_of_t

    def _gain_guess(self, rp, m, n):
        """The gain a co-design starts from: the caller's, or zero."""
        g0 = (rp.feedback_guess if rp.feedback_guess is not None
              else self.feedback_guess)
        return (np.zeros((m, n)) if g0 is None
                else np.atleast_2d(np.asarray(g0, dtype=float)))

    def _gain(self):
        """A callable (t, params, u_column) -> the gain matrix, or None.

        One accessor for all four forms, so nothing downstream has to know which
        was asked for. A constant ignores its arguments; a scheduled gain reads the
        time; a co-designed constant gain reads the static parameters, where the
        solver put it; a co-designed SCHEDULE reads the trailing rows of the
        control, where the transcription put it.
        """
        rp = self._phases[0]
        spec = rp.feedback if rp.feedback is not None else self.feedback
        kind, value = _feedback_kind(spec, rp.ncontrols, rp.nstates)
        n, m, npu = rp.nstates, rp.ncontrols, rp.nparameters
        if kind == "none":
            return None
        if kind == "constant":
            return lambda t, params, uc: value
        if kind == "co-design":
            def from_params(t, params, uc):
                pv = np.asarray(params, dtype=float).ravel()
                if len(pv) < npu + m*n:
                    # Before the first solve there is nothing to read; the guess is
                    # the honest stand-in, and saying zero would silently verify an
                    # open-loop design as though it were the closed-loop one.
                    return self._gain_guess(rp, m, n)
                return pv[npu:npu + m*n].reshape(m, n)

            return from_params
        if kind == "co-design-schedule":
            def from_control(t, params, uc):
                uv = np.asarray(uc, dtype=float).ravel()
                if len(uv) < m + m*n:
                    return self._gain_guess(rp, m, n)
                return uv[m:m + m*n].reshape(m, n)

            return from_control

        # Scheduled: compile the user's expression once, then evaluate numerically.
        tsym = ca.SX.sym("t", 1)
        Kexpr = ca.reshape(ca.vertcat(*[ca.vertcat(*r)
                                        for r in _rows_of(value(tsym))]), m, n)
        f = ca.Function("Kfun", [tsym], [Kexpr])
        return lambda t, params, uc: np.asarray(f(t))

    def _feedback_bounds(self):
        """The bounds a co-designed gain was given, or None if it was not one."""
        rp = self._phases[0]
        spec = rp.feedback if rp.feedback is not None else self.feedback
        kind, _v = _feedback_kind(spec, rp.ncontrols, rp.nstates)
        if kind not in ("co-design", "co-design-schedule"):
            return None
        return (rp.feedback_bounds if rp.feedback_bounds is not None
                else self.feedback_bounds)

    def _designed_gain(self, design, u_nodes):
        """Every value a co-designed gain took, as one flat array.

        A constant gain contributes its m*n entries; a schedule contributes them at
        every node, because a schedule that touches its bound anywhere is touching
        it, and the whole point of looking is to notice that it did.
        """
        rp = self._phases[0]
        spec = rp.feedback if rp.feedback is not None else self.feedback
        kind, _v = _feedback_kind(spec, rp.ncontrols, rp.nstates)
        n, m, npu = rp.nstates, rp.ncontrols, rp.nparameters
        if kind == "co-design":
            if design.parameters is None:
                return None
            pv = np.asarray(design.parameters, dtype=float).ravel()
            return pv[npu:npu + m*n] if len(pv) >= npu + m*n else None
        if kind == "co-design-schedule":
            u = np.asarray(u_nodes, dtype=float)
            # Node-major, so that tiling the entry box over the nodes lines
            # up with it.
            return (u[m:m + m*n, :].T.ravel() if u.shape[0] >= m + m*n
                    else None)
        return None

    def _violation_many(self, thetas, t_nodes, u_nodes, params, nsub=None,
                        want_nodes=False, want_uexcess=False):
        """Violation at every parameter in `thetas`, by one vectorised sweep.

        Integrates all K plants together with fixed-step RK4, `nsub` steps per node
        interval, accumulating the largest path excess as it goes rather than
        storing the trajectories -- the running maximum is all the measure needs
        and it keeps the memory flat in K.

        UNDER FEEDBACK this integrates a CLOSED loop, and it has to, because an
        open-loop verification of a closed-loop design checks a controller nobody
        is going to build. The reference trajectory is carried as an extra column
        of the same sweep: column 0 is the nominal plant under u_bar alone, and
        every other column is driven by u_bar + K (x_k - x_0). That makes the
        reference an integration of the nominal plant rather than a table
        interpolated between nodes, which is an assumption about the
        implementation and is stated as one: the controller is taken to regenerate
        the reference by integrating the nominal model, not to store it at the
        design's node spacing. Storing it instead would add an interpolation error
        of its own, and at this node count that error is not negligible beside the
        tolerances being certified.

        With want_uexcess it also returns, per scenario, how far the REALISED
        control went outside its bounds. The design imposes those bounds at the
        segment boundaries; between them nothing does, and a gain that quietly
        asks for more actuator than exists is a design no plant can execute.
        """
        rp = self._phases[0]
        thetas = np.atleast_2d(np.asarray(thetas, dtype=float))
        Kg = self._gain()
        nreal = len(thetas)
        # Column 0 is the reference when there is feedback; otherwise there is no
        # reference and every column is its own plant.
        if Kg is not None:
            thetas = np.vstack([np.atleast_2d(self._theta_ref), thetas])
        K = len(thetas)
        nsub = self._nsub if nsub is None else nsub
        self._nver += K

        mp = self._maps(K)
        TH = thetas.T                                     # (d, K)
        # A risk measure may have appended static parameters of its own -- CVaR's
        # eta and its slacks -- and the user's equations know nothing about them.
        P = np.tile(np.asarray(params, dtype=float).ravel()[:rp.nparameters]
                    .reshape(-1, 1), (1, K))
        X = self._x0_of(thetas)                           # (n, K)
        N = len(t_nodes)
        nodes = np.empty((rp.nstates, N, K)) if want_nodes else None
        if want_nodes:
            nodes[:, 0, :] = X
        worst = np.zeros(K)
        uex = np.zeros(K)

        ulo = (np.asarray(rp.bounds.lower.controls,
                          dtype=float).ravel()[:rp.ncontrols].reshape(-1, 1)
               if rp.bounds.lower.controls is not None else None)
        uhi = (np.asarray(rp.bounds.upper.controls,
                          dtype=float).ravel()[:rp.ncontrols].reshape(-1, 1)
               if rp.bounds.upper.controls is not None else None)

        mu = rp.ncontrols

        def realise(Ufull, Xc, tt):
            """The control each column actually applies, at time tt.

            Ufull may be wider than the user's control vector: a co-designed gain
            schedule rides in its trailing rows, put there by the transcription.
            """
            ub = Ufull[0:mu]
            if Kg is None:
                return ub
            # Errors are not raised here: a column that has already diverged makes
            # this product non-finite, and that is dealt with once, at the end of
            # the sweep, rather than as a warning per step per node.
            with np.errstate(invalid="ignore", over="ignore"):
                U = ub + Kg(tt, params, Ufull[:, 0]) @ (Xc - Xc[:, 0:1])
            if want_uexcess and ulo is not None:
                np.maximum(uex, np.maximum(np.maximum(ulo - U, U - uhi),
                                           0.0).max(axis=0), out=uex)
            return U

        def path_excess(Xc, Uc, tt):
            if not rp.npath:
                return
            pv = np.asarray(mp["g"](Xc, Uc, P, np.full((1, K), tt), TH))
            np.maximum(worst, _bound_excess_many(pv, rp.bounds.lower.path,
                                                 rp.bounds.upper.path,
                                                 self.path_scale), out=worst)

        for i in range(N - 1):
            ta, tb = t_nodes[i], t_nodes[i + 1]
            h = (tb - ta) / nsub
            ua, ub = u_nodes[:, i:i + 1], u_nodes[:, i + 1:i + 2]
            for sstep in range(nsub):
                w0 = sstep / float(nsub)
                wh = (sstep + 0.5) / float(nsub)
                w1 = (sstep + 1.0) / float(nsub)
                BA = np.tile(ua + w0 * (ub - ua), (1, K))
                BH = np.tile(ua + wh * (ub - ua), (1, K))
                BB = np.tile(ua + w1 * (ub - ua), (1, K))
                tA, tH, tB = ta + sstep * h, ta + (sstep + 0.5) * h, ta + (sstep + 1) * h
                UA = realise(BA, X, tA)
                path_excess(X, UA, tA)
                k1 = np.asarray(mp["f"](X, UA, P, np.full((1, K), tA), TH))
                X2 = X + 0.5 * h * k1
                k2 = np.asarray(mp["f"](X2, realise(BH, X2, tH), P,
                                        np.full((1, K), tH), TH))
                X3 = X + 0.5 * h * k2
                k3 = np.asarray(mp["f"](X3, realise(BH, X3, tH), P,
                                        np.full((1, K), tH), TH))
                X4 = X + h * k3
                k4 = np.asarray(mp["f"](X4, realise(BB, X4, tB), P,
                                        np.full((1, K), tB), TH))
                X = X + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            if want_nodes:
                nodes[:, i + 1, :] = X
        path_excess(X, realise(np.tile(u_nodes[:, -1:], (1, K)), X,
                               t_nodes[-1]), t_nodes[-1])

        if rp.nevents:
            ev = np.asarray(mp["e"](self._x0_of(thetas), X, P,
                                    np.full((1, K), t_nodes[0]),
                                    np.full((1, K), t_nodes[-1]), TH))
            np.maximum(worst, _bound_excess_many(ev, rp.bounds.lower.events,
                                                 rp.bounds.upper.events,
                                                 self.event_scale), out=worst)

        # A column whose trajectory left the reals is not an unknown, it is the
        # worst possible outcome: the closed loop diverged at that parameter. NaN
        # would propagate into the search's ranking and into every max() above it,
        # where it does not compare, so the divergence would be invisible to the
        # very search whose job is to find it. Infinity compares correctly and is
        # the truth. This is reachable under an aggressive ancillary gain, which is
        # exactly where it matters: the design is stable at the scenarios it was
        # built from and unstable a little way outside them.
        bad = ~np.isfinite(worst)
        if bad.any():
            worst = np.where(bad, np.inf, worst)
            uex = np.where(bad, np.inf, uex)

        # Drop the reference column: it is scenario 0 of the design and is scored
        # in its own right when it appears in `thetas`.
        if Kg is not None:
            worst, uex = worst[1:], uex[1:]
            if want_nodes:
                nodes = nodes[:, :, 1:]
        assert len(worst) == nreal
        if want_nodes and want_uexcess:
            return worst, nodes, uex
        if want_nodes:
            return worst, nodes
        if want_uexcess:
            return worst, uex
        return worst

    def _violation(self, theta, t_nodes, u_nodes, params):
        """Violation at a single parameter, by the vectorised path."""
        return float(self._violation_many(np.atleast_2d(theta), t_nodes, u_nodes,
                                          params)[0])

    def _violation_ref(self, theta, t_nodes, u_nodes, params, nsub=8):
        """Violation at a single parameter, by SciPy's adaptive DOP853.

        This is the integrator the reported numbers come from. It is far too slow
        to search with and that is not what it is for.
        """
        from scipy.integrate import solve_ivp
        rp = self._phases[0]
        n = rp.nstates
        th = np.asarray(theta, dtype=float)
        thr = np.asarray(self._theta_ref, dtype=float)
        Kg = self._gain()
        params = np.asarray(params, dtype=float).ravel()[:rp.nparameters]
        # Under feedback the reference is part of the state being integrated, so
        # this carries 2n equations: the nominal plant under u_bar, and the
        # scenario under u_bar + K(x - x_ref). Same closed loop as the vectorised
        # sweep, different integrator -- which is the whole point of having both.
        x = self._x0_of(np.atleast_2d(th))[:, 0]
        if Kg is not None:
            x = np.concatenate([self._x0_of(np.atleast_2d(thr))[:, 0], x])
        x0 = x.copy()
        v = 0.0

        def split(xx):
            return (xx[n:], xx[:n]) if Kg is not None else (xx, None)

        # A co-designed gain SCHEDULE rides in the trailing rows of the control
        # table, so `ufull` here may be wider than the user's control vector. The
        # reference plant is driven by the nominal control alone and must be handed
        # only those rows; handing it the whole column is a shape error from
        # CasADi, and would have been a wrong answer if the shapes had happened to
        # agree.
        def nominal(ufull):
            return np.asarray(ufull).ravel()[0:rp.ncontrols]

        def applied(ufull, xx, tt):
            xk, xr = split(xx)
            ub = nominal(ufull)
            return ub if Kg is None else ub + Kg(tt, params, ufull) @ (xk - xr)

        def excess_path(xx, uu, tt):
            if not rp.npath:
                return 0.0
            pv = np.asarray(self._vg(split(xx)[0], uu, params, tt, th)).ravel()
            return _bound_excess(pv, rp.bounds.lower.path, rp.bounds.upper.path,
                                 self.path_scale)

        for i in range(len(t_nodes) - 1):
            ta, tb = t_nodes[i], t_nodes[i + 1]
            ua, ub = u_nodes[:, i], u_nodes[:, i + 1]

            def rhs(tt, xx, ta=ta, tb=tb, ua=ua, ub=ub):
                ubar = ua + (0.0 if tb == ta else (tt - ta) / (tb - ta)) * (ub - ua)
                xk, xr = split(xx)
                dk = np.asarray(self._vf(xk, applied(ubar, xx, tt), params,
                                         tt, th)).ravel()
                if Kg is None:
                    return dk
                dr = np.asarray(self._vf(xr, nominal(ubar), params,
                                         tt, thr)).ravel()
                return np.concatenate([dr, dk])

            grid = np.linspace(ta, tb, nsub + 1)[1:]
            r = solve_ivp(rhs, (ta, tb), x, t_eval=grid, method="DOP853",
                          rtol=1e-11, atol=1e-13)
            if not r.success:
                raise RuntimeError("robust verifier: integration failed at theta = %s"
                                   % _fmt(th))
            for j, tt in enumerate(grid):
                ww = 0.0 if tb == ta else (tt - ta) / (tb - ta)
                ubar = ua + ww * (ub - ua)
                v = max(v, excess_path(r.y[:, j], applied(ubar, r.y[:, j], tt),
                                       tt))
            x = r.y[:, -1]
        if rp.nevents:
            ev = np.asarray(self._ve(split(x0)[0], split(x)[0], params,
                                     t_nodes[0], t_nodes[-1], th)).ravel()
            v = max(v, _bound_excess(ev, rp.bounds.lower.events,
                                     rp.bounds.upper.events, self.event_scale))
        return v

    def _check_steps(self, thetas, t_nodes, u_nodes, params):
        """Measure the search integrator against the adaptive one.

        Reported rather than asserted. If these two disagree by anything close to
        the slack, the certificate is about the wrong trajectory and the step
        count wants raising -- which is exactly the failure the C++ study of this
        problem found at eight steps per segment.
        """
        worst = 0.0
        for th in thetas:
            a = self._violation(th, t_nodes, u_nodes, params)
            b = self._violation_ref(th, t_nodes, u_nodes, params)
            worst = max(worst, abs(a - b))
        return worst

    def _nodes_from(self, thetas, t_nodes, u_nodes, params):
        """Each scenario's own trajectory through the current control, at the nodes.

        This is the warm start: a state guess that is feasible scenario by scenario
        rather than merely the right shape. One vectorised sweep produces all of
        them.
        """
        _v, nodes = self._violation_many(thetas, t_nodes, u_nodes, params,
                                         want_nodes=True)
        return [nodes[:, :, k] for k in range(len(thetas))]

    # -- the worst-case oracle -----------------------------------------------------

    def _worst_case(self, t_nodes, u_nodes, params, n_seed, n_refine, seed):
        """Search the uncertainty set for the parameter served worst.

        Seeding, then refinement, and both are done as VECTORISED sweeps -- every
        candidate integrated in the same pass -- because a sweep over five hundred
        parameters costs barely more than a sweep over one. That rules out a
        conventional local optimiser, which asks for points one at a time and was
        five times the cost of everything else here when this used one.

        The refinement is therefore a shrinking cloud: sample around each of the
        best few points found so far, evaluate the whole cloud at once, keep the
        best, shrink the radius, repeat. It needs no derivative, which suits a
        violation that is a maximum over constraints and over time and so is not
        differentiable where the active constraint changes -- which is exactly
        where the worst case tends to sit.

        For one parameter the seeding alone is effectively exhaustive. In more
        than one it is not, and nothing here claims otherwise.
        """
        U = self.uncertainty
        rng = np.random.default_rng(seed)
        seeds = np.vstack([U.boundary_points(), U.quasi_random(n_seed, seed=seed)])
        vals = self._violation_many(seeds, t_nodes, u_nodes, params)

        n_refine = max(1, n_refine)
        order = np.argsort(-vals)[:n_refine]
        best = [(float(vals[k]), seeds[k].copy()) for k in order]
        lo, hi = U.box()
        radius = 0.25 * (hi - lo)
        per = max(8, 64 // n_refine)
        for _sweep in range(6):
            cloud = []
            for _v, centre in best:
                q = centre + rng.normal(0.0, 1.0, size=(per, U.dim)) * radius
                cloud.extend(U.clip(z) for z in q)
            cloud = np.array(cloud)
            cv = self._violation_many(cloud, t_nodes, u_nodes, params)
            pool = best + [(float(cv[j]), cloud[j]) for j in range(len(cloud))]
            pool.sort(key=lambda r: -r[0])
            best = pool[:n_refine]
            radius = radius * 0.45
        return best[0][0], np.asarray(best[0][1], dtype=float)

    # -- the solve ------------------------------------------------------------------

    def solve(self, algorithm=None, slack=0.0, risk="nominal",
              cvar_alpha=0.9, mv_lambda=1.0,
              scenarios="sigma-points",
              n_scenarios=None, generate=True, max_iterations=12,
              tighten=0.9, margin=None, n_seed=None, n_refine=3,
              out_of_sample=1000, wait_and_see=0, seed=20260927, nsub=16,
              polish=True, verbose=True, printer=print):
        """Design a control that serves every parameter in the uncertainty set.

        slack          how much violation OF THE BOUNDS THE USER DECLARED a
                       design may leave, in the units of the violation measure.
                       It is not the constraint: a terminal ball of radius 0.02 is
                       declared in the event bounds, and slack is how far outside
                       that ball the design may still stray somewhere in the set.
                       Passing the ball's own radius here would quietly double it.
                       Zero demands the declared bounds hold everywhere; a small
                       positive value allows for the resolution of the search and
                       of the integrator, and is what the margin exists to make
                       attainable.
        risk           "nominal" (the cost of the central scenario), "expectation"
                       (the quadrature rule's estimate of E[J]), "mean-variance"
                       (E[J] + mv_lambda * Var[J]) or "cvar" (the conditional
                       value at risk at level cvar_alpha, by Rockafellar-Uryasev).
                       The last two carry each scenario's cost as a state and need
                       .cost_bounds set.
        cvar_alpha     the CVaR level: 0.9 averages the worst tenth
        mv_lambda      the weight on the variance in mean-variance
        scenarios      "sigma-points", "qmc" or "explicit"
        n_scenarios    how many, for "qmc"
        generate       add the worst parameter found and re-solve, iteratively
        tighten        design against constraints shrunk to this fraction of their
                       two-sided half-width; see margin
        margin         absolute inward margin per event and per path constraint,
                       as dict(events=..., path=...). Overrides tighten.
        n_seed         evaluations used to seed the worst-case search
        out_of_sample  how many parameters to draw for scoring, none of which take
                       any part in the design
        wait_and_see   how many independent solves to average for the lower bound
        nsub           fixed steps per node interval in the search integrator,
                       measured against the adaptive one and reported
        polish         after the loop stops, solve the final scenario set once
                       more from the user's own guess and keep the better answer
        """
        if self.uncertainty is None:
            raise ValueError("RobustProblem: set .uncertainty before solving")
        if self.initial_state is None:
            raise ValueError(
                "RobustProblem: set .initial_state. The verification integrates "
                "the plant from it, and the driver will not guess it from the "
                "event bounds.")
        if not self._phases:
            raise ValueError("RobustProblem: no phase")

        rp = self._phases[0]
        alg = algorithm if algorithm is not None else Algorithm()
        self._check_algorithm(alg)
        (self._vf, self._vg, self._ve, self._vL,
         self._vphi) = self._build_verifier()
        self._map_cache = {}
        self._nsub = nsub
        self._nver = 0
        U = self.uncertainty
        if n_seed is None:
            n_seed = 128 if U.dim == 1 else 64 * 2 ** U.dim
        params = (np.asarray(rp.guess.parameters, dtype=float).ravel()
                  if rp.guess.parameters is not None else np.zeros(rp.nparameters))

        margins = _margins(rp, tighten, margin)
        self._check_equalities(rp, 2 if generate else 1)

        thetas, weights = self._initial_scenarios(scenarios, n_scenarios)
        weights = list(weights)
        # The reference for any ancillary feedback is scenario 0 of the rule, which
        # is its centre. It is fixed here and never changes: the generated
        # scenarios are appended, so scenario 0 stays what it was.
        self._theta_ref = np.asarray(thetas[0], dtype=float)
        out = RobustSolution()

        if verbose:
            printer("Robust design: %s" % U.describe())
            printer("  slack %.4g, %d starting scenarios (%s)"
                    % (slack, len(thetas), scenarios))
            printer("  %3s %4s %11s %12s %s"
                    % ("it", "M", "objective", "worst found", "at theta"))

        guess = None
        design = None
        for it in range(max_iterations):
            prob = self._augment(thetas, margins, guess, weights, risk,
                                 cvar_alpha, mv_lambda,
                                 pguess=(params if it else None))
            sol = prob.solve(alg)
            out.n_solves += 1
            if not sol.status.success:
                # A failed solve with no freedom left is the scenario budget
                # running out, and that is worth saying in as many words: left to
                # IPOPT it arrives as return code -10 with no indication of which
                # remedy is wanted, and neither of them is obvious.
                dof = self._dof(prob._phases[0])[0]
                starved = dof <= 0
                if design is None:
                    raise RuntimeError(
                        "RobustProblem: the first solve failed (%s). Nothing below "
                        "would be meaningful.%s"
                        % (sol.status.error_msg or "NLP return code %d"
                           % sol.status.nlp_return_code,
                           ("\n      " + self._dof_message(prob._phases[0],
                                                           len(thetas), risk))
                           if starved else ""))
                if verbose:
                    if starved:
                        printer("  %3d %4d  cannot carry another scenario.\n      %s"
                                % (it, len(thetas),
                                   self._dof_message(prob._phases[0],
                                                     len(thetas), risk)))
                    else:
                        printer("  %3d %4d   adding theta = %s made the problem "
                                "unsolvable; keeping the previous design"
                                % (it, len(thetas), _fmt(thetas[-1])))
                thetas = thetas[:-1]
                weights = weights[:-1]
                out.budget_exhausted = starved
                break

            design = sol
            # Every static parameter the AUGMENTED problem carries, not only the
            # user's: a co-designed constant gain lives in the trailing entries,
            # and the guard used to read rp.nparameters, so on a problem with no
            # parameters of its own -- which is the common case -- the solved gain
            # was never read back and the verification silently checked the design
            # against the gain it had STARTED from. A closed-loop certificate for
            # the wrong loop.
            if sol.parameters is not None:
                params = np.asarray(sol.parameters, dtype=float).ravel()
            t_nodes = np.asarray(sol.time).ravel()
            u_nodes = np.asarray(sol.controls).reshape(-1, len(t_nodes))
            v, th_worst = self._worst_case(t_nodes, u_nodes, params,
                                           n_seed, n_refine, seed + it)
            out.history.append(dict(iteration=it, M=len(thetas),
                                    objective=float(sol.objective),
                                    violation=v, theta=th_worst.copy()))
            if verbose:
                printer("  %3d %4d %11.6f %12.3e  %s"
                        % (it, len(thetas), sol.objective, v, _fmt(th_worst)))

            if v <= slack:
                out.converged = True
                break
            if not generate or it == max_iterations - 1:
                break

            # Warm start: the previous control, with each scenario's own plant
            # integrated through it so that the state guess is feasible scenario by
            # scenario rather than merely the right shape.
            thetas = list(thetas) + [th_worst]
            weights = weights + [0.0]        # a constraint, not a quadrature node
            guess = (t_nodes, u_nodes,
                     self._nodes_from(np.array(thetas), t_nodes, u_nodes, params))

        # Polish. The loop reaches its final scenario set through a chain of warm
        # starts, and each subproblem is nonconvex, so where it ends up is partly
        # the history of that chain rather than a property of the scenario set. One
        # cold solve of the final set answers the question "was the path to it
        # costing anything?" -- and reports the answer either way.
        if polish:
            cold = self._augment(thetas, margins, None, weights, risk,
                                     cvar_alpha, mv_lambda)
            csol = cold.solve(alg)
            out.n_solves += 1
            better = None
            if csol.status.success:
                ct = np.asarray(csol.time).ravel()
                cu = np.asarray(csol.controls).reshape(-1, len(ct))
                # The cold solve has its OWN static parameters, and under co-design
                # that includes its own gain. Verifying it against the warm chain's
                # parameters would compare two designs by integrating one of them
                # wrongly.
                cparams = (np.asarray(csol.parameters, dtype=float).ravel()
                           if csol.parameters is not None else params)
                cv, _cth = self._worst_case(ct, cu, cparams, n_seed, n_refine,
                                            seed + 4242)
                dt = np.asarray(design.time).ravel()
                du = np.asarray(design.controls).reshape(-1, len(dt))
                dv, _dth = self._worst_case(dt, du, params, n_seed, n_refine,
                                            seed + 4242)
                # Prefer a certified design; among certified ones, the cheaper.
                cert_c, cert_d = cv <= slack, dv <= slack
                if cert_c and (not cert_d or csol.objective < design.objective):
                    better = csol
                elif not cert_c and not cert_d and cv < dv:
                    better = csol
                if verbose:
                    printer("  polish: cold solve of the same %d scenarios gives "
                            "%.6f (worst %.3e)" % (len(thetas), csol.objective, cv))
                    printer("          warm chain gave %.6f (worst %.3e); keeping "
                            "the %s" % (design.objective, dv,
                                        "cold one" if better is not None
                                        else "warm one"))
            if better is not None:
                design = better
                params = cparams

        t_nodes = np.asarray(design.time).ravel()
        u_nodes = np.asarray(design.controls).reshape(-1, len(t_nodes))
        v, th_worst = self._worst_case(t_nodes, u_nodes, params, n_seed,
                                       n_refine, seed + 9999)
        # The certificate is reported from the ADAPTIVE integrator at the worst
        # parameter the search found, and the search integrator is measured against
        # it there and at the design scenarios. A certificate carried by the cheap
        # integrator alone would be a claim about the cheap integrator.
        check_at = np.vstack([np.atleast_2d(th_worst), np.array(thetas)])
        drift = self._check_steps(check_at, t_nodes, u_nodes, params)
        v_ref = self._violation_ref(th_worst, t_nodes, u_nodes, params)
        v = max(v, v_ref)

        # Under feedback, how much actuator does the ancillary gain actually ask
        # for? The design imposes the control bounds on the realised control at the
        # segment boundaries; between them nothing does, and the verification
        # samples far more finely. A design whose corrections saturate is one no
        # plant can execute, whatever its certificate says about the states.
        u_excess = 0.0
        if self._gain() is not None:
            probe = np.vstack([np.array(thetas),
                               self.uncertainty.boundary_points(),
                               self.uncertainty.quasi_random(n_seed, seed=seed + 7)])
            _wv, uex = self._violation_many(probe, t_nodes, u_nodes, params,
                                            want_uexcess=True)
            u_excess = float(uex.max())
        out.design = design
        out.scenarios = np.array(thetas)
        out.time = t_nodes
        # The FULL control table, which a co-designed schedule widens: the
        # verifier needs the gain rows and they belong with the control they ride
        # in. out.gain_schedule pulls them out for a caller who wants to look.
        out.controls = u_nodes
        out.gain_schedule = (u_nodes[rp.ncontrols:]
                             if u_nodes.shape[0] > rp.ncontrols else None)
        out.objective = float(design.objective)
        out.certificate = dict(violation=v, theta=th_worst,
                               evaluations=self._nver, slack=slack,
                               integrator_drift=drift, nsub=nsub,
                               control_excess=u_excess,
                               feedback=self._gain() is not None)
        # A co-designed gain sitting on its bound is the signature of the failure
        # this whole apparatus exists to catch: the optimiser is being paid in cost
        # for a gain it is not being charged for, so it takes all the gain it is
        # allowed. Reported rather than refused, because a gain ON its bound is
        # still a valid design if it certifies -- and if it does not, this says
        # where to look first.
        gb = self._feedback_bounds()
        if gb is not None:
            gv = self._designed_gain(design, u_nodes)
            if gv is not None:
                glo, ghi = _gain_box(gb, rp.ncontrols, rp.nstates)
                # A schedule repeats the box at every node.
                r = gv.size // glo.size
                glo, ghi = np.tile(glo, r), np.tile(ghi, r)
                tol = 1.0e-6 * np.maximum(1.0, ghi - glo)
                out.certificate["gain_entries"] = int(gv.size)
                out.certificate["gain_on_bound"] = int(
                    np.count_nonzero((gv <= glo + tol) | (gv >= ghi - tol)))
                out.certificate["gain_norm"] = float(np.abs(gv).max())
        out.converged = out.converged and v <= slack

        if out_of_sample:
            out.out_of_sample = self._score(t_nodes, u_nodes, params,
                                            out_of_sample, slack, seed)
        if wait_and_see:
            # Against the USER's bounds, not the tightened ones: a lower bound
            # computed from a harder problem than the one being bounded is not a
            # lower bound.
            out.wait_and_see = self._wait_and_see(
                alg, dict(events=None, path=None), wait_and_see, seed)
            out.n_solves += wait_and_see
        out.n_verifications = self._nver
        if verbose:
            printer("")
            out.report(printer)
        return out

    # -- pieces of the solve --------------------------------------------------------

    def _check_algorithm(self, alg):
        ms = getattr(alg, "transcription_method", None) == "multiple-shooting"
        lin = getattr(alg, "ms_control_parameterisation", None) == "linear"
        if not (ms and lin):
            warnings.warn(
                "RobustProblem: the verifier interpolates the returned control "
                "table linearly between nodes, which reproduces the designed "
                "control exactly only for transcription_method='multiple-shooting' "
                "with ms_control_parameterisation='linear'. With any other choice "
                "the verification measures a slightly different control from the "
                "one designed, and the difference is charged to the design.",
                stacklevel=3)

    def _check_equalities(self, rp, M):
        """Warn about a PINNED TERMINAL event, which no robust design can meet.

        An event that depends only on the initial state is a shared initial
        condition and pinning it is exactly right --- every scenario starts from
        the same place. An event that depends on the FINAL state is a terminal
        condition, and pinning that is generically infeasible for more than one
        scenario, because one open-loop control cannot steer several different
        plants to the same point.

        Which is which is decided symbolically rather than guessed: CasADi is
        asked whether each event expression depends on the final state at all.
        Guessing it, as an earlier version of this did, warned about correct
        problems whose only events were initial conditions, and a warning that
        fires on correct code is worse than no warning.
        """
        if not rp.nevents or M <= 1:
            return
        if rp.bounds.lower.events is None or rp.bounds.upper.events is None:
            return
        lo = np.asarray(rp.bounds.lower.events, dtype=float)
        hi = np.asarray(rp.bounds.upper.events, dtype=float)
        n, npar, d = rp.nstates, rp.nparameters, self.uncertainty.dim
        xi = ca.SX.sym("xi", n)
        xf = ca.SX.sym("xf", n)
        p = ca.SX.sym("p", npar)
        t0 = ca.SX.sym("t0", 1)
        tf = ca.SX.sym("tf", 1)
        th = ca.SX.sym("th", d)
        ev = ca.vertcat(rp.events(xi, xf, p, t0, tf, th))
        pinned = []
        for i in range(min(rp.nevents, ev.shape[0], len(lo))):
            if hi[i] - lo[i] <= 0.0 and ca.depends_on(ev[i], xf):
                pinned.append(i)
        if pinned:
            warnings.warn(
                "RobustProblem: events %s depend on the final state and are pinned "
                "to a single value. One open-loop control cannot steer several "
                "different plants to the same point, so the augmented problem is "
                "generically infeasible for more than one scenario. Relax them to "
                "a tolerance." % pinned, stacklevel=3)

    def _initial_scenarios(self, scenarios, n_scenarios):
        U = self.uncertainty
        if scenarios == "sigma-points":
            pts, w = U.sigma_points()
        elif scenarios == "qmc":
            if n_scenarios is None:
                n_scenarios = 2 * U.dim + 1
            # Sobol' balances only on powers of two, so the generator rounds the
            # count up and the surplus is dropped here. Returning more points than
            # were asked for would be a surprise, and with a scenario budget as
            # tight as multiple shooting's it would be an expensive one.
            pts = U.quasi_random(n_scenarios, seed=1)[:n_scenarios]
            w = np.full(len(pts), 1.0 / len(pts))
        elif scenarios == "explicit":
            pts, w = U.sigma_points()
        else:
            raise ValueError("RobustProblem: unknown scenarios=%r" % (scenarios,))
        return [np.asarray(p, dtype=float) for p in pts], w

    def _score(self, t_nodes, u_nodes, params, m, slack, seed):
        """Out-of-sample scoring, split at the edge of the uncertainty set.

        Inside the set is what the design promised. Outside it is what truncating
        the distribution gave up, and it is reported separately rather than
        averaged in, because mixing them measures the truncation and not the
        design.
        """
        rng = np.random.default_rng(seed + 1)
        draws = self.uncertainty.sample(m, rng)
        v = self._violation_many(draws, t_nodes, u_nodes, params)
        mask = np.array([self.uncertainty.contains(th) for th in draws])
        inside, outside = v[mask], v[~mask]
        cost = None
        if self._vL is not None or self._vphi is not None:
            Jd = self._costs_many(draws[mask], t_nodes, u_nodes, params)
            if len(Jd):
                cost = dict(mean=float(Jd.mean()), sd=float(Jd.std()),
                            p90=float(np.percentile(Jd, 90.0)),
                            worst=float(Jd.max()))
        return dict(n=m, n_out=len(outside),
                    mean=float(inside.mean()) if len(inside) else float("nan"),
                    worst=float(inside.max()) if len(inside) else float("nan"),
                    within=(100.0 * float(np.mean(inside <= slack))
                            if len(inside) else float("nan")),
                    worst_out=float(outside.max()) if len(outside) else 0.0,
                    cost=cost)

    def _wait_and_see(self, alg, margins, k, seed):
        """The average of k problems each solved knowing its own parameter.

        E[min] <= min E[.], so this is a lower bound on the best achievable EXPECTED
        cost, and hence on any risk measure that dominates the mean -- which
        mean-variance with a non-negative weight and CVaR both do. The gap to it is
        the value of knowing the uncertainty in advance. It is not a design: the
        controls it produces are all different and their average solves nothing.

        Each subproblem is solved with the plain per-scenario cost, because that is
        what "knowing theta in advance" means. Applying a risk measure to a single
        known parameter would be applying it to a degenerate distribution, where
        CVaR and mean-variance both collapse to the cost itself anyway.

        It is an average over k draws and therefore itself a random quantity; with
        k small it is a noisy bound, and the driver reports it as a diagnostic
        rather than putting it in the certificate.
        """
        rng = np.random.default_rng(seed + 2)
        vals = []
        for th in self.uncertainty.sample(k, rng):
            prob = self._augment([th], margins, None, [1.0], "nominal")
            sol = prob.solve(alg)
            if sol.status.success:
                vals.append(float(sol.objective))
        return float(np.mean(vals)) if vals else None


# ----------------------------------------------------------------------------------
#  helpers
# ----------------------------------------------------------------------------------


def _widen_controls(u, mu, N, g0):
    """Pad a control guess with a co-designed gain schedule's own rows."""
    u = np.asarray(u, dtype=float)
    if u.shape[0] >= mu:
        return u[:mu, :]
    pad = np.tile(np.asarray(g0, dtype=float).reshape(-1, 1), (1, N))
    return np.vstack([u, pad[:mu - u.shape[0], :]])


def _gain_box(bounds, m, n):
    """A co-designed gain's box, flattened, from scalars or from (m, n) arrays.

    Scalars give every entry the same range, which is what one writes first and is
    almost never what one wants: the entries of a gain matrix are not commensurate
    -- one multiplies an angle, the next an angular rate -- so a single number wide
    enough for the largest is far too wide for the rest, and a co-designed gain
    takes every bit of width it is given. Arrays let the box be put around a gain
    already known to stabilise the family, which is the only form of co-design this
    driver has found to be worth anything.
    """
    lo, hi = bounds
    lo = np.broadcast_to(np.asarray(lo, dtype=float), (m, n)).reshape(-1)
    hi = np.broadcast_to(np.asarray(hi, dtype=float), (m, n)).reshape(-1)
    if np.any(hi < lo):
        raise ValueError("RobustProblem: feedback_bounds has an upper bound below "
                         "its lower bound")
    return lo, hi


def _horner(c, t):
    """A polynomial with numpy's coefficient order, evaluated in elementary ops."""
    out = float(c[0])
    for k in range(1, len(c)):
        out = out*t + float(c[k])
    return out


def _feedback_kind(spec, m, n):
    """Classify what was put in .feedback, and check its shape while we are here."""
    if spec is None:
        return ("none", None)
    if isinstance(spec, str):
        if spec not in ("co-design", "co-design-schedule"):
            raise ValueError(
                "RobustProblem: feedback=%r. Give None, a gain matrix, a callable "
                "of t, 'co-design' for a constant gain optimised with the "
                "trajectory, or 'co-design-schedule' for a time-varying one."
                % (spec,))
        return (spec, None)
    if callable(spec):
        return ("scheduled", spec)
    K = np.atleast_2d(np.asarray(spec, dtype=float))
    if K.shape != (m, n):
        raise ValueError(
            "RobustProblem: feedback gain is %s, expected (%d, %d) --- one row "
            "per control, one column per state." % (K.shape, m, n))
    return ("constant", K)


def _rows_of(K):
    """Whatever a scheduled gain returned, as a list of lists."""
    if hasattr(K, "shape") and not isinstance(K, (list, tuple)):
        try:
            return [[K[i, j] for j in range(K.shape[1])] for i in range(K.shape[0])]
        except (TypeError, IndexError):
            pass
    return [list(r) for r in K]


def _tile(v, M):
    if v is None:
        return None
    return list(np.tile(np.asarray(v, dtype=float), M))


def _margins(rp, tighten, margin):
    """Inward margins for the design, one per event and per path constraint.

    Given as absolute numbers they are used as given. Derived from `tighten` they
    are (1 - tighten) times the half-width of each TWO-SIDED bound; a one-sided
    bound has no half-width to take a fraction of and gets no automatic margin,
    which is said out loud because a one-sided constraint left at its bound is
    exactly where the between-scenario overshoot appears.
    """
    if margin is not None:
        return dict(events=margin.get("events"), path=margin.get("path"))
    out = {}
    for key, lo, hi in (("events", rp.bounds.lower.events, rp.bounds.upper.events),
                        ("path", rp.bounds.lower.path, rp.bounds.upper.path)):
        if lo is None or hi is None:
            out[key] = None
            continue
        lo = np.asarray(lo, dtype=float)
        hi = np.asarray(hi, dtype=float)
        half = 0.5 * (hi - lo)
        finite = np.isfinite(half) & (np.abs(half) < 1e29)
        mg = np.where(finite, (1.0 - tighten) * half, 0.0)
        # A pinned bound (half-width zero) is an equality; tightening it is
        # meaningless and would make it infeasible rather than merely tight.
        mg = np.where(half <= 0.0, 0.0, mg)
        if np.any(~finite):
            warnings.warn(
                "RobustProblem: %s constraints %s are one-sided, so `tighten` gives "
                "them no margin. Pass margin=dict(%s=[...]) to set one; without it "
                "the design sits exactly on those bounds and the worst case between "
                "scenarios will exceed them."
                % (key, list(np.where(~finite)[0]), key), stacklevel=4)
        out[key] = mg
    return out


def _shrink(lo, hi, margin):
    if lo is None or hi is None:
        return lo, hi
    lo = np.asarray(lo, dtype=float).copy()
    hi = np.asarray(hi, dtype=float).copy()
    if margin is None:
        return list(lo), list(hi)
    mg = np.asarray(margin, dtype=float)
    wide = hi > lo
    lo = np.where(wide, lo + mg, lo)
    hi = np.where(wide, hi - mg, hi)
    return list(lo), list(hi)


def _bound_excess(value, lo, hi, scale):
    """The largest scaled amount by which `value` falls outside [lo, hi]."""
    if lo is None or hi is None:
        return 0.0
    v = np.asarray(value, dtype=float).ravel()
    lo = np.asarray(lo, dtype=float).ravel()
    hi = np.asarray(hi, dtype=float).ravel()
    s = (np.ones_like(v) if scale is None
         else np.asarray(scale, dtype=float).ravel())
    out = np.maximum(np.maximum(lo - v, v - hi), 0.0) / s
    return float(out.max()) if out.size else 0.0


def _bound_excess_many(values, lo, hi, scale):
    """Column-wise `_bound_excess`: values is (nc, K), the result is (K,)."""
    if lo is None or hi is None:
        return np.zeros(values.shape[1])
    v = np.asarray(values, dtype=float)
    lo = np.asarray(lo, dtype=float).reshape(-1, 1)
    hi = np.asarray(hi, dtype=float).reshape(-1, 1)
    s = (np.ones_like(lo) if scale is None
         else np.asarray(scale, dtype=float).reshape(-1, 1))
    out = np.maximum(np.maximum(lo - v, v - hi), 0.0) / s
    return out.max(axis=0) if out.size else np.zeros(v.shape[1])


