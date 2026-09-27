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
    less probability than the same k does in one --- 3 sigma is 99.73% on a line
    and 97.07% in the plane. The driver reports how much mass is outside, because
    that is what the design is not promising anything about.
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
        if self.wait_and_see is not None:
            printer("  wait-and-see lower bound      : %.6f  (the value of knowing"
                    % self.wait_and_see)
            printer("                                  theta in advance: %.1f%%)"
                    % (100.0 * (self.objective - self.wait_and_see)
                       / abs(self.wait_and_see)))
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

    def _augment(self, thetas, margins, guess=None, weights=None, risk="nominal"):
        """Build the deterministic psopt.Problem for a given scenario set.

        The scenario set does two jobs and the driver keeps them apart. As a
        CONSTRAINT SET every member is a plant the design must serve, and every
        member counts equally --- there is no such thing as a constraint that holds
        with weight one sixth. As a QUADRATURE RULE it estimates the expected cost,
        and there the weights are the rule's and matter.

        So the scenarios added by the generation loop enter the constraints and
        carry weight zero in the objective. They are not quadrature nodes; they are
        the places the design was failing, which is a biased sample of the
        uncertainty by construction. Letting them into the expectation would be
        quietly replacing the risk measure with a worst-case-weighted one.
        """
        rp = self._phases[0]
        M = len(thetas)
        n, m = rp.nstates, rp.ncontrols
        ne, npth = rp.nevents, rp.npath

        prob = Problem(name=self.name)
        ph = prob.add_phase(nstates=n * M, ncontrols=m, nevents=ne * M,
                            npath=npth * M, nparameters=rp.nparameters)
        ph.nodes = list(rp.nodes)

        th_list = [np.asarray(t, dtype=float) for t in thetas]

        def dynamics(x, u, p, t):
            return ca.vertcat(*[rp.dynamics(x[n * k:n * (k + 1)], u, p, t, th_list[k])
                                for k in range(M)])

        ph.dynamics = dynamics
        if npth:
            ph.path = lambda x, u, p, t: ca.vertcat(
                *[rp.path(x[n * k:n * (k + 1)], u, p, t, th_list[k]) for k in range(M)])
        if ne:
            ph.events = lambda xi, xf, p, t0, tf: ca.vertcat(
                *[rp.events(xi[n * k:n * (k + 1)], xf[n * k:n * (k + 1)], p, t0, tf,
                            th_list[k]) for k in range(M)])
        if risk == "nominal":
            # The cost of the first scenario, which is the quadrature rule's centre
            # for every rule the driver offers. Feasibility is robust; the objective
            # is the nominal one, and saying so is the whole of the honesty here.
            if rp.endpoint is not None:
                ph.endpoint = lambda xi, xf, p, t0, tf: rp.endpoint(
                    xi[0:n], xf[0:n], p, t0, tf)
            if rp.integrand is not None:
                ph.integrand = lambda x, u, p, t: rp.integrand(x[0:n], u, p, t)
        elif risk == "expectation":
            w = (np.ones(M) / M if weights is None
                 else np.asarray(weights, dtype=float))
            if len(w) != M:
                raise ValueError("RobustProblem: %d weights for %d scenarios"
                                 % (len(w), M))
            if rp.endpoint is not None:
                ph.endpoint = lambda xi, xf, p, t0, tf: sum(
                    float(w[k]) * rp.endpoint(xi[n * k:n * (k + 1)],
                                              xf[n * k:n * (k + 1)], p, t0, tf)
                    for k in range(M))
            if rp.integrand is not None:
                ph.integrand = lambda x, u, p, t: sum(
                    float(w[k]) * rp.integrand(x[n * k:n * (k + 1)], u, p, t)
                    for k in range(M))
        else:
            raise ValueError(
                "RobustProblem: risk=%r. 'nominal' and 'expectation' are "
                "implemented. Mean-variance and CVaR are not, and they are not "
                "one-liners: the variance of a Lagrange cost across scenarios is "
                "not the integral of anything, so each scenario needs its running "
                "cost carried as an extra state before either can be written."
                % (risk,))

        ph.bounds.lower.states = _tile(rp.bounds.lower.states, M)
        ph.bounds.upper.states = _tile(rp.bounds.upper.states, M)
        ph.bounds.lower.controls = rp.bounds.lower.controls
        ph.bounds.upper.controls = rp.bounds.upper.controls
        ph.bounds.lower.parameters = rp.bounds.lower.parameters
        ph.bounds.upper.parameters = rp.bounds.upper.parameters
        ph.bounds.t0 = rp.bounds.t0
        ph.bounds.tf = rp.bounds.tf
        if ne:
            lo, hi = _shrink(rp.bounds.lower.events, rp.bounds.upper.events,
                             margins["events"])
            ph.bounds.lower.events = _tile(lo, M)
            ph.bounds.upper.events = _tile(hi, M)
        if npth:
            lo, hi = _shrink(rp.bounds.lower.path, rp.bounds.upper.path,
                             margins["path"])
            ph.bounds.lower.path = _tile(lo, M)
            ph.bounds.upper.path = _tile(hi, M)

        N = ph.nodes[-1]
        if guess is not None:
            t_g, u_g, x_per = guess
            ph.guess.time = np.asarray(t_g).reshape(1, N)
            ph.guess.controls = np.asarray(u_g).reshape(m, N)
            ph.guess.states = np.vstack(x_per)
        else:
            ph.guess.time = (np.asarray(rp.guess.time).reshape(1, N)
                             if rp.guess.time is not None
                             else np.linspace(0.0, 1.0, N).reshape(1, N))
            ph.guess.controls = (np.asarray(rp.guess.controls).reshape(m, N)
                                 if rp.guess.controls is not None
                                 else np.zeros((m, N)))
            base = (np.asarray(rp.guess.states) if rp.guess.states is not None
                    else np.zeros((n, N)))
            ph.guess.states = np.vstack([base] * M)
        ph.guess.parameters = rp.guess.parameters
        return prob

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
        return f, g, e

    def _maps(self, K):
        """Mapped versions of the verifier functions, cached by width."""
        if self._map_cache.get("K") != K:
            self._map_cache = dict(
                K=K,
                f=self._vf.map(K),
                g=self._vg.map(K) if self._vg is not None else None,
                e=self._ve.map(K) if self._ve is not None else None)
        return self._map_cache

    def _x0_of(self, thetas):
        """Initial states, one column per parameter vector."""
        if callable(self.initial_state):
            return np.array([np.asarray(self.initial_state(th), dtype=float)
                             for th in thetas]).T
        x0 = np.asarray(self.initial_state, dtype=float)
        return np.tile(x0.reshape(-1, 1), (1, len(thetas)))

    def _violation_many(self, thetas, t_nodes, u_nodes, params, nsub=None,
                        want_nodes=False):
        """Violation at every parameter in `thetas`, by one vectorised sweep.

        Integrates all K plants together with fixed-step RK4, `nsub` steps per node
        interval, accumulating the largest path excess as it goes rather than
        storing the trajectories -- the running maximum is all the measure needs
        and it keeps the memory flat in K.
        """
        rp = self._phases[0]
        thetas = np.atleast_2d(np.asarray(thetas, dtype=float))
        K = len(thetas)
        nsub = self._nsub if nsub is None else nsub
        self._nver += K

        mp = self._maps(K)
        TH = thetas.T                                     # (d, K)
        P = np.tile(np.asarray(params, dtype=float).reshape(-1, 1), (1, K))
        X = self._x0_of(thetas)                           # (n, K)
        N = len(t_nodes)
        nodes = np.empty((rp.nstates, N, K)) if want_nodes else None
        if want_nodes:
            nodes[:, 0, :] = X
        worst = np.zeros(K)

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
                UA = np.tile(ua + w0 * (ub - ua), (1, K))
                UH = np.tile(ua + wh * (ub - ua), (1, K))
                UB = np.tile(ua + w1 * (ub - ua), (1, K))
                tA, tH, tB = ta + sstep * h, ta + (sstep + 0.5) * h, ta + (sstep + 1) * h
                path_excess(X, UA, tA)
                k1 = np.asarray(mp["f"](X, UA, P, np.full((1, K), tA), TH))
                k2 = np.asarray(mp["f"](X + 0.5 * h * k1, UH, P,
                                        np.full((1, K), tH), TH))
                k3 = np.asarray(mp["f"](X + 0.5 * h * k2, UH, P,
                                        np.full((1, K), tH), TH))
                k4 = np.asarray(mp["f"](X + h * k3, UB, P,
                                        np.full((1, K), tB), TH))
                X = X + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            if want_nodes:
                nodes[:, i + 1, :] = X
        path_excess(X, np.tile(u_nodes[:, -1:], (1, K)), t_nodes[-1])

        if rp.nevents:
            ev = np.asarray(mp["e"](self._x0_of(thetas), X, P,
                                    np.full((1, K), t_nodes[0]),
                                    np.full((1, K), t_nodes[-1]), TH))
            np.maximum(worst, _bound_excess_many(ev, rp.bounds.lower.events,
                                                 rp.bounds.upper.events,
                                                 self.event_scale), out=worst)
        return (worst, nodes) if want_nodes else worst

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
        th = np.asarray(theta, dtype=float)
        params = np.asarray(params, dtype=float)
        x = self._x0_of(np.atleast_2d(th))[:, 0]
        x0 = x.copy()
        v = 0.0

        def excess_path(xx, uu, tt):
            if not rp.npath:
                return 0.0
            pv = np.asarray(self._vg(xx, uu, params, tt, th)).ravel()
            return _bound_excess(pv, rp.bounds.lower.path, rp.bounds.upper.path,
                                 self.path_scale)

        for i in range(len(t_nodes) - 1):
            ta, tb = t_nodes[i], t_nodes[i + 1]
            ua, ub = u_nodes[:, i], u_nodes[:, i + 1]

            def rhs(tt, xx, ta=ta, tb=tb, ua=ua, ub=ub):
                w = 0.0 if tb == ta else (tt - ta) / (tb - ta)
                return np.asarray(self._vf(xx, ua + w * (ub - ua), params, tt,
                                           th)).ravel()

            grid = np.linspace(ta, tb, nsub + 1)[1:]
            r = solve_ivp(rhs, (ta, tb), x, t_eval=grid, method="DOP853",
                          rtol=1e-11, atol=1e-13)
            if not r.success:
                raise RuntimeError("robust verifier: integration failed at theta = %s"
                                   % _fmt(th))
            for j, tt in enumerate(grid):
                w = 0.0 if tb == ta else (tt - ta) / (tb - ta)
                v = max(v, excess_path(r.y[:, j], ua + w * (ub - ua), tt))
            x = r.y[:, -1]
        if rp.nevents:
            ev = np.asarray(self._ve(x0, x, params, t_nodes[0], t_nodes[-1],
                                     th)).ravel()
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
        risk           "nominal" (the cost of the central scenario) or
                       "expectation" (the quadrature rule's estimate of E[J])
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
        self._vf, self._vg, self._ve = self._build_verifier()
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
            prob = self._augment(thetas, margins, guess, weights, risk)
            sol = prob.solve(alg)
            out.n_solves += 1
            if not sol.status.success:
                if design is None:
                    raise RuntimeError(
                        "RobustProblem: the first solve failed (%s). Nothing below "
                        "would be meaningful." % (sol.status.error_msg or
                                                  "NLP return code %d"
                                                  % sol.status.nlp_return_code))
                if verbose:
                    printer("  %3d %4d   adding theta = %s made the problem "
                            "unsolvable; keeping the previous design"
                            % (it, len(thetas), _fmt(thetas[-1])))
                thetas = thetas[:-1]
                weights = weights[:-1]
                break

            design = sol
            if sol.parameters is not None and rp.nparameters:
                params = np.asarray(sol.parameters, dtype=float).ravel()
            t_nodes = np.asarray(sol.time).ravel()
            u_nodes = np.asarray(sol.controls).reshape(rp.ncontrols, len(t_nodes))
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
            cold = self._augment(thetas, margins, None, weights, risk)
            csol = cold.solve(alg)
            out.n_solves += 1
            better = None
            if csol.status.success:
                ct = np.asarray(csol.time).ravel()
                cu = np.asarray(csol.controls).reshape(rp.ncontrols, len(ct))
                cv, _cth = self._worst_case(ct, cu, params, n_seed, n_refine,
                                            seed + 4242)
                dt = np.asarray(design.time).ravel()
                du = np.asarray(design.controls).reshape(rp.ncontrols, len(dt))
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

        t_nodes = np.asarray(design.time).ravel()
        u_nodes = np.asarray(design.controls).reshape(rp.ncontrols, len(t_nodes))
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
        out.design = design
        out.scenarios = np.array(thetas)
        out.time = t_nodes
        out.controls = u_nodes
        out.objective = float(design.objective)
        out.certificate = dict(violation=v, theta=th_worst,
                               evaluations=self._nver, slack=slack,
                               integrator_drift=drift, nsub=nsub)
        out.converged = out.converged and v <= slack

        if out_of_sample:
            out.out_of_sample = self._score(t_nodes, u_nodes, params,
                                            out_of_sample, slack, seed)
        if wait_and_see:
            # Against the USER's bounds, not the tightened ones: a lower bound
            # computed from a harder problem than the one being bounded is not a
            # lower bound.
            out.wait_and_see = self._wait_and_see(
                alg, dict(events=None, path=None), wait_and_see, seed, risk)
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
            pts = U.quasi_random(n_scenarios, seed=1)
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
        return dict(n=m, n_out=len(outside),
                    mean=float(inside.mean()) if len(inside) else float("nan"),
                    worst=float(inside.max()) if len(inside) else float("nan"),
                    within=(100.0 * float(np.mean(inside <= slack))
                            if len(inside) else float("nan")),
                    worst_out=float(outside.max()) if len(outside) else 0.0)

    def _wait_and_see(self, alg, margins, k, seed, risk="nominal"):
        """The average of k problems each solved knowing its own parameter.

        E[min] <= min E[.], so this is a lower bound on any implementable design
        and the gap to it is the value of knowing the uncertainty in advance. It is
        not a design: the controls it produces are all different and their average
        solves nothing.

        It is an average over k draws and therefore itself a random quantity; with
        k small it is a noisy bound, and the driver reports it as a diagnostic
        rather than putting it in the certificate.
        """
        rng = np.random.default_rng(seed + 2)
        vals = []
        for th in self.uncertainty.sample(k, rng):
            prob = self._augment([th], margins, None, [1.0], risk)
            sol = prob.solve(alg)
            if sol.status.success:
                vals.append(float(sol.objective))
        return float(np.mean(vals)) if vals else None


# ----------------------------------------------------------------------------------
#  helpers
# ----------------------------------------------------------------------------------


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


