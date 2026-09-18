"""psopt -- Pythonic problem-definition API (B-1: single phase).

Dynamics, costs and events are supplied as Python callables over CasADi symbols;
they are traced, codegen'd and CppAD-taped exactly as a hand-written PSOPT example
(faithful route). The C++ engine is unchanged.
"""
import os
import numpy as np
import casadi as ca
from . import codegen
# Load the extension with RTLD_GLOBAL so PSOPT symbols it carries (e.g. auto_link,
# pulled in via --whole-archive) are visible to the JIT-compiled math .so that the
# driver dlopens for multi-phase linkages.
import sys as _sys, os as _os
_flags = _sys.getdlopenflags()
_sys.setdlopenflags(_flags | _os.RTLD_GLOBAL)
from . import _psopt
_sys.setdlopenflags(_flags)


class PSOPTError(RuntimeError):
    """PSOPT refused the problem, or something went wrong setting it up.

    Raised when the solve came back with a non-zero error_flag: an invalid option,
    a bad dimension, a guess of the wrong shape, an exception thrown inside the
    solve. It is NOT raised when the NLP merely fails to converge -- that is a
    result rather than an error, it comes back through sol.status.success, and a
    script that wants to look at a non-converged trajectory can.

    The exception carries the Solution it came from as .solution, so that
    sol.status and whatever the solve did produce are still reachable:

        try:
            sol = prob.solve(alg)
        except psopt.PSOPTError as e:
            print(e)                      # what PSOPT said
            print(e.solution.status)      # the return codes

    Why this exists. PSOPT's own default is algorithm.on_error = "fail-fast",
    which calls exit(). In C++ that is a defensible default. Inside a Python
    process it terminates the interpreter: no traceback, no exception to catch,
    nothing in a notebook but a dead kernel -- and with print_level = 0 it does
    so in silence, because the diagnostic goes through PSOPT's own printing.
    Three of the examples in python/examples were written to a wrong option value
    during development and each time the script simply stopped with exit status 1
    and no output at all.

    So the Python interface defaults on_error to "fail-soft" and raises this
    instead. Passing on_error="fail-fast" to Algorithm restores PSOPT's own
    behaviour, exit() included.
    """
    def __init__(self, message, solution=None):
        RuntimeError.__init__(self, message)
        self.solution = solution


# ---- guess helpers -------------------------------------------------------------
def ramp(pairs, n):
    """Build an (len(pairs) x n) guess matrix; row i is linspace(start_i, end_i, n).
    e.g. ramp([(1.0, -1.0)], 50) for a single state sweeping 1 -> -1."""
    return np.vstack([np.linspace(a, b, n) for (a, b) in pairs])


def const(values, n):
    """Build a (len(values) x n) guess matrix with each row held constant."""
    return np.vstack([np.full(n, float(v)) for v in values])


def times(t0, tf, n):
    """Row-vector time guess linspace(t0, tf, n) shaped (1 x n)."""
    return np.linspace(t0, tf, n).reshape(1, n)


class _Bound:
    def __init__(self):
        self.states = None
        self.controls = None
        self.events = None
        self.path = None
        self.parameters = None


class _Bounds:
    def __init__(self):
        self.lower = _Bound()
        self.upper = _Bound()
        self.t0 = (0.0, 0.0)
        self.tf = (0.0, 0.0)


class _Guess:
    def __init__(self):
        self.states = None
        self.controls = None
        self.time = None
        self.parameters = None


class Phase:
    def __init__(self, nstates, ncontrols, nevents=0, npath=0, nparameters=0,
                 nobserved=0, nsamples=0, nalgebraic=0):
        self.nstates = nstates
        self.ncontrols = ncontrols
        self.nevents = nevents
        self.npath = npath
        # How many of the controls are the algebraic variables of a semi-explicit index-1
        # DAE: the LAST nalgebraic controls, determined by the FIRST nalgebraic path
        # constraints, which must be equalities. Read only by transcription_method =
        # "multiple-shooting", whose segment integrator then solves them at every stage
        # instead of holding them across a segment.
        self.nalgebraic = nalgebraic
        self.nparameters = nparameters
        self.nobserved = nobserved
        self.nsamples = nsamples
        self.nodes = [20]
        # user maths: callables over CasADi symbols (set by the user)
        self.dynamics = None    # f(x,u,p,t) -> xdot (vector length nstates)
        self.path = None        # g(x,u,p,t) -> path  (vector length npath) or None
        self.integrand = None   # L(x,u,p,t) -> scalar or None (=> 0)
        self.endpoint = None    # phi(xi,xf,p,t0,tf) -> scalar or None (=> 0)
        self.events = None      # e(xi,xf,p,t0,tf) -> vector length nevents or None
        self.observation = None # h(x,u,p,t) -> vector length nobserved (param estimation)
        self.observation_nodes = None  # 1 x nsamples sample times
        self.observations = None       # nobserved x nsamples measured data
        # Weight applied to each scalar residual, nobserved x nsamples. Left as None the
        # estimation is unweighted, which is what PSOPT fills in; set it to weight
        # observations that are not equally trustworthy -- the reciprocal of each
        # measurement's standard deviation is the usual choice.
        self.residual_weights = None
        # Tikhonov regularisation of the parameter vector: this multiple of ||p||^2 is
        # added to the least-squares objective. Zero, the default, is no regularisation.
        self.regularization_factor = None
        self.bounds = _Bounds()
        self.guess = _Guess()
        # Discrete-valued declarations. Each entry is (index, [admissible values]).
        # Integer controls are relaxed by outer convexification over the Cartesian
        # product of their admissible sets and recovered by sum-up rounding; integer
        # parameters cannot chatter and are solved exactly by enumeration, which
        # switches the driver to psopt_solve_integer. See declare_integer_control and
        # declare_integer_parameter below.
        self.integer_controls = []
        self.integer_parameters = []

    def declare_integer_control(self, control_index, values):
        """Restrict a control to a finite set of values.

        The control is written in the dynamics as an ordinary control; PSOPT performs
        the outer convexification internally. Declaring several integer controls in a
        phase convexifies over the product of their admissible sets, so the cost grows
        as that product: two binary controls cost four weights per node.
        """
        vals = [float(v) for v in values]
        if not vals:
            raise ValueError("declare_integer_control: the admissible set is empty")
        if not 0 <= control_index < self.ncontrols:
            raise ValueError("declare_integer_control: control_index %d is outside "
                             "the phase's %d controls" % (control_index, self.ncontrols))
        for k, (idx, _v) in enumerate(self.integer_controls):
            if idx == control_index:
                self.integer_controls[k] = (control_index, vals)
                return
        self.integer_controls.append((control_index, vals))

    def declare_integer_parameter(self, parameter_index, values):
        """Restrict a static parameter to a finite set of values.

        Unlike a control, a static parameter cannot chatter, so no relaxation is
        tight for it. The problem is solved exactly by enumerating the admissible
        combinations across all phases, solving each with the parameters pinned, and
        keeping the best. The number of solves is the product of the set sizes.
        """
        vals = [float(v) for v in values]
        if not vals:
            raise ValueError("declare_integer_parameter: the admissible set is empty")
        if not 0 <= parameter_index < self.nparameters:
            raise ValueError("declare_integer_parameter: parameter_index %d is outside "
                             "the phase's %d parameters" % (parameter_index, self.nparameters))
        for k, (idx, _v) in enumerate(self.integer_parameters):
            if idx == parameter_index:
                self.integer_parameters[k] = (parameter_index, vals)
                return
        self.integer_parameters.append((parameter_index, vals))

    def _build_functions(self):
        x  = ca.SX.sym("x",  self.nstates)
        u  = ca.SX.sym("u",  self.ncontrols)
        p  = ca.SX.sym("p",  self.nparameters)
        t  = ca.SX.sym("t",  1)
        xi = ca.SX.sym("xi", self.nstates)
        xf = ca.SX.sym("xf", self.nstates)
        t0 = ca.SX.sym("t0", 1)
        tf = ca.SX.sym("tf", 1)
        dx   = self.dynamics(x, u, p, t) if self.dynamics else ca.SX.zeros(self.nstates)
        dx   = ca.vertcat(dx) if self.nstates else ca.SX.zeros(0)
        path = self.path(x, u, p, t) if self.path else ca.SX.zeros(0)
        L    = self.integrand(x, u, p, t) if self.integrand else ca.SX(0)
        phi  = self.endpoint(xi, xf, p, t0, tf) if self.endpoint else ca.SX(0)
        ev   = self.events(xi, xf, p, t0, tf) if self.events else ca.SX.zeros(0)
        ev   = ca.vertcat(ev) if self.nevents else ca.SX.zeros(0)
        dae_f = ca.Function("dae", [x, u, p, t], [dx, path])
        L_f   = ca.Function("L",   [x, u, p, t], [L])
        phi_f = ca.Function("phi", [xi, xf, p, t0, tf], [phi])
        ev_f  = ca.Function("ev",  [xi, xf, p, t0, tf], [ev])
        obs_f = None
        if self.observation is not None:
            obs = ca.vertcat(self.observation(x, u, p, t))
            obs_f = ca.Function("obs", [x, u, p, t], [obs])
        return dae_f, L_f, phi_f, ev_f, obs_f


class Algorithm:
    def __init__(self, collocation_method="Legendre", nlp_method="IPOPT",
                 derivatives="automatic", scaling="automatic",
                 nlp_iter_max=1000, nlp_tolerance=1e-6,
                 # general extras (None => leave PSOPT default)
                 print_level=None, hessian=None, ipopt_linear_solver=None,
                 constraint_scaling=None, ipopt_max_cpu_time=None, diagnostic_level=None,
                 # hp-adaptive mesh refinement
                 mesh_refinement=None, mr_max_iterations=None, ode_tolerance=None,
                 mr_max_growth_factor=None, mr_min_order=None, mr_max_order=None,
                 # integrated-residual transcription / Nie-Kerrigan flexible-order
                 transcription_method=None, ir_residual_nodes=None, ir_regularization=None,
                 ir_objective=None, ir_residual_bound=None, ir_dair=None,
                 ir_dair_delta_factor=None, ir_local_order=None,
                 ir_include_path=None, ir_path_weight=None, ir_residual_scaling=None,
                 ir_element_local_controls=None, ir_flexible_mesh=None,
                 ir_min_element_fraction=None,
                 # multiple-shooting transcription
                 ms_steps_per_segment=None, ms_integrator=None,
                 ms_control_parameterisation=None,
                 ms_path_samples=None, ms_flexible_segments=None,
                 ms_min_segment_fraction=None, ms_refine_tolerance=None,
                 ms_algebraic_iterations=None, ms_adaptive_steps=None,
                 ms_max_steps_per_segment=None, ms_implicit_iterations=None,
                 # PSOPT's own SQP solver (nlp_method="SQP")
                 qp_solver=None, qp_restoration=None, sqp_strategy=None,
                 qp_iter_max=None, trust_region=None, trust_region_radius=None,
                 elastic_penalty=None,
                 # parameter estimation
                 parameter_statistics=None, parameter_estimation_norm=None,
                 # transcription and solver detail
                 objective_form=None, defect_scaling=None, diff_matrix=None,
                 jac_sparsity_ratio=None, hess_sparsity_ratio=None,
                 save_sparsity_pattern=None, nsteps_error_integration=None,
                 mr_kappa=None, mr_M1=None, mr_switch_detection=None, switch_order=None,
                 hessian_verify=None, on_error="fail-soft",
                 max_integer_combinations=None):
        self.collocation_method = collocation_method
        self.nlp_method = nlp_method
        self.derivatives = derivatives
        self.scaling = scaling
        self.nlp_iter_max = nlp_iter_max
        self.nlp_tolerance = nlp_tolerance
        self.print_level = print_level
        self.hessian = hessian
        self.ipopt_linear_solver = ipopt_linear_solver
        self.constraint_scaling = constraint_scaling
        self.ipopt_max_cpu_time = ipopt_max_cpu_time
        self.diagnostic_level = diagnostic_level
        self.mesh_refinement = mesh_refinement
        self.mr_max_iterations = mr_max_iterations
        self.ode_tolerance = ode_tolerance
        self.mr_max_growth_factor = mr_max_growth_factor
        self.mr_min_order = mr_min_order
        self.mr_max_order = mr_max_order
        self.transcription_method = transcription_method
        self.ir_residual_nodes = ir_residual_nodes
        self.ir_regularization = ir_regularization
        self.ir_objective = ir_objective
        self.ir_residual_bound = ir_residual_bound
        self.ir_dair = ir_dair
        self.ir_dair_delta_factor = ir_dair_delta_factor
        self.ir_local_order = ir_local_order
        self.ir_include_path = ir_include_path
        self.ir_path_weight = ir_path_weight
        self.ir_residual_scaling = ir_residual_scaling
        self.ir_element_local_controls = ir_element_local_controls
        self.ir_flexible_mesh = ir_flexible_mesh
        self.ir_min_element_fraction = ir_min_element_fraction
        self.ms_steps_per_segment = ms_steps_per_segment
        self.ms_integrator = ms_integrator
        self.ms_control_parameterisation = ms_control_parameterisation
        self.ms_path_samples = ms_path_samples
        self.ms_flexible_segments = ms_flexible_segments
        self.ms_min_segment_fraction = ms_min_segment_fraction
        self.ms_refine_tolerance = ms_refine_tolerance
        self.ms_algebraic_iterations = ms_algebraic_iterations
        self.ms_adaptive_steps = ms_adaptive_steps
        self.ms_max_steps_per_segment = ms_max_steps_per_segment
        self.ms_implicit_iterations = ms_implicit_iterations
        self.qp_solver = qp_solver
        self.qp_restoration = qp_restoration
        self.sqp_strategy = sqp_strategy
        self.qp_iter_max = qp_iter_max
        self.trust_region = trust_region
        self.trust_region_radius = trust_region_radius
        self.elastic_penalty = elastic_penalty
        self.parameter_statistics = parameter_statistics
        self.parameter_estimation_norm = parameter_estimation_norm
        self.objective_form = objective_form
        self.defect_scaling = defect_scaling
        self.diff_matrix = diff_matrix
        self.jac_sparsity_ratio = jac_sparsity_ratio
        self.hess_sparsity_ratio = hess_sparsity_ratio
        self.save_sparsity_pattern = save_sparsity_pattern
        self.nsteps_error_integration = nsteps_error_integration
        self.mr_kappa = mr_kappa
        self.mr_M1 = mr_M1
        self.mr_switch_detection = mr_switch_detection
        self.switch_order = switch_order
        self.hessian_verify = hessian_verify
        self.on_error = on_error
        self.max_integer_combinations = max_integer_combinations


def _col(a):
    return np.ascontiguousarray(np.atleast_2d(np.asarray(a, dtype=float)))


def _phase_dict(ph):
    return {
        "nstates": ph.nstates, "ncontrols": ph.ncontrols, "nparameters": ph.nparameters,
        "nevents": ph.nevents, "npath": ph.npath, "nodes": list(ph.nodes),
        "nalgebraic": getattr(ph, "nalgebraic", 0),
        "states_lower": list(map(float, ph.bounds.lower.states)),
        "states_upper": list(map(float, ph.bounds.upper.states)),
        "controls_lower": list(map(float, ph.bounds.lower.controls or [])),
        "controls_upper": list(map(float, ph.bounds.upper.controls or [])),
        "events_lower": list(map(float, ph.bounds.lower.events or [])),
        "events_upper": list(map(float, ph.bounds.upper.events or [])),
        "path_lower": list(map(float, ph.bounds.lower.path or [])),
        "path_upper": list(map(float, ph.bounds.upper.path or [])),
        "t0_lower": float(ph.bounds.t0[0]), "t0_upper": float(ph.bounds.t0[1]),
        "tf_lower": float(ph.bounds.tf[0]), "tf_upper": float(ph.bounds.tf[1]),
        "nobserved": ph.nobserved, "nsamples": ph.nsamples,
        "parameters_lower": list(map(float, ph.bounds.lower.parameters or [])),
        "parameters_upper": list(map(float, ph.bounds.upper.parameters or [])),
        "guess_states": _col(ph.guess.states),
        "guess_controls": _col(ph.guess.controls) if ph.guess.controls is not None else _col([[]]),
        "guess_time": _col(ph.guess.time),
        "guess_parameters": _col(ph.guess.parameters) if ph.guess.parameters is not None else _col([[]]),
        "observation_nodes": _col(ph.observation_nodes) if ph.observation_nodes is not None else _col([[]]),
        "observations": _col(ph.observations) if ph.observations is not None else _col([[]]),
        "residual_weights": (_col(ph.residual_weights)
                             if getattr(ph, "residual_weights", None) is not None
                             else _col([[]])),
        "regularization_factor": (float(ph.regularization_factor)
                                  if getattr(ph, "regularization_factor", None) is not None
                                  else 0.0),
        "integer_controls": [{"index": i, "values": v} for (i, v) in ph.integer_controls],
        "integer_parameters": [{"index": i, "values": v} for (i, v) in ph.integer_parameters],
    }


def _alg_dict(a):
    d = {"collocation_method": a.collocation_method, "nlp_method": a.nlp_method,
         "derivatives": a.derivatives, "scaling": a.scaling,
         "nlp_iter_max": a.nlp_iter_max, "nlp_tolerance": a.nlp_tolerance}
    optional = ["print_level", "hessian", "ipopt_linear_solver", "constraint_scaling",
                "ipopt_max_cpu_time", "diagnostic_level",
                "mesh_refinement", "mr_max_iterations", "ode_tolerance",
                "mr_max_growth_factor", "mr_min_order", "mr_max_order",
                "transcription_method", "ir_residual_nodes", "ir_regularization",
                "ir_objective", "ir_residual_bound", "ir_dair",
                "ir_dair_delta_factor", "ir_local_order",
                "ir_include_path", "ir_path_weight", "ir_residual_scaling",
                "ir_element_local_controls", "ir_flexible_mesh",
                "ir_min_element_fraction",
                "ms_steps_per_segment", "ms_integrator",
                "ms_control_parameterisation",
                "ms_path_samples", "ms_flexible_segments",
                "ms_min_segment_fraction", "ms_refine_tolerance",
                "ms_algebraic_iterations", "ms_adaptive_steps",
                "ms_max_steps_per_segment", "ms_implicit_iterations",
                "qp_solver", "qp_restoration", "sqp_strategy", "qp_iter_max",
                "trust_region", "trust_region_radius", "elastic_penalty",
                "parameter_statistics", "parameter_estimation_norm",
                "objective_form", "defect_scaling", "diff_matrix",
                "jac_sparsity_ratio", "hess_sparsity_ratio",
                "save_sparsity_pattern", "nsteps_error_integration",
                "mr_kappa", "mr_M1", "mr_switch_detection", "switch_order",
                "hessian_verify", "on_error", "max_integer_combinations"]
    for k in optional:
        v = getattr(a, k, None)
        if v is not None:
            d[k] = v
    return d


class IntegerControlResult:
    """Rounded integer control recovered from the relaxed weights by sum-up rounding.

    control        rounded value on each mesh interval
    mode_index     index into the admissible set on each interval
    interval_widths
    integral_gap   accumulated sum-up-rounding gap
    n_switches     instants at which any declared integer control changes
    """
    def __init__(self, d):
        self.control = np.asarray(d["control"]).ravel()
        self.mode_index = np.asarray(d["mode_index"]).ravel().astype(int)
        self.interval_widths = np.asarray(d["interval_widths"]).ravel()
        self.integral_gap = float(d["integral_gap"])
        self.n_switches = int(d["n_switches"])


class _Status:
    """What the solve did, as opposed to what it found.

    ``objective`` comes back whatever happened, so a script that reads only the
    objective cannot tell a converged solve from one that hit its iteration limit.
    These four say which, and ``success`` is the one-line answer.

    nlp_return_code  the NLP solver's own code. From IPOPT, 0 is "solved" and 1 is
                     "solved to acceptable level"; both are successes. From PSOPT's
                     own SQP, 0 is converged and 1 is the iteration limit reached,
                     which is NOT a success -- SQP_interface reclassifies a
                     budget-exhausted run that is nonetheless acceptable to 0 before
                     returning, so a 1 that survives means the result was not.
    error_flag       non-zero for a set-up failure or a thrown exception, which is a
                     different thing from the NLP not converging.
    error_msg        what that failure was, empty when there was none.
    cpu_time         seconds.
    """
    def __init__(self, d, nlp_method):
        self.nlp_return_code = int(d.get("nlp_return_code", 0))
        self.error_flag = int(d.get("error_flag", 0))
        self.error_msg = d.get("error_msg", "")
        self.cpu_time = float(d.get("cpu_time", 0.0))
        self.mesh_refinement_iterations = int(d.get("mesh_refinement_iterations", 0))
        self.mesh_stats = list(d.get("mesh_stats", []))
        self._nlp_method = nlp_method

    @property
    def success(self):
        if self.error_flag != 0:
            return False
        if self.nlp_return_code == 0:
            return True
        # See the note above: code 1 means different things in the two solvers.
        return self.nlp_return_code == 1 and str(self._nlp_method).upper() == "IPOPT"

    def __repr__(self):
        # Deliberately short and on one line. The error message can be several
        # lines of PSOPT's own banner, and putting it in here made the repr wider
        # than a page -- which is visible in any transcript that prints it, and
        # was visible in the application examples document, where the line ran
        # outside its box. It is still on .error_msg.
        return ("<PSOPT %s rc=%d flag=%d cpu=%.3gs>"
                % ("converged" if self.success else "NOT converged",
                   self.nlp_return_code, self.error_flag, self.cpu_time))


class _PhaseDuals:
    """The multipliers and diagnostics of one phase.

    costates              discrete adjoint, one row per state. This is what a solution
                          is checked against the maximum principle with.
    hamiltonian           the Hamiltonian along the trajectory; constant on an
                          autonomous problem with free final time, which is a useful
                          independent check on a converged solve.
    dual_path             multiplier of each path constraint, when the phase has any.
    dual_events           multiplier of each event constraint, when the phase has any.
    terminal_state        the state at the final node.
    terminal_costate      the costate there, which is the transversality condition.
    relative_local_error  relative local discretisation error per mesh interval, which
                          is what mesh refinement drives down.
    """
    def __init__(self, d):
        def g(k):
            return np.asarray(d[k]) if k in d else None
        self.costates = g("costates")
        self.hamiltonian = g("hamiltonian")
        self.dual_path = g("dual_path")
        self.dual_events = g("dual_events")
        self.terminal_state = g("terminal_state")
        self.terminal_costate = g("terminal_costate")
        self.relative_local_error = g("relative_local_error")


class ParameterStatistics:
    """Covariance and confidence intervals of an estimated parameter vector.

    Present only when algorithm.parameter_statistics = "yes" was asked for and PSOPT
    could form them; otherwise Solution.parameter_statistics is None. A covariance
    that could not be formed and one that happens to be zero are different things, so
    this object's absence is the solver's own verdict rather than an inference.
    """
    def __init__(self, d):
        self.covariance = np.asarray(d["covariance"])
        self.confidence_low = np.asarray(d["confidence_low"]).ravel()
        self.confidence_high = np.asarray(d["confidence_high"]).ravel()
        self.residuals = np.asarray(d["residuals"]).ravel()
        self.sigma_hat = float(d["sigma_hat"])

    @property
    def standard_errors(self):
        """Square roots of the diagonal of the covariance."""
        return np.sqrt(np.clip(np.diag(np.atleast_2d(self.covariance)), 0.0, None))


class Solution:
    def __init__(self, d, nlp_method="IPOPT"):
        self.objective = d["objective"]
        self.states = np.asarray(d["states"])
        self.controls = np.asarray(d["controls"]) if "controls" in d else None
        self.time = np.asarray(d["time"]).ravel()
        self.parameters = np.asarray(d["parameters"]).ravel() if "parameters" in d else None
        self.status = _Status(d, nlp_method)
        self.duals = _PhaseDuals(d)
        self.costates = self.duals.costates
        self.parameter_statistics = (ParameterStatistics(d["parameter_statistics"])
                                     if "parameter_statistics" in d else None)
        # One entry per declared integer control, in declaration order. Note that
        # self.controls is in the weights layout when integer controls are declared:
        # the trailing rows are the product-mode weights, and the rounded controls are
        # here rather than there.
        self.integer_controls = [IntegerControlResult(r)
                                 for r in d.get("integer_controls", [])]
        # Selected value of each declared integer parameter, in declaration order.
        self.integer_parameters = [dict(index=int(r["index"]),
                                        value=float(r["value"]),
                                        mode_index=int(r["mode_index"]))
                                   for r in d.get("integer_parameters", [])]


class MultiSolution:
    def __init__(self, d, nlp_method="IPOPT"):
        self.objective = d["objective"]
        self.states = [np.asarray(s) for s in d["states"]]
        self.controls = [np.asarray(c) for c in d["controls"]]
        self.time = [np.asarray(t).ravel() for t in d["time"]]
        self.parameters = [None if p is None else np.asarray(p).ravel()
                           for p in d.get("parameters", [])]
        self.status = _Status(d, nlp_method)
        # One _PhaseDuals per phase, in phase order.
        self.duals = [_PhaseDuals(x) for x in d.get("duals", [])]
        self.costates = [x.costates for x in self.duals]
        self.dual_linkages = (np.asarray(d["dual_linkages"]).ravel()
                              if "dual_linkages" in d else None)
        self.parameter_statistics = (ParameterStatistics(d["parameter_statistics"])
                                     if "parameter_statistics" in d else None)
        # Per phase, one entry per declared integer control / parameter.
        self.integer_controls = [[IntegerControlResult(r) for r in ph]
                                 for ph in d.get("integer_controls", [])]
        self.integer_parameters = [[dict(index=int(r["index"]), value=float(r["value"]),
                                         mode_index=int(r["mode_index"])) for r in ph]
                                   for ph in d.get("integer_parameters", [])]


class Problem:
    def __init__(self, name="problem"):
        self.name = name
        self._phases = []
        self._links = []

    def add_phase(self, nstates, ncontrols, nevents=0, npath=0, nparameters=0,
                  nobserved=0, nsamples=0, nalgebraic=0):
        ph = Phase(nstates, ncontrols, nevents, npath, nparameters, nobserved, nsamples,
                   nalgebraic)
        self._phases.append(ph)
        return ph

    def link_phases(self, a, b, jumps=None):
        """Auto-link all states/time between phases a and b; optional jumps={state: delta}
        subtract delta from that state's continuity residual (e.g. mass jettison)."""
        self._links.append(dict(a=a, b=b, jumps=dict(jumps or {})))

    def solve(self, algorithm):
        if len(self._phases) == 1 and not self._links:
            return self._solve_single(self._phases[0], algorithm)
        return self._solve_multi(algorithm)

    @staticmethod
    def _check(sol):
        """Raise if the solve reported a set-up failure; otherwise hand it back.

        error_flag is PSOPT's own distinction: non-zero means the problem was
        refused or something threw, and zero means the solve ran -- whether or not
        the NLP converged, which sol.status.success reports separately.
        """
        if sol.status.error_flag != 0:
            full = sol.status.error_msg or ""
            # PSOPT wraps its diagnostic in a banner that tells a C++ user to put a
            # breakpoint on error_message() and take a backtrace. That is not advice a
            # Python caller can act on, and in a traceback it buries the one line that
            # matters, so the exception carries the diagnostic itself and keeps the
            # banner on .details.
            core = full
            if "====>" in full and "<====" in full:
                core = full.split("====>", 1)[1].rsplit("<====", 1)[0]
            core = " ".join(core.split()) or (
                "PSOPT reported error_flag %d with no message" % sol.status.error_flag)
            err = PSOPTError(core, sol)
            err.details = full.strip()
            raise err
        return sol

    def _solve_single(self, ph, algorithm):
        dae_f, L_f, phi_f, ev_f, obs_f = ph._build_functions()
        dims = dict(nx=ph.nstates, nu=ph.ncontrols, npar=ph.nparameters,
                    npath=ph.npath, nevents=ph.nevents, nobs=ph.nobserved)
        so = codegen.compile_single_phase(self.name, dims, dae_f, L_f, phi_f, ev_f, obs_f=obs_f)
        spec = {"outfilename": self.name + ".txt", "so_path": so}
        spec.update(_phase_dict(ph))
        spec["algorithm"] = _alg_dict(algorithm)
        return self._check(Solution(_psopt.solve_single_phase(spec),
                                    algorithm.nlp_method))

    def _solve_multi(self, algorithm):
        phase_funcs, nstates_by_phase = [], []
        for ph in self._phases:
            dae_f, L_f, phi_f, ev_f, _obs = ph._build_functions()
            phase_funcs.append(dict(nx=ph.nstates, npath=ph.npath, nevents=ph.nevents,
                                    dae_f=dae_f, L_f=L_f, phi_f=phi_f, ev_f=ev_f))
            nstates_by_phase.append(ph.nstates)
        so = codegen.compile_multiphase(self.name, phase_funcs, self._links, nstates_by_phase)
        nlinkages = sum(nstates_by_phase[lk["b"] - 1] + 1 for lk in self._links)
        spec = {"outfilename": self.name + ".txt", "so_path": so,
                "nphases": len(self._phases), "nlinkages": nlinkages,
                "phases": [_phase_dict(ph) for ph in self._phases],
                "algorithm": _alg_dict(algorithm)}
        return self._check(MultiSolution(_psopt.solve_multiphase(spec),
                                         algorithm.nlp_method))
