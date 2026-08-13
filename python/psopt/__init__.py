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
                 nobserved=0, nsamples=0):
        self.nstates = nstates
        self.ncontrols = ncontrols
        self.nevents = nevents
        self.npath = npath
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
                 ir_dair_delta_factor=None, ir_local_order=None):
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


def _col(a):
    return np.ascontiguousarray(np.atleast_2d(np.asarray(a, dtype=float)))


def _phase_dict(ph):
    return {
        "nstates": ph.nstates, "ncontrols": ph.ncontrols, "nparameters": ph.nparameters,
        "nevents": ph.nevents, "npath": ph.npath, "nodes": list(ph.nodes),
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
                "ir_dair_delta_factor", "ir_local_order"]
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


class Solution:
    def __init__(self, d):
        self.objective = d["objective"]
        self.states = np.asarray(d["states"])
        self.controls = np.asarray(d["controls"]) if "controls" in d else None
        self.time = np.asarray(d["time"]).ravel()
        self.parameters = np.asarray(d["parameters"]).ravel() if "parameters" in d else None
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
    def __init__(self, d):
        self.objective = d["objective"]
        self.states = [np.asarray(s) for s in d["states"]]
        self.controls = [np.asarray(c) for c in d["controls"]]
        self.time = [np.asarray(t).ravel() for t in d["time"]]
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
                  nobserved=0, nsamples=0):
        ph = Phase(nstates, ncontrols, nevents, npath, nparameters, nobserved, nsamples)
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

    def _solve_single(self, ph, algorithm):
        dae_f, L_f, phi_f, ev_f, obs_f = ph._build_functions()
        dims = dict(nx=ph.nstates, nu=ph.ncontrols, npar=ph.nparameters,
                    npath=ph.npath, nevents=ph.nevents, nobs=ph.nobserved)
        so = codegen.compile_single_phase(self.name, dims, dae_f, L_f, phi_f, ev_f, obs_f=obs_f)
        spec = {"outfilename": self.name + ".txt", "so_path": so}
        spec.update(_phase_dict(ph))
        spec["algorithm"] = _alg_dict(algorithm)
        return Solution(_psopt.solve_single_phase(spec))

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
        return MultiSolution(_psopt.solve_multiphase(spec))
