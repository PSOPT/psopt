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


class Solution:
    def __init__(self, d):
        self.objective = d["objective"]
        self.states = np.asarray(d["states"])
        self.controls = np.asarray(d["controls"]) if "controls" in d else None
        self.time = np.asarray(d["time"]).ravel()
        self.parameters = np.asarray(d["parameters"]).ravel() if "parameters" in d else None


class MultiSolution:
    def __init__(self, d):
        self.objective = d["objective"]
        self.states = [np.asarray(s) for s in d["states"]]
        self.controls = [np.asarray(c) for c in d["controls"]]
        self.time = [np.asarray(t).ravel() for t in d["time"]]


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
