"""bryson_denham with integrated-residual (IR) transcription via the Python API (B-4).

NOTE: IR/NK are work-in-progress engine features with no committed C++ example. This
test validates that the Python binding drives the IR machinery *identically* to a native
driver with the same options (pass-through fidelity), not a published value. With the
default ir_objective="residual" the reported objective is the integrated residual
(a feasibility measure), not the control cost.

Native reference (Hermite-Simpson + integrated-residual): 5.049038e-06."""
import sys, numpy as np, casadi as ca
sys.path.insert(0, "/tmp/psopt_py")
import psopt

prob = psopt.Problem(name="bryson_ir")
ph = prob.add_phase(nstates=3, ncontrols=1, nevents=5, npath=0)
ph.nodes = [50]

ph.dynamics = lambda x, u, p, t: ca.vertcat(x[1], u[0], 0.5 * u[0]**2)
ph.endpoint = lambda xi, xf, p, t0, tf: xf[2]
ph.events   = lambda xi, xf, p, t0, tf: ca.vertcat(xi[0], xi[1], xi[2], xf[0], xf[1])

ph.bounds.lower.states   = [0.0, -10.0, -10.0]
ph.bounds.upper.states   = [1.0/9.0, 10.0, 10.0]
ph.bounds.lower.controls = [-10.0]
ph.bounds.upper.controls = [10.0]
ph.bounds.lower.events   = [0.0, 1.0, 0.0, 0.0, -1.0]
ph.bounds.upper.events   = [0.0, 1.0, 0.0, 0.0, -1.0]
ph.bounds.t0 = (0.0, 0.0)
ph.bounds.tf = (0.0, 50.0)

ph.guess.states   = psopt.ramp([(0.0, 0.0), (1.0, -1.0), (0.0, 0.0)], 50)
ph.guess.controls = psopt.const([0.0], 50)
ph.guess.time     = psopt.times(0.0, 0.5, 50)

alg = psopt.Algorithm(
    collocation_method="Hermite-Simpson",
    transcription_method="integrated-residual",   # IR; ir_objective defaults to "residual"
    print_level=0,
)
sol = prob.solve(alg)

print("\n================ RESULT ================")
print("objective (integrated residual) :", repr(sol.objective))
print("native reference                : 5.049038e-06")
