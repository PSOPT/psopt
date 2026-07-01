"""bryson_denham with hp-adaptive mesh refinement via the Python API (B-4).
Starts from a coarse 10-node mesh; PSOPT refines to meet ode_tolerance.
Native reference (same options): 7 iterations, 10->29 nodes, objective 3.999997e+00."""
import sys, numpy as np, casadi as ca
sys.path.insert(0, "/tmp/psopt_py")
import psopt

prob = psopt.Problem(name="bryson_mesh")
ph = prob.add_phase(nstates=3, ncontrols=1, nevents=5, npath=0)
ph.nodes = [10]                                   # coarse start mesh

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

ph.guess.states   = psopt.ramp([(0.0, 0.0), (1.0, -1.0), (0.0, 0.0)], 10)
ph.guess.controls = psopt.const([0.0], 10)
ph.guess.time     = psopt.times(0.0, 0.5, 10)

alg = psopt.Algorithm(
    collocation_method="Legendre",
    mesh_refinement="automatic", mr_max_iterations=7, ode_tolerance=1.0e-6,
    print_level=0,                                # quiet mode
)
sol = prob.solve(alg)

print("\n================ RESULT ================")
print("objective       :", repr(sol.objective))
print("native ref      : 3.999997e+00")
print("final nodes     :", sol.time.shape[0], "(native ref: 29)")
