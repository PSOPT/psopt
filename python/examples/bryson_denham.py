"""bryson_denham via the PSOPT Python API (B-1 acceptance test).
Target objective (matches the hand-written C++ example): 3.999539e+00."""
import sys, numpy as np, casadi as ca
sys.path.insert(0, "/tmp/psopt_py")
import psopt

prob = psopt.Problem(name="bryson_denham")
ph = prob.add_phase(nstates=3, ncontrols=1, nevents=5, npath=0)
ph.nodes = [50]

# user maths -- pure CasADi, traced and CppAD-taped under the hood
ph.dynamics = lambda x, u, p, t: ca.vertcat(x[1], u[0], 0.5 * u[0]**2)
ph.endpoint = lambda xi, xf, p, t0, tf: xf[2]                       # minimise x3(tf)
ph.events   = lambda xi, xf, p, t0, tf: ca.vertcat(xi[0], xi[1], xi[2], xf[0], xf[1])

ph.bounds.lower.states   = [0.0, -10.0, -10.0]
ph.bounds.upper.states   = [1.0/9.0, 10.0, 10.0]
ph.bounds.lower.controls = [-10.0]
ph.bounds.upper.controls = [10.0]
ph.bounds.lower.events   = [0.0, 1.0, 0.0, 0.0, -1.0]
ph.bounds.upper.events   = [0.0, 1.0, 0.0, 0.0, -1.0]
ph.bounds.t0 = (0.0, 0.0)
ph.bounds.tf = (0.0, 50.0)

x0 = np.zeros((3, 50))
x0[1, :] = np.linspace(1.0, -1.0, 50)
ph.guess.states   = x0
ph.guess.controls = np.zeros((1, 50))
ph.guess.time     = np.linspace(0.0, 0.5, 50).reshape(1, 50)

alg = psopt.Algorithm(collocation_method="Legendre", nlp_iter_max=1000, nlp_tolerance=1e-6)
sol = prob.solve(alg)

print("\n================ RESULT ================")
print("objective         :", repr(sol.objective))
print("target (C++)      : 3.999539e+00")
print("states shape      :", sol.states.shape)
print("controls shape    :", sol.controls.shape)
print("time samples      :", sol.time.shape)
