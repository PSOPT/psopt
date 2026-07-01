"""cracking (parameter estimation) via the PSOPT Python API (B-3 acceptance test).
3 parameters in the dynamics, 0 controls, an observation function fitting 21 noisy
samples of both states. Target objective (least-squares residual): 4.319519e-03."""
import sys, numpy as np, casadi as ca
sys.path.insert(0, "/tmp/psopt_py")
import psopt

# ---- measured data (from cracking.cxx) ----
y1meas = [1.0,0.8105,0.6208,0.5258,0.4345,0.3903,0.3342,0.3034,0.2735,0.2405,0.2283,
          0.2071,0.1669,0.153,0.1339,0.1265,0.12,0.099,0.087,0.077,0.069]
y2meas = [0.0,0.2,0.2886,0.301,0.3215,0.3123,0.2716,0.2551,0.2258,0.1959,0.1789,
          0.1457,0.1198,0.0909,0.0719,0.0561,0.046,0.028,0.019,0.014,0.01]
tmeas  = [0.0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,
          0.3,0.35,0.4,0.45,0.5,0.55,0.65,0.75,0.85,0.95]

prob = psopt.Problem(name="cracking")
ph = prob.add_phase(nstates=2, ncontrols=0, nparameters=3, nobserved=2, nsamples=21)
ph.nodes = [80]

# dynamics: parameters theta1..theta3 drive the kinetics
ph.dynamics    = lambda x, u, p, t: ca.vertcat(-(p[0] + p[2]) * x[0]**2,
                                               p[0] * x[0]**2 - p[1] * x[1])
ph.observation = lambda x, u, p, t: ca.vertcat(x[0], x[1])     # observe both states

ph.bounds.lower.states = [0.0, 0.0]
ph.bounds.upper.states = [2.0, 2.0]
ph.bounds.lower.parameters = [0.0, 0.0, 0.0]
ph.bounds.upper.parameters = [20.0, 20.0, 20.0]
ph.bounds.t0 = (0.0, 0.0)
ph.bounds.tf = (0.95, 0.95)

ph.observation_nodes = np.array(tmeas).reshape(1, 21)
ph.observations      = np.vstack([y1meas, y2meas])

sg = np.zeros((2, 40))
sg[0, :] = np.linspace(1.0, 0.069, 40)
sg[1, :] = np.linspace(0.30, 0.01, 40)
ph.guess.states     = sg
ph.guess.time       = np.linspace(0.0, 0.95, 40).reshape(1, 40)
ph.guess.parameters = np.zeros((3, 1))

alg = psopt.Algorithm(collocation_method="Hermite-Simpson", nlp_iter_max=1000, nlp_tolerance=1e-6)
sol = prob.solve(alg)

print("\n================ RESULT ================")
print("objective    :", repr(sol.objective))
print("target (C++) : 4.319519e-03")
print("parameters   :", sol.parameters)
