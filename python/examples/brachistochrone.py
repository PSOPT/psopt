"""Brachistochrone -- the minimum-time descent curve, and what the costates say.

The bead slides without friction from (0,0) to x = 2, y = 2 under gravity; the
control is the path angle theta. Minimising the final time gives the cycloid
Johann Bernoulli asked for in 1696.

This is the Python port of examples/brac1 in the C++ distribution, whose answer
is tf = 8.247591e-01 s. It is also the example to read for what a solve returns
beyond the trajectory, because a minimum-time problem has a sharp check on it:
with the objective written as the Mayer term tf, the Hamiltonian along an optimal
trajectory of an autonomous free-final-time problem is identically -1. If the
Hamiltonian PSOPT reports is not flat and not -1, the solve has not converged
however plausible the objective looks.
"""
import numpy as np
import casadi as ca
from _common import psopt, report

N = 40
prob = psopt.Problem(name="brachistochrone")
ph = prob.add_phase(nstates=3, ncontrols=1, nevents=5)
ph.nodes = [N]

# x' = v sin(theta),  y' = v cos(theta),  v' = g cos(theta)
ph.dynamics = lambda x, u, p, t: ca.vertcat(x[2] * ca.sin(u[0]),
                                            x[2] * ca.cos(u[0]),
                                            9.8 * ca.cos(u[0]))
ph.endpoint = lambda xi, xf, p, t0, tf: tf                  # minimise the final time
ph.events = lambda xi, xf, p, t0, tf: ca.vertcat(xi[0], xi[1], xi[2], xf[0], xf[1])

ph.bounds.lower.states = [0.0, 0.0, 0.0]
ph.bounds.upper.states = [20.0, 20.0, 20.0]
ph.bounds.lower.controls = [0.0]
ph.bounds.upper.controls = [2.0 * np.pi]
ph.bounds.lower.events = [0.0, 0.0, 0.0, 2.0, 2.0]
ph.bounds.upper.events = [0.0, 0.0, 0.0, 2.0, 2.0]
ph.bounds.t0 = (0.0, 0.0)
ph.bounds.tf = (0.0, 10.0)                                   # the final time is free

ph.guess.states = np.vstack([np.linspace(0.0, 2.0, N),
                             np.linspace(0.0, 2.0, N),
                             np.linspace(0.0, 2.0, N)])
ph.guess.controls = np.ones((1, N))
ph.guess.time = np.linspace(0.0, 2.0, N).reshape(1, N)

alg = psopt.Algorithm(collocation_method="Legendre", nlp_tolerance=1.0e-6,
                      nlp_iter_max=1000, print_level=0)
sol = prob.solve(alg)

ok = report(sol, "Brachistochrone (C++ example brac1)", reference=8.247591e-01, tol=1.0e-5)

# ---- what the solve produced beyond the trajectory ------------------------------
print("final time     : %.9g s" % sol.time[-1])
print("costates       : %s   (one row per state)" % (sol.costates.shape,))

# lambda_x and lambda_y are constant here: neither x nor y appears in the dynamics
# or the cost, so their adjoint equations are lambda' = 0. That is a statement about
# the problem, not about PSOPT, which makes it a fair check on the solver.
lam = sol.costates
print("lambda_x       : %+.6e ... %+.6e" % (lam[0, 0], lam[0, -1]))
print("lambda_y       : %+.6e ... %+.6e" % (lam[1, 0], lam[1, -1]))
print("                 both constant: neither x nor y appears in f or L")

H = sol.duals.hamiltonian.ravel()
print("Hamiltonian    : mean %+.6f, spread %.2e" % (H.mean(), H.max() - H.min()))
print("                 theory: identically -1")
# The gate is 1e-3 rather than the NLP tolerance. The Hamiltonian is reconstructed
# from the discrete adjoint on 40 nodes, so what it agrees with -1 to is set by the
# discretisation, not by how tightly the NLP was solved. On this mesh it comes out
# at about 1.4e-4, and refining the mesh is what moves it.
ok = ok and abs(H.mean() + 1.0) < 1.0e-3 and (H.max() - H.min()) < 5.0e-3

# The cycloid through (0,0) and (2,2) can be solved for directly, which gives an
# answer that owes nothing to PSOPT. With x = a(theta - sin theta),
# y = a(1 - cos theta), the endpoint condition x = y fixes theta and then a, and
# tf = sqrt(a/g) * theta.
from scipy.optimize import brentq                                   # noqa: E402
th = brentq(lambda s: (s - np.sin(s)) - (1.0 - np.cos(s)), 1.0, 2.0 * np.pi - 1.0e-9)
a = 2.0 / (1.0 - np.cos(th))
tf_exact = np.sqrt(a / 9.8) * th
print("closed form    : tf = %.9g s   (cycloid, (0,0) to (2,2))" % tf_exact)
print("difference     : %.2e" % abs(sol.time[-1] - tf_exact))
ok = ok and abs(sol.time[-1] - tf_exact) < 1.0e-6

print("\nRESULT: %s" % ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
