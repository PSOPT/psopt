"""The hypersensitive problem -- and what mesh refinement does to it.

    minimise  1/2 * integral_0^50 (x^2 + u^2) dt
    subject to  x' = -x^3 + u,  x(0) = 1.5,  x(50) = 1.

The name is Rao and Mease's. The solution has boundary layers at each end that
decay over an O(1) interval and a long, almost constant interior arc, so a uniform
mesh spends nearly all of its nodes where nothing happens and resolves neither end.
It is the standard test of an adaptive mesh.

This is the Python port of examples/hypersensitive in the C++ distribution, whose
answer is 1.330826e+00. It is the example to read for hp-adaptive mesh refinement
from Python, and for sol.status.mesh_stats, which is the per-iteration table PSOPT
prints at the end of a run.
"""
import numpy as np
import casadi as ca
from _common import psopt, report

prob = psopt.Problem(name="hypersensitive")
ph = prob.add_phase(nstates=1, ncontrols=1, nevents=2)
# Two entries: manual refinement solves on 25 nodes, then on 50, starting the
# second solve from the first solution.
ph.nodes = [25, 50]

ph.dynamics = lambda x, u, p, t: ca.vertcat(-x[0] ** 3 + u[0])
ph.integrand = lambda x, u, p, t: 0.5 * (x[0] ** 2 + u[0] ** 2)
ph.events = lambda xi, xf, p, t0, tf: ca.vertcat(xi[0], xf[0])

ph.bounds.lower.states = [-50.0]
ph.bounds.upper.states = [50.0]
ph.bounds.lower.controls = [-50.0]
ph.bounds.upper.controls = [50.0]
ph.bounds.lower.events = [1.5, 1.0]
ph.bounds.upper.events = [1.5, 1.0]
ph.bounds.t0 = (0.0, 0.0)
ph.bounds.tf = (50.0, 50.0)

N0 = 25
ph.guess.states = np.linspace(1.5, 1.0, N0).reshape(1, N0)
ph.guess.controls = np.zeros((1, N0))
ph.guess.time = np.linspace(0.0, 50.0, N0).reshape(1, N0)

alg = psopt.Algorithm(collocation_method="Legendre", nlp_tolerance=1.0e-6,
                      nlp_iter_max=1000, print_level=0)
sol = prob.solve(alg)

ok = report(sol, "Hypersensitive (C++ example hypersensitive)",
            reference=1.330826e+00, tol=1.0e-5)

print("mesh iterations: %d" % sol.status.mesh_refinement_iterations)
print("%-4s %-10s %6s %8s %12s %10s" % ("iter", "method", "nodes", "NLP its", "ODE error", "CPU (s)"))
for k, m in enumerate(sol.status.mesh_stats):
    print("%-4d %-10s %6d %8d %12.3e %10.3f"
          % (k + 1, m["method"], m["nnodes"], m["n_jacobian_evals"],
             m["epsilon_max"], m["cpu_time"]))

# Where the nodes ended up is the whole point: the interior arc needs almost none
# of them and the two boundary layers need nearly all.
t = sol.time
print("nodes in t<1   : %d" % int(np.sum(t < 1.0)))
print("nodes in 1..49 : %d" % int(np.sum((t >= 1.0) & (t <= 49.0))))
print("nodes in t>49  : %d" % int(np.sum(t > 49.0)))
print("x at the ends  : %.6f -> %.6f, interior minimum %.3e"
      % (sol.states[0, 0], sol.states[0, -1], np.abs(sol.states[0]).min()))

print("\nRESULT: %s" % ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
