"""Breakwell's problem -- a state constraint with an exact answer.

    minimise  1/2 * integral_0^1 u^2 dt
    subject to  x' = v,  v' = u,
                x(0) = 0, v(0) = 1,  x(1) = 0, v(1) = -1,
                x(t) <= l  with  l = 0.1.

The constraint bites: the unconstrained optimum reaches x = 1/4. With l < 1/6 the
solution has a boundary arc and the optimal cost is 4/(9l) exactly -- here 40/9 =
4.444444..., which is the reference this example checks against and which owes
nothing to PSOPT.

This is the Python port of examples/breakwell in the C++ distribution. It is the
one to read for a pure state bound (imposed through the state bounds rather than
as a path constraint) and for what the costates do on a constrained arc: lambda_x
is piecewise constant, jumping where the trajectory enters and leaves the boundary.
"""
import numpy as np
import casadi as ca
from _common import psopt, report

N = 200
L = 0.1
prob = psopt.Problem(name="breakwell")
ph = prob.add_phase(nstates=2, ncontrols=1, nevents=4)
ph.nodes = [N]

ph.dynamics = lambda x, u, p, t: ca.vertcat(x[1], u[0])
ph.integrand = lambda x, u, p, t: 0.5 * u[0] ** 2
ph.events = lambda xi, xf, p, t0, tf: ca.vertcat(xi[0], xi[1], xf[0], xf[1])

ph.bounds.lower.states = [-2.0, -2.0]
ph.bounds.upper.states = [L, 2.0]                 # the state constraint x <= 0.1
ph.bounds.lower.controls = [-10.0]
ph.bounds.upper.controls = [10.0]
ph.bounds.lower.events = [0.0, 1.0, 0.0, -1.0]
ph.bounds.upper.events = [0.0, 1.0, 0.0, -1.0]
ph.bounds.t0 = (0.0, 0.0)
ph.bounds.tf = (1.0, 1.0)

ph.guess.states = np.vstack([np.zeros(N), np.ones(N)])
ph.guess.controls = np.zeros((1, N))
ph.guess.time = np.linspace(0.0, 1.0, N).reshape(1, N)

alg = psopt.Algorithm(collocation_method="Legendre", nlp_tolerance=1.0e-6,
                      nlp_iter_max=1000, print_level=0)
sol = prob.solve(alg)

exact = 4.0 / (9.0 * L)
ok = report(sol, "Breakwell (C++ example breakwell)", reference=exact, tol=1.0e-5)

x = sol.states[0]
print("max x(t)       : %.9f   (the bound is %.3f)" % (x.max(), L))
print("time on bound  : %.3f of the horizon" % (np.mean(x > L - 1.0e-6),))

# On a boundary arc the state constraint is active and its multiplier is non-zero,
# so lambda_x is constant on each free arc and changes between them. Reading the
# three levels off the costate is the clearest illustration of that structure the
# problem offers.
lam_x = sol.costates[0]
print("lambda_x       : first %+.4f, middle %+.4f, last %+.4f"
      % (lam_x[0], lam_x[N // 2], lam_x[-1]))
print("               : constant on each arc; it changes where the bound is met")

print("\nRESULT: %s" % ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
