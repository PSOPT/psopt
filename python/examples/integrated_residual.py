"""A singular arc, and why integrated-residual transcription exists.

    minimise  integral_0^tf (1 + x2^2)/2 dt
    subject to  x1' = x2,  x2' = u,  |u| <= 1,
                (x1,x2)(0) = (2,0),  (x1,x2)(tf) = (5,0),  tf free.

Locatelli, *Optimal Control of a Double Integrator*, problem 11.2-1a. The
Hamiltonian is linear in u, so the control is bang-bang wherever the switching
function is non-zero, and singular where it vanishes on an interval. The exact
solution is u = +1 on [0,1], u = 0 on the singular arc [1,3] where x2 is held at
1, and u = -1 on [3,4]; tf* = 4 and J* = 10/3.

Collocation approximates the control by a polynomial through the mesh, which a
bang-bang control is not, so the answer is right and the control oscillates
around the switches. Integrated-residual transcription minimises the integral of
the squared dynamics residual rather than forcing it to zero at collocation
points, which does not ask the control to be smooth, and the arcs come out flat.

Selecting it from Python is transcription_method="integrated-residual"; the ir_*
family configures it. Two of those settings are worth knowing before the first
attempt. It is built on the Hermite-Simpson mesh and PSOPT refuses any other
collocation_method. And ir_objective chooses what is minimised: the default
"residual" reports the integrated residual, a feasibility measure rather than the
control cost, while "cost" minimises the problem's own objective and then needs
the dynamics enforced some other way -- either ir_residual_bound >= 0, the robust
constraint form used here, or a large enough ir_regularization for the penalty
form. Asking for "cost" with neither leaves the dynamics essentially unenforced,
and the answer is then whatever the bounds allow. This example solves the problem both ways and compares the
control against the exact bang-singular-bang solution.
"""
import numpy as np
import casadi as ca
from _common import psopt, report

EXACT_J = 10.0 / 3.0
EXACT_TF = 4.0
N = 60


def build():
    prob = psopt.Problem(name="singular_di")
    ph = prob.add_phase(nstates=2, ncontrols=1, nevents=4)
    ph.nodes = [N]
    ph.dynamics = lambda x, u, p, t: ca.vertcat(x[1], u[0])
    ph.integrand = lambda x, u, p, t: 0.5 * (1.0 + x[1] ** 2)
    ph.events = lambda xi, xf, p, t0, tf: ca.vertcat(xi[0], xi[1], xf[0], xf[1])
    ph.bounds.lower.states = [-10.0, -5.0]
    ph.bounds.upper.states = [10.0, 5.0]
    ph.bounds.lower.controls = [-1.0]
    ph.bounds.upper.controls = [1.0]
    ph.bounds.lower.events = [2.0, 0.0, 5.0, 0.0]
    ph.bounds.upper.events = [2.0, 0.0, 5.0, 0.0]
    ph.bounds.t0 = (0.0, 0.0)
    ph.bounds.tf = (2.0, 8.0)
    ph.guess.states = np.vstack([np.linspace(2.0, 5.0, N), np.ones(N)])
    ph.guess.controls = np.zeros((1, N))
    ph.guess.time = np.linspace(0.0, 4.0, N).reshape(1, N)
    return prob


def exact_u(t):
    return np.where(t < 1.0, 1.0, np.where(t < 3.0, 0.0, -1.0))


print("Singular-arc double integrator (Locatelli 11.2-1a).")
print("Exact: u = +1 on [0,1], 0 on [1,3], -1 on [3,4]; tf* = 4, J* = 10/3.\n")

results = {}
for label, alg in (
    ("Legendre collocation",
     psopt.Algorithm(collocation_method="Legendre", nlp_tolerance=1.0e-6,
                     nlp_iter_max=2000, print_level=0)),
    # Integrated-residual transcription is built on the Hermite-Simpson mesh and
    # PSOPT says so if it is asked for anything else.
    ("Integrated residual",
     psopt.Algorithm(collocation_method="Hermite-Simpson", nlp_tolerance=1.0e-6,
                     nlp_iter_max=2000, print_level=0,
                     transcription_method="integrated-residual",
                     ir_residual_nodes=3, ir_objective="cost",
                     ir_residual_bound=1.0e-9)),
):
    sol = build().solve(alg)
    results[label] = sol

ok = True
for label, sol in results.items():
    ok = report(sol, label, reference=EXACT_J, tol=2.0e-3) and ok
    t, u = sol.time, sol.controls[0]
    # How far the computed control is from the exact one, away from the two
    # switching instants where any discretisation is entitled to disagree.
    away = (np.abs(t - 1.0) > 0.15) & (np.abs(t - 3.0) > 0.15)
    err = np.max(np.abs(u[away] - exact_u(t[away])))
    # Flatness of the singular arc, which is what the two methods differ on.
    sing = (t > 1.2) & (t < 2.8)
    print("final time     : %.6f   (exact %.1f)" % (t[-1], EXACT_TF))
    print("max |u - u*|   : %.3e   (away from the two switches)" % err)
    print("singular arc   : u ranges over %.3e on 1.2 < t < 2.8" % np.ptp(u[sing]))

print("\nThe integrated-residual control is flatter on the singular arc; both")
print("objectives are within a fraction of a per cent of 10/3.")
print("\nRESULT: %s" % ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
