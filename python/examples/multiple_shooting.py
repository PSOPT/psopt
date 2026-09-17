"""Multiple shooting from Python, against an answer that is exact.

    minimise  1/2 * integral_0^1 u^2 dt
    subject to  x1' = x2,  x2' = u,  (x1,x2)(0) = (0,0),  (x1,x2)(1) = (1,0).

The minimum-energy double integrator. Pontryagin gives u*(t) = 6 - 12t and
J* = 6 exactly, so this example checks against arithmetic rather than against a
previous run of anything.

It solves the same problem twice, by collocation and by multiple shooting, to
show what selecting the second costs and what it changes. In multiple shooting
the mesh interval is not a collocation point but a *segment*: the NLP variables
are the state at the segment boundaries and the control parameters inside them,
an embedded integrator marches each segment, and the constraints are the
mismatches at the joins. The trajectory is therefore an integrated one by
construction, which is the property the method is chosen for.

The settings worth knowing are all here:

    transcription_method="multiple-shooting"   selects it
    ms_integrator                              the marching scheme (RK4 here)
    ms_steps_per_segment                       integrator steps between joins
    ms_control_parameterisation                "constant" or "linear" in a segment
    ms_path_samples                            where path constraints are imposed
"""
import numpy as np
import casadi as ca
from _common import psopt, report

EXACT = 6.0


def build(nodes):
    prob = psopt.Problem(name="mineng_di_ms")
    ph = prob.add_phase(nstates=2, ncontrols=1, nevents=4)
    ph.nodes = [nodes]
    ph.dynamics = lambda x, u, p, t: ca.vertcat(x[1], u[0])
    ph.integrand = lambda x, u, p, t: 0.5 * u[0] ** 2
    ph.events = lambda xi, xf, p, t0, tf: ca.vertcat(xi[0], xi[1], xf[0], xf[1])
    ph.bounds.lower.states = [-5.0, -5.0]
    ph.bounds.upper.states = [5.0, 5.0]
    ph.bounds.lower.controls = [-50.0]
    ph.bounds.upper.controls = [50.0]
    ph.bounds.lower.events = [0.0, 0.0, 1.0, 0.0]
    ph.bounds.upper.events = [0.0, 0.0, 1.0, 0.0]
    ph.bounds.t0 = (0.0, 0.0)
    ph.bounds.tf = (1.0, 1.0)
    ph.guess.states = np.zeros((2, nodes))
    ph.guess.controls = np.zeros((1, nodes))
    ph.guess.time = np.linspace(0.0, 1.0, nodes).reshape(1, nodes)
    return prob


print("Minimum-energy double integrator.  Exact optimum J* = 6, u*(t) = 6 - 12t.")

# ---- collocation, for comparison -------------------------------------------------
prob = build(40)
sol_c = prob.solve(psopt.Algorithm(collocation_method="Legendre", nlp_tolerance=1.0e-8,
                                   print_level=0))
ok = report(sol_c, "Legendre collocation, 40 nodes", reference=EXACT, tol=1.0e-8)
print("CPU            : %.3f s" % sol_c.status.cpu_time)

# ---- multiple shooting -----------------------------------------------------------
# 20 segments, four RK4 steps in each, with the control linear inside a segment.
# Linear rather than constant because the optimal control is linear in t here, and a
# piecewise-constant parameterisation would be the thing limiting the accuracy
# rather than the integrator.
prob = build(21)
sol_m = prob.solve(psopt.Algorithm(
    transcription_method="multiple-shooting",
    ms_integrator="RK4",
    ms_steps_per_segment=4,
    ms_control_parameterisation="linear",
    nlp_tolerance=1.0e-8, nlp_iter_max=1000, print_level=0))
ok = report(sol_m, "Multiple shooting, 20 segments x 4 RK4 steps",
            reference=EXACT, tol=1.0e-6) and ok
print("CPU            : %.3f s" % sol_m.status.cpu_time)

# ---- the control, against the closed form ----------------------------------------
for name, sol in (("collocation", sol_c), ("shooting", sol_m)):
    t = sol.time
    u = sol.controls[0]
    err = np.max(np.abs(u - (6.0 - 12.0 * t)))
    print("%-12s : max |u(t) - (6 - 12t)| = %.3e" % (name, err))
    ok = ok and err < 1.0e-3

print("\nRESULT: %s" % ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
