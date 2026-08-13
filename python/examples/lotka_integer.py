"""Lotka-Volterra fishing with an integer control, through the PSOPT Python API.

Mirrors the native examples/lotka: a binary fishing control switches the fleet on
and off, and the objective is to steer the prey-predator system to its coexistence
steady state (1,1) over a fixed horizon.

The dynamics are written once with the fishing control as an ordinary control.
declare_integer_control tells PSOPT that it takes values in {0,1}; the outer
convexification (weight controls plus an SOS1 constraint) happens inside the
solver, and the rounded integer control comes back on the solution through
sum-up rounding.

Reference: the relaxed optimum is the infimum of the integer problem,
Phi* = 1.34408 (mintOC). No optimal integer control exists; the integer objective
approaches Phi* from above as the mesh is refined.
"""
import numpy as np
import casadi as ca
import psopt

C0, C1 = 0.4, 0.2
PHI_STAR = 1.34408
T_END = 12.0


def solve(nnodes):
    prob = psopt.Problem(name="lotka_integer")
    ph = prob.add_phase(nstates=3, ncontrols=1, nevents=3, npath=0)
    ph.nodes = [nnodes]

    # x0 prey, x1 predator, x2 accumulated cost. w is the fishing control, written
    # here as an ordinary control; its integrality is declared below.
    ph.dynamics = lambda x, u, p, t: ca.vertcat(
        x[0] - x[0] * x[1] - C0 * x[0] * u[0],
        -x[1] + x[0] * x[1] - C1 * x[1] * u[0],
        (x[0] - 1.0) ** 2 + (x[1] - 1.0) ** 2,
    )
    ph.endpoint = lambda xi, xf, p, t0, tf: xf[2]
    ph.events = lambda xi, xf, p, t0, tf: ca.vertcat(xi[0], xi[1], xi[2])

    ph.bounds.lower.states = [0.0, 0.0, 0.0]
    ph.bounds.upper.states = [2.0, 2.0, 100.0]
    # The bounds of an integer control are ignored; written for readability.
    ph.bounds.lower.controls = [0.0]
    ph.bounds.upper.controls = [1.0]
    ph.bounds.lower.events = [0.5, 0.7, 0.0]
    ph.bounds.upper.events = [0.5, 0.7, 0.0]
    ph.bounds.t0 = (0.0, 0.0)
    ph.bounds.tf = (T_END, T_END)

    # The one line that makes the control discrete.
    ph.declare_integer_control(0, [0.0, 1.0])

    tg = np.linspace(0.0, T_END, nnodes)
    ph.guess.time = tg.reshape(1, nnodes)
    ph.guess.controls = 0.5 * np.ones((1, nnodes))
    xs = np.zeros((3, nnodes))
    xs[0, :], xs[1, :] = 0.5, 0.7
    ph.guess.states = xs

    alg = psopt.Algorithm(collocation_method="trapezoidal", nlp_iter_max=3000,
                          nlp_tolerance=1e-7, print_level=0)
    return prob.solve(alg), ph


def simulate(control, widths):
    """Forward-simulate the rounded binary control to get its true objective."""
    x = np.array([0.5, 0.7, 0.0])
    sub = 40

    def rhs(s, w):
        return np.array([s[0] - s[0] * s[1] - C0 * s[0] * w,
                         -s[1] + s[0] * s[1] - C1 * s[1] * w,
                         (s[0] - 1.0) ** 2 + (s[1] - 1.0) ** 2])

    for w, h in zip(control, widths):
        dt = h / sub
        for _ in range(sub):
            k1 = rhs(x, w)
            k2 = rhs(x + 0.5 * dt * k1, w)
            k3 = rhs(x + 0.5 * dt * k2, w)
            k4 = rhs(x + dt * k3, w)
            x = x + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return x[2]


if __name__ == "__main__":
    print("\n  Lotka-Volterra fishing with a declared integer control")
    print("  reference relaxed optimum Phi* = %.5f\n" % PHI_STAR)
    print("  %6s  %14s  %14s  %12s  %10s" %
          ("nodes", "relaxed Phi", "integer Phi", "gap to Phi*", "switches"))
    print("  " + "-" * 66)
    for n in (40, 80, 160, 320):
        sol, _ph = solve(n)
        ic = sol.integer_controls[0]          # one entry per declared integer control
        integer_phi = simulate(ic.control, ic.interval_widths)
        print("  %6d  %14.6f  %14.6f  %12.3e  %10d" %
              (n, sol.states[2, -1], integer_phi, integer_phi - PHI_STAR, ic.n_switches))
    print("  " + "-" * 66)
    print("  sol.controls is in the weights layout when an integer control is")
    print("  declared; the rounded control is in sol.integer_controls[k].control.\n")
