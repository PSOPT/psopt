"""A discrete-valued static parameter, through the PSOPT Python API.

Mirrors the native examples/int_static_linear. A double integrator is driven by a
constant acceleration p that may only take the values {0,1,2,3}, and the cost
penalises the terminal velocity away from a target of 2.3. Since v(1) = p, the
closed-form answer is p = 2 with J = (2 - 2.3)^2 = 0.09.

The point of the example is why a static parameter is not treated like an integer
control. Outer convexification is tight for a control because a control may switch
arbitrarily fast, so a fractional weight is realisable in the limit by chattering
between admissible values. A static parameter cannot chatter: a convex blend of two
admissible values is not a realisable design, no mesh refinement closes the gap, and
the relaxed optimum is a bound that may be attained by no admissible design. Here
the relaxed problem would answer p = 2.3 exactly, with J = 0, which is not a design
that can be built.

PSOPT therefore solves these exactly, by enumerating the admissible combinations,
solving each with the parameters pinned, and keeping the best. Declaring an integer
parameter switches the driver to psopt_solve_integer automatically.
"""
import numpy as np
import casadi as ca
import psopt

TARGET = 2.3
ADMISSIBLE = [0.0, 1.0, 2.0, 3.0]


def build(nnodes=20):
    prob = psopt.Problem(name="integer_parameter")
    ph = prob.add_phase(nstates=2, ncontrols=0, nevents=2, npath=0, nparameters=1)
    ph.nodes = [nnodes]

    # s_dot = v, v_dot = p  with p a static parameter
    ph.dynamics = lambda x, u, p, t: ca.vertcat(x[1], p[0])
    ph.endpoint = lambda xi, xf, p, t0, tf: (xf[1] - TARGET) ** 2
    ph.events = lambda xi, xf, p, t0, tf: ca.vertcat(xi[0], xi[1])

    ph.bounds.lower.states = [-10.0, -10.0]
    ph.bounds.upper.states = [10.0, 10.0]
    ph.bounds.lower.events = [0.0, 0.0]
    ph.bounds.upper.events = [0.0, 0.0]
    # A sensible relaxed range; the enumeration pins the parameter to each
    # admissible value in turn, so these bounds only have to contain them.
    ph.bounds.lower.parameters = [0.0]
    ph.bounds.upper.parameters = [3.0]
    ph.bounds.t0 = (0.0, 0.0)
    ph.bounds.tf = (1.0, 1.0)

    tg = np.linspace(0.0, 1.0, nnodes)
    ph.guess.time = tg.reshape(1, nnodes)
    ph.guess.states = np.vstack([np.zeros(nnodes), np.linspace(0.0, 2.0, nnodes)])
    ph.guess.parameters = np.array([[1.0]])
    return prob, ph


if __name__ == "__main__":
    alg = psopt.Algorithm(collocation_method="Legendre", nlp_iter_max=1000,
                          nlp_tolerance=1e-8, print_level=0)

    prob, ph = build()
    ph.declare_integer_parameter(0, ADMISSIBLE)
    sol = prob.solve(alg)

    sel = sol.integer_parameters[0]
    p_expected, J_expected = 2.0, (2.0 - TARGET) ** 2

    print("\n  Integer static parameter, admissible set %s, target c = %.1f\n"
          % (ADMISSIBLE, TARGET))
    print("  enumerated  ->  p = %.6f,  J = %.6f" % (sel["value"], sol.objective))
    print("  closed form ->  p = %.6f,  J = %.6f" % (p_expected, J_expected))

    ok = (abs(sel["value"] - p_expected) < 1e-9
          and abs(sol.objective - J_expected) < 1e-6)
    print("\n  %s\n" % ("PASS" if ok else "FAIL"))

    # For contrast: the same problem with the parameter left continuous. Its optimum
    # is p = 2.3 with J = 0, a design that cannot be built, which is why the relaxed
    # answer is not usable for a static parameter.
    prob2, _ph2 = build()
    sol2 = prob2.solve(alg)
    print("  relaxed (parameter left continuous) ->  p = %.6f,  J = %.6f"
          % (sol2.parameters[0], sol2.objective))
    print("  -- not realisable: the design must come from the admissible set.\n")
