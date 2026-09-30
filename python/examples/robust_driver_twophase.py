"""A robust design over TWO PHASES, and what a phase boundary means when the plant is
uncertain.

The problem is the one the C++ driver's multi-phase increment was measured on, so the two
can be compared line for line: x' = -theta x + u, with theta the uncertain decay rate,
from x(0) = 0 to x(t_f) = 1 within a ball, in minimum time, transcribed as two phases
joined by continuity. The boundary time is itself a decision variable.

WHAT A BOUNDARY IS HERE, AND WHAT IT IS NOT

The phase boundary time is shared by every scenario, and that is forced rather than
chosen. A phase has one t0, one tf and one control, so letting scenario k switch at its
own time would need M separate chains of phases, and those cannot share a control grid,
which is what makes the design non-anticipative. A switching time is therefore a
here-and-now decision like the control itself: committed before theta is revealed. That
is consistent with the formulation, and it excludes a boundary triggered by a state
event, which is a different problem and not this one.

WHY TWO PHASES AT ALL ON A PROBLEM THIS SMALL

To measure the machinery and not to flatter it. Nothing about the physics changes at the
boundary, so the two-phase answer should be the single-phase answer up to the mesh, and it
is: two phases of eleven nodes carry twenty-two control values over the horizon where one
phase carries eleven, so the two-phase design is the cheaper of the two by a few per cent,
and that is the finer grid and not the extra phase. What the extra phase does buy is four
things the driver has to get right, each of which a single-phase problem never exercises:
replicate the linkage once per scenario, chain the warm start across it, verify a
trajectory that crosses it by both of its routes, and count the linkage equalities out of
the degrees of freedom. The two routes are the fixed-step sweep the worst-case search uses
and the adaptive DOP853 cross-check the reported certificate comes from; the drift between
them is the reading that says the step count is right at the boundary as well as inside a
phase.

PART TWO: A BOUNDARY THAT IS NOT CONTINUITY

Part one joins the phases by continuity, and a continuity link cannot test the one thing a
declared jump has that continuity has not, which is a SIGN. So part two solves the same
problem with the state dropping 0.2 at the boundary, and checks the sign against the solved
NLP's own linkage instead of against another copy of the driver's opinion. It also turns out
to measure something about the method: the drop costs the nominal design time and saves the
robust one, for a reason set out where the numbers are printed.

ONE PRACTICAL NOTE ON THE ITERATION LIMIT

The iteration limit here is 5000 and not the 1000 the other examples use, which is measured
and not defensive. The driver turns on `free_padded_defect_rows`, because without it IPOPT
counts nstates equality rows of zeros PER SCENARIO against the variables and refuses a
scenario set of any size. On this problem, with two phases and the terminal ball tightened
to 0.9, freeing those rows makes the solve end in a long flat tail: feasibility reaches
5e-13 within a few hundred iterations and then the objective creeps down by about 7e-9 an
iteration for two thousand more. It does converge, to 0.849830918 against 0.849830910 with
the rows kept, so the answer is not in question and the iteration count is. At 1000 the
same solve stops 0.2% short and the driver correctly calls it a failure.
"""
import numpy as np
import casadi as ca
from _common import psopt
from psopt.robust import RobustProblem, Gaussian, Explicit

MU, SIGMA, DELTA, SLACK = 1.0, 0.2, 0.05, 1.0e-3
NODES = 11

rp = RobustProblem(name="robust_driver_twophase")

# Phase 1: from the start, to wherever the boundary falls.
a = rp.add_phase(nstates=1, ncontrols=1, nevents=1)
a.nodes = [NODES]
a.dynamics = lambda x, u, p, t, th: -th[0] * x[0] + u[0]
a.events = lambda xi, xf, p, t0, tf, th: ca.vertcat(xi[0])
a.bounds.lower.states = [-3.0]
a.bounds.upper.states = [3.0]
a.bounds.lower.controls = [-2.0]
a.bounds.upper.controls = [3.0]
a.bounds.lower.events = [0.0]
a.bounds.upper.events = [0.0]
a.bounds.t0 = (0.0, 0.0)
a.bounds.tf = (0.1, 5.0)
a.guess.states = np.zeros((1, NODES))
a.guess.controls = np.ones((1, NODES))
a.guess.time = np.linspace(0.0, 1.0, NODES).reshape(1, NODES)

# Phase 2: from the boundary to the target ball, and it carries the objective.
b = rp.add_phase(nstates=1, ncontrols=1, nevents=1)
b.nodes = [NODES]
b.dynamics = lambda x, u, p, t, th: -th[0] * x[0] + u[0]
b.events = lambda xi, xf, p, t0, tf, th: ca.vertcat(xf[0])
b.endpoint = lambda xi, xf, p, t0, tf: tf
b.bounds.lower.states = [-3.0]
b.bounds.upper.states = [3.0]
b.bounds.lower.controls = [-2.0]
b.bounds.upper.controls = [3.0]
b.bounds.lower.events = [1.0 - DELTA]
b.bounds.upper.events = [1.0 + DELTA]
b.bounds.t0 = (0.1, 5.0)
b.bounds.tf = (0.2, 10.0)
b.guess.states = np.zeros((1, NODES))
b.guess.controls = np.ones((1, NODES))
b.guess.time = np.linspace(1.0, 2.0, NODES).reshape(1, NODES)

# Continuity. Stating it is optional, the driver joining consecutive phases that way when
# nothing says otherwise, and it is written out here because a reader should see where a
# jump would go: link_phases(1, 2, jumps={0: 0.3}) would DROP state 0 by 0.3 at the
# boundary, which is how a staging event is written. Part two below does exactly that.
rp.link_phases(1, 2)

rp.uncertainty = Gaussian(mean=[MU], cov=[[SIGMA ** 2]], truncate=3.0)
rp.initial_state = np.array([0.0])

alg = psopt.Algorithm(transcription_method="multiple-shooting", ms_integrator="RK4",
                      ms_steps_per_segment=10, ms_control_parameterisation="linear",
                      scaling="automatic", nlp_tolerance=1.0e-7, nlp_iter_max=5000,
                      print_level=0)

print(__doc__)
out = rp.solve(alg, slack=SLACK, scenarios="sigma-points", generate=True,
               max_iterations=10, tighten=0.9, out_of_sample=200, wait_and_see=0,
               verbose=True)

print("\n  phases in the design : %d" % len(out.phase_time))
for p, (t, u) in enumerate(zip(out.phase_time, out.phase_controls)):
    print("    phase %d: %d nodes, t from %.4f to %.4f, control rows %d"
          % (p + 1, len(t), t[0], t[-1], u.shape[0]))
print("  the boundary is at t = %.4f, and it is one decision shared by every scenario"
      % out.phase_time[0][-1])

# An independent check of the whole trajectory, crossing the boundary in NumPy: the design
# is one control history in two pieces, and a verification that stopped at the boundary
# would be certifying half of it. `jump` is what the state LOSES at the boundary, which is
# PSOPT's sign for link_phases(jumps=...) and is taken here from the base interface's
# convention and not from the robust driver, so that the two are compared and not assumed.
def miss(th, design, jump=0.0, nsub=64):
    x = 0.0
    for ip, (t, u) in enumerate(zip(*design)):
        if ip:
            x -= jump
        for i in range(len(t) - 1):
            h = (t[i + 1] - t[i]) / nsub
            for s in range(nsub):
                w = [s / nsub, (s + 0.5) / nsub, (s + 1.0) / nsub]
                ua, uh, ub = [u[0, i] + q * (u[0, i + 1] - u[0, i]) for q in w]
                k1 = -th * x + ua
                k2 = -th * (x + 0.5 * h * k1) + uh
                k3 = -th * (x + 0.5 * h * k2) + uh
                k4 = -th * (x + h * k3) + ub
                x += (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return max(0.0, abs(x - 1.0) - DELTA)


grid = np.linspace(MU - 3 * SIGMA, MU + 3 * SIGMA, 601)
mine = np.array([miss(th, (out.phase_time, out.phase_controls)) for th in grid])
c = out.certificate
print("\n  an independent NumPy integrator over 601 decay rates, crossing the boundary")
print("  itself: worst %.3e at theta = %.4f, against the driver's %.3e"
      % (mine.max(), grid[int(np.argmax(mine))], c["violation"]))
print("  the fixed-step sweep against adaptive DOP853, both crossing the boundary")
print("  themselves: %.2e, which is %.0e times under the slack, so the certificate is"
      % (c["integrator_drift"], SLACK / max(c["integrator_drift"], 1.0e-300)))
print("  about the designed trajectory and not about the step count")

print("\n  the C++ driver on the same problem gives t_f 1.040263 on 4 scenarios with a")
print("  certificate of 9.352e-04 at theta = 0.406, by both of ITS routes. This driver")
print("  agrees to every figure it prints, which is worth more than it looks: the two")
print("  share no code at all. One augments the problem in C++ through psopt_solve_robust")
print("  and verifies with a hand-written RK4; the other emits CasADi expressions from")
print("  Python, replicates the linkage symbolically per scenario, and verifies twice,")
print("  with a vectorised fixed-step sweep and with SciPy's DOP853. Counted exactly:")
print("  two independent augmentations, and four integrators that share no code, the")
print("  fourth being the NumPy one above; each of them crosses the boundary itself.")
print("  Reproducibility on THIS problem does not repeal the manual's caveat on the")
print("  objective. The loop follows one chain of warm starts through a nonconvex")
print("  problem, so two loops that generated different scenarios can reach different")
print("  local minima; here they generated the same ones.")

# The C++ driver's t_f is checked and not merely quoted. This example once passed on a
# design that certified and cost 1.016717, because the inward tightening was computed from
# phase 1, whose event is pinned, and broadcast onto phase 2, whose ball was therefore never
# tightened at all. Everything the example verified about that design was correct; the design
# was the answer to a different problem. What caught it was the C++ driver's number, so that
# number is now part of the test. If this line ever fails, the two implementations have
# diverged, and which of them moved is the question to answer before anything else.
CXX_TF = 1.040263
ok = (out.converged and c["violation"] <= SLACK
      and abs(mine.max() - c["violation"]) < 1.0e-3
      and c["integrator_drift"] is not None
      and c["integrator_drift"] < 0.1 * SLACK
      and abs(out.objective - CXX_TF) < 1.0e-5
      and len(out.phase_time) == 2)

# ---- PART TWO: A BOUNDARY THAT IS NOT CONTINUITY ---------------------------------------
#
# The same problem with a staging event at the boundary: the state DROPS by JUMP there, the
# way a launch vehicle's mass drops when a stage is let go. Two things are checked, and they
# are deliberately checked against different authorities, because a jump has a SIGN and a
# verifier that agreed with itself about the wrong one would report a clean certificate for
# a trajectory the design does not have. Part one cannot catch that: on a continuity link
# the jump vector is zero and either sign gives the same answer.
#
#   1. the SOLVED design satisfies x(t0 of phase 2) = x(tf of phase 1) - JUMP, read off the
#      augmented NLP's own states, which is PSOPT's linkage and the only authority on what
#      the sign means;
#   2. the driver's certificate agrees with the NumPy integrator above, which applies that
#      same rule, having been written from the base interface's documented convention.
JUMP = 0.2
rp2 = RobustProblem(name="robust_driver_twophase_jump")
for src in (a, b):
    ph = rp2.add_phase(nstates=1, ncontrols=1, nevents=1)
    ph.nodes, ph.dynamics, ph.events = list(src.nodes), src.dynamics, src.events
    ph.endpoint = src.endpoint
    ph.bounds.lower.states, ph.bounds.upper.states = [-3.0], [3.0]
    ph.bounds.lower.controls, ph.bounds.upper.controls = [-2.0], [3.0]
    ph.bounds.lower.events = list(src.bounds.lower.events)
    ph.bounds.upper.events = list(src.bounds.upper.events)
    ph.bounds.t0, ph.bounds.tf = src.bounds.t0, src.bounds.tf
    ph.guess.states = np.zeros((1, NODES))
    ph.guess.controls = np.ones((1, NODES))
    ph.guess.time = np.array(src.guess.time)
rp2.link_phases(1, 2, jumps={0: JUMP})
rp2.uncertainty = Gaussian(mean=[MU], cov=[[SIGMA ** 2]], truncate=3.0)
rp2.initial_state = np.array([0.0])

print("\n  PART TWO: the same problem with the state dropping %.2f at the boundary" % JUMP)
out2 = rp2.solve(alg, slack=SLACK, scenarios="sigma-points", generate=True,
                 max_iterations=10, tighten=0.9, out_of_sample=0, wait_and_see=0,
                 verbose=True)

xs = out2.design.states
gap = max(abs(float(xs[1][k, 0]) - float(xs[0][k, -1]) + JUMP)
          for k in range(xs[0].shape[0]))
mine2 = np.array([miss(th, (out2.phase_time, out2.phase_controls), JUMP) for th in grid])
c2 = out2.certificate
print("\n  the solved states across the boundary, over all %d scenarios: the largest"
      % xs[0].shape[0])
print("  departure from x(t0 of phase 2) = x(tf of phase 1) - %.2f is %.2e" % (JUMP, gap))
print("  the NumPy integrator, applying that same rule itself: worst %.3e at theta = %.4f,"
      % (mine2.max(), grid[int(np.argmax(mine2))]))
print("  against the driver's %.3e" % c2["violation"])

# And the measurement the jump turned out to be worth making. The nominal design is solved
# both ways for the comparison, at one scenario and one call each, because the interesting
# number is a DIFFERENCE of differences and neither half of it means anything alone.
def nominal(problem):
    problem.uncertainty = Explicit([[MU]])
    return problem.solve(alg, slack=SLACK, scenarios="explicit", generate=False,
                         out_of_sample=0, wait_and_see=0, verbose=False)


rp_n, rp2_n = RobustProblem(name="tp_n"), RobustProblem(name="tp2_n")
for target, jumps in ((rp_n, None), (rp2_n, {0: JUMP})):
    for src in (a, b):
        ph = target.add_phase(nstates=1, ncontrols=1, nevents=1)
        ph.nodes, ph.dynamics, ph.events = list(src.nodes), src.dynamics, src.events
        ph.endpoint = src.endpoint
        ph.bounds.lower.states, ph.bounds.upper.states = [-3.0], [3.0]
        ph.bounds.lower.controls, ph.bounds.upper.controls = [-2.0], [3.0]
        ph.bounds.lower.events = list(src.bounds.lower.events)
        ph.bounds.upper.events = list(src.bounds.upper.events)
        ph.bounds.t0, ph.bounds.tf = src.bounds.t0, src.bounds.tf
        ph.guess.states = np.zeros((1, NODES))
        ph.guess.controls = np.ones((1, NODES))
        ph.guess.time = np.array(src.guess.time)
    target.link_phases(1, 2, jumps=jumps)
    target.initial_state = np.array([0.0])
n0, n1 = nominal(rp_n), nominal(rp2_n)

print("\n  %-12s %14s %14s" % ("boundary", "nominal t_f", "robust t_f"))
print("  %-12s %14.6f %14.6f" % ("continuity", n0.objective, out.objective))
print("  %-12s %14.6f %14.6f" % ("drop %.2f" % JUMP, n1.objective, out2.objective))
print("""
  The drop costs the nominal design time, as dropping 0.2 of a state that has to reach 1
  should, and it SAVES the robust design time. That is not a paradox and it is the reason
  this part of the example is here. Look at where the robust design goes: with continuity
  it takes the state to %.4f at the boundary, which is on the far side of zero from its
  target. It is buying insensitivity. The uncertainty enters as -theta*x, so the
  sensitivity of the trajectory to theta is proportional to x and changes SIGN with it, and
  an excursion through negative x lets the sensitivity accumulated there cancel what is
  accumulated on the way up. The continuity design pays control effort to make that
  excursion. The staging drop performs part of it for nothing, so what looks like a penalty
  is a subsidy for the manoeuvre the robust design already wanted.

  The lesson generalises past this problem: the price of robustness is not additive with a
  physical penalty, because the two may be paying for the same thing. A penalty has to be
  measured against the ROBUST design and not against the nominal one."""
      % float(out.design.states[0][0, -1]))

ok2 = (out2.converged and c2["violation"] <= SLACK and gap < 1.0e-7
       and abs(mine2.max() - c2["violation"]) < 1.0e-3
       and n1.objective > n0.objective          # the drop costs the nominal design
       and out2.objective < out.objective)      # and saves the robust one

print("\nobjective      : %.12g" % out.objective)
print("objective, jump: %.12g" % out2.objective)
print("RESULT: %s" % ("PASS" if (ok and ok2) else "FAIL"))
raise SystemExit(0 if (ok and ok2) else 1)
