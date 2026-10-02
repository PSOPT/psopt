// include/robust.h
//
// Robust optimal control by scenario augmentation, v1.
//
// A robust optimal control problem asks for ONE control, committed before an
// uncertain parameter is revealed, that serves every plant in an uncertainty set.
// Replacing that set by a finite list of SCENARIOS turns it into one ordinary
// deterministic optimal control problem -- M copies of the state against one copy
// of the control, which is what makes the design non-anticipative -- and psopt()
// solves it as it stands.
//
// This header holds the parts of that construction which are the library's rather
// than the user's:
//
//   * the uncertainty, which is three things and keeps them apart: a SET to be
//     robust over, a DISTRIBUTION to sample from, and a rule for placing
//     scenarios;
//   * pure, unit-testable building blocks -- unscented scenario sets, a
//     low-discrepancy sequence over the set, containment, projection, boundary
//     points, the coverage of a truncated set, and the worst-case oracle;
//   * the driver psopt_solve_robust, which runs the outer loop: solve, find the
//     parameter the design serves worst, add it, re-solve.
//
// psopt() itself is never modified. The driver sits in the same architectural
// slot as psopt_solve_integer (include/integer_parameters.h): an outer loop over
// an unmodified solver.
//
// WHAT THE USER WRITES, AND WHY THERE ARE TWO WAYS TO WRITE IT
//
// In the Python interface the driver builds the augmented problem itself, because
// the user's equations arrive as CasADi expressions it can replicate. In C++ they
// arrive as dae() and events() with fixed signatures, taped by CppAD, and nothing
// can rewrite an existing pair. What the library can do is call a nominal pair that
// takes the uncertain parameter as an extra argument, once per scenario, from inside
// an augmented dae of its own. That is RobustModel, and a caller who fills one
// writes the nominal problem and nothing else: the augmentation, the bounds with
// their inward tightening, the warm start, the risk measure, the ancillary feedback
// and the verification integrator are all the library's.
//
// The other way remains. A caller who fills spec.setup writes the augmented problem
// themselves and supplies spec.violation to verify it; python/examples/robust_arm.py is
// that route written out, in Python, and reads as C++ would. It is the route to take when
// the nominal equations cannot be evaluated
// numerically -- a dae reading get_delayed_state or get_interpolated_state, for
// instance -- and it is the stronger claim about a certificate, the verification then
// sharing no code with the library at all. Both routes use the same driver: the
// scenario rule, the oracle, the generation loop and the polish step are the same for
// every problem either way.
//
// Copyright (c) Victor M. Becerra, 2026. Part of the PSOPT library (LGPL).

#ifndef PSOPT_ROBUST_H
#define PSOPT_ROBUST_H

#include "psopt.h"

#include <vector>
#include <cmath>
#include <cstddef>
#include <limits>

//////////////////////////////////////////////////////////////////////////
///////////////////  The uncertainty  ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// A Gaussian's SET is the ellipsoid at `truncate` standard deviations; a
// Uniform's set is its box; an Explicit uncertainty is its own list. The set is
// what the design is robust over and what the worst case is sought in; the
// distribution is what out-of-sample checks draw from. Conflating them is how a
// design ends up being scored against its own truncation.
struct RobustUncertainty {
    enum Kind { GAUSSIAN, UNIFORM, EXPLICIT };

    Kind        kind;
    RowVectorXd mean;       // GAUSSIAN
    MatrixXd    cov;        // GAUSSIAN
    double      truncate;   // GAUSSIAN: the set's radius, in standard deviations
    RowVectorXd lo, hi;     // UNIFORM
    std::vector<RowVectorXd> points;    // EXPLICIT
    RowVectorXd              weights;   // EXPLICIT

    RobustUncertainty() : kind(EXPLICIT), truncate(3.0) {}
};

RobustUncertainty robust_gaussian(const RowVectorXd& mean, const MatrixXd& cov,
                                  double truncate);
RobustUncertainty robust_uniform(const RowVectorXd& lo, const RowVectorXd& hi);
RobustUncertainty robust_explicit(const std::vector<RowVectorXd>& points);

int robust_dimension(const RobustUncertainty& U);

//////////////////////////////////////////////////////////////////////////
///////////////////  Pure building blocks  ///////////////////////////////
//////////////////////////////////////////////////////////////////////////

// The inverse standard normal CDF, by Acklam's rational approximation refined by
// one Halley step. Relative error below 1e-15 over the open unit interval, which
// is far more than the low-discrepancy seeding below needs; it is refined anyway
// because the same function is the natural one to expose and a user may want it
// for something sharper.
double robust_norm_ppf(double p);

// The regularized lower incomplete gamma function P(a, x), by the series for
// x < a+1 and the continued fraction otherwise.
double robust_lower_gamma(double a, double x);

// What fraction of an n-dimensional standard normal lies inside the ellipsoid at
// `k` standard deviations. This is chi-square with n degrees of freedom at k^2,
// and it falls away with dimension in a way that catches people out: at k = 3 it
// is 99.73% on a line, 98.89% in the plane and 97.07% in three dimensions. A habit
// formed on scalar uncertainty truncates a multi-parameter set more tightly than
// it means to, so the number is worth being able to print.
double robust_set_coverage(double k, int n);

// The unscented scenario set: 2n+1 points reproducing the mean and the covariance
// exactly, with kappa = 3 - n. For one Gaussian parameter these are the three
// points of the Gauss-Hermite rule, {mu, mu +- sqrt(3) sigma} with weights
// {2/3, 1/6, 1/6}; the sigma points of the estimation literature and the
// Gauss-Hermite nodes of the quadrature literature are the same three numbers.
//
// Returns false, leaving the outputs untouched, when n > 3 would make the central
// weight negative. A negative weight is admissible in a quadrature and not in a
// scenario set, every member of which is a constraint the design must satisfy;
// use robust_low_discrepancy instead.
bool robust_sigma_points(const RobustUncertainty& U,
                         std::vector<RowVectorXd>& points, RowVectorXd& weights);

// m low-discrepancy points spread over the SET, by a Halton sequence. For a
// Gaussian the radii are stratified over (0, 1] with 1 attained, because the
// worst parameter very often sits on the boundary and a sequence that reaches it
// only by accident would miss exactly the point a certificate is about.
void robust_low_discrepancy(const RobustUncertainty& U, int m, int skip,
                            std::vector<RowVectorXd>& points);

// Extreme points of the set: the 2n ends of a Gaussian's principal axes, the
// corners of a box up to eight dimensions and its face centres beyond, or an
// Explicit uncertainty's own list. Added by hand to every worst-case search for
// the reason above.
void robust_boundary_points(const RobustUncertainty& U,
                            std::vector<RowVectorXd>& points);

bool        robust_contains(const RobustUncertainty& U, const RowVectorXd& theta);
RowVectorXd robust_clip(const RobustUncertainty& U, const RowVectorXd& theta);

// Distance from the mean in units of the covariance, for a Gaussian; the largest
// normalised excursion beyond the box centre for a Uniform. The quantity a
// marginal interval throws away.
double robust_mahalanobis(const RobustUncertainty& U, const RowVectorXd& theta);

//////////////////////////////////////////////////////////////////////////
///////////////////  The worst-case oracle  //////////////////////////////
//////////////////////////////////////////////////////////////////////////

// What a solve produced, in a form that can be kept.
//
// Sol owns raw arrays and frees them in its destructor but defines no copy
// constructor or copy assignment -- a documented rule-of-three gap that
// psopt_solve_integer works around in its own way -- so a Sol cannot be held from
// one iteration to the next. This carries the parts of one that the outer loop and
// the user actually need: the trajectory to warm-start the next solve from, and
// the trajectory to verify.
// One phase of a solved design. A single-phase problem has one of these and the fields of
// RobustDesign below are it; a multi-phase problem has one per phase, in order.
struct RobustPhaseTrajectory {
    MatrixXd time;          // 1 x N
    MatrixXd controls;      // ncontrols x N, the nodal table
    MatrixXd states;        // nstates x N   (all scenario copies, as solved)
    MatrixXd controls_full; // ncontrols x (2N-1), or empty; see RobustDesign
    MatrixXd time_full;     // 1 x (2N-1), or empty
};

struct RobustDesign {
    // PHASE 1. These five fields are the first phase and nothing else, which is the whole
    // design when there is one phase and is why they are spelled without a phase index: a
    // single-phase problem, which is most of them, reads exactly as it did before the
    // driver learned about phases. A multi-phase design is in `phase` below, whose first
    // entry holds these same five.
    MatrixXd time;          // 1 x N
    MatrixXd controls;      // ncontrols x N, the nodal table
    MatrixXd states;        // nstates x N   (all scenario copies, as solved)
    MatrixXd parameters;    // static parameters, if any, shared by every phase
    double   objective;
    bool     valid;

    // Every phase, in order, the first of them the one above. A violation function written
    // for a multi-phase problem reads this; one written for a single-phase problem need
    // never know it exists.
    std::vector<RobustPhaseTrajectory> phase;

    // The COMPLETE control history, node and midpoint values interleaved, and the times
    // that go with it. Filled for the two transcriptions that carry a control at the
    // midpoint of every interval, Hermite-Simpson and multiple shooting under the
    // "quadratic" parameterisation, and empty for every other one, where the nodal table
    // IS the history.
    //
    // Anything that has to reproduce the designed control must prefer these, a
    // verification integrator above all. Reading the nodal values alone under those two
    // parameterisations gives two thirds of the control variables and none of the
    // curvature, and the error is not in the safe direction: a parabola read as a chord
    // flatters a design as readily as it damns one, because the parabola lies outside the
    // chord on one side. Measured on the two-link arm, the same Hermite-Simpson design is
    // charged 8.3e-02 read correctly and 2.6e+00 read as chords, and the same multiple
    // shooting quadratic design is charged 2.9 correctly and 1.8 as chords.
    MatrixXd controls_full; // ncontrols x (2N-1), or empty
    MatrixXd time_full;     // 1 x (2N-1), or empty

    RobustDesign() : objective(0.0), valid(false) {}

    // How many phases this design has, which is one whenever the caller never built more.
    int nphases() const { return phase.empty() ? (valid ? 1 : 0) : (int) phase.size(); }
};

// How badly a given design serves one parameter vector: zero when every declared
// constraint holds along an independently integrated trajectory, and otherwise the
// amount by which the worst of them is exceeded. Supplied by the user, because the
// independent integrator is the user's -- and it must BE independent: checking a
// design against the integrator that produced it checks nothing.
typedef double (*RobustViolationFn)(const RowVectorXd& theta, void* user_data);

// Search the set for the parameter served worst. Seeding by boundary points and a
// low-discrepancy sequence, then a shrinking cloud around the best few, which
// needs no derivative -- and must not want one, because the violation is a maximum
// over constraints and over time and is not differentiable where the active
// constraint changes, which is exactly where the worst case tends to sit.
//
// For one parameter the seeding alone is effectively exhaustive. In more than one
// it is not, and nothing here claims otherwise: the result is the worst that was
// FOUND, over robust_worst_case_evaluations(n_seed, n_refine, U) evaluations, and
// the inner problem is not concave.
double robust_worst_case(const RobustUncertainty& U,
                         RobustViolationFn violation, void* user_data,
                         int n_seed, int n_refine, unsigned seed,
                         RowVectorXd& worst_at);

// How many calls to the violation function robust_worst_case will make. Exact, so
// a caller can budget before starting rather than discover afterwards.
long robust_worst_case_evaluations(const RobustUncertainty& U, int n_seed,
                                   int n_refine);

//////////////////////////////////////////////////////////////////////////
///////////////////  The scenario budget  ////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// Degrees of freedom the transcription has left, as a LOWER bound.
//
//     dof = free initial states + ncontrols * nodes + free time endpoints
//           + free static parameters - equality events
//
// A scenario brings nstates values at t0 that the defect equations do not determine,
// and nstates pinned initial conditions that determine them. THOSE CANCEL, so a
// scenario is free of itself, and what costs one is an equality event beyond its
// initial condition -- a pinned TERMINAL state above all, which one open-loop control
// cannot meet for several different plants anyway. Relaxing such an event to a
// tolerance costs nothing at all, an inequality event taking no degree of freedom.
//
// The exception is a transcription that collocates every stored node, which the global
// Lobatto schemes (Legendre, Chebyshev) do and nothing else does. There the defect
// block holds one condition per stored value per state and the initial condition is
// one more, so the over-determination is real, it costs nstates per scenario, and no
// node count rescues it. Measured on the two-link arm: refused above twelve scenarios
// whatever the mesh. Multiple shooting and trapezoidal collocation carry eighty on the
// same problem.
//
// WHAT THIS SAID BEFORE, AND WHY IT WAS WRONG. The free initial states were missing
// from the sum, and the formula nevertheless predicted the observed wall exactly -- the
// arm at 25 nodes refused above twelve scenarios, a two-state problem at 31 nodes above
// fifteen. Both were phantom. PSOPT's defect block holds nstates*(norder+1) rows and a
// scheme with norder intervals fills only nstates*norder of them; the rest were written
// as zeros with bounds [0,0], and IPOPT counts equality rows against variables. The arm
// had 51 degrees of freedom at every scenario count and was refused for having -9.
// NLP_bounds now frees those rows, and this is the honest arithmetic.
//
// It remains a lower bound because an equality event can be linearly dependent on the
// defect equations, and a dependent constraint removes no freedom. So a non-positive
// count does not prove the problem is over-determined: use it to EXPLAIN a failure,
// never to refuse in advance. Left to IPOPT the failure arrives as
// Not_Enough_Degrees_Of_Freedom, return code -10, after the whole problem has been
// assembled and taped, and the message names no remedy at all.
int robust_degrees_of_freedom(Prob& problem, int iphase, Alg& algorithm);

// The same count for a problem of SEVERAL phases: the per-phase counts summed, less the
// linkage equalities, which belong to the problem and to no phase. A sum that had not
// subtracted them would report freedom the problem does not have, which is the direction
// that lets a starved solve through. `robust_linkage_equalities` counts them, reading the
// linkage bounds when the caller set them and taking all of them as equalities when the
// caller left them unset, which is what auto_link builds.
//
// `robust_tightest_phase` names the phase with the least freedom, which is the one whose
// arithmetic a failure diagnosis should print: the total says a multi-phase problem is
// starved, and this says where.
int robust_problem_degrees_of_freedom(Prob& problem, Alg& algorithm);
int robust_linkage_equalities(Prob& problem);
int robust_tightest_phase(Prob& problem, Alg& algorithm);

// Is every entry of the initial guess a finite number?
//
// Worth its own check, because of how the alternative presents itself. A guess
// built by integrating a previous design can diverge -- the horizon it integrates
// over is a decision variable, and a robust final time grows as scenarios are
// added, so a step size that was ample at the start need not be by the end. The
// infinities then reach the constraint evaluation as NaN, some rows are never
// written, and PSOPT's constraint-coverage guard reports a defect in the library.
// The information is there, in a NaN warning printed just above it, but the fatal
// message names the wrong culprit. Checking first turns that into one sentence
// naming the array.
//
// On failure `offender` is set to "states", "controls", "time" or "parameters".
bool robust_guess_is_finite(Prob& problem, int iphase, const char** offender);

// The same arithmetic written out, with the remedy that fits, into a caller's buffer.
// `nscenarios` is used only to report the cost per scenario.
void robust_dof_message(Prob& problem, int iphase, int nscenarios, Alg& algorithm,
                        char* buffer, size_t buffer_size);

//////////////////////////////////////////////////////////////////////////
///////////////////  The model interface  ////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// WHY THIS EXISTS
//
// Writing the augmentation by hand was 173 of the 317 code lines of the arm example
// before this interface existed: eight for the dae loop, thirteen for the events loop,
// fifty-one to size the phase and replicate the bounds, forty-nine for the warm start, and
// forty-five for the verification integrator. Only the first eight are what people expect
// the work to be.
//
// A RobustModel is the NOMINAL problem, stated once, with the uncertain parameter as one
// extra argument. psopt_solve_robust then builds the augmented problem itself: it
// installs a dae and an events function that loop the user's nominal pair over the
// scenario list, sizes the phase for M copies of the state, replicates the bounds with
// the inward tightening applied, and builds the warm start by integrating each
// scenario's own plant through the previous control.
//
// THE ONE DESIGN DECISION THAT MAKES THIS POSSIBLE
//
// theta arrives as `const double*`, not as adouble. A scenario is DATA: it is not a
// decision variable and it must not reach the derivative tape. Keeping it a plain number
// has a consequence worth more than the tidiness: the same nominal dae can be called
// numerically, outside any taping context, which is what lets the library build the warm
// start (and, in due course, a verification integrator) out of the user's own equations
// instead of asking for a second hand-written copy of them. See robust_dae_value.
//
// WHAT THE NOMINAL FUNCTIONS MAY NOT DO
//
// They must be pure arithmetic on their arguments. A dae that reaches into the tape
// through get_delayed_state, get_interpolated_state or auto_link cannot be evaluated
// numerically, so such a problem keeps its own spec.setup and writes the augmentation by
// hand. The Python driver has the same restriction, for the same reason.

// What is minimised when the cost is a random variable. One number per realisation
// makes "minimise the cost" an incomplete problem statement, and which summary of the
// distribution is minimised is a modelling decision with consequences.
//
//   ROBUST_NOMINAL         the cost of the central scenario. Feasibility is robust and
//                          the objective is the nominal one, which is the right choice
//                          when the cost is shared, as it is in a minimum-time problem.
//   ROBUST_EXPECTATION     the scenario rule's estimate of E[J].
//   ROBUST_MEAN_VARIANCE   E[J] + mv_lambda * Var[J].
//   ROBUST_CVAR            the conditional value at risk at cvar_alpha, the mean of the
//                          worst 1 - alpha of the distribution, by the device of
//                          Rockafellar and Uryasev.
//
// The first two are weighted sums of per-scenario costs, so they go straight into the
// integrand and the endpoint cost and need nothing else. The last two need each
// scenario's cost as a quantity in its own right, because the variance of a Lagrange
// cost across scenarios is not the integral of anything, so the augmentation carries one
// extra state per scenario whose derivative is that scenario's integrand and whose
// initial value is pinned to zero. A state needs bounds, which is what
// RobustModel::cost_lower and cost_upper are for.
// Where the ancillary gain comes from, if there is one.
//
// A design with no gain is OPEN LOOP: one control history, committed before the
// uncertainty is revealed, serving every plant in the set unaided. That is honest and over
// a long horizon it is expensive, costing a factor of about three in final time on the
// two-link arm. With a gain the design becomes a tube,
//
//     u_k(t) = u_bar(t) + K ( x_k(t) - x_ref(t) )
//
// in which u_bar remains the only control decision variable and the reference is scenario
// 0 of the rule, the centre of the set, which therefore runs open loop by construction.
// The design stays a here-and-now design: u_bar and K are both fixed before the parameter
// is revealed and nothing adapts to it. What changes is that one history no longer has to
// serve every plant by itself.
//
//   ROBUST_FEEDBACK_NONE       open loop.
//   ROBUST_FEEDBACK_CONSTANT   a gain the caller gives, in model.feedback_gain.
//   ROBUST_FEEDBACK_SCHEDULED  a gain the caller gives as a function of time, in
//                              model.feedback_schedule. Unlike the Python driver, which
//                              must emit its schedule for CppAD to tape and so is limited
//                              to a polynomial in t, a C++ schedule is an ordinary adouble
//                              function: an interpolant through the conditional helpers of
//                              psopt.h, a table lookup, anything tape-safe.
//   ROBUST_FEEDBACK_CODESIGN   the gain's entries become static parameters and are
//                              optimised with the trajectory.
//   ROBUST_FEEDBACK_CODESIGN_SCHEDULE
//                              a time-varying gain optimised as extra CONTROLS, which is
//                              the representation to reach for when a co-designed schedule
//                              is wanted: the transcription already gives a control a time
//                              profile at the resolution the trajectory has, so the gain
//                              inherits one with no basis to choose.
//
// CO-DESIGN IS OFFERED AND IS NOT RELIABLE, and what that means is worth stating
// carefully, because the measurements point both ways.
//
// Most co-designed variants tried on the two-link arm, constant and scheduled, with bounds
// from wide to a box around a gain known to certify, on scenario sets of three, nine and
// sixteen points, came out CHEAPER on the objective than a given LQR gain and FAILED their
// certificate, by between one and five orders of magnitude against a slack of 1e-3. In one
// run of the C++ driver a co-designed constant gain in a box of plus or minus five, started
// from that LQR gain, certified: t_f 3.4210 against 3.5583 for the given gain, with no
// violation anywhere in the set and a realised control 6.4e-04 outside its bounds, and
// three integrators written independently of each other agree on that verdict. The Python
// driver on the same problem, the same box and the same guess reached a different local
// minimum at nearly the same objective which misses by 6.1e-02, and Python's verifier
// passes the C++ design when handed it, so the disagreement is between the two SOLVES and
// not between the two verifications.
//
// So a co-designed gain sometimes certifies, and whether it does is a property of the local
// minimum the solve happens to reach and not of the formulation. The reason it cannot be
// relied on is structural. A scenario set enters the design as a constraint set, so the
// gain is rewarded for making those M plants cheap and charged nothing for what it does to
// the rest of the family; and a gain multiplies a deviation that is itself a function of
// the uncertain parameter, so its leverage on an unsampled plant is bounded by nothing the
// design can see. Neither a denser set nor a tighter bound repairs that. A gain chosen by a
// criterion that quantifies over the whole family, which is what solving a Riccati equation
// does, carries no such hazard.
//
// The practical reading: co-design is worth trying, its result is worth nothing until
// verified, and the verification is where the decision is made. The driver reports the
// largest entry of the gain it designed and how many entries sit on their bound, because a
// gain held at its bound is being chosen by the bound and not by the problem.
enum RobustFeedbackKind {
    ROBUST_FEEDBACK_NONE = 0,
    ROBUST_FEEDBACK_CONSTANT,
    ROBUST_FEEDBACK_SCHEDULED,
    ROBUST_FEEDBACK_CODESIGN,
    ROBUST_FEEDBACK_CODESIGN_SCHEDULE
};

// A gain schedule the caller gives. `gain` receives ncontrols*nstates entries in row-major
// order and is written in adouble arithmetic, because it is evaluated inside the taped
// dae where the time may itself depend on a free final time.
typedef void (*RobustGainFn)(adouble& time, adouble* gain, void* user_data);

enum RobustRisk {
    ROBUST_NOMINAL = 0,
    ROBUST_EXPECTATION,
    ROBUST_MEAN_VARIANCE,
    ROBUST_CVAR
};

typedef void (*RobustDaeFn)(adouble* derivatives, adouble* path, adouble* states,
                            adouble* controls, adouble* parameters, adouble& time,
                            const double* theta, int ntheta,
                            adouble* xad, int iphase, Workspace* workspace);

// How one phase joins the next, for one scenario: `nlink` values the driver holds at zero.
// The states and times are those of the two phases either side of the boundary.
typedef void (*RobustLinkFn)(adouble* link,
                             adouble* final_states_prev, adouble& tf_prev,
                             adouble* initial_states_next, adouble& t0_next,
                             adouble* parameters, const double* theta, int ntheta,
                             void* user_data);

// The same in numbers, for the verifier: the state one phase ends at, and the state the
// next phase starts from. Continuity is x0_next = xf_prev, which is what the driver assumes
// when neither this nor RobustModel::link is given.
typedef void (*RobustJumpFn)(const double* theta, int ntheta, int from_phase,
                             const double* final_states_prev, double t,
                             double* initial_states_next, void* user_data);

// The state a scenario starts from, when it depends on the uncertain parameter. Writes
// nstates values. See RobustModel::initial_state_fn for what it is for and what it does
// not do.
typedef void (*RobustInitialStateFn)(const double* theta, int ntheta, double* x0,
                                     void* user_data);

typedef void (*RobustEventsFn)(adouble* e, adouble* initial_states,
                               adouble* final_states, adouble* parameters,
                               adouble& t0, adouble& tf,
                               const double* theta, int ntheta,
                               adouble* xad, int iphase, Workspace* workspace);

// ---- one phase of a nominal problem ---------------------------------------------------
//
// A model IS its first phase, by inheritance below, so a single-phase problem is written
// exactly as it always was and nothing here needs to be known about. A problem of several
// phases fills RobustModel::later_phases with the rest, in order.
struct RobustPhase {
    // ---- the nominal sizes, per scenario ----------------------------------
    int nstates, ncontrols, nevents, npath;
    RowVectorXi nodes;

    // ---- the nominal maths ------------------------------------------------
    // dae and events take theta; the costs do not, being shared across the scenarios.
    // A cost that differed per scenario would be a risk measure, which needs a cost
    // state per scenario and is built by the driver instead.
    RobustDaeFn    dae;
    RobustEventsFn events;
    adouble (*endpoint_cost)(adouble* initial_states, adouble* final_states,
                             adouble* parameters, adouble& t0, adouble& tf,
                             adouble* xad, int iphase, Workspace* workspace);
    adouble (*integrand_cost)(adouble* states, adouble* controls, adouble* parameters,
                              adouble& time, adouble* xad, int iphase,
                              Workspace* workspace);

    // ---- the nominal bounds, stated once ----------------------------------
    RowVectorXd states_lower, states_upper;
    RowVectorXd controls_lower, controls_upper;
    RowVectorXd events_lower, events_upper;
    RowVectorXd path_lower, path_upper;
    double t0_lower, t0_upper, tf_lower, tf_upper;

    // ---- the nominal guess ------------------------------------------------
    MatrixXd guess_states;      // nstates x N
    MatrixXd guess_controls;    // ncontrols x N
    MatrixXd guess_time;        // 1 x N

    // Per-constraint scale factors for the violation measure, so that a metre and a
    // radian are not added together. One entry per nominal event and per nominal path
    // constraint of THIS phase. Left empty they are all one and the measure is the raw
    // infinity norm, which is right only when the constraints share units.
    RowVectorXd event_scale, path_scale;

    // Inward margins, one per constraint and per END, overriding what `tighten` would
    // derive for this phase's events and path constraints. Empty, which is the default,
    // every margin is derived as before. A single entry is broadcast to every constraint
    // in the vector; otherwise there must be one entry per nominal event or per nominal
    // path constraint of THIS phase, for the same reason the scales are per phase: they
    // are the length of the phase's own vectors.
    //
    // They exist because `tighten` is necessarily symmetric while a bound's two ends
    // often do not mean the same thing. A heating rate declared on [0, 1.0e6] has an
    // upper end that is a requirement and a lower end that is a statement of sign: the
    // quantity cannot be negative, and the design is not being asked to hold it away
    // from zero. A margin derived from the half width gives that lower end 5.0e4, which
    // a vehicle at an entry interface misses by a factor of four and can do nothing
    // about, its initial state being pinned by its own events. The design is then
    // infeasible for a reason that has nothing to do with robustness, and the solver
    // reports it as an infeasible problem rather than as a tightening that overreached.
    // Setting the lower margins to zero leaves the statement of sign alone and tightens
    // the requirement, which is what was meant.
    //
    // A margin must be finite and must not be negative: a negative one would loosen a
    // bound the user declared, which is not a margin but a different problem. A pair
    // that closes a bound's box is refused, an empty box being exactly the silent
    // infeasibility these exist to prevent. A pinned bound is an equality and takes no
    // margin from either source. An explicitly given margin IS applied to a one-sided
    // bound, which is how such a bound gets the margin the symmetric derivation has no
    // half width to produce.
    RowVectorXd events_margin_lower, events_margin_upper;
    RowVectorXd path_margin_lower,   path_margin_upper;

    RobustPhase()
        : nstates(0), ncontrols(0), nevents(0), npath(0),
          dae(0), events(0), endpoint_cost(0), integrand_cost(0),
          t0_lower(0.0), t0_upper(0.0), tf_lower(0.0), tf_upper(0.0) {}
};

struct RobustModel : public RobustPhase {
    // ---- the phases after the first ---------------------------------------
    //
    // Empty for a single-phase problem, which is what the model was until it learned about
    // phases and is still how most are written. Filled, the phases run in the order given,
    // the model itself being the first, and the driver joins them.
    //
    // THE PHASE BOUNDARY TIMES ARE SHARED BY EVERY SCENARIO, and that is forced rather than
    // chosen. PSOPT gives a phase one t0, one tf and one control, so letting scenario k
    // switch at its own time would need M separate chains of phases, and those cannot share
    // a control grid, which is what makes the design non-anticipative. A switching time is
    // therefore a here-and-now decision like the control itself. That is consistent with
    // the formulation and it excludes a boundary triggered by a state event.
    std::vector<RobustPhase> later_phases;

    // How one phase joins the next, for ONE scenario, as adouble maths. Left null the
    // driver joins them by CONTINUITY, which is what most multi-phase problems want: every
    // state of the scenario carried across unchanged, and the time carried across once,
    // being shared. Set, it writes `nlink` values the driver holds at zero, once per
    // scenario per boundary, and `jump` below has to say the same thing in numbers.
    RobustLinkFn link;
    int          nlink;        // values `link` writes, per scenario per boundary
    void*        link_data;

    // The same statement for the VERIFIER, which integrates one scenario across the
    // boundary and cannot solve an implicit linkage. Left null with `link` also null the
    // verifier carries the state across unchanged, which is continuity. Left null with
    // `link` set, the verifier says so and refuses, because carrying the state across a
    // linkage that is not continuity would verify a trajectory the design does not have.
    RobustJumpFn jump;
    void*        jump_data;

    // ---- what belongs to the whole problem and not to a phase --------------
    int         nparameters;
    RowVectorXd parameters_lower, parameters_upper;
    RowVectorXd guess_parameters;
    // The state every scenario starts from, which the warm start and the verification
    // integrate out of. It is asked for rather than read off the event bounds, because the
    // driver cannot in general tell which events are the initial conditions, and guessing
    // would put a silent error in the one place they both rest on.
    RowVectorXd initial_state;

    // The same thing when the initial state DEPENDS ON THE PARAMETER, which is the case
    // whenever the uncertainty is in where the plant starts rather than in how it behaves.
    // Set this and `initial_state` is not read.
    //
    // What it changes and what it does not. The design already knows about a parameterised
    // initial condition, because the nominal `events` receive theta and a caller writes
    // e = x(t0) - x0(theta) there as they would write any other event. What the library
    // cannot infer from that is the VALUE, so without this function the warm start and the
    // verifier integrate every scenario out of the fixed vector above while the design
    // holds each one to its own starting point, and the certificate is then about a
    // different plant from the one that was designed. That failure is silent, which is why
    // this exists.
    //
    // The function is called with a scenario and must write nstates values. It has to agree
    // with the events: the library does not check that it does, having no way to, and a
    // disagreement between them is the same silent error in a new place.
    //
    // Measured on x' = -x + u with x(0) = theta over a 3 sigma set of 0.3, brought to x = 1
    // within 0.02. The DESIGN is right either way, at t_f = 2.8134, which is the closed form
    // ln(0.3/(0.9*0.02)) this problem happens to have. What differs is the picture the
    // library forms of it: with the function set the loop stops at five scenarios with a
    // certificate of zero, which an integrator written independently of the library
    // confirms; without it the loop runs to twelve scenarios chasing a violation of 3.0e-01
    // that the same independent integrator puts at zero. The error runs pessimistic there
    // because the initial condition is itself an event, and would run optimistic on a
    // problem where it is not. What is common to both is that the certificate is not about
    // the design.
    RobustInitialStateFn initial_state_fn;
    void*                initial_state_data;   // passed to it unchanged

    // Inward tightening of the design's constraints, as a fraction of each two-sided
    // half width, without which the generation loop cannot terminate: the scenarios are
    // satisfied to the stated tolerance exactly, so the violation between two
    // neighbouring scenarios is necessarily a little larger. 0.9 keeps nine tenths.
    double tighten;

    // ---- the ancillary gain, if any ---------------------------------------
    RobustFeedbackKind feedback_kind;
    MatrixXd           feedback_gain;      // ncontrols x nstates, for CONSTANT
    RobustGainFn       feedback_schedule;  // for SCHEDULED
    void*              feedback_data;      // passed to feedback_schedule
    // Bounds on a co-designed gain's entries, each either 1x1 and broadcast or
    // ncontrols x nstates. A co-designed gain is a decision variable and nothing else
    // bounds it; an unbounded one runs away into saturation, where the realised control is
    // not the control that was designed. Required by the two co-designed kinds.
    MatrixXd           feedback_lower, feedback_upper;
    // A gain to start a co-design from. Starting at zero starts it at the open-loop
    // design, which is the expensive local minimum the gain exists to escape. Once the
    // generation loop has a gain of its own the setup starts from THAT instead, the
    // scenario being added being one the previous gain nearly served.
    MatrixXd           feedback_guess;

    // The range each scenario's cost is certainly inside, needed by the two risk
    // measures that carry that cost as a state. Asked for rather than guessed: a bound
    // that turned out to be active would silently change the risk measure into
    // something else. Left as they are, both not-a-number, ROBUST_MEAN_VARIANCE and
    // ROBUST_CVAR are refused.
    double cost_lower, cost_upper;

    // Substeps per node interval in the DEFAULT verification integrator. Sixteen was
    // measured on the arm: at eight the transcription reported a design landing exactly
    // on its tolerance ball that an independent integrator found outside it, at sixteen
    // the two agree to 2e-05, and at thirty-two to the printed precision. A setting
    // validated on the nominal problem has to be validated again on the robust one,
    // whose trajectory is longer and gentler.
    int verify_substeps;

    // Substeps per node interval in the warm start's integration, doubled on demand.
    // The horizon being integrated over is a decision variable and a robust final time
    // grows as scenarios are added, so a count that was ample at the start need not be
    // by the end; when the integration returns a non-finite value the count is doubled,
    // up to robust_warm_substeps_max, and failing that the nominal guess is used.
    int warm_substeps;

    // Where the algorithm options go, and this is not optional decoration.
    //
    // psopt_level2_setup RESETS every field of the Alg it is given to that field's
    // default. A hand-written setup therefore assigns its options AFTER calling it, and
    // the generated setup has to give the caller the same moment. This callback is that
    // moment: it runs immediately after every psopt_level2_setup, once per scenario
    // count, so whatever it sets is what the solve uses.
    //
    // Options set on the Alg before psopt_solve_robust is called are NOT honoured; they
    // are overwritten by the level 2 setup of the first iteration. The first version of
    // this interface tried to snapshot them and put them back, which silently reinstated
    // the empty strings of every field the caller had never touched and left the solve
    // with no collocation method at all. A callback has no such list to get wrong.
    void (*configure)(Alg& algorithm, void* user_data);
    void*  configure_data;

    // The initialiser list is in DECLARATION ORDER, which is the order the members are
    // actually initialised in whatever order this list is written. GCC says so with
    // -Wreorder, which CI turns into an error through PSOPT_STRICT_WARNINGS, and it is
    // worth keeping that way: a list whose order is a fiction is one place a reader can
    // be misled about what depends on what.
    RobustModel()
        : RobustPhase(),
          link(0), nlink(0), link_data(0), jump(0), jump_data(0),
          nparameters(0),
          initial_state_fn(0), initial_state_data(0),
          tighten(0.9),
          feedback_kind(ROBUST_FEEDBACK_NONE), feedback_schedule(0), feedback_data(0),
          cost_lower(std::numeric_limits<double>::quiet_NaN()),
          cost_upper(std::numeric_limits<double>::quiet_NaN()),
          verify_substeps(16), warm_substeps(8),
          configure(0), configure_data(0) {}
};

const int robust_warm_substeps_max = 256;

// Inward margins for a set of bounds: (1 - tighten) times the half width of each
// TWO-SIDED bound, and zero for a pinned or one-sided one. A pinned bound is an
// equality, and tightening it would make it infeasible instead of merely tight; a
// one-sided bound has no half width to take a fraction of, and a design left sitting on
// such a bound is exactly where the between-scenario overshoot appears, so the driver
// says so rather than silently leaving it alone.
//
// `nonefinite` is set to the number of one-sided entries found.
void robust_margins(const RowVectorXd& lower, const RowVectorXd& upper,
                    double tighten, RowVectorXd& margin, int& n_one_sided);

// The same with per-end overrides, which is what the design is actually tightened with.
// `override_lower` and `override_upper` are each empty (derive from `tighten`), a single
// entry (broadcast) or one entry per constraint; any other length is a caller error that
// robust_model_is_usable refuses before a solve starts, and that this function treats as
// absent rather than guessing which constraint a short vector meant.
//
// An explicitly given margin is applied whether the bound is two-sided or not, a pinned
// bound excepted. `n_one_sided` counts only the ends left WITHOUT a margin, those being
// the ends the driver has something to tell the user about.
void robust_margins(const RowVectorXd& lower, const RowVectorXd& upper, double tighten,
                    const RowVectorXd& override_lower, const RowVectorXd& override_upper,
                    RowVectorXd& margin_lower, RowVectorXd& margin_upper,
                    int& n_one_sided);

// Whether one override vector is a length this interface accepts for `n` constraints:
// zero (absent), one (broadcast) or n. Exposed so that a caller validating a model of
// its own asks the same question the driver asks.
bool robust_margin_override_is_sized(const RowVectorXd& override_vector, int n);

// Evaluate a RobustModel's nominal dae numerically at one scenario, on plain doubles.
//
// This is the mechanism the model interface rests on, and it is exposed because it is
// also what a user needs in order to write a verification integrator or a warm start of
// their own out of the equations they have already written once. `derivatives` and
// `path` receive the values; `path` may be null when the model declares no path
// constraints. xad and workspace are passed as null, so the nominal dae must not use
// them, as the note above requires.
void robust_dae_value(const RobustModel& model, const double* theta, int ntheta,
                      const double* states, const double* controls,
                      const double* parameters, double time,
                      double* derivatives, double* path);

//////////////////////////////////////////////////////////////////////////
///////////////////  The driver  /////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

struct RobustSpec {
    // ---- what to solve ----------------------------------------------------
    RobustUncertainty uncertainty;

    // Violation of the user's declared constraints that a design may leave. It is
    // NOT the constraint: a terminal ball of radius 0.02 belongs in the event
    // bounds, and slack is how far outside that ball the design may still stray
    // somewhere in the set. Passing the ball's own radius here would quietly
    // double it. Zero demands the declared bounds hold everywhere, which an inward
    // margin on the design is what makes attainable.
    double slack;

    // What is minimised. See RobustRisk. The scenario set does two jobs and the driver
    // keeps them apart: as a CONSTRAINT SET every member counts equally, and as a
    // QUADRATURE RULE for the risk measure the weights are the rule's. So a scenario the
    // generation loop adds is imposed in full and carries weight ZERO in the objective,
    // being the place the design was failing and therefore a biased sample of the
    // uncertainty by construction.
    RobustRisk risk;
    double     cvar_alpha;     // the CVaR level; 0.9 averages the worst tenth
    double     mv_lambda;      // the weight on the variance in mean-variance

    // The quadrature weights of the starting scenario set. Left empty, the driver takes
    // them from the unscented rule when it builds that set, or makes them uniform when
    // the caller has filled spec.scenarios. Generated scenarios append a zero.
    RowVectorXd weights;

    int      max_iterations;   // scenario-generation iterations
    int      n_seed;           // low-discrepancy seeds per worst-case search
    int      n_refine;         // cloud centres kept during refinement
    unsigned seed;             // reproducibility
    bool     verbose;

    // Decline the driver's one change to the caller's Alg. psopt_solve_robust sets
    // algorithm.free_padded_defect_rows, because the rows it frees are nstates per
    // SCENARIO of counted equalities holding no dynamics and they are the whole of what
    // caps the scenario count. Set this to keep them, and read the note on
    // Alg::free_padded_defect_rows first: keeping them means IPOPT refuses the augmented
    // problem with Not_Enough_Degrees_Of_Freedom at about
    // (ncontrols * nodes) / nstates scenarios, whatever freedom the design still has.
    bool     keep_padded_defect_rows;

    // Solve the final scenario set once more from the caller's own guess, and keep
    // that design if it certifies and is cheaper. One extra call to psopt(), two when
    // the cold answer loses and the warm one has to be put back into the caller's Sol.
    //
    // It is on by default because the warm chain the loop needs is also what limits
    // what it returns. The chain is required while the scenario set is small and
    // growing: each subproblem is nonconvex, and a guess that is merely the right
    // shape leaves the augmented problem with M copies of an infeasible arc. But it
    // follows ONE homotopy, and it finishes at a local minimum that the final scenario
    // set does not require. Measured on examples/robust_driver, the two-link arm at 25
    // nodes with a terminal ball of 0.03: the chain finishes at t_f = 8.9713 with a
    // worst violation of 6.7e-05, and a cold solve of its own twelve scenarios reaches
    // t_f = 7.6387 with NO violation anywhere in the set -- fifteen per cent faster and
    // a cleaner certificate, confirmed by a scan of 4001 payloads through an integrator
    // independent of the transcription. One extra solve buys that.
    //
    // A reported objective is therefore an upper bound on what the method can deliver,
    // and is not the value of the robust problem. The cold answer is kept only when the
    // caller's own violation function certifies it, so this cannot trade a certificate
    // for a cheaper number.
    bool     polish;
    // The rule itself is robust_prefer_cold below.

    // The nominal problem, when the caller wants the driver to build the augmented one.
    // Left null, `setup` below is required and the augmentation is the caller's. Set,
    // `setup` is not called at all and may be left null.
    // psopt_solve_robust's five-argument overload sets this.
    RobustModel* model;

    // ---- what the user supplies -------------------------------------------

    // Build `problem` and `algorithm` for this scenario list. Called before every
    // solve, with the previous solution when there is one so that a warm start can
    // be built from it -- and it should be: each subproblem is nonconvex, and a
    // guess that is merely the right shape leaves the augmented problem with M
    // copies of an infeasible arc and the solver free to wander to a distant local
    // minimum.
    //
    // The user does the level 1 and level 2 setup here, sizes the phase for
    // scenarios.size() copies of the state, sets the bounds, and points
    // problem.user_data at whatever the augmented dae and events need to read the
    // scenario list.
    void (*setup)(Prob& problem, Alg& algorithm,
                  const std::vector<RowVectorXd>& scenarios,
                  const RobustDesign& previous, void* user_data);

    // The violation of the CURRENT design at one parameter vector, by the user's own
    // independent integrator. Required when the augmentation is the caller's; optional
    // when spec.model is set, in which case leaving it null selects
    // robust_model_violation above and setting it overrides that with the caller's, which
    // is the stronger claim. Read the note on robust_model_violation before relying on
    // the default.
    double (*violation)(const RowVectorXd& theta, const RobustDesign& design,
                        void* user_data);

    void* user_data;

    // ---- what the driver reports ------------------------------------------
    std::vector<RowVectorXd> scenarios;       // the final scenario set
    RobustDesign             design;          // its trajectory, kept
    double                   certificate;     // worst violation found over the set
    RowVectorXd              certificate_at;  // and where
    long                     evaluations;     // violation calls that went into it
    int                      n_solves;        // calls to psopt()
    bool                     converged;       // nothing worse than slack was found
    // The largest amount by which the realised control of the final design leaves the
    // user's control bounds, sampled between the nodes where nothing constrains it. Zero
    // without feedback. A design whose corrections quietly ask for more actuator than
    // exists is one no plant can execute, whatever its certificate says about the states.
    double                   control_excess;
    RowVectorXd              reference_parameter;  // scenario 0 of the starting rule
    // A co-designed gain sitting on its bound is the signature of the failure the note on
    // RobustFeedbackKind records: the optimiser is paid in cost for a gain it is not
    // charged for, so it takes all the gain it is allowed. Reported and not refused,
    // because a gain on its bound is still a valid design if it certifies, and if it does
    // not this says where to look first. Zero for every kind of gain the caller gives.
    int                      gain_entries, gain_on_bound;
    double                   gain_norm;       // the largest |K| entry of the design
    bool                     own_verifier;    // true when spec.violation produced the
                                              // certificate, false when the library's
                                              // robust_model_violation did
    bool                     budget_exhausted;// a solve starved of freedom

    RobustSpec()
        : slack(0.0), risk(ROBUST_NOMINAL), cvar_alpha(0.9), mv_lambda(1.0),
          max_iterations(12), n_seed(128), n_refine(3), seed(20260927u),
          verbose(true), keep_padded_defect_rows(false), polish(true), model(0),
          setup(0), violation(0), user_data(0),
          certificate(0.0), evaluations(0), n_solves(0), converged(false),
          control_excess(0.0), gain_entries(0), gain_on_bound(0), gain_norm(0.0),
          own_verifier(false), budget_exhausted(false) {}
};

// The sizes of the augmented problem for M scenarios under a given risk measure. The
// generated setup uses this rather than computing the sizes inline, so that a caller who
// wants to check what a scenario count will cost gets the same answer the driver will
// use, and so that the arithmetic exists in one place.
//
// The layout it describes, which three separate places read back:
//
//   states      nstates per scenario, then one cost state per scenario for the two
//               measures that carry one
//   controls    the nominal controls
//   parameters  the user's, then CVaR's eta and one slack per scenario
//   events      nevents per scenario, then one pinned zero per cost state, then one
//               Rockafellar-Uryasev row per scenario under CVaR
//   path        npath per scenario
void robust_augmented_sizes(const RobustModel& model, int M, RobustRisk risk,
                            int& nstates, int& ncontrols, int& nevents,
                            int& npath, int& nparameters);

// The same for ONE PHASE of a model of several, the sizes of a phase being its own. The
// call above is this one on the model's first phase, which is the whole of a single-phase
// problem.
void robust_augmented_phase_sizes(const RobustModel& model, const RobustPhase& phase,
                                  int M, RobustRisk risk,
                                  int& nstates, int& ncontrols, int& nevents,
                                  int& npath, int& nparameters);

// Evaluate a RobustModel's nominal integrand numerically at one scenario. The companion
// of robust_dae_value, and what the generated setup uses to seed the cost states of a
// risk measure that carries them: starting them at zero would start the objective of a
// mean-variance or CVaR design at a value its own trajectory contradicts.
double robust_integrand_value(const RobustModel& model,
                              const double* states, const double* controls,
                              const double* parameters, double time);

// The state one scenario starts from: model.initial_state_fn if the caller set one, and
// model.initial_state otherwise. Exposed because a caller who writes their own verifier
// has to start it where the library's would, and because a caller who has a parameterised
// initial condition needs one place that answers the question.
void robust_initial_state(const RobustModel& model, const double* theta, int ntheta,
                          double* x0);

// Evaluate a RobustModel's nominal events numerically at one scenario. The companion of
// robust_dae_value, and needed for the same reason.
// One phase's dae evaluated numerically, for a model of several phases. The call above is
// this one on the model's first phase.
void robust_phase_dae_value(const RobustPhase& phase, const RobustModel& model,
                            const double* theta, int ntheta,
                            const double* states, const double* controls,
                            const double* parameters, double time,
                            double* derivatives, double* path);

void robust_events_value(const RobustModel& model, const double* theta, int ntheta,
                         const double* initial_states, const double* final_states,
                         const double* parameters, double t0, double tf,
                         double* e);

// How a design's control is to be read BETWEEN the stored values. The transcription
// decides it, and everything that reproduces the designed control has to use the same
// reading or it integrates a different controller: the library's verifier, the warm start,
// and a verification integrator a caller writes of their own.
//
// Exposed for the same reason robust_phase_dae_value is. RobustDesign's note on
// controls_full tells a caller that anything reproducing the designed control must prefer
// the complete table, and gives the measured cost of not doing so; until these were public
// the header asked for that and supplied nothing to do it with, and the library's own warm
// start was reading the nodal table as a chord while its verifier read the parabola. Two
// routes that must agree about one controller should not each carry their own arithmetic.
//
//   ROBUST_HELD      the stored value, held across the interval
//   ROBUST_LINEAR    the chord between the two stored values
//   ROBUST_PARABOLA  the parabola through node, midpoint and node, from controls_full
//
// `robust_control_shape` asks the algorithm and the design together: the parabola needs a
// complete table to read, so a design that reports none is read as a chord whatever the
// transcription was. A caller that wants the parabola should check for it, because a
// silent fall back to the chord is charged to the design.
enum RobustControlShape { ROBUST_HELD, ROBUST_LINEAR, ROBUST_PARABOLA };

RobustControlShape robust_control_shape(Alg& algorithm, const MatrixXd& controls_full,
                                        const MatrixXd& time_full);

// The control on interval i at fraction w of it, w running from zero at node i to one at
// node i+1. `nc` is the number of ROWS of the design's control table, which under a
// co-designed gain schedule is wider than the user's control vector: the gain rides in the
// trailing rows, put there by the transcription, and a caller reading only the user's
// controls would read a design whose gain never varies.
void robust_control_at(const MatrixXd& controls, const MatrixXd& controls_full,
                       RobustControlShape shape, int i, double w, int nc, double* u);
void robust_control_at(const RobustPhaseTrajectory& traj, RobustControlShape shape,
                       int i, double w, int nc, double* u);

// How badly a design serves one parameter vector, by an integrator built from the model's
// own equations. This is what psopt_solve_robust uses when the caller supplies a model and
// leaves spec.violation null, and it is exposed so that a caller can call it directly,
// compare it against a verifier of their own, or use it as the starting point for one.
//
// It integrates the plant from model.initial_state with a fixed-step RK4 at
// model.verify_substeps steps per node interval, reading the control the way `algorithm`
// means it: held under a constant parameterisation, the ramp under a linear one, and the
// parabola through node, midpoint and node for the two that carry a midpoint control,
// taken from design.controls_full. It returns the largest scaled amount by which the
// trajectory falls outside the DECLARED event and path bounds, and zero when it stays
// inside all of them. The declared bounds are the model's own, not the tightened ones the
// design was solved against: the tightening is a margin for the design and not a change
// to what is being verified.
//
// WHAT INDEPENDENCE THIS DOES AND DOES NOT GIVE. It is independent of the TRANSCRIPTION:
// a different integrator, a different step control, its own code path, and it shares with
// the solve only the user's equations. It is not independent of the LIBRARY. A design
// checked against the integrator that produced it checks nothing, and this is not that,
// but a caller who wants the stronger claim writes their own and sets spec.violation, as
// python/examples/robust_arm.py does on the Python side. The driver says which of the two
// produced a certificate.
// Under an ancillary gain two further arguments matter, and both default so that code
// written against the open-loop version keeps compiling. `theta_ref` is the parameter of
// the reference trajectory, scenario 0 of the rule, which the verifier integrates beside
// the plant because the correction is a deviation from it: passing null with a gain in
// place makes the reference the plant itself, which is not the controller that was
// designed. `control_excess` receives the largest amount by which the REALISED control
// leaves the user's control bounds, which is reported separately rather than added to the
// violation, because a design whose corrections ask for more actuator than exists is a
// different fault from one that misses a constraint, and the two want different remedies.
// The history the verifier walks, handed back so that a caller can show what was checked
// without integrating the design a second time.
//
// WHY THIS IS IN THE LIBRARY AND NOT IN THE CALLER. A certificate is one number and a
// picture of a certificate is a history, so a caller who wants to display what the number is
// about has two routes: re-integrate the design, or be given what the verifier already
// computed. The first route is a second integrator, and a second integrator is a second
// CONTROLLER unless every detail agrees: the substep count, the control shape, the realised
// control under an ancillary gain, the initial-state hook, and the normalisation of the
// bound excess. Several of those are not public, so a caller writing their own would
// reproduce them by eye, and a figure that disagrees with the certificate printed beside it
// is worse than no figure at all. Recording costs the verifier nothing it was not already
// doing.
//
// One entry per phase, in order, appended in phase order. Each column but the last is one
// RK4 STEP START, so a phase with N nodes and `verify_substeps` steps per node interval
// gives (N-1)*verify_substeps + 1 columns. `controls` is the REALISED control, the one the
// plant saw, which under an ancillary gain is not the designed one. `path` is empty when the
// phase declares no path constraints.
//
// THE LAST COLUMN IS RECORDED AND NOT CHECKED, and the distinction matters. The verifier
// evaluates the path constraints at every step start, so the last point it checks is one
// substep before the final node; the final node's own path values are therefore outside the
// certificate. That is the behaviour this library has always had, and adding the final node
// to the check would move every certificate already in print, so it is not changed here. The
// final column is filled so that a caller can see where the trajectory ended and can see
// that value for itself.
//
// A verification that returns a non-finite violation abandons the integration, so the
// history is then whatever had been walked when it gave up and should not be read as a
// trajectory.
struct RobustVerifyHistory {
    MatrixXd time;      // 1 x M
    MatrixXd states;    // nstates x M
    MatrixXd controls;  // ncontrols x M, the REALISED control
    MatrixXd path;      // npath x M, or empty when the phase declares none
};

double robust_model_violation(const RobustModel& model, Alg& algorithm,
                              const RowVectorXd& theta, const RobustDesign& design,
                              const RowVectorXd* theta_ref = 0,
                              double* control_excess = 0,
                              std::vector<RobustVerifyHistory>* history = 0);

// Which of two candidate designs the polish step keeps: prefer one the caller's own
// violation function certifies, and among certified designs the cheaper; if neither
// certifies, the one that comes closer.
//
// Exposed and named because it is the single decision in this driver that could
// silently trade a certificate for a cheaper objective, and a rule that can be tested
// on its own is one that can be relied on. A cold design that undercuts the warm one on
// the objective while violating the constraints somewhere in the set is not a design,
// and this returns false for it however large the saving.
bool robust_prefer_cold(double cold_objective, double cold_worst,
                        double warm_objective, double warm_worst, double slack);

// Design a control that serves every parameter in the uncertainty set.
//
// Starts from the unscented scenario set, or from spec.scenarios if the caller has
// filled it. Then: solve; find the parameter the design serves worst over the
// whole set; if that is worse than the slack, add it as a scenario and solve
// again. That is a cutting plane on a semi-infinite constraint, and it is what
// lets the loop finish with a certificate rather than a statistic.
//
// Returns the psopt() return code of the design in `solution`.
[[nodiscard]] int psopt_solve_robust(Sol& solution, RobustSpec& spec,
                                     Prob& problem, Alg& algorithm);

// The same driver, with the augmentation built from a nominal model instead of by the
// caller's own setup. Equivalent to setting spec.model and calling the overload above.
// spec.violation is still the caller's: the verification has to be independent of the
// transcription, and in this version the library does not offer one.
[[nodiscard]] int psopt_solve_robust(Sol& solution, RobustSpec& spec,
                                     RobustModel& model,
                                     Prob& problem, Alg& algorithm);

#endif // PSOPT_ROBUST_H
