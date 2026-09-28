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
// WHAT THE USER STILL WRITES, AND WHY
//
// In the Python interface the driver builds the augmented problem itself, because
// the user's equations arrive as CasADi expressions it can replicate. In C++ they
// arrive as dae() and events() with fixed signatures, taped by CppAD, and nothing
// can rewrite them. So in C++ the AUGMENTATION is the user's -- the dae loops over
// the scenarios it finds in problem.user_data, exactly as examples/robust_arm.cxx
// does -- and the DRIVER is the library's. That division is the honest one for
// this language, and it is where the work actually is: the scenario rule, the
// oracle, the warm start and the loop are the same for every problem, while the
// augmented dae is three lines of the user's own.
//
// Copyright (c) Victor M. Becerra, 2026. Part of the PSOPT library (LGPL).

#ifndef PSOPT_ROBUST_H
#define PSOPT_ROBUST_H

#include "psopt.h"

#include <vector>
#include <cmath>
#include <cstddef>

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
struct RobustDesign {
    MatrixXd time;          // 1 x N
    MatrixXd controls;      // ncontrols x N
    MatrixXd states;        // nstates x N   (all scenario copies, as solved)
    MatrixXd parameters;    // static parameters, if any
    double   objective;
    bool     valid;

    RobustDesign() : objective(0.0), valid(false) {}
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

    // The violation of the CURRENT design at one parameter vector, by the user's
    // own independent integrator.
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
    bool                     budget_exhausted;// a solve starved of freedom

    RobustSpec()
        : slack(0.0), max_iterations(12), n_seed(128), n_refine(3), seed(20260927u),
          verbose(true), keep_padded_defect_rows(false),
          setup(0), violation(0), user_data(0),
          certificate(0.0), evaluations(0), n_solves(0), converged(false),
          budget_exhausted(false) {}
};

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

#endif // PSOPT_ROBUST_H
