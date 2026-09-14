//////////////////////////////////////////////////////////////////////////////
// test_multiple_shooting.cpp
//
// The multiple-shooting transcription: segment-start states and controls are
// the decision variables, the trajectory between them is produced by a
// fixed-step RK4 recorded on the same tape as everything else, and the rows
// that carry the collocation defect under every other transcription carry the
// matching condition x_{k+1} - phi(x_k, u_k) = 0 instead.
//
// The first test is the one worth having. On the minimum-energy double
// integrator the transcription's OWN discrete problem has a closed-form answer:
// with the control constrained to be piecewise constant on M equal segments,
//
//     min (1/2) int_0^1 u^2 dt   s.t.  xddot = u,  (0,0) -> (1,0)
//
// is a least-norm problem in the M-dimensional space of piecewise-constant
// controls, and projecting the two linear constraints onto that space gives
//
//     J*(M) = 6 M^2 / (M^2 - 1),
//
// which is 6.25 at M=5 and tends to the continuous optimum 6 as M grows. A
// transcription that returns this number to machine precision is solving the
// problem it defines exactly, and one that returns something close to 6 instead
// is solving a different problem and getting a plausible answer -- which is the
// failure mode worth separating out, because it looks like success.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>

#include <cmath>
#include <string>

namespace msh {

// case 0: minimum energy, fixed tf, closed-form discrete optimum.
// case 1: harmonic oscillator, for which RK4 is not exact, so the integrator's
//         own order can be measured.
// case 2: minimum time, free tf, |u| <= 1, (0,0) -> (1,0), tf* = 2.
// case 3: Bryson-Denham, min (1/2)int u^2 with x <= 1/9 and J* = 4, whose active arc is what
//         shows whether a path constraint imposed at the segment boundaries is the constraint
//         the user wrote.
// case 4: minimum time with ASYMMETRIC control bounds, u in [-1,2], so the single switch
//         falls at tf/3 rather than tf/2. tf* = sqrt(3).
static int g_case = 0;

const double TF_SQRT3 = 1.7320508075688772;

adouble endpoint_cost(adouble*, adouble*, adouble*, adouble&, adouble& tf, adouble*,
                      int, Workspace*)
{ return ( g_case == 2 || g_case == 4 ) ? tf : (adouble) 0.0; }

adouble integrand_cost(adouble*, adouble* u, adouble*, adouble&, adouble*, int, Workspace*)
{ return ( g_case == 2 || g_case == 4 ) ? (adouble) 0.0 : 0.5*u[0]*u[0]; }

void dae(adouble* d, adouble* path, adouble* s, adouble* c, adouble*, adouble&,
         adouble*, int, Workspace*)
{
    d[0] = s[1];
    d[1] = ( g_case == 1 ) ? ( -9.0*s[0] + c[0] ) : c[0];
    if ( g_case == 3 ) path[0] = s[0];
}

void events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&,
            adouble*, int, Workspace*)
{ e[0] = i[0]; e[1] = i[1]; e[2] = f[0]; e[3] = f[1]; }

void linkages(adouble*, adouble*, Workspace*) {}

// The same problem as case 0, split into two phases joined by a linkage, so that the answer
// can be asked to be independent of where the split is drawn.
static int g_two_phase_events = 0;

void events2(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&,
             adouble*, int iphase, Workspace*)
{
    (void) g_two_phase_events;
    if ( iphase == 1 ) { e[0] = i[0]; e[1] = i[1]; }
    else               { e[0] = f[0]; e[1] = f[1]; }
}

void linkages2(adouble* l, adouble* xad, Workspace* w)
{
    adouble xf[2], x0[2];
    get_final_states(xf, xad, 1, w);
    get_initial_states(x0, xad, 2, w);
    l[0] = xf[0] - x0[0];
    l[1] = xf[1] - x0[1];
    l[2] = get_final_time(xad,1,w) - get_initial_time(xad,2,w);
}

struct Run {
    int flag; double J; double tf; double err_est;
    int nvars;                  // mesh_stats[0].nvars: the decision vector's actual length
    MatrixXd u_nodes;           // controls at the segment boundaries
    MatrixXd u_hs, t_hs;        // node and midpoint controls interleaved, and their times
    MatrixXd x_nodes, t_nodes;  // the states at the segment boundaries, and their times
    int segments;               // how many segments the phase ended with
};

static Run solve(int which, int segments, int steps,
                 const std::string& upar = "constant", int path_samples = 0,
                 bool flexible_segments = false,
                 const std::string& integrator = "RK4",
                 bool automatic = false, double refine_tol = 1.0e-3)
{
    g_case = which;

    Alg algorithm; Sol solution; Prob problem;
    Run out; out.flag = -1; out.J = 0.0; out.tf = 0.0; out.err_est = 0.0; out.nvars = 0;
    out.segments = 0;

    const int nodes = segments + 1;

    problem.name        = "multiple shooting";
    problem.outfilename = "test_multiple_shooting.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = ( which == 3 ) ? 1 : 0;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -5.0, -5.0;
    problem.phases(1).bounds.upper.states   <<  5.0,  5.0;
    problem.phases(1).bounds.lower.controls(0) = ( which == 2 || which == 4 ) ? -1.0 : -30.0;
    problem.phases(1).bounds.upper.controls(0) = ( which == 2 ) ?  1.0
                                               : ( which == 4 ) ?  2.0 : 30.0;
    if ( which == 3 ) {
        problem.phases(1).bounds.lower.path(0) = -5.0;
        problem.phases(1).bounds.upper.path(0) =  1.0/9.0;
        problem.phases(1).bounds.lower.events << 0.0, 1.0, 0.0, -1.0;
        problem.phases(1).bounds.upper.events << 0.0, 1.0, 0.0, -1.0;
    }
    else {
        problem.phases(1).bounds.lower.events   << 0.0, 0.0, 1.0, 0.0;
        problem.phases(1).bounds.upper.events   << 0.0, 0.0, 1.0, 0.0;
    }
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = ( which == 2 || which == 4 ) ? 0.5 : 1.0;
    problem.phases(1).bounds.upper.EndTime   = ( which == 2 || which == 4 ) ? 8.0 : 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(2, nodes);
    if ( which == 3 ) problem.phases(1).guess.states.row(1) = linspace( 1.0, -1.0, nodes);
    else              problem.phases(1).guess.states.row(0) = linspace( 0.0,  1.0, nodes);
    problem.phases(1).guess.controls = zeros(1, nodes);
    problem.phases(1).guess.time     = linspace(0.0, ( which == 2 ) ? 2.0
                                                   : ( which == 4 ) ? 1.73 : 1.0, nodes);

    algorithm.nlp_method            = "IPOPT";
    algorithm.scaling               = "automatic";
    algorithm.derivatives           = "automatic";
    algorithm.nlp_iter_max          = 2000;
    algorithm.nlp_tolerance         = 1.0e-10;
    algorithm.print_level           = 0;
    algorithm.mesh_refinement       = automatic ? "automatic" : "manual";
    algorithm.mr_max_iterations     = 7;
    algorithm.collocation_method    = "Hermite-Simpson";
    algorithm.transcription_method  = "multiple-shooting";
    algorithm.ms_steps_per_segment  = steps;
    algorithm.ms_control_parameterisation = upar;
    algorithm.ms_path_samples             = path_samples;
    algorithm.ms_flexible_segments        = flexible_segments;
    algorithm.ms_integrator               = integrator;
    algorithm.ms_refine_tolerance         = refine_tol;

    out.flag = psopt(solution, problem, algorithm);
    if (out.flag == 0) {
        out.J = solution.cost;
        MatrixXd t = solution.get_time_in_phase(1);
        out.tf = t(0, (int) t.cols() - 1);
        MatrixXd e = solution.get_relative_local_error_in_phase(1);
        for (int q = 0; q < e.size(); q++) out.err_est = std::max(out.err_est, std::fabs(e(q)));
        out.nvars   = solution.mesh_stats[0].nvars;
        out.u_nodes = solution.get_controls_in_phase(1);
        out.u_hs    = solution.get_hs_controls_in_phase(1);
        out.t_hs    = solution.get_hs_time_in_phase(1);
        out.x_nodes = solution.get_states_in_phase(1);
        out.t_nodes = solution.get_time_in_phase(1);
        out.segments = (int) out.t_nodes.cols() - 1;
    }
    return out;
}

// ---------------------------------------------------------------------------
// Case 1's segment map in closed form. xddot = -9x + u with u held over the
// segment is linear with a constant forcing, so writing y = x1 - u/9 gives
// ybar'' = -9 ybar and the flow over a span D is
//
//     y(D)  =  y_k cos 3D + (x2_k/3) sin 3D,     x2(D) = -3 y_k sin 3D + x2_k cos 3D.
//
// This is what makes case 1 the right problem for an integrator test: the true
// segment map is available without a reference integration, so the scheme's own
// error can be measured rather than estimated.
static void exact_segment_map(double x1, double x2, double u, double D,
                              double& x1e, double& x2e)
{
    const double y = x1 - u/9.0;
    const double s3 = std::sin(3.0*D), c3 = std::cos(3.0*D);
    x1e = y*c3 + (x2/3.0)*s3 + u/9.0;
    x2e = -3.0*y*s3 + x2*c3;
}

// The largest per-segment error of the scheme, in the SAME normalisation PSOPT
// reports: divided by w_i = max( max_k |x_i(t_k)|, max_k |xdot_i(t_k)| ) + 1.
// Comparing an estimate with a differently normalised truth measures nothing.
static double true_relative_local_error(const Run& r, int M)
{
    double w1 = 0.0, w2 = 0.0;
    for (int k = 0; k <= M; k++) {
        const double uk = r.u_nodes(0, (k < M) ? k : M-1);
        w1 = std::max(w1, std::max(std::fabs(r.x_nodes(0,k)), std::fabs(r.x_nodes(1,k))));
        w2 = std::max(w2, std::max(std::fabs(r.x_nodes(1,k)),
                                   std::fabs(-9.0*r.x_nodes(0,k) + uk)));
    }
    w1 += 1.0; w2 += 1.0;
    double worst = 0.0;
    for (int k = 0; k < M; k++) {
        double x1e, x2e;
        exact_segment_map(r.x_nodes(0,k), r.x_nodes(1,k), r.u_nodes(0,k),
                          r.t_nodes(0,k+1) - r.t_nodes(0,k), x1e, x2e);
        worst = std::max(worst, std::max(std::fabs(x1e - r.x_nodes(0,k+1))/w1,
                                         std::fabs(x2e - r.x_nodes(1,k+1))/w2));
    }
    return worst;
}

static double discrete_optimum(double M) { return 6.0*M*M/(M*M - 1.0); }

// Case 1's continuous optimum, from the controllability Gramian: for
// xddot = -w^2 x + u driven from rest to (1,0) on [0,1] the minimum-energy cost
// is (1/2) x_T' W^-1 x_T. With w = 3 this is 8.600..., which is NOT case 0's 6 --
// a distinction that cost one wrong assertion before it was noticed.
static double oscillator_optimum()
{
    const double w = 3.0, T = 1.0;
    const double s2 = T/2.0 - std::sin(2*w*T)/(4*w);
    const double c2 = T/2.0 + std::sin(2*w*T)/(4*w);
    const double sc = (1.0 - std::cos(2*w*T))/(4*w);
    const double W11 = s2/(w*w), W12 = sc/w, W22 = c2;
    return 0.5*W22/(W11*W22 - W12*W12);
}

// The parabola through (u_k, ubar_k, u_{k+1}) at local coordinate x in [0,1] of segment k,
// read from the interleaved arrays the solution reports. This is the control the segment
// integrator was handed, so anything asked of the control is asked of this.
static double quad_control(const Run& r, int k, double x)
{
    const double L0 =  2.0*(x-0.5)*(x-1.0);
    const double Lm = -4.0*x*(x-1.0);
    const double L1 =  2.0*x*(x-0.5);
    return L0*r.u_hs(0,2*k) + Lm*r.u_hs(0,2*k+1) + L1*r.u_hs(0,2*k+2);
}

// ---------------------------------------------------------------------------
// The minimum-energy TRIPLE integrator, whose optimal control is exactly a
// parabola. The costates satisfy lam1' = 0, lam2' = -lam1, lam3' = -lam2 and
// u = -lam3, so u* is quadratic in t; equivalently, the minimiser of int u^2
// subject to three linear functionals of u lies in the span of their Riesz
// representers, which are {1, 1-t, (1-t)^2/2} and no higher. A parabola meeting
// the three terminal conditions is unique, so the two coincide and the
// continuous optimum is reachable by a piecewise-quadratic control at ANY
// number of segments.
// ---------------------------------------------------------------------------

adouble endpoint_cost3(adouble*, adouble*, adouble*, adouble&, adouble&, adouble*,
                       int, Workspace*) { return (adouble) 0.0; }
adouble integrand_cost3(adouble*, adouble* u, adouble*, adouble&, adouble*, int, Workspace*)
{ return 0.5*u[0]*u[0]; }
void dae3(adouble* d, adouble*, adouble* s, adouble* c, adouble*, adouble&,
          adouble*, int, Workspace*)
{ d[0] = s[1]; d[1] = s[2]; d[2] = c[0]; }
void events3(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&,
             adouble*, int, Workspace*)
{ e[0]=i[0]; e[1]=i[1]; e[2]=i[2]; e[3]=f[0]; e[4]=f[1]; e[5]=f[2]; }

// (0,0,0) -> (1,0,0) on [0,1]: u*(t) = 60 - 360 t + 360 t^2 and J* = 360 exactly.
const double TRIPLE_JSTAR = 360.0;

static Run solve_triple(int segments, int steps, const std::string& upar)
{
    Alg algorithm; Sol solution; Prob problem;
    Run out; out.flag = -1; out.J = 0.0; out.tf = 0.0; out.err_est = 0.0; out.nvars = 0;

    const int nodes = segments + 1;

    problem.name        = "multiple shooting, triple integrator";
    problem.outfilename = "test_multiple_shooting_3.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 3;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 6;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -50.0, -500.0, -5000.0;
    problem.phases(1).bounds.upper.states   <<  50.0,  500.0,  5000.0;
    problem.phases(1).bounds.lower.controls(0) = -1.0e5;
    problem.phases(1).bounds.upper.controls(0) =  1.0e5;
    problem.phases(1).bounds.lower.events   << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
    problem.phases(1).bounds.upper.events   << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0; problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &integrand_cost3;
    problem.endpoint_cost  = &endpoint_cost3;
    problem.dae            = &dae3;
    problem.events         = &events3;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(3, nodes);
    problem.phases(1).guess.controls = zeros(1, nodes);
    problem.phases(1).guess.time     = linspace(0.0, 1.0, nodes);

    algorithm.nlp_method            = "IPOPT";
    algorithm.scaling               = "automatic";
    algorithm.derivatives           = "automatic";
    algorithm.nlp_iter_max          = 2000;
    algorithm.nlp_tolerance         = 1.0e-12;
    algorithm.print_level           = 0;
    algorithm.mesh_refinement       = "manual";
    algorithm.collocation_method    = "Hermite-Simpson";
    algorithm.transcription_method  = "multiple-shooting";
    algorithm.ms_steps_per_segment  = steps;
    algorithm.ms_control_parameterisation = upar;

    out.flag = psopt(solution, problem, algorithm);
    if (out.flag == 0) {
        out.J     = solution.cost;
        out.nvars = solution.mesh_stats[0].nvars;
        out.u_hs  = solution.get_hs_controls_in_phase(1);
        out.t_hs  = solution.get_hs_time_in_phase(1);
    }
    return out;
}

// The same minimum-energy problem in two phases of `segments` each, joined by a linkage.
static Run solve_two_phase(int segments, int steps)
{
    g_case = 0;

    Alg algorithm; Sol solution; Prob problem;
    Run out; out.flag = -1; out.J = 0.0; out.tf = 0.0; out.err_est = 0.0;

    const int nodes = segments + 1;

    problem.name        = "multiple shooting, two phases";
    problem.outfilename = "test_multiple_shooting_2p.txt";
    problem.nphases     = 2;
    problem.nlinkages   = 3;
    psopt_level1_setup(problem);

    for (int p = 1; p <= 2; p++) {
        problem.phases(p).nstates   = 2;
        problem.phases(p).ncontrols = 1;
        problem.phases(p).nevents   = 2;
        problem.phases(p).npath     = 0;
        problem.phases(p).nodes     << nodes;
    }
    psopt_level2_setup(problem, algorithm);

    for (int p = 1; p <= 2; p++) {
        problem.phases(p).bounds.lower.states   << -5.0, -5.0;
        problem.phases(p).bounds.upper.states   <<  5.0,  5.0;
        problem.phases(p).bounds.lower.controls(0) = -30.0;
        problem.phases(p).bounds.upper.controls(0) =  30.0;
        problem.phases(p).guess.states   = zeros(2, nodes);
        problem.phases(p).guess.controls = zeros(1, nodes);
    }
    problem.phases(1).bounds.lower.events << 0.0, 0.0;
    problem.phases(1).bounds.upper.events << 0.0, 0.0;
    problem.phases(2).bounds.lower.events << 1.0, 0.0;
    problem.phases(2).bounds.upper.events << 1.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 0.5; problem.phases(1).bounds.upper.EndTime   = 0.5;
    problem.phases(2).bounds.lower.StartTime = 0.5; problem.phases(2).bounds.upper.StartTime = 0.5;
    problem.phases(2).bounds.lower.EndTime   = 1.0; problem.phases(2).bounds.upper.EndTime   = 1.0;
    problem.phases(1).guess.states.row(0) = linspace(0.0, 0.5, nodes);
    problem.phases(1).guess.time          = linspace(0.0, 0.5, nodes);
    problem.phases(2).guess.states.row(0) = linspace(0.5, 1.0, nodes);
    problem.phases(2).guess.time          = linspace(0.5, 1.0, nodes);

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events2;
    problem.linkages       = &linkages2;

    algorithm.nlp_method            = "IPOPT";
    algorithm.scaling               = "automatic";
    algorithm.derivatives           = "automatic";
    algorithm.nlp_iter_max          = 2000;
    algorithm.nlp_tolerance         = 1.0e-10;
    algorithm.print_level           = 0;
    algorithm.mesh_refinement       = "manual";
    algorithm.collocation_method    = "Hermite-Simpson";
    algorithm.transcription_method  = "multiple-shooting";
    algorithm.ms_steps_per_segment  = steps;

    out.flag = psopt(solution, problem, algorithm);
    if (out.flag == 0) out.J = solution.cost;
    return out;
}

} // namespace msh


// ---------------------------------------------------------------------------
// The transcription attains the exact optimum of the problem it defines.
// ---------------------------------------------------------------------------

TEST(MultipleShooting, ItAttainsTheExactOptimumOfItsOwnDiscreteProblem)
{
    const int segs[4] = { 5, 10, 20, 40 };
    for (int q = 0; q < 4; q++) {
        const msh::Run r = msh::solve(0, segs[q], 10);
        ASSERT_EQ(r.flag, 0) << "failed at " << segs[q] << " segments";
        EXPECT_NEAR(r.J, msh::discrete_optimum(segs[q]), 1.0e-8)
            << segs[q] << " segments: got " << r.J
            << " against " << msh::discrete_optimum(segs[q]);
    }
}


// ---------------------------------------------------------------------------
// And the answer does not depend on the step count HERE, which is not a
// tautology but a property of this problem: the dynamics are a chain of
// integrators driven by a control that is constant across a segment, so the
// exact solution over a segment is a cubic in time and RK4 is exact for it.
// A transcription whose answer moved with the step count on this problem would
// be integrating something other than what it was given.
// ---------------------------------------------------------------------------

TEST(MultipleShooting, TheStepCountDoesNotMatterWhereTheSchemeIsExact)
{
    const msh::Run coarse = msh::solve(0, 10, 1);
    const msh::Run fine   = msh::solve(0, 10, 40);

    ASSERT_EQ(coarse.flag, 0);
    ASSERT_EQ(fine.flag,   0);

    EXPECT_NEAR(coarse.J, fine.J, 1.0e-9);
    EXPECT_NEAR(coarse.J, msh::discrete_optimum(10.0), 1.0e-8);
}


// ---------------------------------------------------------------------------
// Where the scheme is NOT exact, the reported discretisation error must fall
// like the fourth power of the step, because that is the order of the scheme and
// because the error of a shooting transcription IS the integrator's error: the
// matching conditions hold the segment ends exactly, so there is no defect
// between the nodes to contribute anything else.
//
// The estimate is Richardson's -- the same segment run again at half the step,
// the difference divided by 2^4 - 1 -- so this test is also what pins the
// estimator itself to the scheme it is estimating.
// ---------------------------------------------------------------------------

TEST(MultipleShooting, TheSegmentIntegratorConvergesAtFourthOrder)
{
    const msh::Run s2 = msh::solve(1, 10, 2);
    const msh::Run s4 = msh::solve(1, 10, 4);
    const msh::Run s8 = msh::solve(1, 10, 8);

    ASSERT_EQ(s2.flag, 0);
    ASSERT_EQ(s4.flag, 0);
    ASSERT_EQ(s8.flag, 0);

    ASSERT_GT(s8.err_est, 0.0) << "the error estimate is identically zero";

    const double r1 = s2.err_est/s4.err_est;
    const double r2 = s4.err_est/s8.err_est;

    EXPECT_GT(r1, 8.0)  << "halving the step reduced the error by only " << r1;
    EXPECT_LT(r1, 32.0) << "halving the step reduced the error by " << r1;
    EXPECT_GT(r2, 12.0) << "halving the step reduced the error by only " << r2;
    EXPECT_LT(r2, 20.0) << "halving the step reduced the error by " << r2;
}


// ---------------------------------------------------------------------------
// A free final time. The segment duration enters the integrator, so this is the
// first thing that needs the derivative path to be more than a formality: the
// end state of every segment depends on tf through its own step length.
// ---------------------------------------------------------------------------

TEST(MultipleShooting, TheFinalTimeCanBeFree)
{
    const msh::Run r = msh::solve(2, 20, 10);

    ASSERT_EQ(r.flag, 0);
    EXPECT_NEAR(r.tf, 2.0, 1.0e-6) << "tf = " << r.tf;
}


// ---------------------------------------------------------------------------
// A piecewise-LINEAR control, and the cleanest demonstration there is of why the
// control parameterisation is not a detail. The optimal control of the
// minimum-energy double integrator is u*(t) = 6 - 12t, exactly linear, so a
// linear parameterisation CONTAINS the answer: the transcription attains the
// continuous optimum 6 at five segments, to machine precision, where the
// piecewise-constant form needs infinitely many (6.25, 6.061, 6.015, ...).
//
// A shooting method compared against collocation with a piecewise-constant
// control loses, and loses for a reason that has nothing to do with shooting.
// ---------------------------------------------------------------------------

TEST(MultipleShooting, ALinearControlIsExactWhereTheOptimalControlIsLinear)
{
    const int segs[3] = { 5, 10, 20 };
    for (int q = 0; q < 3; q++) {
        const msh::Run r = msh::solve(0, segs[q], 10, "linear");
        ASSERT_EQ(r.flag, 0) << "failed at " << segs[q] << " segments";
        EXPECT_NEAR(r.J, 6.0, 1.0e-9)
            << segs[q] << " segments: got " << r.J << " for a control the "
            << "parameterisation can represent exactly";
    }
}


// ---------------------------------------------------------------------------
// A path constraint imposed only at the segment boundaries is not the path
// constraint the user wrote, and the way it fails is the way that matters: it
// fails in the direction that looks like success.
//
// Bryson and Denham's problem has J* = 4 with x <= 1/9, and no feasible
// trajectory can cost less. Enforced at the boundaries only, multiple shooting
// returns 3.993 at ten segments and 3.998 at twenty -- BELOW the optimum, which
// is the constraint leaking between the boundaries and nothing else. Sampling
// inside the segments puts the answer back above J*, where a restricted control
// parameterisation must leave it.
// ---------------------------------------------------------------------------

TEST(MultipleShooting, APathConstraintLeaksBetweenTheSegmentBoundaries)
{
    const msh::Run edges   = msh::solve(3, 10, 10, "linear", 0);
    const msh::Run sampled = msh::solve(3, 10, 10, "linear", 2);

    ASSERT_EQ(edges.flag,   0);
    ASSERT_EQ(sampled.flag, 0);

    EXPECT_LT(edges.J, 4.0)
        << "the boundary-only form returned " << edges.J << ", which is not below the "
        << "optimum -- this test is then measuring something else";
    EXPECT_GT(sampled.J, 4.0)
        << "sampling inside the segments left the cost below the optimum: " << sampled.J;
    EXPECT_LT(sampled.J, 4.01) << "J = " << sampled.J;
}


// ---------------------------------------------------------------------------
// Multiple phases, joined by linkages. The invariant worth testing is not that
// two phases solve, but that the answer does not depend on where the phase
// boundary was drawn: one phase of twenty segments and two phases of ten are the
// same discretisation of the same problem, and both must return the closed-form
// discrete optimum.
// ---------------------------------------------------------------------------

TEST(MultipleShooting, TwoPhasesAgreeWithOne)
{
    const msh::Run one = msh::solve(0, 20, 10);
    const msh::Run two = msh::solve_two_phase(10, 10);

    ASSERT_EQ(one.flag, 0);
    ASSERT_EQ(two.flag, 0) << "the two-phase problem failed to solve";

    EXPECT_NEAR(two.J, one.J, 1.0e-9)
        << "two phases " << two.J << " against one phase " << one.J;
    EXPECT_NEAR(two.J, msh::discrete_optimum(20.0), 1.0e-8);
}


// ---------------------------------------------------------------------------
// The costates. They are not the multipliers of the matching conditions: those
// are a discrete adjoint, but the covector mapping PSOPT applies to collocation
// defects has a quadrature weight and a differentiation matrix in it and is not
// written for them -- applied to them it returned +96 where the answer is -12.
// They come instead from integrating the adjoint equation backwards along the
// converged primal, which is defined for any transcription that produces a
// trajectory, and which PSOPT already does for the residual-box solves.
//
// The closed form for this problem is l1 = -12 and l2 = 12t - 6.
// ---------------------------------------------------------------------------

TEST(MultipleShooting, TheCostatesAreRecoveredFromThePrimal)
{
    msh::g_case = 0;

    Alg algorithm; Sol solution; Prob problem;
    const int segments = 20, nodes = segments + 1;

    problem.name        = "multiple shooting costates";
    problem.outfilename = "test_multiple_shooting_costates.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -5.0, -5.0;
    problem.phases(1).bounds.upper.states   <<  5.0,  5.0;
    problem.phases(1).bounds.lower.controls(0) = -30.0;
    problem.phases(1).bounds.upper.controls(0) =  30.0;
    problem.phases(1).bounds.lower.events   << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.upper.events   << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0;
    problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &msh::integrand_cost;
    problem.endpoint_cost  = &msh::endpoint_cost;
    problem.dae            = &msh::dae;
    problem.events         = &msh::events;
    problem.linkages       = &msh::linkages;

    problem.phases(1).guess.states   = zeros(2, nodes);
    problem.phases(1).guess.states.row(0) = linspace(0.0, 1.0, nodes);
    problem.phases(1).guess.controls = zeros(1, nodes);
    problem.phases(1).guess.time     = linspace(0.0, 1.0, nodes);

    algorithm.nlp_method            = "IPOPT";
    algorithm.scaling               = "automatic";
    algorithm.derivatives           = "automatic";
    algorithm.nlp_iter_max          = 2000;
    algorithm.nlp_tolerance         = 1.0e-10;
    algorithm.print_level           = 0;
    algorithm.mesh_refinement       = "manual";
    algorithm.collocation_method    = "Hermite-Simpson";
    algorithm.transcription_method  = "multiple-shooting";
    algorithm.ms_steps_per_segment  = 10;
    algorithm.ms_control_parameterisation = "linear";

    ASSERT_EQ(psopt(solution, problem, algorithm), 0);

    MatrixXd L = solution.get_dual_costates_in_phase(1);
    MatrixXd T = solution.get_time_in_phase(1);
    ASSERT_EQ(L.cols(), T.cols());

    double e1 = 0.0, e2 = 0.0;
    for (int q = 0; q < T.cols(); q++) {
        e1 = std::max( e1, std::fabs( L(0,q) + 12.0 ) );
        e2 = std::max( e2, std::fabs( L(1,q) - (12.0*T(0,q) - 6.0) ) );
    }
    EXPECT_LT(e1, 1.0e-6) << "max |lambda_1 + 12| = " << e1;
    EXPECT_LT(e2, 1.0e-6) << "max |lambda_2 - (12t-6)| = " << e2;
}


// ---------------------------------------------------------------------------
// Flexible segment boundaries, and the thing they remove: a dependence of the
// answer on an arithmetic coincidence.
//
// Minimum time with u in [-1,2] switches at tf/3 and tf is free, so a UNIFORM
// partition has a boundary on the switch exactly when the segment count is
// divisible by three -- and the answer is then correct to eight digits, and
// wrong by parts in a thousand when it is not. Nothing about the problem
// changes between M = 9 and M = 10.
//
// With the boundaries free, every segment count gives the same answer, because
// the optimisation puts a boundary where the solution needs one rather than
// where the partition happened to put it.
// ---------------------------------------------------------------------------

TEST(MultipleShooting, AUniformPartitionIsRightOnlyByCoincidence)
{
    const msh::Run lucky   = msh::solve(4,  9, 10);   // 9 divisible by 3
    const msh::Run unlucky = msh::solve(4, 10, 10);   // 10 is not

    ASSERT_EQ(lucky.flag,   0);
    ASSERT_EQ(unlucky.flag, 0);

    const double e_lucky   = std::fabs(lucky.tf   - msh::TF_SQRT3)/msh::TF_SQRT3;
    const double e_unlucky = std::fabs(unlucky.tf - msh::TF_SQRT3)/msh::TF_SQRT3;

    EXPECT_LT(e_lucky,   1.0e-7) << "M = 9, tf = "  << lucky.tf;
    EXPECT_GT(e_unlucky, 1.0e-4) << "M = 10, tf = " << unlucky.tf
        << " -- if this is small the coincidence is not being exercised";
}


TEST(MultipleShooting, FlexibleSegmentsRemoveTheCoincidence)
{
    const int segs[4] = { 5, 7, 10, 20 };     // none divisible by three
    for (int q = 0; q < 4; q++) {
        const msh::Run fixed = msh::solve(4, segs[q], 10, "constant", 0, false);
        const msh::Run flex  = msh::solve(4, segs[q], 10, "constant", 0, true);

        ASSERT_EQ(fixed.flag, 0) << "M = " << segs[q];
        ASSERT_EQ(flex.flag,  0) << "the flexible partition failed at M = " << segs[q];

        const double e_fixed = std::fabs(fixed.tf - msh::TF_SQRT3)/msh::TF_SQRT3;
        const double e_flex  = std::fabs(flex.tf  - msh::TF_SQRT3)/msh::TF_SQRT3;

        EXPECT_GT(e_fixed, 1.0e-4) << "M = " << segs[q] << ": the uniform partition is "
            << "unexpectedly accurate, so this test is measuring nothing";
        EXPECT_LT(e_flex,  1.0e-7) << "M = " << segs[q] << ": flexible tf = " << flex.tf;
    }
}


// ---------------------------------------------------------------------------
// The quadratic control parameterisation.
//
// With the dynamics integrated to whatever ms_steps_per_segment buys, the
// accuracy of the answer is capped by the control parameterisation and by
// nothing else. The order of that approximation is therefore the order of the
// method, and the test that says so is an exactness test rather than a
// convergence table: on a problem whose optimal control IS a parabola, a
// parameterisation that carries parabolas must return the CONTINUOUS optimum at
// any number of segments, and one that carries ramps cannot.
//
// The residue left by the quadratic form is the INTEGRATOR's, not the
// parameterisation's, and the two are told apart by holding the segment count
// and refining the step: RK4 on a quadrature is Simpson's rule, exact through
// cubics, and here x1' = x2 is quartic in t. Refining the step drives the
// quadratic form's error to round-off and leaves the linear form's exactly
// where it was.
// ---------------------------------------------------------------------------

TEST(MultipleShooting, AQuadraticControlIsExactWhereTheOptimalControlIsQuadratic)
{
    const int segs[3] = { 3, 5, 10 };
    for (int q = 0; q < 3; q++) {
        const msh::Run r = msh::solve_triple(segs[q], 80, "quadratic");
        ASSERT_EQ(r.flag, 0) << "failed at " << segs[q] << " segments";
        EXPECT_NEAR(r.J, msh::TRIPLE_JSTAR, 1.0e-6)
            << segs[q] << " segments: got " << r.J << " for a control the "
            << "parameterisation can represent exactly";
    }
}


TEST(MultipleShooting, ThreeQuadraticSegmentsBeatTwentyLinearOnes)
{
    const msh::Run quad = msh::solve_triple( 3, 80, "quadratic");
    const msh::Run lin  = msh::solve_triple(20, 80, "linear");

    ASSERT_EQ(quad.flag, 0);
    ASSERT_EQ(lin.flag,  0);

    const double e_quad = std::fabs(quad.J - msh::TRIPLE_JSTAR)/msh::TRIPLE_JSTAR;
    const double e_lin  = std::fabs(lin.J  - msh::TRIPLE_JSTAR)/msh::TRIPLE_JSTAR;

    EXPECT_GT(e_lin, 1.0e-6) << "the ramp is unexpectedly accurate here (" << e_lin
        << "), so this test is measuring nothing";
    EXPECT_LT(e_quad, e_lin/100.0)
        << "quadratic at 3 segments: " << e_quad << ";  linear at 20: " << e_lin;
    EXPECT_LT(quad.nvars, lin.nvars)
        << "and it should be doing it with fewer variables: " << quad.nvars
        << " against " << lin.nvars;
}


TEST(MultipleShooting, RefiningTheStepRemovesTheQuadraticFormsErrorAndNotTheRampsError)
{
    const msh::Run q_coarse = msh::solve_triple(5,  10, "quadratic");
    const msh::Run q_fine   = msh::solve_triple(5, 320, "quadratic");
    const msh::Run l_coarse = msh::solve_triple(5,  10, "linear");
    const msh::Run l_fine   = msh::solve_triple(5, 320, "linear");

    ASSERT_EQ(q_coarse.flag, 0);  ASSERT_EQ(q_fine.flag, 0);
    ASSERT_EQ(l_coarse.flag, 0);  ASSERT_EQ(l_fine.flag, 0);

    const double eq_c = std::fabs(q_coarse.J - msh::TRIPLE_JSTAR);
    const double eq_f = std::fabs(q_fine.J   - msh::TRIPLE_JSTAR);
    const double el_c = std::fabs(l_coarse.J - msh::TRIPLE_JSTAR);
    const double el_f = std::fabs(l_fine.J   - msh::TRIPLE_JSTAR);

    EXPECT_LT(eq_f, eq_c/1000.0)
        << "the quadratic form's error did not fall with the step: "
        << eq_c << " -> " << eq_f << ", so it is not the integrator's error";
    EXPECT_GT(el_f, el_c/2.0)
        << "the ramp's error fell with the step: " << el_c << " -> " << el_f
        << ", so it is not the parameterisation's error";
}


// The layout written in three places. Adding a block to one of them and not the
// others does not fail a check; it corrupts the heap. The quadratic form adds
// exactly one control variable per SEGMENT, so the length of the decision
// vector has to move by exactly that and by nothing else.
TEST(MultipleShooting, TheQuadraticFormAddsExactlyOneVariablePerSegment)
{
    const int segs[3] = { 5, 10, 20 };
    for (int q = 0; q < 3; q++) {
        const msh::Run lin  = msh::solve(0, segs[q], 10, "linear");
        const msh::Run quad = msh::solve(0, segs[q], 10, "quadratic");
        ASSERT_EQ(lin.flag,  0);
        ASSERT_EQ(quad.flag, 0);
        EXPECT_EQ(quad.nvars - lin.nvars, segs[q])
            << "M = " << segs[q] << ": linear has " << lin.nvars
            << " variables and quadratic " << quad.nvars;
    }
}


// A caller who reads solution.get_controls_in_phase alone under this
// parameterisation is reading two thirds of the control variables and calling
// it the control history. The midpoint values are reported through the same
// interleaved pair Hermite-Simpson uses, and they are genuinely free: on a
// problem whose optimal control has curvature they differ from the average of
// the two nodes they sit between, which is what the ramp would have given.
TEST(MultipleShooting, TheMidpointControlsAreReportedAndAreNotTheNodeAverage)
{
    const int M = 6;
    const msh::Run r = msh::solve_triple(M, 40, "quadratic");
    ASSERT_EQ(r.flag, 0);
    ASSERT_EQ(r.u_hs.cols(), 2*M + 1)
        << "the interleaved control history has " << r.u_hs.cols()
        << " columns, not the 2M+1 a midpoint control implies";
    ASSERT_EQ(r.t_hs.cols(), 2*M + 1);

    double gap = 0.0, scale = 0.0;
    for (int k = 0; k < M; k++) {
        const double avg = 0.5*( r.u_hs(0,2*k) + r.u_hs(0,2*k+2) );
        gap   = std::max(gap,   std::fabs(r.u_hs(0,2*k+1) - avg));
        scale = std::max(scale, std::fabs(r.u_hs(0,2*k+1)));
    }
    EXPECT_GT(gap, 0.01*scale)
        << "every midpoint control is the average of its neighbours (gap " << gap
        << " against scale " << scale << "), which is what a ramp would give";

    // and the midpoints sit where they say they do
    for (int k = 0; k < M; k++)
        EXPECT_NEAR(r.t_hs(0,2*k+1), 0.5*(r.t_hs(0,2*k) + r.t_hs(0,2*k+2)), 1.0e-12);
}


// A limitation, pinned so that it cannot change quietly. The parabola through
// three values inside the control's box need not stay inside it: it overshoots
// by a quarter of the second difference. The sharpest second difference a
// solution can present is a jump -- and putting a segment boundary exactly on a
// jump is what ms_flexible_segments is FOR, so the two features collide
// precisely where each is doing its job. On the minimum-time problem with
// u in [-1,2] the reported control reaches about 2.375, and that is the control
// the segment integrator is handed, not an artefact of plotting.
TEST(MultipleShooting, TheParabolaCanLeaveTheControlBounds)
{
    const msh::Run r = msh::solve(4, 10, 10, "quadratic", 0, true);
    ASSERT_EQ(r.flag, 0);
    ASSERT_EQ(r.u_hs.cols(), 2*10 + 1);

    double worst = 0.0;
    for (int k = 0; k < 10; k++)
        for (int q = 0; q <= 100; q++) {
            const double u = msh::quad_control(r, k, ((double) q)/100.0);
            worst = std::max(worst, std::max(u - 2.0, -1.0 - u));
        }

    EXPECT_GT(worst, 0.1)
        << "the parabola stayed within the control bounds (worst excursion "
        << worst << ") -- if this holds, the caution in validate.cxx and in the "
        << "Alg comment is no longer describing the code";

    // and the values it interpolates are themselves inside the bounds, so this
    // is the interpolant leaving the box and not the solver ignoring it
    for (int c = 0; c < r.u_hs.cols(); c++) {
        EXPECT_LE(r.u_hs(0,c),  2.0 + 1.0e-7);
        EXPECT_GE(r.u_hs(0,c), -1.0 - 1.0e-7);
    }
}


// ---------------------------------------------------------------------------
// A choice of explicit scheme.
//
// algorithm.ms_integrator = "RK8" is Cooper and Verner's eleven-stage
// eighth-order formula. What it is for is not accuracy for its own sake but
// accuracy per RIGHT-HAND-SIDE EVALUATION: the tape a shooting transcription
// builds is stages x steps x segments long, and both the memory it occupies and
// the time to evaluate the constraints are linear in that. Four orders of
// convergence for 2.75 times the stages is a large trade in the right direction
// whenever the integrator error matters at all.
//
// The scheme's own error is measured here against the closed-form segment map of
// case 1 rather than against PSOPT's estimate of it, so that the scheme and the
// estimator are tested separately and neither can excuse the other.
// ---------------------------------------------------------------------------

// Each scheme has to be measured in the window where its error is above
// round-off and below saturation, and those windows do not overlap: on this
// problem RK8 is already at 2e-15 with four steps across a tenth of the horizon,
// where RK4 still has five useful decades. Measuring both on one configuration
// would mean measuring one of them on noise.
static void check_order(const char* scheme, int M, const int* steps, int nsteps,
                        double lo, double hi)
{
    double prev = 0.0;
    int checked = 0;
    for (int q = 0; q < nsteps; q++) {
        const msh::Run r = msh::solve(1, M, steps[q], "constant", 0, false, scheme);
        ASSERT_EQ(r.flag, 0) << scheme << " failed at " << steps[q] << " steps";
        const double e = msh::true_relative_local_error(r, M);
        if ( prev > 0.0 && e > 1.0e-12 ) {
            const double ratio = prev/e;
            EXPECT_GT(ratio, lo) << scheme << ", " << steps[q] << " steps: ratio " << ratio;
            EXPECT_LT(ratio, hi) << scheme << ", " << steps[q] << " steps: ratio " << ratio;
            checked++;
        }
        prev = e;
    }
    EXPECT_GT(checked, 0) << scheme << ": every ratio was at the round-off floor, so this "
        << "measured nothing -- the configuration, not the scheme, is what failed";
}

TEST(MultipleShooting, TheEighthOrderSchemeConvergesAtEighthOrder)
{
    const int s4[3] = { 2, 4, 8 };
    check_order("RK4", 10, s4, 3, 12.0, 20.0);          // 2^4 = 16

    // A mistyped tableau coefficient does not fail: it gives a scheme of lower
    // order that still converges to the right answer. This is where it shows.
    // Two steps across a third of the horizon already leaves RK8 within a factor
    // of a hundred of round-off, so the window here is one halving wide. One
    // ratio in [150, 400] separates order eight from every lower order.
    const int s8[2] = { 1, 2 };
    check_order("RK8", 3, s8, 2, 150.0, 400.0);         // 2^8 = 256
}


TEST(MultipleShooting, TheEighthOrderSchemeBuysMoreAccuracyPerEvaluation)
{
    const int M = 10;
    // 8 RK4 steps is 32 evaluations per segment; 2 RK8 steps is 22. The cheaper
    // one has to be the more accurate one, or the scheme is not worth having.
    const msh::Run r4 = msh::solve(1, M, 8, "constant", 0, false, "RK4");
    const msh::Run r8 = msh::solve(1, M, 2, "constant", 0, false, "RK8");
    ASSERT_EQ(r4.flag, 0);
    ASSERT_EQ(r8.flag, 0);

    const double e4 = msh::true_relative_local_error(r4, M);
    const double e8 = msh::true_relative_local_error(r8, M);

    EXPECT_LT(e8, e4/10.0)
        << "RK8 at 22 evaluations per segment gave " << e8
        << " against RK4 at 32 evaluations giving " << e4;
}


// The estimate has to be the error of the trajectory the user HAS, not of one
// that was computed to produce it. Writing C h^p for the leading error of a
// scheme of order p, the difference between the run at h and the run at h/2 is
// C h^p (2^-p - 1), so the solved run's error is that difference over
// (1 - 2^-p); dividing by (2^p - 1) instead gives the HALF-STEP run's error,
// which is smaller by exactly 2^p.
//
// That factor was in the code from the day the branch was written and was
// invisible, because a constant factor leaves every convergence ratio right. It
// is caught here by comparing against a truth computed outside PSOPT, and the
// test covers both schemes because the factor is 2^p and therefore differs
// between them -- 16 and 256.
TEST(MultipleShooting, TheErrorEstimateIsTheErrorOfTheRunThatWasSolved)
{
    struct Cfg { const char* scheme; int M; int steps; };
    const Cfg cfg[4] = { {"RK4", 10, 2}, {"RK4", 10, 4}, {"RK8", 3, 1}, {"RK8", 3, 2} };

    for (int q = 0; q < 4; q++) {
        const msh::Run r = msh::solve(1, cfg[q].M, cfg[q].steps, "constant", 0, false,
                                      cfg[q].scheme);
        ASSERT_EQ(r.flag, 0) << cfg[q].scheme << " failed at " << cfg[q].steps << " steps";
        const double truth = msh::true_relative_local_error(r, cfg[q].M);
        ASSERT_GT(truth, 1.0e-13) << cfg[q].scheme << " at " << cfg[q].steps
            << " steps is at the round-off floor, so this comparison measures nothing";
        const double ratio = r.err_est/truth;
        EXPECT_GT(ratio, 0.8) << cfg[q].scheme << " at " << cfg[q].steps << " steps: reported "
            << r.err_est << " against a true " << truth;
        EXPECT_LT(ratio, 1.25) << cfg[q].scheme << " at " << cfg[q].steps << " steps: reported "
            << r.err_est << " against a true " << truth;
    }
}


// And the caution that belongs with the option. Raising the order of the
// integrator improves the ANSWER only where the integrator was the binding
// error, and under a piecewise-constant control it never is: the control
// parameterisation caps the cost at O(h^2) and no scheme can lift that.
TEST(MultipleShooting, AHigherOrderSchemeCannotLiftTheControlParameterisationsCap)
{
    const msh::Run c4 = msh::solve(0, 20, 10, "constant", 0, false, "RK4");
    const msh::Run c8 = msh::solve(0, 20, 10, "constant", 0, false, "RK8");
    ASSERT_EQ(c4.flag, 0);
    ASSERT_EQ(c8.flag, 0);

    // The same discrete problem, so the same answer: the segment map of case 0
    // is a polynomial both schemes integrate exactly.
    EXPECT_NEAR(c4.J, c8.J, 1.0e-9)
        << "constant control: RK4 gave " << c4.J << " and RK8 " << c8.J;
    EXPECT_NEAR(c8.J, msh::discrete_optimum(20), 1.0e-8);
}


// ---------------------------------------------------------------------------
// Automatic segment refinement.
//
// The question a shooting mesh asks is not the one the other drivers answer.
// More segments does not mean a better approximation of the DYNAMICS here --
// those are integrated to whatever ms_steps_per_segment and ms_integrator buy,
// however many segments there are. What the segment count controls is the
// resolution of the CONTROL PARAMETERISATION and the coverage of the PATH
// CONSTRAINTS, and the indicator measures those.
//
// The counts below are path-dependent -- the refinement makes a discrete choice
// each iteration from an estimate computed on the previous solve -- so the
// thresholds are set from measurement with room, and the claims are about what
// the refinement achieves rather than about which mesh it took to get there.
// The element-refinement tests learned that the hard way.
// ---------------------------------------------------------------------------

TEST(MultipleShooting, RefinementMakesTheAnswerIndependentOfTheStartingMesh)
{
    // Minimum time with u in [-1,2]: the switch falls at tf/3, so a uniform
    // partition resolves it only when the segment count is divisible by three.
    // None of these three is, and every one of them is wrong by parts in a
    // thousand before refinement.
    const int segs[3] = { 7, 10, 13 };
    for (int q = 0; q < 3; q++) {
        const msh::Run fixed = msh::solve(4, segs[q], 10, "constant", 0, false, "RK4", false);
        const msh::Run autom = msh::solve(4, segs[q], 10, "constant", 0, false, "RK4", true);

        ASSERT_EQ(fixed.flag, 0) << "M = " << segs[q];
        ASSERT_EQ(autom.flag, 0) << "automatic refinement failed at M = " << segs[q];

        const double e_fixed = std::fabs(fixed.tf - msh::TF_SQRT3)/msh::TF_SQRT3;
        const double e_auto  = std::fabs(autom.tf - msh::TF_SQRT3)/msh::TF_SQRT3;

        EXPECT_GT(e_fixed, 1.0e-4) << "M = " << segs[q] << ": the uniform partition is "
            << "unexpectedly accurate, so this test is measuring nothing";
        EXPECT_LT(e_auto, 1.0e-6) << "M = " << segs[q] << " refined to " << autom.segments
            << " segments and gave tf = " << autom.tf;
        EXPECT_GT(autom.segments, segs[q]) << "the mesh did not grow";
    }
}


// And the other half of the same claim: a mesh that ALREADY resolves the corner
// is left alone. This is the patch-177 trap in a new place -- an estimator built
// for smooth solutions, pointed at a discontinuity the mesh has already
// resolved, flags it every iteration however fine the mesh becomes. The
// indicator avoids it by forming its departure from a window extended to the
// left and one extended to the right and taking the SMALLER: a corner sitting on
// a boundary spoils exactly one of the two.
TEST(MultipleShooting, AMeshThatAlreadyResolvesTheCornerIsNotShattered)
{
    // 12 is divisible by three, so a uniform boundary sits on the switch.
    const msh::Run resolved = msh::solve(4, 12, 10, "constant", 0, false, "RK4", true);
    ASSERT_EQ(resolved.flag, 0);

    EXPECT_LT(std::fabs(resolved.tf - msh::TF_SQRT3)/msh::TF_SQRT3, 1.0e-6);
    EXPECT_LE(resolved.segments, 24)
        << "the mesh grew from 12 to " << resolved.segments
        << " on a partition that already had a boundary on the switch -- the "
        << "indicator is flagging a corner it should be ignoring";

    // whereas one segment more, and the switch is inside a segment again
    const msh::Run inside = msh::solve(4, 13, 10, "constant", 0, false, "RK4", true);
    ASSERT_EQ(inside.flag, 0);
    EXPECT_GT(inside.segments, resolved.segments)
        << "a mesh with the switch inside a segment (" << inside.segments
        << ") was refined no more than one with the switch on a boundary ("
        << resolved.segments << ")";
}


TEST(MultipleShooting, RefinementImprovesASmoothProblem)
{
    // The oscillator, whose optimal control is smooth and reachable by none of
    // the three parameterisations, so refinement is answering the question it
    // was designed for and nothing else is confusing the measurement.
    const msh::Run coarse = msh::solve(1, 5, 10, "constant", 0, false, "RK4", false);
    const msh::Run autom  = msh::solve(1, 5, 10, "constant", 0, false, "RK4", true);

    ASSERT_EQ(coarse.flag, 0);
    ASSERT_EQ(autom.flag,  0) << "automatic refinement failed on a smooth problem";

    EXPECT_GT(autom.segments, 5) << "the mesh did not grow";
    EXPECT_LT(autom.J, coarse.J)
        << "refined " << autom.J << " against coarse " << coarse.J
        << " -- this is a minimisation, so more segments cannot cost more";

    const double Jstar = msh::oscillator_optimum();
    const double e_coarse = std::fabs(coarse.J - Jstar);
    const double e_auto   = std::fabs(autom.J  - Jstar);
    EXPECT_LT(e_auto, 0.1*e_coarse)
        << "the refinement was worth less than a factor of ten: coarse "
        << coarse.J << " (error " << e_coarse << "), refined " << autom.J
        << " (error " << e_auto << "), against J* = " << Jstar;
}


// The path constraints are the other thing the segment count controls, and they
// are the half that produces a WRONG answer rather than an inaccurate one: the
// state between two boundaries is not a decision variable, so a constraint
// imposed only at boundaries leaks, and the optimiser returns a cost BELOW the
// true optimum. Nothing in the NLP reports that, because every constraint the
// NLP was given is satisfied; the indicator has to go and look.
TEST(MultipleShooting, RefinementShrinksAPathConstraintLeak)
{
    const msh::Run fixed = msh::solve(3, 10, 10, "linear", 0, false, "RK4", false);
    const msh::Run autom = msh::solve(3, 10, 10, "linear", 0, false, "RK4", true);

    ASSERT_EQ(fixed.flag, 0);
    ASSERT_EQ(autom.flag, 0) << "automatic refinement failed with a path constraint";

    // J* = 4 and no feasible trajectory costs less, so 4 - J is the leak.
    const double leak_fixed = 4.0 - fixed.J;
    const double leak_auto  = std::fabs(4.0 - autom.J);

    EXPECT_GT(leak_fixed, 1.0e-3) << "the coarse mesh is not leaking (J = " << fixed.J
        << "), so this test is measuring nothing";
    EXPECT_LT(leak_auto, leak_fixed/20.0)
        << "the leak went from " << leak_fixed << " to " << leak_auto
        << " over " << fixed.segments << " -> " << autom.segments << " segments";
}


// ---------------------------------------------------------------------------
// Equality path constraints.
//
// Bryson's maximum-range problem: max x(tf) with xdot = v u1, ydot = v u2,
// vdot = a - g u2 and u1^2 + u2^2 = 1, which is an equality in the CONTROLS
// ALONE. That distinction decides everything about how it behaves here.
//
// Between two segment boundaries the state is not a decision variable, so an
// equality involving the state cannot be made to hold there by anything but the
// dynamics. An equality in the controls alone is a different matter: under a
// piecewise-CONSTANT parameterisation a control that satisfies it at the start
// of a segment satisfies it everywhere on that segment, exactly, and the
// transcription represents the constraint with no error at all.
//
// And such a component must NOT be sampled inside the segments. Sampling an
// equality adds an equation and no unknown; IPOPT says so in as many words --
// Not_Enough_Degrees_Of_Freedom, return code -10 -- and PSOPT then returned the
// initial guess with a success flag. Every run of this problem with
// ms_path_samples > 0 failed that way before ms_samplable_path_indices.
// ---------------------------------------------------------------------------

namespace mspath {

adouble endpoint_cost(adouble*, adouble* f, adouble*, adouble&, adouble&, adouble*,
                      int, Workspace*) { return -f[0]; }
adouble integrand_cost(adouble*, adouble*, adouble*, adouble&, adouble*, int, Workspace*)
{ return (adouble) 0.0; }

void dae(adouble* d, adouble* path, adouble* s, adouble* c, adouble*, adouble&,
         adouble*, int, Workspace*)
{
    const double g = 1.0, a = 0.5*g;
    d[0] = s[2]*c[0];
    d[1] = s[2]*c[1];
    d[2] = a - g*c[1];
    path[0] = c[0]*c[0] + c[1]*c[1];
}

void events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&,
            adouble*, int, Workspace*)
{ e[0]=i[0]; e[1]=i[1]; e[2]=i[2]; e[3]=f[1]; }

struct Run { int flag; int rc; double J; MatrixXd u; int M; };

static Run solve(int segments, const std::string& upar, int path_samples,
                 int diagnostic_level = 0)
{
    Alg algorithm; Sol solution; Prob problem;
    Run out; out.flag = -1; out.rc = -99; out.J = 0.0; out.M = 0;

    const int nodes = segments + 1;
    problem.name        = "multiple shooting, equality path";
    problem.outfilename = "test_multiple_shooting_eq.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 3;
    problem.phases(1).ncontrols = 2;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 1;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -10.0, -10.0, -10.0;
    problem.phases(1).bounds.upper.states   <<  10.0,  10.0,  10.0;
    problem.phases(1).bounds.lower.controls << -10.0, -10.0;
    problem.phases(1).bounds.upper.controls <<  10.0,  10.0;
    problem.phases(1).bounds.lower.path(0)  = 1.0;      // an EQUALITY
    problem.phases(1).bounds.upper.path(0)  = 1.0;
    problem.phases(1).bounds.lower.events   << 0.0, 0.0, 0.0, 0.0;
    problem.phases(1).bounds.upper.events   << 0.0, 0.0, 0.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 2.0; problem.phases(1).bounds.upper.EndTime   = 2.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &msh::linkages;

    problem.phases(1).guess.states   = zeros(3, nodes);
    problem.phases(1).guess.controls = zeros(2, nodes);
    problem.phases(1).guess.controls.row(0) = ones(1, nodes);
    problem.phases(1).guess.time     = linspace(0.0, 2.0, nodes);

    algorithm.nlp_method            = "IPOPT";
    algorithm.scaling               = "automatic";
    algorithm.derivatives           = "automatic";
    algorithm.nlp_iter_max          = 2000;
    algorithm.nlp_tolerance         = 1.0e-9;
    algorithm.print_level           = 0;
    algorithm.mesh_refinement       = "manual";
    algorithm.collocation_method    = "Hermite-Simpson";
    algorithm.transcription_method  = "multiple-shooting";
    algorithm.ms_steps_per_segment  = 20;
    algorithm.ms_control_parameterisation = upar;
    algorithm.ms_path_samples             = path_samples;
    algorithm.diagnostic_level            = diagnostic_level;

    out.flag = psopt(solution, problem, algorithm);
    out.rc   = solution.nlp_return_code;
    if (out.flag == 0) {
        out.J = solution.cost;
        out.u = solution.get_controls_in_phase(1);
        out.M = (int) out.u.cols() - 1;
    }
    return out;
}

// The worst |u|^2 - 1 anywhere along a segment, using the control representation
// the transcription actually integrated -- which is the only place the question
// can be asked, the nodal values satisfying it by construction.
static double worst_violation(const Run& r, const std::string& upar)
{
    double worst = 0.0;
    for (int k = 0; k < r.M; k++)
        for (int q = 0; q <= 40; q++) {
            const double s = ((double) q)/40.0;
            double u1, u2;
            if ( upar == "constant" ) { u1 = r.u(0,k); u2 = r.u(1,k); }
            else { u1 = (1-s)*r.u(0,k) + s*r.u(0,k+1);
                   u2 = (1-s)*r.u(1,k) + s*r.u(1,k+1); }
            worst = std::max( worst, std::fabs(u1*u1 + u2*u2 - 1.0) );
        }
    return worst;
}

} // namespace mspath


TEST(MultipleShooting, AControlOnlyEqualityPathConstraintIsExactUnderAHeldControl)
{
    const mspath::Run r = mspath::solve(10, "constant", 0);
    ASSERT_EQ(r.flag, 0) << "IPOPT return code " << r.rc;

    EXPECT_LT(mspath::worst_violation(r, "constant"), 1.0e-12)
        << "a constant control satisfying the equality at the start of a segment "
        << "satisfies it everywhere on that segment, so this should be round-off";
    EXPECT_NEAR(r.J, -1.7944638, 0.02)
        << "against the Hermite-Simpson collocation answer on fifty nodes";
}


TEST(MultipleShooting, AContinuousControlLeavesAControlOnlyEqualityBetweenTheNodes)
{
    // The chord between two points on the unit circle lies inside it, so the
    // linear form satisfies the equality at the nodes and nowhere else. This is
    // a limitation rather than a defect, and it is pinned so that it cannot stop
    // being true quietly.
    const mspath::Run coarse = mspath::solve(10, "linear", 0);
    const mspath::Run fine   = mspath::solve(40, "linear", 0);
    ASSERT_EQ(coarse.flag, 0) << "IPOPT return code " << coarse.rc;
    ASSERT_EQ(fine.flag,   0) << "IPOPT return code " << fine.rc;

    const double v_coarse = mspath::worst_violation(coarse, "linear");
    const double v_fine   = mspath::worst_violation(fine,   "linear");

    EXPECT_GT(v_coarse, 1.0e-3) << "the ramp is unexpectedly satisfying the equality "
        << "between the nodes (" << v_coarse << "), so this test measures nothing";
    EXPECT_LT(v_fine, v_coarse/4.0)
        << "four times the segments should buy about sixteen: " << v_coarse
        << " -> " << v_fine;
}


// The regression this all turns on: asking for interior samples must not make
// the problem infeasible. Before ms_samplable_path_indices every one of these
// returned the initial guess with a success flag and an IPOPT return code of
// -10, Not_Enough_Degrees_Of_Freedom.
TEST(MultipleShooting, SamplingDoesNotOverDetermineAnEqualityPathConstraint)
{
    const mspath::Run none = mspath::solve(10, "constant", 0);
    ASSERT_EQ(none.flag, 0) << "IPOPT return code " << none.rc;

    for (int samples : { 1, 3 }) {
        const mspath::Run r = mspath::solve(10, "constant", samples);
        ASSERT_EQ(r.flag, 0) << samples << " interior samples: IPOPT return code " << r.rc
            << " (-10 is Not_Enough_Degrees_Of_Freedom, which is what sampling an "
            << "equality used to produce)";
        EXPECT_NE(r.rc, -10);
        EXPECT_NEAR(r.J, none.J, 1.0e-6)
            << "an equality component is not sampled, so asking for samples must not "
            << "change the answer: " << none.J << " against " << r.J;
    }
}


// ---------------------------------------------------------------------------
// The solution diagnostics run with this transcription.
//
// They were refused while it had no costates. It has them now, and the part
// that matters most -- the rank and conditioning of the constraint Jacobian --
// never depended on the transcription at all: it re-tapes the constraints at
// the final iterate and factorises them, which is the same question whatever
// wrote the rows.
//
// What the report says on the problem above is the point of enabling it. With
// the equality depending on the controls ALONE the terminal control slot is
// pinned to its neighbour, so the path row at the final node duplicates the row
// before it and the report reads "RANK DEFICIENT BY 5 (1 beyond the empty
// rows)". With the same constraint made state-dependent it reads "deficient by
// 4, which is exactly the empty rows: the constraints that are actually written
// are independent". The four are structural -- three unused matching rows from
// the shared layout and the duration row, whose two times are both fixed here --
// and counting them separately is what makes the check usable at all, since
// otherwise it can never read as full rank under this transcription.
//
// The report prints rather than returns, so what is asserted here is that
// asking for it is allowed and changes nothing about the answer. The numbers
// themselves are demonstrated in the scratch record and the manual.
// ---------------------------------------------------------------------------

TEST(MultipleShooting, TheSolutionDiagnosticsRunAndChangeNothing)
{
    const mspath::Run quiet = mspath::solve(10, "constant", 0, 0);
    const mspath::Run loud  = mspath::solve(10, "constant", 0, 2);

    ASSERT_EQ(quiet.flag, 0) << "IPOPT return code " << quiet.rc;
    ASSERT_EQ(loud.flag,  0) << "asking for diagnostics made the solve fail: IPOPT return code "
                             << loud.rc;
    EXPECT_EQ(quiet.M, loud.M);
    EXPECT_NEAR(quiet.J, loud.J, 1.0e-12)
        << "a diagnostic that changes the answer is not a diagnostic: "
        << quiet.J << " against " << loud.J;
}
