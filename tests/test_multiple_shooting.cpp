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
static int g_case = 0;

adouble endpoint_cost(adouble*, adouble*, adouble*, adouble&, adouble& tf, adouble*,
                      int, Workspace*)
{ return ( g_case == 2 ) ? tf : (adouble) 0.0; }

adouble integrand_cost(adouble*, adouble* u, adouble*, adouble&, adouble*, int, Workspace*)
{ return ( g_case == 2 ) ? (adouble) 0.0 : 0.5*u[0]*u[0]; }

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

struct Run { int flag; double J; double tf; double err_est; };

static Run solve(int which, int segments, int steps,
                 const std::string& upar = "constant", int path_samples = 0)
{
    g_case = which;

    Alg algorithm; Sol solution; Prob problem;
    Run out; out.flag = -1; out.J = 0.0; out.tf = 0.0; out.err_est = 0.0;

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
    problem.phases(1).bounds.lower.controls(0) = ( which == 2 ) ?  -1.0 : -30.0;
    problem.phases(1).bounds.upper.controls(0) = ( which == 2 ) ?   1.0 :  30.0;
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
    problem.phases(1).bounds.lower.EndTime   = ( which == 2 ) ? 0.5 : 1.0;
    problem.phases(1).bounds.upper.EndTime   = ( which == 2 ) ? 8.0 : 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(2, nodes);
    if ( which == 3 ) problem.phases(1).guess.states.row(1) = linspace( 1.0, -1.0, nodes);
    else              problem.phases(1).guess.states.row(0) = linspace( 0.0,  1.0, nodes);
    problem.phases(1).guess.controls = zeros(1, nodes);
    problem.phases(1).guess.time     = linspace(0.0, ( which == 2 ) ? 2.0 : 1.0, nodes);

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
    algorithm.ms_control_parameterisation = upar;
    algorithm.ms_path_samples             = path_samples;

    out.flag = psopt(solution, problem, algorithm);
    if (out.flag == 0) {
        out.J = solution.cost;
        MatrixXd t = solution.get_time_in_phase(1);
        out.tf = t(0, (int) t.cols() - 1);
        MatrixXd e = solution.get_relative_local_error_in_phase(1);
        for (int q = 0; q < e.size(); q++) out.err_est = std::max(out.err_est, std::fabs(e(q)));
    }
    return out;
}

static double discrete_optimum(double M) { return 6.0*M*M/(M*M - 1.0); }

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
