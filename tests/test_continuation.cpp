//////////////////////////////////////////////////////////////////////////////
// test_continuation.cpp
//
// set_guess_from_solution() copies a solution into the initial guess of a
// problem, so that a sequence of related solves -- a parameter scan, a profile
// likelihood, a homotopy in a bound or a penalty weight -- behaves as one
// calculation rather than as many independent ones.
//
// The failure it prevents is not a failure to converge, which is why it is
// worth a test of its own. Started from a common cold guess, most points of a
// scan reach the same solution branch and a few do not: they converge, report
// success, and return a different local minimum. Nothing in the output
// distinguishes them, and on a plot of the scan they are indistinguishable from
// real structure in whatever is being scanned. Two of seventeen points did this
// while the profile likelihood of examples/cracking was being built, on a
// profile that is otherwise flat to two parts in 100,000.
//
// The test problem is chosen so that the whole scan is known in closed form:
//
//     minimise  int_0^1 (1/2) u^2 dt   subject to   xdot = u/p,
//     x(0) = 0,  x(1) = 1,  and p pinned by an event constraint.
//
// The dynamics integrate to 1 = (1/p) int u dt, so the minimum-energy control
// is the constant u = p and J*(p) = p^2/2 for every pinned value of p. A scan
// over p therefore has an exactly known answer at every point, which is what
// lets the test assert that continuation has not moved the answer.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>
#include <cmath>

namespace continuation_test {

adouble endpoint_cost(adouble* i, adouble* f, adouble* p, adouble& t0,
                      adouble& tf, adouble* xad, int iphase, Workspace* w)
{ return 0.0; }

adouble integrand_cost(adouble* states, adouble* controls, adouble* p,
                       adouble& t, adouble* xad, int iphase, Workspace* w)
{ return 0.5*controls[0]*controls[0]; }

void dae(adouble* d, adouble* path, adouble* states, adouble* controls,
         adouble* p, adouble& time, adouble* xad, int iphase, Workspace* w)
{ d[0] = controls[0]/p[0]; }

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{ e[0] = i[0]; e[1] = f[0]; e[2] = p[0]; }

void linkages(adouble* linkages, adouble* xad, Workspace* w) { }

// Everything that does not change from one solve to the next. The scan below
// then alters exactly two numbers per point: the bounds of the event that pins
// the parameter.
void setup(Prob& problem, Alg& algorithm)
{
    problem.name        = "Continuation test";
    problem.outfilename = "test_continuation.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates     = 1;
    problem.phases(1).ncontrols   = 1;
    problem.phases(1).nevents     = 3;
    problem.phases(1).npath       = 0;
    problem.phases(1).nparameters = 1;
    problem.phases(1).nodes       << 24;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states     << -10.0;
    problem.phases(1).bounds.upper.states     <<  10.0;
    problem.phases(1).bounds.lower.controls   << -20.0;
    problem.phases(1).bounds.upper.controls   <<  20.0;
    problem.phases(1).bounds.lower.parameters <<   0.5;
    problem.phases(1).bounds.upper.parameters <<  10.0;
    problem.phases(1).bounds.lower.events     << 0.0, 1.0, 2.5;
    problem.phases(1).bounds.upper.events     << 0.0, 1.0, 2.5;
    problem.phases(1).bounds.lower.StartTime  = 0.0;
    problem.phases(1).bounds.upper.StartTime  = 0.0;
    problem.phases(1).bounds.lower.EndTime    = 1.0;
    problem.phases(1).bounds.upper.EndTime    = 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states     = zeros(1,10);
    problem.phases(1).guess.controls   = ones(1,10);
    problem.phases(1).guess.time       = linspace(0.0, 1.0, 10);
    problem.phases(1).guess.parameters = ones(1,1);

    algorithm.nlp_method        = "IPOPT";
    algorithm.scaling           = "automatic";
    algorithm.derivatives       = "automatic";
    algorithm.collocation_method = "Hermite-Simpson";
    algorithm.nlp_iter_max      = 500;
    algorithm.nlp_tolerance     = 1.0e-8;
    algorithm.print_level       = 0;
}

void pin_parameter(Prob& problem, double v)
{
    problem.phases(1).bounds.lower.events(2) = v;
    problem.phases(1).bounds.upper.events(2) = v;
}

} // namespace continuation_test


TEST(Continuation, TheGuessIsTheSolution)
{
    using namespace continuation_test;

    Alg algorithm; Prob problem; Sol solution;
    setup(problem, algorithm);

    ASSERT_EQ(psopt(solution, problem, algorithm), 0);
    ASSERT_EQ(solution.error_flag, 0);
    EXPECT_NEAR(solution.cost, 0.5*2.5*2.5, 1.0e-6);

    // Before: the guess is the ten-point cold one supplied in setup().
    ASSERT_EQ(problem.phases(1).guess.states.cols(), 10);

    set_guess_from_solution(problem, solution);

    const MatrixXd& x = solution.get_states_in_phase(1);
    const MatrixXd& u = solution.get_controls_in_phase(1);
    const MatrixXd& t = solution.get_time_in_phase(1);
    const MatrixXd& p = solution.get_parameters_in_phase(1);

    ASSERT_EQ(problem.phases(1).guess.states.cols(),   x.cols());
    ASSERT_EQ(problem.phases(1).guess.controls.cols(), u.cols());
    ASSERT_EQ(problem.phases(1).guess.time.cols(),     t.cols());
    ASSERT_EQ(problem.phases(1).guess.parameters.rows(), p.rows());

    for (int j = 0; j < x.cols(); j++) {
        EXPECT_DOUBLE_EQ(problem.phases(1).guess.states(0,j),   x(0,j));
        EXPECT_DOUBLE_EQ(problem.phases(1).guess.controls(0,j), u(0,j));
        EXPECT_DOUBLE_EQ(problem.phases(1).guess.time(0,j),     t(0,j));
    }
    EXPECT_DOUBLE_EQ(problem.phases(1).guess.parameters(0,0), p(0,0));

    // The guess so produced must be usable: re-solving from it returns the
    // same answer rather than tripping over its own output.
    ASSERT_EQ(psopt(solution, problem, algorithm), 0);
    EXPECT_NEAR(solution.cost, 0.5*2.5*2.5, 1.0e-6);
}


TEST(Continuation, AScanTracedByContinuationKeepsTheClosedForm)
{
    using namespace continuation_test;

    Alg algorithm; Prob problem; Sol solution;
    setup(problem, algorithm);

    // Trace outwards from the value the problem is set up with, as a profile
    // likelihood would, each solve starting from the one before it.
    const double pv[] = { 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };
    const int    n    = (int)(sizeof pv / sizeof pv[0]);

    double J[n];
    for (int k = 0; k < n; k++) {
        pin_parameter(problem, pv[k]);
        ASSERT_EQ(psopt(solution, problem, algorithm), 0) << "at p = " << pv[k];
        ASSERT_EQ(solution.error_flag, 0)                 << "at p = " << pv[k];
        J[k] = solution.cost;
        set_guess_from_solution(problem, solution);

        // J*(p) = p^2/2 exactly, so every point of the scan can be checked
        // against the closed form rather than against its neighbours.
        EXPECT_NEAR(J[k], 0.5*pv[k]*pv[k], 1.0e-5) << "at p = " << pv[k];

        // The parameter really is pinned where it was asked to be.
        MatrixXd pp = solution.get_parameters_in_phase(1);
        EXPECT_NEAR(pp(0,0), pv[k], 1.0e-6);
    }

    // And the scan is monotone, which is the shape a plot of it would show.
    for (int k = 1; k < n; k++) EXPECT_GT(J[k], J[k-1]);
}


TEST(Continuation, AnEmptySolutionLeavesTheGuessAlone)
{
    using namespace continuation_test;

    // A Sol that has never been through psopt() carries nothing usable. The
    // helper must leave the user's guess in place rather than replacing it
    // with an empty matrix, which would fail later and far from the cause.
    Alg algorithm; Prob problem; Sol fresh;
    setup(problem, algorithm);

    MatrixXd states_before   = problem.phases(1).guess.states;
    MatrixXd controls_before = problem.phases(1).guess.controls;

    set_guess_from_solution(problem, fresh);

    ASSERT_EQ(problem.phases(1).guess.states.cols(),   states_before.cols());
    ASSERT_EQ(problem.phases(1).guess.controls.cols(), controls_before.cols());
    for (int j = 0; j < states_before.cols(); j++) {
        EXPECT_DOUBLE_EQ(problem.phases(1).guess.states(0,j),   states_before(0,j));
        EXPECT_DOUBLE_EQ(problem.phases(1).guess.controls(0,j), controls_before(0,j));
    }
}
