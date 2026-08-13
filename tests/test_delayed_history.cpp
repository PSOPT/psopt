//////////////////////////////////////////////////////////////////////////////
// test_delayed_history.cpp
//
// The initial history of a delay differential equation. For t in [t0-d, t0) the
// right hand side asks for the trajectory before the phase begins, so a history
//
//     x(t) = phi(t),   u(t) = psi(t),   t < t0
//
// is data of the problem in the same sense as the initial condition is. PSOPT
// uses the constant history phi(t) = x(t0), psi(t) = u(t0) unless
// problem.initial_history is set, and both routes have closed-form answers here.
//
//     x1' = -x1(t-1),   x1(0) = 1
//     x2' =  u(t-1),    x2(0) = 0,   u == 1 pinned by its bounds,   t in [0,3]
//
// By the method of steps, one delay interval at a time:
//
//   phi(t) = 1+t, psi(t) = 1-t   ->   x1(3) = -3/8,   x2(3) = 7/2
//   default constant history     ->   x1(3) = -1/6,   x2(3) = 3
//
// Both histories meet the trajectory at t = 0, so neither introduces a jump; a
// history that did would make the right hand side discontinuous at t0+d, which no
// collocation scheme resolves to better than the interval containing it.
//
// 31 nodes place a node at t = 1 and t = 2, the points at which the delayed
// arguments cross into the history and where the solution carries its kinks; a kink
// inside an interval rather than on a node costs an order. The exact solution is
// piecewise polynomial of low degree, so this coarse mesh is enough to identify the
// terminal values to eight figures, and the test stays cheap.
//
// The same problem is distributed as examples/delay_history.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>
#include <cmath>

namespace delayed_history {

adouble endpoint_cost(adouble* i, adouble* f, adouble* p, adouble& t0,
                      adouble& tf, adouble* xad, int iphase, Workspace* w)
{ return 0.0; }

adouble integrand_cost(adouble* s, adouble* c, adouble* p, adouble& t,
                       adouble* xad, int iphase, Workspace* w)
{ return 0.0; }

void dae(adouble* d, adouble* path, adouble* states, adouble* controls,
         adouble* p, adouble& time, adouble* xad, int iphase, Workspace* w)
{
    adouble x1_delayed, u_delayed;
    const double delay = 1.0;

    // Both indices count from zero, as every other index of the interface does.
    get_delayed_state  ( &x1_delayed, 0, iphase, time, delay, xad, w);
    get_delayed_control( &u_delayed,  0, iphase, time, delay, xad, w);

    d[0] = -x1_delayed;
    d[1] =  u_delayed;
}

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{ e[0] = i[0]; e[1] = i[1]; }

void linkages(adouble* linkages, adouble* xad, Workspace* w) { }

void history(adouble* states, adouble* controls, adouble* parameters,
             adouble& time, adouble* xad, int iphase, Workspace* w)
{
    states[0]   = 1.0 + time;    // phi(t) = 1 + t
    states[1]   = 0.0;           // never evaluated at a delayed time
    controls[0] = 1.0 - time;    // psi(t) = 1 - t
}

// Solve once, with or without the user history, and return the terminal states.
static void solve_it(bool use_history, double& x1f, double& x2f)
{
    Alg algorithm; Sol solution; Prob problem;
    const int N = 31;

    problem.name        = "Delayed initial history";
    problem.outfilename = "test_delayed_history.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 2;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << N;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -10.0, -10.0;
    problem.phases(1).bounds.upper.states   <<  10.0,  10.0;
    problem.phases(1).bounds.lower.controls << 1.0;      // u is pinned
    problem.phases(1).bounds.upper.controls << 1.0;
    problem.phases(1).bounds.lower.events   << 1.0, 0.0;
    problem.phases(1).bounds.upper.events   << 1.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 3.0;
    problem.phases(1).bounds.upper.EndTime   = 3.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    if (use_history) problem.initial_history = &history;

    MatrixXd x_guess = zeros(2,N);
    x_guess.row(0)   = ones(1,N);
    problem.phases(1).guess.states   = x_guess;
    problem.phases(1).guess.controls = ones(1,N);
    problem.phases(1).guess.time     = linspace(0.0,3.0,N);

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.nlp_iter_max       = 1000;
    algorithm.nlp_tolerance      = 1.e-8;
    algorithm.collocation_method = "Hermite-Simpson";
    algorithm.print_level        = 0;

    psopt(solution, problem, algorithm);

    MatrixXd x = solution.get_states_in_phase(1);
    x1f = x(0, x.cols()-1);
    x2f = x(1, x.cols()-1);
}

} // namespace delayed_history

// With problem.initial_history unset, PSOPT clamps the delayed time to t0, which is
// the constant history phi(t) = x1(0) = 1, psi(t) = u(0) = 1.
TEST(DelayedVariables, DefaultConstantHistory)
{
    double x1f, x2f;
    delayed_history::solve_it(false, x1f, x2f);

    EXPECT_NEAR(x1f, -1.0/6.0, 1.0e-5);
    EXPECT_NEAR(x2f,  3.0,     1.0e-5);
}

// With a user history the answers change, and by the amount the method of steps
// predicts: this is the check that the history is really being used, and used at
// the right times.
TEST(DelayedVariables, UserSuppliedHistory)
{
    double x1f, x2f;
    delayed_history::solve_it(true, x1f, x2f);

    EXPECT_NEAR(x1f, -0.375, 1.0e-5);
    EXPECT_NEAR(x2f,  3.5,   1.0e-5);

    // The two histories must give different answers, or the test above would pass
    // on an implementation that ignored problem.initial_history altogether.
    double x1_default, x2_default;
    delayed_history::solve_it(false, x1_default, x2_default);
    EXPECT_GT(std::fabs(x1f - x1_default), 0.2);
    EXPECT_GT(std::fabs(x2f - x2_default), 0.4);
}
