//////////////////////////////////////////////////////////////////////////
//////////////////      delay_history.cxx           //////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:  User-supplied initial history for a delay problem  //////
//////// Last modified:  13 August 2026                            ///////
//////// Reference:      method of steps, any text on delay        ///////
////////                 differential equations                    ///////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2026        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

// A delay differential equation is not well posed without an initial history: for
// t in [t0-d, t0) its right hand side asks for the trajectory at times before the
// phase begins. PSOPT assumes the constant history phi(t) = x(t0), psi(t) = u(t0)
// unless problem.initial_history is set. This example exists to check that both
// routes are right, because both have closed-form answers here.
//
//     x1' = -x1(t-1),          x1(0) = 1
//     x2' =  u(t-1),           x2(0) = 0,     u(t) = 1 pinned by its bounds
//
// on t in [0,3], with the history
//
//     phi(t) = 1 + t,     psi(t) = 1 - t,     t < 0.
//
// Both meet the trajectory at t = 0, so neither introduces a jump: phi(0-) = x1(0) = 1
// and psi(0-) = u(0) = 1. A history that does not meet it is admissible, and PSOPT will
// use it, but the resulting discontinuity in the right hand side cannot be resolved by
// any collocation scheme to better than one interval, so it makes a poor test.
//
// By the method of steps, integrating one delay interval at a time:
//
//     x1(3) = -3/8 = -0.375              with this history
//     x1(3) = -1/6 = -0.166667           with the constant history phi(t) = x1(0) = 1
//     x2(3) =  7/2 =  3.5                with this history
//     x2(3) =  3                         with the constant history psi(t) = u(0) = 1
//
// The two histories therefore give plainly different answers, and the example runs
// the problem twice, once each way, and reports both against the closed form.
//
// There is nothing to optimize here: the control is pinned and the trajectory is
// determined by the dynamics. That is deliberate, so that any discrepancy is the
// treatment of the history and nothing else.

#include "psopt.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////
///////////////////  Initial history function  ///////////////////////////
//////////////////////////////////////////////////////////////////////////

// Called only with time < t0. It must fill every state and every control, using
// adouble arithmetic, exactly as the other user functions do.

void initial_history(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    states[0]   = 1.0 + time;    // phi(t) = 1 + t
    states[1]   = 0.0;           // x2 is never evaluated at a delayed time
    controls[0] = 1.0 - time;    // psi(t) = 1 - t
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return 0.0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{
    return 0.0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
    adouble x1_delayed, u_delayed;
    double  d = 1.0;

    // The indices count from zero, like every other index of the interface.
    get_delayed_state  ( &x1_delayed, 0, iphase, time, d, xad, workspace);
    get_delayed_control( &u_delayed,  0, iphase, time, d, xad, workspace);

    derivatives[0] = -x1_delayed;
    derivatives[1] =  u_delayed;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    e[0] = initial_states[0];
    e[1] = initial_states[1];
}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
    // single phase problem
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Solve the problem once  ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

static void run_case(bool use_history, double& x1f, double& x2f)
{
    Alg  algorithm;
    Sol  solution;
    Prob problem;

    problem.name        = use_history ? "Delay problem with a user history"
                                      : "Delay problem with the default history";
    problem.outfilename = use_history ? "delay_history_user.txt"
                                      : "delay_history_default.txt";

    problem.nphases   = 1;
    problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 2;
    problem.phases(1).npath     = 0;
    // 121 nodes put a node exactly at t = 1 and t = 2, the points at which the delayed
    // arguments cross into the history and where the solution of a delay equation carries
    // its kinks. A kink inside an interval, rather than on a node, costs an order.
    problem.phases(1).nodes     << 121;

    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -10.0, -10.0;
    problem.phases(1).bounds.upper.states   <<  10.0,  10.0;

    // The control is pinned, so the trajectory is determined by the dynamics alone.
    problem.phases(1).bounds.lower.controls << 1.0;
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

    if ( use_history ) problem.initial_history = &initial_history;

    int nnodes = 121;
    MatrixXd x_guess = zeros(2,nnodes);
    x_guess.row(0)   = ones(1,nnodes);
    problem.phases(1).guess.states   = x_guess;
    problem.phases(1).guess.controls = ones(1,nnodes);
    problem.phases(1).guess.time     = linspace(0.0,3.0,nnodes);

    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-8;
    algorithm.collocation_method          = "Hermite-Simpson";
    algorithm.print_level                 = 0;

    psopt(solution, problem, algorithm);

    MatrixXd x = solution.get_states_in_phase(1);
    x1f = x(0, x.cols()-1);
    x2f = x(1, x.cols()-1);
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{
    double x1_user, x2_user, x1_def, x2_def;

    run_case(true,  x1_user, x2_user);
    run_case(false, x1_def,  x2_def);

    const double x1_user_exact = -0.375,      x2_user_exact = 3.5;
    const double x1_def_exact  = -1.0/6.0,    x2_def_exact  = 3.0;

    printf("\n  Initial history of a delay differential equation\n");
    printf("  x1' = -x1(t-1), x1(0)=1;  x2' = u(t-1), x2(0)=0, u == 1;  t in [0,3]\n\n");
    printf("  %-34s %14s %14s %12s\n", "history", "x1(3)", "closed form", "error");
    printf("  %s\n", "-----------------------------------------------------------------------------");
    printf("  %-34s %14.8f %14.8f %12.2e\n", "phi(t) = 1+t, psi(t) = 1-t (user)",
           x1_user, x1_user_exact, fabs(x1_user-x1_user_exact));
    printf("  %-34s %14.8f %14.8f %12.2e\n", "phi(t) = x1(0), psi(t) = u(0)",
           x1_def, x1_def_exact, fabs(x1_def-x1_def_exact));
    printf("\n  %-34s %14s %14s %12s\n", "history", "x2(3)", "closed form", "error");
    printf("  %s\n", "-----------------------------------------------------------------------------");
    printf("  %-34s %14.8f %14.8f %12.2e\n", "phi(t) = 1+t, psi(t) = 1-t (user)",
           x2_user, x2_user_exact, fabs(x2_user-x2_user_exact));
    printf("  %-34s %14.8f %14.8f %12.2e\n", "phi(t) = x1(0), psi(t) = u(0)",
           x2_def, x2_def_exact, fabs(x2_def-x2_def_exact));

    bool ok =    fabs(x1_user-x1_user_exact) < 1.0e-4
              && fabs(x2_user-x2_user_exact) < 1.0e-4
              && fabs(x1_def -x1_def_exact ) < 1.0e-4
              && fabs(x2_def -x2_def_exact ) < 1.0e-4;

    printf("\n  %s\n\n", ok ? "PASS" : "FAIL");

    printf("  The history is data of the problem, in the same sense as the initial\n");
    printf("  condition is: the two runs above solve different problems, and neither\n");
    printf("  answer is a discretization artefact.\n\n");

    return ok ? 0 : 1;
}
