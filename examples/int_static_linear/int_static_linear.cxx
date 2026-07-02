//////////////////////////////////////////////////////////////////////////
//////////////////      int_static_linear.cxx      ///////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Integer static parameter (minimal example)   //////////////
//////// Last modified: 02 July 2026                        ///////////////
//////// Reference:     none (constructed analytical test)  ///////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2026        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

// Purpose
// -------
// The smallest self-contained demonstration of an integer (discrete-valued)
// static parameter. A single static parameter p may take only the admissible
// values {0,1,2,3}; psopt_solve_integer enumerates them, solves the fixed-
// parameter OCP for each, and returns the best.
//
// System (a linear double integrator), t in [0,1]:
//
//     s_dot = v          (position)
//     v_dot = p          (velocity;  p a static parameter = constant accel.)
//     s(0) = 0,  v(0) = 0     ==>   v(1) = p.
//
// Cost (Mayer term on the terminal velocity):
//
//     J = ( v(1) - c )^2 = (p - c)^2,   c = 2.3.
//
// The discrete optimum is the admissible value nearest c: p = 2, J = 0.09.
// The example declares p integer, solves with psopt_solve_integer, reads the
// selected value back with reconstruct_integer_parameter, and checks it
// against the closed form.
//
// Declaring several integer parameters (on one or more phases) is the same
// pattern: psopt_solve_integer enumerates their full Cartesian product, up to
// algorithm.max_integer_combinations. See include/integer_parameters.h.

#include "integer_parameters.h"

#include <iomanip>

using namespace std;

//////////////////////////////////////////////////////////////////////////
///////////////////  Problem constants  //////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static const double TARGET = 2.3;   // c in J = (v(1) - c)^2

//////////////////////////////////////////////////////////////////////////
///////////////////  End point (Mayer) cost function  ////////////////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    adouble vf = final_states[1];      // terminal velocity v(1) = p
    adouble d  = vf - TARGET;
    return d*d;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Integrand (Lagrange) cost function  /////////////////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{
    return 0.0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's  ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
    adouble v = states[1];
    adouble p = parameters[0];

    derivatives[0] = v;          // s_dot = v
    derivatives[1] = p;          // v_dot = p
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function  /////////////////////////
//////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    e[0] = initial_states[0];    // s(0)
    e[1] = initial_states[1];    // v(0)
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function  /////////////////
//////////////////////////////////////////////////////////////////////////

void linkages(adouble* linkages, adouble* xad, Workspace* workspace)
{
    // Single phase problem - no linkages.
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Main routine  ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main(void)
{
    Alg  algorithm;
    Sol  solution;
    Prob problem;

    problem.name        = "Integer static parameter (minimal example)";
    problem.outfilename = "int_static_linear.txt";

    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates     = 2;      // position, velocity
    problem.phases(1).ncontrols   = 0;      // no controls
    problem.phases(1).nparameters = 1;      // the static parameter p
    problem.phases(1).nevents     = 2;      // s(0) = 0, v(0) = 0
    problem.phases(1).npath       = 0;
    problem.phases(1).nodes       = (RowVectorXi(1) << 20).finished();

    psopt_level2_setup(problem, algorithm);

    // State bounds (generous; s(1) = p/2, v(1) = p, with p in [0,3]).
    problem.phases(1).bounds.lower.states(0) = -20.0;
    problem.phases(1).bounds.upper.states(0) =  20.0;
    problem.phases(1).bounds.lower.states(1) = -20.0;
    problem.phases(1).bounds.upper.states(1) =  20.0;

    // Static parameter bounds. The convex hull of the admissible set is a
    // sensible relaxed range; psopt_solve_integer pins the parameter to each
    // admissible value in turn and restores these bounds afterwards.
    problem.phases(1).bounds.lower.parameters(0) = 0.0;
    problem.phases(1).bounds.upper.parameters(0) = 3.0;

    // Events: s(0) = 0, v(0) = 0.
    problem.phases(1).bounds.lower.events(0) = 0.0;
    problem.phases(1).bounds.upper.events(0) = 0.0;
    problem.phases(1).bounds.lower.events(1) = 0.0;
    problem.phases(1).bounds.upper.events(1) = 0.0;

    // Fixed time horizon [0,1].
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0;
    problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    // Initial guess. Midpoint of the parameter box is a neutral start.
    int nnodes = problem.phases(1).nodes(0);
    problem.phases(1).guess.states        = zeros(2, nnodes);
    problem.phases(1).guess.time          = linspace(0.0, 1.0, nnodes);
    problem.phases(1).guess.parameters(0) = 1.5;

    algorithm.nlp_method       = "IPOPT";
    algorithm.scaling          = "automatic";
    algorithm.derivatives      = "automatic";
    algorithm.nlp_iter_max     = 1000;
    algorithm.nlp_tolerance    = 1.e-8;
    algorithm.mesh_refinement  = "manual";
    algorithm.print_level      = 0;         // quiet per-combination solves

    //--------------------------------------------------------------------
    // Declare the static parameter (index 0 of phase 1) as integer, giving
    // its admissible values, then solve by enumeration.
    //--------------------------------------------------------------------
    RowVectorXd values(4);
    values << 0.0, 1.0, 2.0, 3.0;
    declare_integer_parameter(problem, 1, 0, values);

    int rc = psopt_solve_integer(solution, problem, algorithm);

    // Read the selected admissible value back from the solved problem.
    IntegerParameterReconstruction sel =
        reconstruct_integer_parameter(solution, problem, 1, 0);

    //--------------------------------------------------------------------
    // Report and check against the closed form (p = 2, J = 0.09).
    //--------------------------------------------------------------------
    cout.setf(ios::fixed);
    cout << setprecision(6);

    cout << "\n=====================================================================\n";
    cout << " Integer static parameter (target c = " << TARGET << ")\n";
    cout << "=====================================================================\n\n";
    cout << " Admissible values : {0, 1, 2, 3}\n";
    cout << " Solver return code: " << rc << "\n\n";
    cout << " Selected p        = " << sel.value
         << "   (admissible index " << sel.index << ")\n";
    cout << " Optimal cost J    = " << solution.cost << "\n\n";

    const double p_expected = 2.0;
    const double J_expected = (p_expected - TARGET)*(p_expected - TARGET);   // 0.09
    const double tol        = 1.e-5;

    bool ok = (solution.error_flag == 0)
           && (fabs(sel.value    - p_expected) < 1.e-9)
           && (fabs(solution.cost - J_expected) < tol);

    cout << " Closed form       ->  p = " << p_expected
         << ",  J = " << J_expected << "\n";
    cout << "=====================================================================\n";
    cout << " Result: " << (ok ? "PASS" : "FAIL") << "\n";
    cout << "=====================================================================\n";

    return ok ? 0 : 1;
}

//////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     /////////////////////////////
//////////////////////////////////////////////////////////////////////////
