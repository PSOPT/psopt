//////////////////////////////////////////////////////////////////////////
//////////////////           mixed_ic.cxx           //////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:  Mixed continuous + integer control (remap check)    //////
//////// Method: integer-control declaration API vs manual weights   //////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2025        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////
//
// A small, essentially convex regulation problem with BOTH a continuous control
// u and a binary control v in {0,1} in the same phase. Its purpose is to exercise
// the continuous+integer control remap that the single-control benchmarks (F-8,
// Lotka) do not: it is solved twice and the objectives are compared.
//
//   states:  x0, x1 and cost accumulator x2
//   dynamics: x0' = x1
//             x1' = -x0 - 0.2 x1 + u + 0.5 v
//             x2' = x0^2 + x1^2 + u^2
//   x(0) = (1, 0, 0),  T = 2,  u in [-3,3],  v in {0,1}
//   minimise x2(T)
//
//   Formulation A (declaration API): two controls, u at index 0 and the integer
//   control v at index 1; declare_integer_control does the convexification. The
//   continuous control at user index 0 must be remapped to internal index 0 with
//   the two mode weights appended after it (the k<c case in the wrapper).
//
//   Formulation B (manual): three controls u, w0, w1 with an explicit SOS1 path
//   constraint and the convex combination written by hand.
//
// The problem is convex in the relaxed weights, so both formulations reach the
// same optimum; matching objectives validate the mixed remap end to end.
//
// AI provenance: drafted with AI assistance.

#include "psopt.h"
#include "integer_controls.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////
///////////////////  Shared cost / events / linkages //////////////////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return final_states[2];
}

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{
    return 0.0;
}

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    e[0] = initial_states[0];
    e[1] = initial_states[1];
    e[2] = initial_states[2];
}

void linkages(adouble* linkages, adouble* xad, Workspace* workspace) {}

//////////////////////////////////////////////////////////////////////////
///////////////////  Formulation A: declaration API dae ///////////////////
//////////////////////////////////////////////////////////////////////////

void dae_decl(adouble* d, adouble* path, adouble* states, adouble* controls,
              adouble* parameters, adouble& time, adouble* xad,
              int iphase, Workspace* workspace)
{
    adouble x0 = states[0], x1 = states[1];
    adouble u  = controls[0];   // continuous control (user index 0)
    adouble v  = controls[1];   // integer control    (user index 1)

    d[0] = x1;
    d[1] = -x0 - 0.2*x1 + u + 0.5*v;
    d[2] = x0*x0 + x1*x1 + u*u;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Formulation B: manual weights dae ////////////////////
//////////////////////////////////////////////////////////////////////////

void dae_manual(adouble* d, adouble* path, adouble* states, adouble* controls,
                adouble* parameters, adouble& time, adouble* xad,
                int iphase, Workspace* workspace)
{
    adouble x0 = states[0], x1 = states[1];
    adouble u  = controls[0];
    adouble w0 = controls[1], w1 = controls[2];

    d[0] = x1;
    d[1] = -x0 - 0.2*x1 + u + 0.5*(w0*0.0 + w1*1.0);
    d[2] = x0*x0 + x1*x1 + u*u;
    path[0] = w0 + w1;   // SOS1
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Solvers //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static void common_bounds_and_alg(Prob& problem, Alg& algorithm, int nnodes)
{
    problem.phases(1).bounds.lower.states << -10.0, -10.0, 0.0;
    problem.phases(1).bounds.upper.states <<  10.0,  10.0, 1.0e4;

    problem.phases(1).bounds.lower.events << 1.0, 0.0, 0.0;
    problem.phases(1).bounds.upper.events << 1.0, 0.0, 0.0;

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 2.0;
    problem.phases(1).bounds.upper.EndTime   = 2.0;

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.nlp_iter_max       = 2000;
    algorithm.nlp_tolerance      = 1.e-8;
    algorithm.collocation_method = "trapezoidal";
    algorithm.mesh_refinement    = "manual";
}

static double solve_declaration(int nnodes)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "mixed_ic (declaration API)"; problem.outfilename = "mixed_ic_A.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 3;
    problem.phases(1).ncontrols = 2;   // u, v
    problem.phases(1).nevents   = 3;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << nnodes;
    psopt_level2_setup(problem, algorithm);

    RowVectorXd values(2); values << 0.0, 1.0;
    declare_integer_control(problem, 1, 1, values);   // v is control index 1

    problem.phases(1).bounds.lower.controls << -3.0, 0.0;   // u in [-3,3]; v slot ignored
    problem.phases(1).bounds.upper.controls <<  3.0, 1.0;
    common_bounds_and_alg(problem, algorithm, nnodes);

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae_decl;
    problem.events         = &events;
    problem.linkages       = &linkages;

    RowVectorXd tg = linspace(0.0, 2.0, nnodes);
    DMatrix sg(3, nnodes);
    sg.row(0) = linspace(1.0, 0.0, nnodes);
    sg.row(1) = zeros(1, nnodes);
    sg.row(2) = zeros(1, nnodes);
    DMatrix cg(2, nnodes);
    cg.row(0) = zeros(1, nnodes);   // u guess
    cg.row(1) = zeros(1, nnodes);   // v guess (value 0 -> mode 0)
    problem.phases(1).guess.states   = sg;
    problem.phases(1).guess.controls = cg;
    problem.phases(1).guess.time     = tg;

    psopt(solution, problem, algorithm);
    DMatrix x = solution.get_states_in_phase(1);
    return x(2, x.cols() - 1);
}

static double solve_manual(int nnodes)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "mixed_ic (manual weights)"; problem.outfilename = "mixed_ic_B.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 3;
    problem.phases(1).ncontrols = 3;   // u, w0, w1
    problem.phases(1).nevents   = 3;
    problem.phases(1).npath     = 1;   // SOS1
    problem.phases(1).nodes     << nnodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.controls << -3.0, 0.0, 0.0;
    problem.phases(1).bounds.upper.controls <<  3.0, 1.0, 1.0;
    problem.phases(1).bounds.lower.path << 1.0;
    problem.phases(1).bounds.upper.path << 1.0;
    common_bounds_and_alg(problem, algorithm, nnodes);

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae_manual;
    problem.events         = &events;
    problem.linkages       = &linkages;

    RowVectorXd tg = linspace(0.0, 2.0, nnodes);
    DMatrix sg(3, nnodes);
    sg.row(0) = linspace(1.0, 0.0, nnodes);
    sg.row(1) = zeros(1, nnodes);
    sg.row(2) = zeros(1, nnodes);
    DMatrix cg(3, nnodes);
    cg.row(0) = zeros(1, nnodes);          // u
    cg.row(1) = ones(1, nnodes);           // w0 = 1
    cg.row(2) = zeros(1, nnodes);          // w1 = 0
    problem.phases(1).guess.states   = sg;
    problem.phases(1).guess.controls = cg;
    problem.phases(1).guess.time     = tg;

    psopt(solution, problem, algorithm);
    DMatrix x = solution.get_states_in_phase(1);
    return x(2, x.cols() - 1);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Main /////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main(void)
{
    const int nnodes = 40;

    double JA = solve_declaration(nnodes);
    double JB = solve_manual(nnodes);

    printf("\n============= mixed continuous + integer remap check =============\n");
    printf("Formulation A (declaration API) : objective = %.8f\n", JA);
    printf("Formulation B (manual weights)  : objective = %.8f\n", JB);
    printf("Absolute difference             : %.3e\n", fabs(JA - JB));
    printf("=================================================================\n");
    if (fabs(JA - JB) < 1.e-5)
        printf("PASS: the declaration API and the manual weights agree; the mixed\n"
               "continuous+integer control remap is validated end to end.\n\n");
    else
        printf("MISMATCH: objectives differ by more than the tolerance.\n\n");

    return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     //////////////////////////////
//////////////////////////////////////////////////////////////////////////
