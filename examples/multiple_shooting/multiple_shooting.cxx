//////////////////////////////////////////////////////////////////////////
//////////////       multiple_shooting.cxx          //////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////            PSOPT  Example              //////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Multiple shooting, and how to read its answers        ////
////////                                                              ////
//////// In the multiple-shooting transcription the decision variables ///
//////// are the states at the SEGMENT BOUNDARIES and the controls;    ///
//////// the trajectory inside a segment is produced by a fixed-step   ///
//////// RK4 recorded on the same tape as everything else, and the     ///
//////// constraints that tie the segments together are the matching   ///
//////// conditions                                                    ///
////////                                                               ///
////////     x_{k+1} - phi( x_k, u_k, p, h_k ) = 0.                    ///
////////                                                               ///
//////// It is selected with                                           ///
////////                                                               ///
////////     algorithm.transcription_method = "multiple-shooting";     ///
////////                                                               ///
//////// and configured with ms_steps_per_segment, which sets how many ///
//////// RK4 steps cross a segment; ms_control_parameterisation, which ///
//////// is "constant" or "linear"; and ms_path_samples, which sets how ///
//////// many interior points of a segment the path constraints are    ///
//////// also enforced at.                                             ///
////////                                                               ///
//////// This example makes four points, each with a number attached,  ///
//////// and two of them are cautions rather than selling points.      ///
////////                                                               ///
//////// Reference for the method: H. G. Bock and K. J. Plitt, "A      ///
//////// multiple shooting algorithm for direct solution of optimal    ///
//////// control problems", IFAC World Congress, 1984.                 ///
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2026         ///////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"
#include <cstdio>
#include <cmath>

using namespace PSOPT;

//////////////////////////////////////////////////////////////////////////
///////////////////  Problem functions  //////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// 0: minimum energy, (0,0) -> (1,0) on [0,1].  J* = 6, u*(t) = 6 - 12t,
//    costates l1 = -12 and l2 = 12t - 6.
// 1: Bryson and Denham's problem, the same dynamics with x <= 1/9 and
//    x(0)=0, v(0)=1, x(1)=0, v(1)=-1.  J* = 4.
static int problem_case = 0;

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{ return 0.0; }

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase, Workspace* workspace)
{ return 0.5*controls[0]*controls[0]; }

void dae(adouble* derivatives, adouble* path, adouble* states, adouble* controls,
         adouble* parameters, adouble& time, adouble* xad, int iphase,
         Workspace* workspace)
{
    derivatives[0] = states[1];
    derivatives[1] = controls[0];
    if ( problem_case == 1 ) path[0] = states[0];
}

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    e[0] = initial_states[0];  e[1] = initial_states[1];
    e[2] = final_states[0];    e[3] = final_states[1];
}

void linkages(adouble* linkages, adouble* xad, Workspace* workspace) {}

//////////////////////////////////////////////////////////////////////////
///////////////////  One solve  ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

struct Row { int flag; double J; double l1_err; double l2_err; };

static Row solve_it(int which, const char* transcription, int segments, int steps,
                    const char* upar, int path_samples, bool costates)
{
    problem_case = which;

    Alg algorithm; Sol solution; Prob problem;
    Row out; out.flag = -1; out.J = 0.0; out.l1_err = 0.0; out.l2_err = 0.0;

    const int nodes = segments + 1;

    problem.name        = "Multiple shooting";
    problem.outfilename = "multiple_shooting.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = ( which == 1 ) ? 1 : 0;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -5.0, -5.0;
    problem.phases(1).bounds.upper.states   <<  5.0,  5.0;
    problem.phases(1).bounds.lower.controls(0) = -30.0;
    problem.phases(1).bounds.upper.controls(0) =  30.0;

    if ( which == 1 ) {
        problem.phases(1).bounds.lower.path(0) = -5.0;
        problem.phases(1).bounds.upper.path(0) =  1.0/9.0;
        problem.phases(1).bounds.lower.events << 0.0, 1.0, 0.0, -1.0;
        problem.phases(1).bounds.upper.events << 0.0, 1.0, 0.0, -1.0;
    }
    else {
        problem.phases(1).bounds.lower.events << 0.0, 0.0, 1.0, 0.0;
        problem.phases(1).bounds.upper.events << 0.0, 0.0, 1.0, 0.0;
    }
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0;
    problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(2, nodes);
    if ( which == 1 ) problem.phases(1).guess.states.row(1) = linspace( 1.0, -1.0, nodes);
    else              problem.phases(1).guess.states.row(0) = linspace( 0.0,  1.0, nodes);
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
    algorithm.transcription_method  = transcription;
    algorithm.ms_steps_per_segment  = steps;
    algorithm.ms_control_parameterisation = upar;
    algorithm.ms_path_samples             = path_samples;

    out.flag = psopt(solution, problem, algorithm);
    if (out.flag != 0) return out;

    out.J = solution.cost;

    if (costates) {
        DMatrix L = solution.get_dual_costates_in_phase(1);
        DMatrix T = solution.get_time_in_phase(1);
        for (int q = 0; q < T.cols(); q++) {
            out.l1_err = fmax( out.l1_err, fabs( L(0,q) - (-12.0) ) );
            out.l2_err = fmax( out.l2_err, fabs( L(1,q) - (12.0*T(0,q) - 6.0) ) );
        }
    }
    return out;
}

// The optimum of the transcription's OWN problem when the control is held constant on M
// equal segments. The continuous problem is least-norm in u; restricting u to the
// M-dimensional space of piecewise-constant functions and projecting the two linear
// constraints onto it gives this in closed form, and it is what a correct implementation
// must return -- not the continuous optimum 6, which it approaches from above as M grows.
static double discrete_optimum(double M) { return 6.0*M*M/(M*M - 1.0); }

//////////////////////////////////////////////////////////////////////////
///////////////////  Main  ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main(void)
{
    printf("================================================================================\n");
    printf("  Multiple shooting, and how to read its answers\n");
    printf("================================================================================\n");

    printf("\n  1. It answers its own problem exactly, which is not the same as answering\n");
    printf("     the continuous one. With the control held constant on M segments the\n");
    printf("     discrete optimum is 6M^2/(M^2-1) in closed form:\n\n");
    printf("        M      computed J        6M^2/(M^2-1)     difference\n");
    const int segs[4] = { 5, 10, 20, 40 };
    for (int q = 0; q < 4; q++) {
        Row r = solve_it(0, "multiple-shooting", segs[q], 10, "constant", 0, false);
        printf("      %3d    %.10f     %.10f     %.1e\n",
               segs[q], r.J, discrete_optimum(segs[q]),
               fabs(r.J - discrete_optimum(segs[q])));
    }
    printf("\n     A result close to 6 instead would mean the transcription was solving a\n");
    printf("     different problem and getting a plausible answer, which is the failure mode\n");
    printf("     worth being able to tell apart, because it looks like success.\n");

    printf("\n  2. The control parameterisation is not a detail. This problem's optimal\n");
    printf("     control is u*(t) = 6 - 12t, exactly linear, so a linear parameterisation\n");
    printf("     contains the answer and attains the continuous optimum at five segments:\n\n");
    printf("        M      constant u        linear u\n");
    for (int q = 0; q < 3; q++) {
        Row a = solve_it(0, "multiple-shooting", segs[q], 10, "constant", 0, false);
        Row b = solve_it(0, "multiple-shooting", segs[q], 10, "linear",   0, false);
        printf("      %3d    %.10f     %.10f\n", segs[q], a.J, b.J);
    }
    printf("\n     Compare a shooting method carrying a piecewise-constant control against\n");
    printf("     collocation's piecewise polynomials and it loses, for a reason that has\n");
    printf("     nothing to do with shooting.\n");

    printf("\n  3. A path constraint imposed only at the segment boundaries is not the path\n");
    printf("     constraint that was written down. Bryson and Denham's problem has J* = 4\n");
    printf("     with x <= 1/9, and NO feasible trajectory can cost less than that:\n\n");
    printf("        M    boundaries only   2 interior samples   4 interior samples\n");
    const int bds[2] = { 10, 20 };
    for (int q = 0; q < 2; q++) {
        Row a = solve_it(1, "multiple-shooting", bds[q], 10, "linear", 0, false);
        Row b = solve_it(1, "multiple-shooting", bds[q], 10, "linear", 2, false);
        Row c = solve_it(1, "multiple-shooting", bds[q], 10, "linear", 4, false);
        printf("      %3d    %.7f         %.7f            %.7f\n", bds[q], a.J, b.J, c.J);
    }
    printf("\n     Below the optimum is not a better answer. Between two segment boundaries\n");
    printf("     there is an integrator and nothing at all constraining what it does, so the\n");
    printf("     constraint leaks; algorithm.ms_path_samples closes the gap.\n");

    printf("\n  4. The costates come back, and they are not the multipliers of the matching\n");
    printf("     conditions. Those are a discrete adjoint, but the covector mapping PSOPT\n");
    printf("     applies to collocation defects is not written for them. The costates below\n");
    printf("     are recovered by integrating the adjoint equation backwards along the\n");
    printf("     converged primal, which is defined for any transcription that produces a\n");
    printf("     trajectory. Against the closed form l1 = -12 and l2 = 12t - 6:\n\n");
    {
        Row a = solve_it(0, "collocation",       20, 10, "linear",   0, true);
        Row b = solve_it(0, "multiple-shooting", 20, 10, "linear",   0, true);
        Row c = solve_it(0, "multiple-shooting", 20, 10, "constant", 0, true);
        printf("        collocation                  max|l1+12| = %.2e   max|l2-(12t-6)| = %.2e\n", a.l1_err, a.l2_err);
        printf("        multiple shooting, linear    max|l1+12| = %.2e   max|l2-(12t-6)| = %.2e\n", b.l1_err, b.l2_err);
        printf("        multiple shooting, constant  max|l1+12| = %.2e   max|l2-(12t-6)| = %.2e\n", c.l1_err, c.l2_err);
    }
    printf("\n     The last line is larger because the constant-control solution really is\n");
    printf("     1.5 per cent away from the continuous optimum, so its costates are too.\n");
    printf("     That is the transcription being consistent with itself, not an error.\n");

    printf("\n--------------------------------------------------------------------------------\n");
    printf("  What this transcription does not yet have: automatic mesh refinement, which\n");
    printf("  asks a different question here (how many segments is a question about the\n");
    printf("  control parameterisation, while the integration error is set by\n");
    printf("  ms_steps_per_segment), and an implicit integrator, without which an index-1\n");
    printf("  DAE cannot be propagated. Refine by giving a sequence of segment counts in\n");
    printf("  problem.phases(i).nodes; each mesh is hot-started from the one before it.\n");

    return 0;
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
