//////////////////////////////////////////////////////////////////////////
//////////////////            flexmesh.cxx          //////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////            PSOPT  Example              //////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: A switching time a fixed mesh cannot see             /////
////////                                                              ////
//////// The integrated-residual transcription represents the state   ////
//////// and control of an element by polynomials of degree d through ////
//////// that element's own nodes.  A discontinuity in the INTERIOR   ////
//////// of an element is something those polynomials cannot be, so   ////
//////// a bang-bang switching time that falls inside an element is   ////
//////// a wall: the solve converges, it converges to the wrong       ////
//////// answer, and tightening the residual box does not move it,    ////
//////// because the residual being driven down is the residual of a  ////
//////// problem the mesh cannot state.                               ////
////////                                                              ////
//////// The remedy of Nie and Kerrigan is a FLEXIBLE MESH: the       ////
//////// element boundaries become decision variables, so the         ////
//////// optimisation puts a boundary ON the switching time instead   ////
//////// of refining around it.  In PSOPT that is one line,           ////
////////                                                              ////
////////     algorithm.ir_flexible_mesh = true;                       ////
////////                                                              ////
//////// with algorithm.ir_min_element_fraction setting the floor on  ////
//////// an element's width as a fraction of the uniform width.       ////
////////                                                              ////
//////// The problem is chosen so that the fixed mesh cannot win by   ////
//////// accident.  Minimum time for a double integrator from (0,0)   ////
//////// to (1,0) with SYMMETRIC control bounds switches at exactly   ////
//////// tf/2, which a uniform partition into an even number of       ////
//////// elements already has a boundary at -- a test the fixed mesh  ////
//////// passes for a reason that has nothing to do with the mesh.    ////
//////// Asymmetric bounds, u in [-1, 2], move the switch to tf/3:    ////
////////                                                              ////
////////     xdot1 = x2,  xdot2 = u,   (x1,x2): (0,0) -> (1,0)        ////
////////     u in [-1,2],  minimise tf                                ////
////////                                                              ////
////////     tf*  = sqrt(3) = 1.7320508,   switch at t = tf/3,        ////
////////     x1*(t) = t^2 and x2*(t) = 2t on the first arc,           ////
////////                                                              ////
//////// and no uniform partition into four elements has a boundary   ////
//////// at one third.                                                ////
////////                                                              ////
//////// Reference: Y. Nie and E. C. Kerrigan, "Solving optimal       ////
//////// control problems with non-smooth solutions using an          ////
//////// integrated residual method and flexible mesh", 2022 IEEE     ////
//////// 61st Conference on Decision and Control.                     ////
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

static const double TF_EXACT = 1.7320508075688772;   // sqrt(3)

//////////////////////////////////////////////////////////////////////////
///////////////////  Problem functions  //////////////////////////////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return tf;                                   // minimum time
}

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    return 0.0;
}

void dae(adouble* derivatives, adouble* path, adouble* states, adouble* controls,
         adouble* parameters, adouble& time, adouble* xad, int iphase,
         Workspace* workspace)
{
    derivatives[0] = states[1];
    derivatives[1] = controls[0];
}

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    e[0] = initial_states[0];        // x1(0) = 0
    e[1] = initial_states[1];        // x2(0) = 0
    e[2] = final_states[0];          // x1(tf) = 1
    e[3] = final_states[1];          // x2(tf) = 0
}

void linkages(adouble* linkages, adouble* xad, Workspace* workspace) {}

//////////////////////////////////////////////////////////////////////////
///////////////////  One solve  ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

struct Row { int flag; double tf; double rel_err; double arc1; };

static Row solve_it(bool flexible, double residual_bound, int nodes, int d,
                    bool verbose)
{
    Alg algorithm; Sol solution; Prob problem;
    Row out; out.flag = -1; out.tf = 0.0; out.rel_err = 1.0; out.arc1 = 1.0;

    problem.name        = "A switch a fixed mesh cannot see";
    problem.outfilename = "flexmesh.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -2.0, -2.0;
    problem.phases(1).bounds.upper.states   <<  2.0,  2.0;
    problem.phases(1).bounds.lower.controls << -1.0;
    problem.phases(1).bounds.upper.controls <<  2.0;
    problem.phases(1).bounds.lower.events   << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.upper.events   << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 0.5;
    problem.phases(1).bounds.upper.EndTime   = 6.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(2, nodes);
    problem.phases(1).guess.states.row(0) = linspace(0.0, 1.0, nodes);
    problem.phases(1).guess.controls = zeros(1, nodes);
    problem.phases(1).guess.time     = linspace(0.0, 1.73, nodes);

    algorithm.nlp_method            = "IPOPT";
    algorithm.scaling               = "automatic";
    algorithm.derivatives           = "automatic";
    algorithm.nlp_iter_max          = 2000;
    algorithm.nlp_tolerance         = 1.0e-8;
    algorithm.print_level           = 0;
    algorithm.mesh_refinement       = "manual";
    algorithm.collocation_method    = "Hermite-Simpson";
    algorithm.transcription_method  = "integrated-residual";
    algorithm.ir_local_order        = d;
    algorithm.ir_residual_nodes     = d + 2;

    // Minimum time needs the residual BOX, not the residual objective. Asked to
    // minimise the residual on a free horizon, the transcription answers with the
    // smoothest trajectory it can find, which is the longest one, and tf simply
    // runs to its upper bound. That is a correct answer to a question nobody meant
    // to ask. Here the cost is minimised subject to |r| <= delta instead.
    algorithm.ir_objective          = "cost";
    algorithm.ir_residual_bound     = residual_bound;

    algorithm.ir_flexible_mesh      = flexible;

    out.flag = psopt(solution, problem, algorithm);
    if (out.flag != 0) return out;

    MatrixXd t = solution.get_time_in_phase(1);
    MatrixXd x = solution.get_states_in_phase(1);
    const int n = (int) t.cols();

    out.tf      = t(0, n-1);
    out.rel_err = fabs(out.tf - TF_EXACT)/TF_EXACT;

    // Do the reported node TIMES belong to the reported node STATES? On the first
    // arc the control sits at its upper bound from rest at the origin, so
    // x1(t) = t^2 exactly, and the state played no part in deciding where the
    // nodes went. Stop short of the switch at tf/3 = 0.577.
    out.arc1 = 0.0;
    for (int k = 0; k < n; k++) {
        if ( t(0,k) > 0.55 ) break;
        out.arc1 = fmax( out.arc1, fabs( x(0,k) - t(0,k)*t(0,k) ) );
    }

    if (verbose) {
        printf("      element boundaries (t):");
        for (int e = 0; e*d < n; e++) printf(" %8.5f", t(0, e*d));
        printf("\n      the switch is at tf/3 = %8.5f\n", out.tf/3.0);
    }
    return out;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Main  ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main(void)
{
    const int nodes = 17;       // 16 intervals = 4 elements of degree 4
    const int d     = 4;

    printf("================================================================================\n");
    printf("  A switching time a fixed mesh cannot see\n");
    printf("  min tf,  xdot1 = x2, xdot2 = u,  u in [-1,2],  (0,0) -> (1,0)\n");
    printf("  exact:  tf* = sqrt(3) = %.7f,  switch at tf/3 = %.7f\n", TF_EXACT, TF_EXACT/3.0);
    printf("  %d nodes, local order %d, so %d elements of degree %d\n",
           nodes, d, (nodes-1)/d, d);
    printf("================================================================================\n");
    printf("  residual box      mesh        tf             rel err     max|x1-t^2| on arc 1\n");
    printf("--------------------------------------------------------------------------------\n");

    const double boxes[4] = { 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7 };
    for (int i = 0; i < 4; i++) {
        Row f = solve_it(false, boxes[i], nodes, d, false);
        Row g = solve_it(true,  boxes[i], nodes, d, false);
        printf("    %7.0e       fixed     %.9f    %.3e   %.2e\n",
               boxes[i], f.tf, f.rel_err, f.arc1);
        printf("    %7.0e       flexible  %.9f    %.3e   %.2e\n",
               boxes[i], g.tf, g.rel_err, g.arc1);
    }

    printf("--------------------------------------------------------------------------------\n");
    printf("  The fixed mesh is stuck at eight parts in a thousand at every tolerance: the\n");
    printf("  residual box is not what limits it, the mesh is. The flexible mesh clears it by\n");
    printf("  three orders of magnitude on the same nodes and the same element degree:\n\n");
    (void) solve_it(true, 1.0e-6, nodes, d, true);
    printf("\n  Note what it does with the freedom. Rather than place one boundary exactly on\n");
    printf("  the switch it BRACKETS the switch inside a thin element about one per cent of\n");
    printf("  the horizon wide, which isolates the discontinuity: the elements on either side\n");
    printf("  are then smooth, and a degree-4 polynomial is an excellent representation of a\n");
    printf("  parabola. algorithm.ir_min_element_fraction is the floor on how thin that\n");
    printf("  element may become; a width free to reach zero gives a singular local problem.\n");
    printf("\n  The last column reads differently in the two rows. For the flexible mesh it is\n");
    printf("  a check that the mesh PSOPT REPORTS is the mesh it solved on -- x1(t) = t^2 on\n");
    printf("  the first arc, and x1 took no part in deciding where the boundaries went, so it\n");
    printf("  is an independent test of the reported node times. For the fixed mesh the same\n");
    printf("  number is simply the error in the trajectory, whose times were never in doubt.\n");

    return 0;
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
