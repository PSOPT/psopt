//////////////////////////////////////////////////////////////////////////
//////////////////            fuller.cxx           ////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////            PSOPT  Example              //////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Fuller's problem - a chattering optimal control  /////////
////////                                                          ////////
//////// This example showcases the integrated-residual family on a  ////
//////// problem whose optimal control chatters (infinitely many     ////
//////// bang-bang switches of decreasing duration).  It contrasts   ////
//////// three things on the same problem:                           ////
////////   A. direct collocation - bang-bang, but with poor and      ////
////////      uncontrolled accuracy between the mesh nodes;          ////
////////   B. the integrated-residual optimality step with an        ////
////////      explicit residual box (min cost s.t. |r| <= delta):    ////
////////      the cost becomes a certified lower bound that rises     ////
////////      towards the true optimum as delta is tightened - the   ////
////////      accuracy/cost Pareto trade-off collocation cannot      ////
////////      offer;                                                 ////
////////   C. the Nie-Kerrigan flexible-order local representation:   ////
////////      a higher local order lowers the residual error for a   ////
////////      given box (order 0 -> 2 gives ~5x here).               ////
////////                                                              ////
//////// Reference: L. Nita, E. C. Kerrigan, E. M. G. Vila and       ////
//////// Y. Nie, "Solving optimal control problems with non-smooth   ////
//////// solutions using an integrated residual method and flexible  ////
//////// mesh", 2022 IEEE 61st CDC.  Their eq. (11):                 ////
////////   min integral_0^T p^2 dt   s.t.  pddot = u,                ////
////////   u in [-0.01, 0.01], p(0)=p(T)=pdot(T)=0, pdot(0)=1, T=300 ////
//////// State form: x1 = p, x2 = pdot; xdot1 = x2, xdot2 = u.       ////
////////                                                              ////
//////// Notes on method usage (honest limitations):                 ////
////////  * The cold-start integrated-residual solve of a chattering ////
////////    problem is fragile, so the optimality step is warm-      ////
////////    started from a collocation solution (the feasibility ->  ////
////////    optimality workflow of the reference).                   ////
////////  * The explicit residual box (algorithm.ir_residual_bound)  ////
////////    is used here so that each solve targets a directly        ////
////////    specified accuracy - the reference's own interface - and ////
////////    so the accuracy/cost Pareto front can be traced out by    ////
////////    varying delta.  The automatic ir_dair schedule            ////
////////    delta = K*(h/T)^2 (horizon-normalised, hence scale        ////
////////    invariant) is also viable on this long horizon (T=300);   ////
////////    the explicit box is chosen only for the clearer control   ////
////////    it gives over the delta sequence.                         ////
////////  * Raising the local order beyond 2 does not further reduce ////
////////    the error here: on a fixed mesh the switch discontinui-  ////
////////    ties induce Gibbs-like local error that higher-degree    ////
////////    polynomials do not remove.  The reference addresses this ////
////////    with a flexible (moving) mesh, which PSOPT does not have.////
////////                                                              ////
//////// This example was prepared with AI assistance (Claude) and   ////
//////// validated against the behaviour reported in the reference.  ////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2025        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"
#include <cstdio>
#include <string>

using namespace PSOPT;

//////////////////////////////////////////////////////////////////////////
///////////////////  Problem functions  //////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static const double U_MAX = 0.01;   // control bound
static const double T_HOR = 300.0;  // fixed time horizon

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return 0.0;
}

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{
    adouble p = states[0];
    return p*p;                       // integral of p^2 dt
}

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
    derivatives[0] = states[1];       // pdot  = x2
    derivatives[1] = controls[0];     // pddot = u
}

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    e[0] = initial_states[0];         // p(0)     = 0
    e[1] = initial_states[1];         // pdot(0)  = 1
    e[2] = final_states[0];           // p(T)     = 0
    e[3] = final_states[1];           // pdot(T)  = 0
}

void linkages(adouble* linkages, adouble* xad, Workspace* workspace) {}

//////////////////////////////////////////////////////////////////////////
///////////////////  Small result record  ////////////////////////////////
//////////////////////////////////////////////////////////////////////////

struct Row { double J; int switches; double umax, R2, rle; };

// Build and solve one Fuller instance.  If guess_* are non-empty they are used
// as the initial guess (warm start); otherwise a plain flat guess is used.
static Row solve_fuller(int nodes, const std::string& trans,
                        const std::string& ir_obj, int ir_order,
                        double ir_bound, int ir_resnodes, double nlptol,
                        const DMatrix& gx, const DMatrix& gu, const DMatrix& gt)
{
    Alg  algorithm; Sol solution; Prob problem;
    problem.name = "Fuller chattering problem";
    problem.outfilename = "fuller.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -200.0, -5.0;
    problem.phases(1).bounds.upper.states   <<  200.0,  5.0;
    problem.phases(1).bounds.lower.controls << -U_MAX;
    problem.phases(1).bounds.upper.controls <<  U_MAX;
    problem.phases(1).bounds.lower.events   << 0.0, 1.0, 0.0, 0.0;
    problem.phases(1).bounds.upper.events   << 0.0, 1.0, 0.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = T_HOR;
    problem.phases(1).bounds.upper.EndTime   = T_HOR;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    if (gx.cols() > 0) {
        problem.phases(1).guess.states   = gx;
        problem.phases(1).guess.controls = gu;
        problem.phases(1).guess.time     = gt;
    } else {
        DMatrix sg(2, nodes);
        sg << linspace(0.0, 0.0, nodes), linspace(1.0, 0.0, nodes);
        problem.phases(1).guess.states   = sg;
        problem.phases(1).guess.controls = -0.001*ones(1, nodes);
        problem.phases(1).guess.time     = linspace(0.0, T_HOR, nodes);
    }

    algorithm.nlp_method     = "IPOPT";
    algorithm.scaling        = "automatic";
    algorithm.derivatives    = "automatic";
    algorithm.nlp_iter_max   = 5000;
    algorithm.nlp_tolerance  = nlptol;
    algorithm.transcription_method = trans;
    algorithm.collocation_method   = "Hermite-Simpson";
    algorithm.print_level    = 0;
    if (trans == "integrated-residual") {
        algorithm.ir_residual_nodes = ir_resnodes;
        if (ir_order >= 2) algorithm.ir_local_order = ir_order;
        if (ir_bound >= 0.0) {           // explicit residual box: min cost s.t. |r|<=delta
            algorithm.ir_objective     = "cost";
            algorithm.ir_residual_bound = ir_bound;
        } else {
            algorithm.ir_objective     = ir_obj;
        }
    }

    psopt(solution, problem, algorithm);

    DMatrix x = solution.get_states_in_phase(1);
    DMatrix u = solution.get_controls_in_phase(1);
    DMatrix errm = solution.get_relative_local_error_in_phase(1);
    int n = (int)u.cols();

    Row r; r.J = solution.cost; r.switches = 0; r.umax = 0.0; r.R2 = 0.0;
    for (int k = 0; k < n; ++k) r.umax = fmax(r.umax, fabs(u(0,k)));
    for (int k = 0; k+1 < n; ++k) if (u(0,k)*u(0,k+1) < 0.0) r.switches++;
    for (int k = 1; k+1 < n; ++k) r.R2 += fabs(u(0,k+1) - 2.0*u(0,k) + u(0,k-1));
    r.rle = 0.0;
    for (int i = 0; i < errm.rows(); ++i)
        for (int j = 0; j < errm.cols(); ++j) r.rle = fmax(r.rle, fabs(errm(i,j)));
    return r;
}

// Robust collocation solve, returned as a warm-start trajectory.
static void warmstart(int nodes, DMatrix& x, DMatrix& u, DMatrix& t)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Fuller warm start"; problem.outfilename = "fuller_warm.txt";
    problem.nphases = 1; problem.nlinkages = 0; psopt_level1_setup(problem);
    problem.phases(1).nstates = 2; problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents = 4; problem.phases(1).npath = 0;
    problem.phases(1).nodes << nodes; psopt_level2_setup(problem, algorithm);
    problem.phases(1).bounds.lower.states   << -200.0, -5.0;
    problem.phases(1).bounds.upper.states   <<  200.0,  5.0;
    problem.phases(1).bounds.lower.controls << -U_MAX;
    problem.phases(1).bounds.upper.controls <<  U_MAX;
    problem.phases(1).bounds.lower.events   << 0.0, 1.0, 0.0, 0.0;
    problem.phases(1).bounds.upper.events   << 0.0, 1.0, 0.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = T_HOR; problem.phases(1).bounds.upper.EndTime = T_HOR;
    problem.integrand_cost = &integrand_cost; problem.endpoint_cost = &endpoint_cost;
    problem.dae = &dae; problem.events = &events; problem.linkages = &linkages;
    DMatrix sg(2, nodes);
    sg << linspace(0.0, 0.0, nodes), linspace(1.0, 0.0, nodes);
    problem.phases(1).guess.states = sg;
    problem.phases(1).guess.controls = -0.001*ones(1, nodes);
    problem.phases(1).guess.time = linspace(0.0, T_HOR, nodes);
    algorithm.nlp_method = "IPOPT"; algorithm.scaling = "automatic";
    algorithm.derivatives = "automatic"; algorithm.nlp_iter_max = 5000;
    algorithm.nlp_tolerance = 1.0e-7; algorithm.collocation_method = "Hermite-Simpson";
    algorithm.print_level = 0;
    psopt(solution, problem, algorithm);
    x = solution.get_states_in_phase(1);
    u = solution.get_controls_in_phase(1);
    t = solution.get_time_in_phase(1);
}

static void print_row(const char* label, int nodes, const Row& r)
{
    printf("  %-30s N=%-4d  J=%.6e  switches=%-3d  umax=%.4f  R2=%.4f  rle=%.2e\n",
           label, nodes, r.J, r.switches, r.umax, r.R2, r.rle);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Main  ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main(void)
{
    DMatrix none;   // empty -> flat guess

    printf("================================================================================\n");
    printf("  Fuller's chattering problem: collocation vs integrated-residual / Nie-Kerrigan\n");
    printf("================================================================================\n");

    // ---- Part A: direct collocation baseline (bang-bang, but uncontrolled accuracy) ----
    // One robust collocation solve (N=80) also serves as the warm start for B and C.
    printf("  A. direct collocation (Hermite-Simpson)\n");
    DMatrix wx, wu, wt; warmstart(80, wx, wu, wt);
    Row colA = solve_fuller(80, "collocation", "residual", 0, -1.0, 4, 1.0e-7, none, none, none);
    print_row("collocation", 80, colA);

    // ---- Part B: integrated-residual optimality step, explicit residual box ----
    // Warm-started from the collocation solution; the cost rises as a certified lower
    // bound as the box delta is tightened.
    printf("  B. integrated-residual, min cost s.t. |r| <= delta (warm-started, cubic-Hermite)\n");
    const double deltas[4] = {1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6};
    for (int i = 0; i < 4; ++i) {
        Row r = solve_fuller(80, "integrated-residual", "cost", 0, deltas[i], 4, 1.0e-7, wx, wu, wt);
        char lab[64]; snprintf(lab, sizeof(lab), "  box delta = %.0e", deltas[i]);
        print_row(lab, 80, r);
    }

    // ---- Part C: Nie-Kerrigan flexible-order improves accuracy for a given box ----
    // Local order 2 needs an even element count (here N=81 -> 40 elements of degree 2);
    // the N=80 collocation trajectory is reused as the (interpolated) warm start.
    printf("  C. Nie-Kerrigan local order at fixed box delta = 1e-3 (warm-started)\n");
    Row c0 = solve_fuller(81, "integrated-residual", "cost", 0, 1.0e-3, 4, 1.0e-6, wx, wu, wt);
    print_row("  local order 0 (cubic-Hermite)", 81, c0);
    Row c2 = solve_fuller(81, "integrated-residual", "cost", 2, 1.0e-3, 6, 1.0e-6, wx, wu, wt);
    print_row("  local order 2 (Nie-Kerrigan)", 81, c2);

    printf("--------------------------------------------------------------------------------\n");
    printf("  Collocation attains only rle ~ 1e-2 with no way to control it; the residual\n");
    printf("  box drives the cost up to a specified accuracy, and local order 2 lowers the\n");
    printf("  residual error further for the same box.  See the header for the fixed-mesh\n");
    printf("  high-order limitation on this discontinuous solution.\n");

    return 0;
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
