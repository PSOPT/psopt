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
//////// and configured with ms_integrator, which chooses the explicit  ///
//////// scheme ("RK4" or "RK8"); ms_steps_per_segment, which sets how  ///
//////// many of its steps cross a segment; ms_control_parameterisation,///
//////// which is "constant", "linear" or "quadratic"; ms_path_samples, ///
//////// which sets how many interior points of a segment the path      ///
//////// constraints are also enforced at; and ms_flexible_segments,    ///
//////// which lets the segment boundaries move. Setting                ///
//////// mesh_refinement = "automatic" turns on segment refinement,     ///
//////// governed by ms_refine_tolerance.                               ///
////////                                                               ///
//////// This example makes eight points, each with a number attached, ///
//////// and four of them are cautions rather than selling points.     ///
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

#include <cstring>
using namespace PSOPT;

//////////////////////////////////////////////////////////////////////////
///////////////////  Problem functions  //////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// 0: minimum energy, (0,0) -> (1,0) on [0,1].  J* = 6, u*(t) = 6 - 12t,
//    costates l1 = -12 and l2 = 12t - 6.
// 1: Bryson and Denham's problem, the same dynamics with x <= 1/9 and
//    x(0)=0, v(0)=1, x(1)=0, v(1)=-1.  J* = 4.
// 2: minimum time with u in [-1,2], so the single switch falls at tf/3 and tf is
//    free.  tf* = sqrt(3).
// 3: the same minimum-energy problem on an OSCILLATOR, xddot = -w^2 x + u, whose
//    optimal control is a sinusoid: smooth, and not a polynomial, so no
//    parameterisation here represents it exactly and the three are measured
//    against the same unreachable answer.  J* comes from the controllability
//    Gramian in closed form.
static int problem_case = 0;

static const double TF_SQRT3 = 1.7320508075688772;
static const double W_OSC    = 10.0;

// J* = (1/2) x_T' W^-1 x_T for xddot = -w^2 x + u on [0,1] from rest to (1,0).
static double oscillator_optimum(void)
{
    const double w = W_OSC, T = 1.0;
    const double s2 = T/2.0 - sin(2*w*T)/(4*w);
    const double c2 = T/2.0 + sin(2*w*T)/(4*w);
    const double sc = (1.0 - cos(2*w*T))/(4*w);
    const double W11 = s2/(w*w), W12 = sc/w, W22 = c2;
    return 0.5*W22/(W11*W22 - W12*W12);
}

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{ return ( problem_case == 2 ) ? tf : (adouble) 0.0; }

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase, Workspace* workspace)
{ return ( problem_case == 2 ) ? (adouble) 0.0 : 0.5*controls[0]*controls[0]; }

void dae(adouble* derivatives, adouble* path, adouble* states, adouble* controls,
         adouble* parameters, adouble& time, adouble* xad, int iphase,
         Workspace* workspace)
{
    derivatives[0] = states[1];
    derivatives[1] = ( problem_case == 3 ) ? ( -W_OSC*W_OSC*states[0] + controls[0] )
                                           : controls[0];
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

struct Row { int flag; double J; double l1_err; double l2_err; double u_out_of_bounds;
             double err_est; int segments; double t_switch_gap; double cpu; };

static Row solve_it(int which, const char* transcription, int segments, int steps,
                    const char* upar, int path_samples, bool costates,
                    bool flexible_segments = false, const char* integrator = "RK4",
                    bool automatic = false)
{
    problem_case = which;

    Alg algorithm; Sol solution; Prob problem;
    Row out; out.flag = -1; out.J = 0.0; out.l1_err = 0.0; out.l2_err = 0.0;
    out.u_out_of_bounds = 0.0; out.err_est = 0.0; out.segments = 0; out.cpu = 0.0;
    out.t_switch_gap = -1.0;

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

    if ( which == 3 ) {
        problem.phases(1).bounds.lower.states   << -50.0, -500.0;
        problem.phases(1).bounds.upper.states   <<  50.0,  500.0;
    }
    else {
        problem.phases(1).bounds.lower.states   << -5.0, -5.0;
        problem.phases(1).bounds.upper.states   <<  5.0,  5.0;
    }
    problem.phases(1).bounds.lower.controls(0) = ( which == 2 ) ? -1.0
                                               : ( which == 3 ) ? -2000.0 : -30.0;
    problem.phases(1).bounds.upper.controls(0) = ( which == 2 ) ?  2.0
                                               : ( which == 3 ) ?  2000.0 :  30.0;

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
    problem.phases(1).bounds.lower.EndTime   = ( which == 2 ) ? 0.5 : 1.0;
    problem.phases(1).bounds.upper.EndTime   = ( which == 2 ) ? 8.0 : 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(2, nodes);
    if ( which == 1 )      problem.phases(1).guess.states.row(1) = linspace( 1.0, -1.0, nodes);
    else if ( which != 3 ) problem.phases(1).guess.states.row(0) = linspace( 0.0,  1.0, nodes);
    problem.phases(1).guess.controls = zeros(1, nodes);
    problem.phases(1).guess.time     = linspace(0.0, ( which == 2 ) ? 1.73 : 1.0, nodes);

    algorithm.nlp_method            = "IPOPT";
    algorithm.scaling               = "automatic";
    algorithm.derivatives           = "automatic";
    algorithm.nlp_iter_max          = 2000;
    algorithm.nlp_tolerance         = ( which == 3 ) ? 1.0e-12 : 1.0e-10;
    algorithm.print_level           = 0;
    algorithm.mesh_refinement       = automatic ? "automatic" : "manual";
    algorithm.mr_max_iterations     = 7;
    algorithm.ms_refine_tolerance   = 1.0e-3;
    algorithm.collocation_method    = "Hermite-Simpson";
    algorithm.transcription_method  = transcription;
    algorithm.ms_steps_per_segment  = steps;
    algorithm.ms_control_parameterisation = upar;
    algorithm.ms_path_samples             = path_samples;
    algorithm.ms_flexible_segments        = flexible_segments;
    algorithm.ms_integrator               = integrator;

    out.flag = psopt(solution, problem, algorithm);
    if (out.flag != 0) return out;

    out.J = solution.cost;
    {
        DMatrix E = solution.get_relative_local_error_in_phase(1);
        for (int q = 0; q < E.size(); q++) out.err_est = fmax(out.err_est, fabs(E(q)));
    }
    // Summed over the mesh iterations, which for a manual single-count run is one.
    for (int q = 0; q < 64; q++) {
        if ( solution.mesh_stats[q].nnodes <= 0 ) break;
        out.cpu += solution.mesh_stats[q].CPU_time;
    }
    {
        DMatrix T = solution.get_time_in_phase(1);
        out.segments = (int) T.cols() - 1;
        if ( which == 2 ) {
            out.J = T(0, T.cols()-1);      // the final time is the answer here
            // how close did any boundary come to the switch, which is at tf/3?
            const double sw = out.J/3.0;
            out.t_switch_gap = 1.0e30;
            for (int q = 0; q < T.cols(); q++)
                out.t_switch_gap = fmin( out.t_switch_gap, fabs(T(0,q) - sw) );
        }
    }

    // How far outside its own bounds does the control the integrator was handed go? The
    // constant and linear forms cannot leave the box their values lie in; the parabola can.
    // The interleaved arrays are the whole control history, node and midpoint together, and
    // they are what get_hs_controls_in_phase reports under this parameterisation too.
    {
        DMatrix U  = solution.get_controls_in_phase(1);
        DMatrix Uh = solution.get_hs_controls_in_phase(1);
        const double ulo = problem.phases(1).bounds.lower.controls(0);
        const double uup = problem.phases(1).bounds.upper.controls(0);
        const int M = (int) U.cols() - 1;
        if ( Uh.cols() == 2*M + 1 ) {
            for (int k = 0; k < M; k++)
                for (int q = 0; q <= 100; q++) {
                    const double x = ((double) q)/100.0;
                    const double uu =  2.0*(x-0.5)*(x-1.0)*Uh(0,2*k)
                                     - 4.0*x*(x-1.0)      *Uh(0,2*k+1)
                                     + 2.0*x*(x-0.5)      *Uh(0,2*k+2);
                    out.u_out_of_bounds = fmax(out.u_out_of_bounds,
                                               fmax(uu - uup, ulo - uu));
                }
        }
    }

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
/////////  A semi-explicit index-1 DAE, and the two ways to pose it  //////
//////////////////////////////////////////////////////////////////////////
//
//   xdot1 = x2,   xdot2 = u - z,   0 = z^3 + z - x1,
//   min (1/2) int_0^1 u^2 dt,   (0,0) -> (1,0).
//
// Index 1 because dg/dz = 3z^2 + 1 never vanishes; and z has no closed form in
// terms of x1, so it cannot be eliminated by hand and the example is not
// secretly an ODE.
//
// dae_form selects which of three problems is being written:
//   0  z is an ordinary control, g = 0 is an equality PATH CONSTRAINT
//   1  INDEX-REDUCED: z is a third state with zdot = x2/(3z^2+1), g = 0 at t0
//   2  the same user code as 0, with nalgebraic = 1 declared
// 0 and 2 are the SAME model. What differs is what PSOPT is told about it.
static int dae_form = 0;

adouble dae_endpoint(adouble*, adouble*, adouble*, adouble&, adouble&, adouble*,
                     int, Workspace*) { return (adouble) 0.0; }
adouble dae_integrand(adouble*, adouble* c, adouble*, adouble&, adouble*, int, Workspace*)
{ return 0.5*c[0]*c[0]; }

void dae_dynamics(adouble* d, adouble* path, adouble* s, adouble* c, adouble*, adouble&,
                  adouble*, int, Workspace*)
{
    if ( dae_form == 1 ) {
        d[0] = s[1];
        d[1] = c[0] - s[2];
        d[2] = s[1]/(3.0*s[2]*s[2] + 1.0);
    }
    else {
        d[0] = s[1];
        d[1] = c[0] - c[1];
        path[0] = c[1]*c[1]*c[1] + c[1] - s[0];
    }
}

void dae_events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&,
                adouble*, int, Workspace*)
{
    e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; e[3]=f[1];
    if ( dae_form == 1 ) e[4] = i[2]*i[2]*i[2] + i[2] - i[0];
}

struct DaeRow { int flag; double J; double residual; };

//////////////////////////////////////////////////////////////////////////
/////////  A problem whose integrator error is LOCALISED  /////////////////
//////////////////////////////////////////////////////////////////////////
//
//   xdot1 = x2,   xdot2 = u + A exp(-((t-1/2)/sigma)^2/2),
//   min (1/2) int_0^1 u^2 dt,   (0,0) -> (1,0).
//
// Perfectly smooth, but its SCALE is sigma, so the integrator's error lives in
// the two or three segments covering the bump and is negligible in the rest. A
// uniform step count has to be set by the worst segment and is then paid for by
// all of them.
static const double BUMP_SIG = 0.02, BUMP_AMP = 40.0;

adouble bump_endpoint(adouble*, adouble*, adouble*, adouble&, adouble&, adouble*,
                      int, Workspace*) { return (adouble) 0.0; }
adouble bump_integrand(adouble*, adouble* c, adouble*, adouble&, adouble*, int, Workspace*)
{ return 0.5*c[0]*c[0]; }
void bump_dynamics(adouble* d, adouble*, adouble* s, adouble* c, adouble*, adouble& t,
                   adouble*, int, Workspace*)
{
    adouble z = (t - 0.5)/BUMP_SIG;
    d[0] = s[1];
    d[1] = c[0] + BUMP_AMP*exp(-0.5*z*z);
}
void bump_events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&,
                 adouble*, int, Workspace*)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; e[3]=f[1]; }

struct BumpRow { int flag; double J; double eps; int iters; };

//////////////////////////////////////////////////////////////////////////
/////////  A STIFF problem, with and without an algebraic relation  ///////
//////////////////////////////////////////////////////////////////////////
//
//  stiff_case 0:  xdot1 = -lam x1 + u,   xdot2 = x1               (stiff, linear)
//  stiff_case 1:  the same with lam = 1                           (not stiff)
//  stiff_case 2:  xdot1 = -lam(x1 - 0.3 sin x2) + u,  xdot2 = x1^2 + x2/2   (nonlinear)
//  stiff_case 3:  xdot1 = x2,  xdot2 = u - z - lam(x2 - 1),  0 = z^3 + z - x1
//
// min (1/2) int_0^1 u^2 with x(0) = 0 and x2(1) fixed.
static double STIFF_LAM = 1000.0;
static int    stiff_case = 0;

adouble stiff_endpoint(adouble*, adouble*, adouble*, adouble&, adouble&, adouble*,
                       int, Workspace*) { return (adouble) 0.0; }
adouble stiff_integrand(adouble*, adouble* c, adouble*, adouble&, adouble*, int, Workspace*)
{ return 0.5*c[0]*c[0]; }
void stiff_dynamics(adouble* d, adouble* path, adouble* s, adouble* c, adouble*, adouble&,
                    adouble*, int, Workspace*)
{
    if ( stiff_case == 2 ) {
        d[0] = -STIFF_LAM*(s[0] - 0.3*sin(s[1])) + c[0];
        d[1] = s[0]*s[0] + 0.5*s[1];
    }
    else if ( stiff_case == 3 ) {
        d[0] = s[1];
        d[1] = c[0] - c[1] - STIFF_LAM*(s[1] - 1.0);
        path[0] = c[1]*c[1]*c[1] + c[1] - s[0];
    }
    else {
        d[0] = -STIFF_LAM*s[0] + c[0];
        d[1] = s[0];
    }
}
void stiff_events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&,
                  adouble*, int, Workspace*)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[1]; }

struct StiffRow { int flag; double J; double gmax; };

static StiffRow solve_stiff(int segments, int steps, const char* integrator, int mimp = 4)
{
    Alg algorithm; Sol solution; Prob problem;
    StiffRow out; out.flag = -1; out.J = 0.0; out.gmax = -1.0;
    const int nodes = segments + 1;
    const int nc    = ( stiff_case == 3 ) ? 2 : 1;

    problem.name        = "Multiple shooting, stiff";
    problem.outfilename = "multiple_shooting_stiff.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = nc;
    problem.phases(1).nevents   = 3;
    problem.phases(1).npath     = ( stiff_case == 3 ) ? 1 : 0;
    if ( stiff_case == 3 ) problem.phases(1).nalgebraic = 1;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states << -100.0, -100.0;
    problem.phases(1).bounds.upper.states <<  100.0,  100.0;
    if ( stiff_case == 3 ) {
        problem.phases(1).bounds.lower.controls << -1.0e5, -20.0;
        problem.phases(1).bounds.upper.controls <<  1.0e5,  20.0;
        problem.phases(1).bounds.lower.path(0) = 0.0;
        problem.phases(1).bounds.upper.path(0) = 0.0;
    }
    else {
        problem.phases(1).bounds.lower.controls(0) = -1.0e5;
        problem.phases(1).bounds.upper.controls(0) =  1.0e5;
    }
    const double xT = ( stiff_case == 3 ) ? 3.0 : 1.0;
    problem.phases(1).bounds.lower.events << 0.0, 0.0, xT;
    problem.phases(1).bounds.upper.events << 0.0, 0.0, xT;
    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0; problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &stiff_integrand;
    problem.endpoint_cost  = &stiff_endpoint;
    problem.dae            = &stiff_dynamics;
    problem.events         = &stiff_events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states        = zeros(2, nodes);
    problem.phases(1).guess.states.row(1) = linspace(0.0, xT, nodes);
    problem.phases(1).guess.controls      = zeros(nc, nodes);
    problem.phases(1).guess.time          = linspace(0.0, 1.0, nodes);

    algorithm.nlp_method            = "IPOPT";
    algorithm.scaling               = "automatic";
    algorithm.derivatives           = "automatic";
    algorithm.nlp_iter_max          = 3000;
    algorithm.nlp_tolerance         = 1.0e-10;
    algorithm.print_level           = 0;
    algorithm.mesh_refinement       = "manual";
    algorithm.collocation_method    = "Hermite-Simpson";
    algorithm.transcription_method  = "multiple-shooting";
    algorithm.ms_control_parameterisation = "constant";
    algorithm.ms_steps_per_segment  = steps;
    algorithm.ms_integrator         = integrator;
    algorithm.ms_implicit_iterations = mimp;

    out.flag = psopt(solution, problem, algorithm);
    if ( out.flag != 0 ) return out;
    out.J = solution.cost;
    if ( stiff_case == 3 ) {
        const MatrixXd u = solution.get_controls_in_phase(1);
        const MatrixXd x = solution.get_states_in_phase(1);
        out.gmax = 0.0;
        for (int q = 0; q < (int) u.cols(); q++) {
            const double z = u(1,q), x1 = x(0,q);
            out.gmax = fmax( out.gmax, fabs(z*z*z + z - x1) );
        }
    }
    return out;
}


static BumpRow solve_bump(int segments, int steps, bool adaptive, const char* integrator,
                          int print_level = 0)
{
    Alg algorithm; Sol solution; Prob problem;
    BumpRow out; out.flag = -1; out.J = 0.0; out.eps = -1.0; out.iters = 0;
    const int max_iter = 8;
    const int nodes = segments + 1;

    problem.name        = "Multiple shooting, localised forcing";
    problem.outfilename = "multiple_shooting_steps.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -50.0, -50.0;
    problem.phases(1).bounds.upper.states   <<  50.0,  50.0;
    problem.phases(1).bounds.lower.controls(0) = -200.0;
    problem.phases(1).bounds.upper.controls(0) =  200.0;
    problem.phases(1).bounds.lower.events   << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.upper.events   << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0; problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &bump_integrand;
    problem.endpoint_cost  = &bump_endpoint;
    problem.dae            = &bump_dynamics;
    problem.events         = &bump_events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states        = zeros(2, nodes);
    problem.phases(1).guess.states.row(0) = linspace(0.0, 1.0, nodes);
    problem.phases(1).guess.controls      = zeros(1, nodes);
    problem.phases(1).guess.time          = linspace(0.0, 1.0, nodes);

    algorithm.nlp_method            = "IPOPT";
    algorithm.scaling               = "automatic";
    algorithm.derivatives           = "automatic";
    algorithm.nlp_iter_max          = 3000;
    algorithm.nlp_tolerance         = 1.0e-10;
    algorithm.print_level           = print_level;
    algorithm.collocation_method    = "Hermite-Simpson";
    algorithm.transcription_method  = "multiple-shooting";
    algorithm.ms_control_parameterisation = "constant";
    algorithm.ms_steps_per_segment  = steps;
    algorithm.ms_integrator         = integrator;
    algorithm.ode_tolerance         = 1.0e-8;
    algorithm.ms_refine_tolerance   = 1.0e9;   // segments held fixed: only steps under test
    if ( adaptive ) {
        algorithm.mesh_refinement   = "automatic";
        algorithm.mr_max_iterations = max_iter;
        algorithm.ms_adaptive_steps = true;
    }
    else {
        algorithm.mesh_refinement   = "manual";
    }

    out.flag = psopt(solution, problem, algorithm);
    if ( out.flag != 0 ) return out;
    out.J = solution.cost;
    const int nrows = adaptive ? max_iter : 1;   // mesh_stats is sized by the iteration count
    int last = 0;
    for (int q = 0; q < nrows; q++) {
        if ( solution.mesh_stats[q].nnodes <= 0 ) break;
        last = q;
    }
    out.iters = last + 1;
    out.eps   = solution.mesh_stats[last].epsilon_max;
    return out;
}


// z solving z^3 + z = x, outside any tape.
static double dae_z_of(double x)
{
    double z = x;
    for (int it = 0; it < 200; it++) {
        const double r = z*z*z + z - x, d = 3.0*z*z + 1.0, dz = r/d;
        z -= dz;
        if ( fabs(dz) < 1.0e-16 ) break;
    }
    return z;
}

static DaeRow solve_dae(int form, const char* transcription, int segments, int steps,
                        const char* integrator, int malg = 4)
{
    dae_form = form;
    Alg algorithm; Sol solution; Prob problem;
    DaeRow out; out.flag = -1; out.J = 0.0; out.residual = -1.0;

    const bool colloc = ( strcmp(transcription, "collocation") == 0 );
    const int nodes = colloc ? 41 : segments + 1;
    const int ns = ( form == 1 ) ? 3 : 2;
    const int nc = ( form == 1 ) ? 1 : 2;
    const int ne = ( form == 1 ) ? 5 : 4;

    problem.name        = "Multiple shooting, index-1 DAE";
    problem.outfilename = "multiple_shooting_dae.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = ns;
    problem.phases(1).ncontrols = nc;
    problem.phases(1).nevents   = ne;
    problem.phases(1).npath     = ( form == 1 ) ? 0 : 1;
    if ( form == 2 ) problem.phases(1).nalgebraic = 1;   // <<< the whole of the declaration
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    if ( form == 1 ) {
        problem.phases(1).bounds.lower.states   << -10.0, -10.0, -10.0;
        problem.phases(1).bounds.upper.states   <<  10.0,  10.0,  10.0;
        problem.phases(1).bounds.lower.controls(0) = -100.0;
        problem.phases(1).bounds.upper.controls(0) =  100.0;
        problem.phases(1).bounds.lower.events   << 0.0, 0.0, 1.0, 0.0, 0.0;
        problem.phases(1).bounds.upper.events   << 0.0, 0.0, 1.0, 0.0, 0.0;
    }
    else {
        problem.phases(1).bounds.lower.states   << -10.0, -10.0;
        problem.phases(1).bounds.upper.states   <<  10.0,  10.0;
        problem.phases(1).bounds.lower.controls << -100.0, -10.0;
        problem.phases(1).bounds.upper.controls <<  100.0,  10.0;
        problem.phases(1).bounds.lower.path(0)  = 0.0;      // an EQUALITY: the algebraic equation
        problem.phases(1).bounds.upper.path(0)  = 0.0;
        problem.phases(1).bounds.lower.events   << 0.0, 0.0, 1.0, 0.0;
        problem.phases(1).bounds.upper.events   << 0.0, 0.0, 1.0, 0.0;
    }
    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0; problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &dae_integrand;
    problem.endpoint_cost  = &dae_endpoint;
    problem.dae            = &dae_dynamics;
    problem.events         = &dae_events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states        = zeros(ns, nodes);
    problem.phases(1).guess.states.row(0) = linspace(0.0, 1.0, nodes);
    problem.phases(1).guess.controls      = zeros(nc, nodes);
    problem.phases(1).guess.time          = linspace(0.0, 1.0, nodes);

    algorithm.nlp_method            = "IPOPT";
    algorithm.scaling               = "automatic";
    algorithm.derivatives           = "automatic";
    algorithm.nlp_iter_max          = 3000;
    algorithm.nlp_tolerance         = 1.0e-11;
    algorithm.print_level           = 0;
    algorithm.mesh_refinement       = "manual";
    algorithm.collocation_method    = "Hermite-Simpson";
    algorithm.transcription_method  = transcription;
    algorithm.ms_steps_per_segment  = steps;
    algorithm.ms_integrator         = integrator;
    algorithm.ms_control_parameterisation = "constant";
    algorithm.ms_algebraic_iterations     = malg;

    out.flag = psopt(solution, problem, algorithm);
    if ( out.flag != 0 ) return out;
    out.J = solution.cost;
    if ( form == 1 || colloc ) return out;

    // The worst |z^3 + z - x1| along the trajectory, reconstructed far more finely
    // than the transcription integrated it. For form 0 the algebraic variable is a
    // control held across the segment while x1 moves under it; for form 2 it is
    // solved wherever it is needed, so the reconstruction solves it too -- which is
    // exactly what the two formulations respectively claim.
    const MatrixXd t = solution.get_time_in_phase(1);
    const MatrixXd u = solution.get_controls_in_phase(1);
    const MatrixXd x = solution.get_states_in_phase(1);
    const int M = (int) t.cols() - 1;
    const int NSUB = 400;
    double worst = 0.0, x1 = x(0,0), x2 = x(1,0);
    for (int k = 0; k < M; k++) {
        const double a = t(0,k), b = t(0,k+1), h = (b-a)/NSUB;
        const double uk = u(0,k), zk = u(1,k);
        for (int q = 0; q < NSUB; q++) {
            const double zz = ( form == 0 ) ? zk : dae_z_of(x1);
            worst = fmax( worst, fabs(zz*zz*zz + zz - x1) );
            auto f = [&](double p, double r, double& dp, double& dr) {
                dp = r; dr = uk - ( (form == 0) ? zk : dae_z_of(p) ); };
            double k1a,k1b,k2a,k2b,k3a,k3b,k4a,k4b;
            f(x1,x2,k1a,k1b);                    f(x1+h/2*k1a, x2+h/2*k1b, k2a,k2b);
            f(x1+h/2*k2a, x2+h/2*k2b, k3a,k3b);  f(x1+h*k3a,   x2+h*k3b,   k4a,k4b);
            x1 += h/6*(k1a+2*k2a+2*k3a+k4a);
            x2 += h/6*(k1b+2*k2b+2*k3b+k4b);
        }
    }
    out.residual = worst;
    return out;
}

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

    printf("\n     More generally, the control parameterisation is what caps the accuracy of\n");
    printf("     the answer, because the dynamics are integrated to whatever\n");
    printf("     ms_steps_per_segment buys and that is independent of the mesh. A control\n");
    printf("     error of O(h^p) gives a cost error of O(h^2p), the cost being stationary and\n");
    printf("     quadratic at the optimum, and the three forms are p = 1, 2 and 3. On an\n");
    printf("     oscillator whose optimal control is a sinusoid -- smooth, and reachable by\n");
    printf("     none of them -- with a fine segment integrator so that only the\n");
    printf("     parameterisation is being measured:\n\n");
    printf("        M      constant u          linear u            quadratic u\n");
    {
        const double Jo = oscillator_optimum();
        double pa = 0.0, pb = 0.0, pc = 0.0;
        const int osegs[4] = { 10, 20, 40, 80 };
        for (int q = 0; q < 4; q++) {
            Row a = solve_it(3, "multiple-shooting", osegs[q], 200, "constant",  0, false);
            Row b = solve_it(3, "multiple-shooting", osegs[q], 200, "linear",    0, false);
            Row c = solve_it(3, "multiple-shooting", osegs[q], 200, "quadratic", 0, false);
            const double ea = fabs(a.J - Jo)/Jo;
            const double eb = fabs(b.J - Jo)/Jo;
            const double ec = fabs(c.J - Jo)/Jo;
            printf("      %3d    %.3e", osegs[q], ea);
            if (pa > 0.0) printf(" /%5.1f", pa/ea); else printf("       ");
            printf("   %.3e", eb);
            if (pb > 0.0) printf(" /%5.1f", pb/eb); else printf("       ");
            printf("   %.3e", ec);
            if (pc > 0.0) printf(" /%5.1f", pc/ec); else printf("       ");
            printf("\n");
            pa = ea; pb = eb; pc = ec;
        }
    }
    printf("\n     The ratio after each doubling is the order: 4 = 2^2, 16 = 2^4 and 64 = 2^6.\n");
    printf("\n     The two error sources are SEPARABLE, and only one of them binds. Holding the\n");
    printf("     segment count and the control form fixed and varying only the integrator, the\n");
    printf("     dynamics can be driven to round-off while the cost does not move at all:\n\n");
    printf("        steps/segment   |J-J*|/J*    reported error   CPU\n");
    {
        const double Jo = oscillator_optimum();
        for (int st : { 2, 4, 10, 40, 160 }) {
            Row r = solve_it(3, "multiple-shooting", 40, st, "linear", 0, false);
            printf("        %13d   %.3e    %.3e        %.2f s\n", st,
                   fabs(r.J - Jo)/Jo, r.err_est, r.cpu);
        }
    }
    printf("\n     Eighty times the integration work buys nothing, because at forty segments with\n");
    printf("     a ramped control the answer is limited by the control and not by the dynamics.\n");
    printf("     That is the single most useful thing to know about this transcription: the\n");
    printf("     integrator's error is bought off at linear cost in CPU and NONE in decision\n");
    printf("     variables, and what is left is the control parameterisation's.\n");
    printf("     The quadratic form carries one extra control variable per segment, in the\n");
    printf("     slot Hermite-Simpson uses for its midpoint control, and reports it through\n");
    printf("     solution.get_hs_controls_in_phase. Reading get_controls_in_phase alone there\n");
    printf("     gives two thirds of the control variables and none of the curvature.\n");

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

    printf("\n  5. The segment boundaries can move, and on a problem with a switch that is\n");
    printf("     worth six orders of magnitude. With u in [-1,2] the single switch falls at\n");
    printf("     tf/3 and tf is free, so a UNIFORM partition has a boundary on the switch\n");
    printf("     exactly when the segment count is divisible by three -- and nothing else\n");
    printf("     about the problem changes between M = 9 and M = 10:\n\n");
    printf("        M      uniform partition   flexible (ms_flexible_segments = true)\n");
    const int msegs[6] = { 5, 6, 7, 9, 10, 20 };
    for (int q = 0; q < 6; q++) {
        Row a = solve_it(2, "multiple-shooting", msegs[q], 10, "constant", 0, false, false);
        Row b = solve_it(2, "multiple-shooting", msegs[q], 10, "constant", 0, false, true);
        printf("      %3d      %.3e%s        %.3e\n", msegs[q],
               fabs(a.J - TF_SQRT3)/TF_SQRT3, (msegs[q] % 3 == 0) ? "  *" : "   ",
               fabs(b.J - TF_SQRT3)/TF_SQRT3);
    }
    printf("      (* the segment count is divisible by three, so a boundary lands on the\n");
    printf("       switch by arithmetic rather than by design)\n");
    printf("\n     The flexible column does not depend on M at all, which is the point. What\n");
    printf("     it costs is a few times the CPU of the uniform solve, and at M = 5 it is\n");
    printf("     both faster and six orders of magnitude better than the uniform partition\n");
    printf("     at M = 40.\n");
    printf("\n     algorithm.ms_min_segment_fraction is the floor on a segment width, and its\n");
    printf("     default of 0.8 is not a safety device -- it is what makes the problem well\n");
    printf("     posed. A moving partition is a free-knot approximation problem, and those\n");
    printf("     are degenerate: moving a boundary and adjusting the controls either side of\n");
    printf("     it compensates, so the objective is flat along a manifold and there is no\n");
    printf("     isolated solution to converge to. A tight floor removes most of that\n");
    printf("     manifold. At a floor of 0.5 or below these solves stop converging.\n");
    printf("\n     For the same reason, do NOT turn this on for a problem whose optimal\n");
    printf("     control is smooth. There is then no corner for a boundary to sit on, the\n");
    printf("     degeneracy is complete, and the solve does not converge at any floor --\n");
    printf("     measured, on the problem of point 1 above with the dynamics made\n");
    printf("     oscillatory. This facility is for solutions with corners.\n");

    printf("\n  6. The parabola can leave the control bounds, and it does so worst exactly\n");
    printf("     where the flexible partition is most useful. A quadratic through three\n");
    printf("     values inside the box need not stay inside it -- it overshoots by a quarter\n");
    printf("     of the second difference -- and the sharpest second difference a solution\n");
    printf("     can present is a jump, which is what a moving boundary is there to sit on.\n");
    printf("     On the minimum-time problem above, with u in [-1,2]:\n\n");
    printf("        parameterisation   partition    tf rel. error   u outside [-1,2] by\n");
    {
        const char* forms[3] = { "constant", "linear", "quadratic" };
        for (int f = 0; f < 3; f++)
            for (int g = 0; g < 2; g++) {
                Row r = solve_it(2, "multiple-shooting", 10, 10, forms[f], 0, false, g == 1);
                printf("        %-18s %-12s %.3e       %.2e\n",
                       forms[f], (g == 1) ? "flexible" : "uniform",
                       fabs(r.J - TF_SQRT3)/TF_SQRT3, r.u_out_of_bounds);
            }
    }
    printf("\n     The control the segment integrator was handed reached 2.375 in the last\n");
    printf("     row, nineteen per cent above its own upper bound, and ms_path_samples does\n");
    printf("     not help: that samples the PATH constraints, and this is a variable bound.\n");
    printf("     Note also that the two continuous forms are far worse than the constant one\n");
    printf("     on this problem whatever the partition, because a continuous control cannot\n");
    printf("     represent a jump. Use \"constant\" where the optimal control has corners or\n");
    printf("     rides its bounds, and the higher forms where it is smooth and interior.\n");

    printf("\n  7. The segment integrator has a choice of explicit scheme, and what the choice\n");
    printf("     buys is accuracy per RIGHT-HAND-SIDE EVALUATION rather than accuracy for its\n");
    printf("     own sake. The tape is stages x steps x segments long, and both the memory it\n");
    printf("     takes and the time to evaluate the constraints are linear in that, so four\n");
    printf("     extra orders of convergence for 2.75 times the stages is a large trade in the\n");
    printf("     right direction. On the oscillator of point 2, with a held control:\n\n");
    printf("        scheme   steps   RHS evals/segment   discretisation error\n");
    {
        const int s4[4] = { 2, 4, 8, 16 };
        const int s8[3] = { 1, 2, 3 };
        for (int q = 0; q < 4; q++) {
            Row r = solve_it(3, "multiple-shooting", 10, s4[q], "constant", 0, false, false, "RK4");
            printf("        RK4     %5d   %17d   %.3e\n", s4[q], 4*s4[q], r.err_est);
        }
        for (int q = 0; q < 3; q++) {
            Row r = solve_it(3, "multiple-shooting", 10, s8[q], "constant", 0, false, false, "RK8");
            printf("        RK8     %5d   %17d   %.3e\n", s8[q], 11*s8[q], r.err_est);
        }
    }
    printf("\n     Twenty-two evaluations of RK8 are twenty times more accurate than thirty-two\n");
    printf("     of RK4, and thirty-three are thirty-four times more accurate than sixty-four.\n");
    printf("     The reported figure is the error of the trajectory that was returned: a\n");
    printf("     Richardson difference scaled by 1/(1 - 2^-p), which is why the scheme's order\n");
    printf("     p has to be known and not assumed.\n");

    printf("\n     What a better integrator does NOT do is make the ANSWER better once the\n");
    printf("     integrator is no longer the binding error, which on most problems it is not.\n");
    printf("     The cost of the oscillator problem against its closed form, at twenty\n");
    printf("     segments:\n\n");
    printf("        control       steps    RK4            RK8\n");
    {
        const double Jo = oscillator_optimum();
        for (const char* form : { "constant", "quadratic" })
            for (int st : { 2, 10 }) {
                Row a = solve_it(3, "multiple-shooting", 20, st, form, 0, false, false, "RK4");
                Row b = solve_it(3, "multiple-shooting", 20, st, form, 0, false, false, "RK8");
                printf("        %-12s %5d    %.6e   %.6e\n", form, st,
                       fabs(a.J-Jo)/Jo, fabs(b.J-Jo)/Jo);
            }
    }
    printf("\n     Under a held control the two rows at ten steps agree to five figures and the\n");
    printf("     answer is 2.3 per cent from the optimum either way: what limits it is the\n");
    printf("     control, and no scheme can lift that. Under the quadratic control the two\n");
    printf("     differ by three orders of magnitude at two steps, because there the\n");
    printf("     integrator WAS the binding error -- and RK8 reaches in two steps, at\n");
    printf("     twenty-two evaluations, the floor RK4 needs ten steps and forty to reach.\n");
    printf("     Raise the order of the integrator to stop the dynamics being the limit;\n");
    printf("     raise the order of the control to move the limit itself.\n");

    printf("\n  8. The segment count can be chosen automatically, and the question it answers\n");
    printf("     is not the one the other mesh-refinement drivers answer. On a collocation\n");
    printf("     mesh, more nodes means a better approximation of the DYNAMICS; here the\n");
    printf("     dynamics are integrated to whatever ms_steps_per_segment and ms_integrator\n");
    printf("     buy, however many segments there are. What the segment count controls is the\n");
    printf("     resolution of the CONTROL PARAMETERISATION and the coverage of the PATH\n");
    printf("     CONSTRAINTS, so that is what the indicator measures, and\n");
    printf("     algorithm.ms_refine_tolerance is the tolerance on it -- not ode_tolerance,\n");
    printf("     which still reports the integrator's error and is still fixed by steps.\n");
    printf("\n     On the minimum-time problem above, from uniform partitions that do not\n");
    printf("     resolve the switch (M not divisible by three):\n\n");
    printf("        start M   uniform tf rel.err   refined M   refined tf rel.err   boundary\n");
    for (int m : { 7, 10, 13 }) {
        Row a = solve_it(2, "multiple-shooting", m, 10, "constant", 0, false, false, "RK4", false);
        Row b = solve_it(2, "multiple-shooting", m, 10, "constant", 0, false, false, "RK4", true);
        printf("        %5d     %.3e           %5d       %.3e            %.1e\n", m,
               fabs(a.J - TF_SQRT3)/TF_SQRT3, b.segments,
               fabs(b.J - TF_SQRT3)/TF_SQRT3, b.t_switch_gap);
    }
    printf("\n     The last column is how far the nearest segment boundary ended up from the\n");
    printf("     switch. The refinement finds it, and the answer stops depending on where\n");
    printf("     the partition started -- which is what ms_flexible_segments buys too, by a\n");
    printf("     different route and without needing a floor to stay well posed.\n");
    printf("\n     The estimator is built so that a corner ALREADY sitting on a boundary is\n");
    printf("     not flagged, which is not automatic: an estimator written for smooth\n");
    printf("     solutions sees the jump and refines the same place for ever. The departure\n");
    printf("     of the control representation from a richer one is formed twice, from a\n");
    printf("     window extended to the left and one extended to the right, and the SMALLER\n");
    printf("     is taken -- a corner on a boundary spoils exactly one of the two, a corner\n");
    printf("     inside a segment spoils both. From a partition that already resolves the\n");
    printf("     switch:\n\n");
    {
        Row c = solve_it(2, "multiple-shooting", 12, 10, "constant", 0, false, false, "RK4", true);
        printf("        started at 12 segments (divisible by three), ended at %d,"
               " tf rel.err %.2e\n", c.segments, fabs(c.J - TF_SQRT3)/TF_SQRT3);
    }
    printf("\n     And on the path-constraint leak of point 3, which is the other thing the\n");
    printf("     segment count controls. J* = 4 and no feasible trajectory costs less, so\n");
    printf("     4 - J is the leak:\n\n");
    printf("        mode                      M      J            leak\n");
    {
        Row a = solve_it(1, "multiple-shooting", 10, 10, "linear", 0, false, false, "RK4", false);
        Row b = solve_it(1, "multiple-shooting", 10, 10, "linear", 0, false, false, "RK4", true);
        Row c = solve_it(1, "multiple-shooting", 10, 10, "linear", 2, false, false, "RK4", true);
        printf("        fixed,      0 samples  %3d    %.7f    %+.2e\n", a.segments, a.J, 4.0-a.J);
        printf("        automatic,  0 samples  %3d    %.7f    %+.2e\n", b.segments, b.J, 4.0-b.J);
        printf("        automatic,  2 samples  %3d    %.7f    %+.2e\n", c.segments, c.J, 4.0-c.J);
    }
    printf("\n     A positive number in the last column is the constraint leaking. Refinement\n");
    printf("     alone takes the leak from parts in a hundred to parts in a million, because\n");
    printf("     narrower segments leave less room between the points where the constraint\n");
    printf("     is imposed; adding two interior samples puts the answer ABOVE the optimum,\n");
    printf("     which is where a restricted control parameterisation has to leave it. The\n");
    printf("     two are complementary and neither replaces the other: ms_path_samples says\n");
    printf("     WHERE inside a segment the constraint is enforced, and the segment count\n");
    printf("     says how far apart those places can be.\n");

    printf("\n  9. A semi-explicit index-1 DAE can be propagated, and it takes one number to\n");
    printf("     say so. The problem below is\n\n");
    printf("        xdot1 = x2,   xdot2 = u - z,   0 = z^3 + z - x1,   min (1/2) int u^2,\n\n");
    printf("     index 1 because dg/dz = 3z^2 + 1 never vanishes, and with no closed form for\n");
    printf("     z in terms of x1, so it cannot be eliminated by hand.\n");
    {
        const DaeRow ref = solve_dae(0, "collocation", 0, 0, "RK4");
        const DaeRow a   = solve_dae(0, "multiple-shooting", 10, 20, "RK4");
        const DaeRow b   = solve_dae(1, "multiple-shooting", 10, 20, "RK4");
        const DaeRow c   = solve_dae(2, "multiple-shooting", 10, 20, "RK4");
        printf("\n     Collocation on 41 nodes gives the reference, J = %.9f. At ten segments\n",
               ref.J);
        printf("     and twenty RK4 steps:\n\n");
        printf("        how the algebraic relation is posed        J             worst |g|\n");
        printf("        z a control, g = 0 a path constraint       %.9f   %.2e\n",
               a.J, a.residual);
        printf("        index-reduced (z a state, g = 0 at t0)     %.9f   (an invariant)\n",
               b.J);
        printf("        phases(1).nalgebraic = 1                   %.9f   %.2e\n",
               c.J, c.residual);
    }
    printf("\n     The first row is what the problem looks like without the declaration: g = 0\n");
    printf("     is imposed where the trajectory is a decision variable, which is the segment\n");
    printf("     boundaries, and between them z is HELD while x1 moves under it. The relation\n");
    printf("     drifts at first order in the segment width and the cost is wrong in its\n");
    printf("     first figure. That is not a bug to be tuned away; it is what the formulation\n");
    printf("     says.\n");
    printf("\n     The last row solves g = 0 for z at every stage of every step, which is the\n");
    printf("     half-explicit scheme. The residual is the inner iteration's and not a\n");
    printf("     quantity that accumulates, and the answer agrees with index reduction -- an\n");
    printf("     entirely separate route to the same class, sharing no code with it -- to ten\n");
    printf("     figures. The convention is a count, like every other size: the LAST\n");
    printf("     nalgebraic controls are the algebraic variables and the FIRST nalgebraic\n");
    printf("     path constraints are their equations, which must be equalities. The user\n");
    printf("     code for the first and last rows above is character for character the same.\n");
    printf("\n     It keeps the order of ms_integrator, which is a theorem and not a hope: for\n");
    printf("     index 1 the algebraic relation defines z = G(x) locally, so a half-explicit\n");
    printf("     method IS the explicit method applied to the reduced ordinary system. At ten\n");
    printf("     segments, against the same formulation at sixty RK8 steps:\n\n");
    printf("        scheme   steps    |J - J_fine|    ratio\n");
    {
        const DaeRow fine = solve_dae(2, "multiple-shooting", 10, 60, "RK8");
        double prev = -1.0;
        for (int st : { 2, 4, 8, 16 }) {
            const DaeRow r = solve_dae(2, "multiple-shooting", 10, st, "RK4");
            const double e = fabs(r.J - fine.J);
            printf("        RK4     %5d    %.3e", st, e);
            if ( prev > 0.0 ) printf("       %5.1f", prev/e);
            printf("\n");
            prev = e;
        }
        const DaeRow r8 = solve_dae(2, "multiple-shooting", 10, 2, "RK8");
        printf("        RK8     %5d    %.3e\n", 2, fabs(r8.J - fine.J));
    }
    printf("\n     Sixteen is 2^4, so RK4 is still fourth order on a DAE; and two steps of RK8\n");
    printf("     are already at the floor, which is what an eighth-order table is for.\n");
    printf("\n     algorithm.ms_algebraic_iterations is how many iterations each stage solve\n");
    printf("     spends, and its default of 4 is derived rather than tuned. The iteration is\n");
    printf("     Broyden's, so nothing in it has to be differentiated and the nested\n");
    printf("     automatic differentiation a DAE capability is usually said to need never\n");
    printf("     arises; the count is FIXED and unrolled, because a loop whose length depends\n");
    printf("     on the values could not be taped, and would make the constraint function\n");
    printf("     non-smooth in the decision variables. Warm-started from the previous stage,\n");
    printf("     whose algebraic variables differ by O(h), the secant recurrence gives\n");
    printf("     exponents 1, 2, 3, 5, 8 -- so m iterations support a scheme of order 2, 3, 5\n");
    printf("     and 8. Measured, at ten segments and four RK8 steps:\n\n");
    printf("        iterations   |J - J_fine|\n");
    {
        const DaeRow fine = solve_dae(2, "multiple-shooting", 10, 60, "RK8");
        for (int m : { 1, 2, 3, 4, 5 }) {
            const DaeRow r = solve_dae(2, "multiple-shooting", 10, 4, "RK8", m);
            printf("        %10d   %.3e\n", m, fabs(r.J - fine.J));
        }
    }
    printf("\n     The fifth iteration buys nothing because there is nothing left to buy.\n");
    printf("\n     What this does NOT do is make a stiff problem tractable. The differential\n");
    printf("     part of the scheme is explicit and inherits its stability restriction\n");
    printf("     exactly: this adds a class of problem, not a stability region.\n");

    printf("\n 10. The step count can be chosen automatically too, and per SEGMENT. It is a\n");
    printf("     different question from the segment count and it has its own driver, because\n");
    printf("     the two chase the separable error sources this transcription is built around:\n");
    printf("     segments control the control parameterisation and the path coverage, steps\n");
    printf("     control the integrator, and neither can fix the other's error. Set\n");
    printf("     algorithm.ms_adaptive_steps = true and ms_steps_per_segment stops being a\n");
    printf("     number to guess and becomes a starting point; ode_tolerance becomes the thing\n");
    printf("     that is met.\n");
    printf("\n     The adaptivity is BETWEEN solves, and that is not a limitation of the\n");
    printf("     implementation. A step count that varied with the decision variables would\n");
    printf("     make the constraint function non-smooth in them -- a step-acceptance test\n");
    printf("     flipping as the iterate moves changes the discrete map the matching condition\n");
    printf("     is written on -- so Newton would be given derivatives that do not describe\n");
    printf("     its own residual. That is Bock's reason for freezing the discretisation, and\n");
    printf("     it is why this needs mesh_refinement = \"automatic\": there has to be another\n");
    printf("     solve for a new step count to be used in.\n");
    printf("\n     On a problem whose integrator error is LOCALISED --\n\n");
    printf("        xdot1 = x2,  xdot2 = u + %.0f exp(-((t-1/2)/%.2f)^2/2)\n\n", BUMP_AMP, BUMP_SIG);
    printf("     -- the forcing is smooth but its scale is %.2f, so the error lives in the two\n", BUMP_SIG);
    printf("     or three of twenty segments that cover it. ode_tolerance = 1e-08:\n\n");
    printf("        uniform steps   reported error   stage evaluations per solve\n");
    for (int st : { 5, 10, 20, 40 }) {
        const BumpRow r = solve_bump(20, st, false, "RK4");
        printf("        %13d   %.3e        %d\n", st, r.eps, 20*st*4);
    }
    {
        const BumpRow a = solve_bump(20, 7, true, "RK4");
        printf("\n     A uniform count has to be 20 to get under the tolerance, which is 1600\n");
        printf("     stage evaluations on every segment whether it needs them or not. Started\n");
        printf("     at 7 and adapted, the same problem converges in %d solves to a reported\n",
               a.iters);
        printf("     error of %.3e with a table running from 2 steps to 15 -- 87 steps over\n", a.eps);
        printf("     twenty segments, 348 stage evaluations. That is a factor of 4.6 in tape\n");
        printf("     length for the same accuracy, and a factor of 3.4 against the best a\n");
        printf("     single per-phase count could do, since a per-phase count is set by the\n");
        printf("     worst segment and 15 x 20 is 1200.\n");
    }
    printf("\n     PSOPT prints the table's min, max and total after each adaptation. A FLAT\n");
    printf("     table is the useful negative result: it says the problem did not need this\n");
    printf("     and a uniform count would have done as well.\n");
    {
        const BumpRow r8 = solve_bump(20, 2, true, "RK8");
        printf("\n     And it composes with the scheme rather than duplicating it. The driver\n");
        printf("     asks how many steps of whatever table is in force are needed, so RK8 does\n");
        printf("     not want a different tolerance -- it wants fewer steps, and the driver\n");
        printf("     finds out how many: started at 2, it converges in %d solves to %.3e on a\n",
               r8.iters, r8.eps);
        printf("     table of 24 steps over twenty segments, 264 stage evaluations.\n");
    }
    printf("\n     algorithm.ms_max_steps_per_segment is the ceiling, and reaching it is worth\n");
    printf("     reading as a diagnosis rather than as a limit: an explicit scheme that cannot\n");
    printf("     resolve a segment in two hundred steps is usually meeting STIFFNESS, which no\n");
    printf("     step count fixes cheaply and which this transcription does not serve.\n");

    printf("\n 11. And for STIFF dynamics the segment integrator has two implicit schemes.\n");
    printf("     ms_integrator = \"TRBDF2\" and \"ESDIRK3\" are ESDIRKs: explicit first stage,\n");
    printf("     one repeated diagonal, stiffly accurate and L-stable, of classical order 2\n");
    printf("     and 3. This is a different KIND of thing from RK4 and RK8, not two more\n");
    printf("     tables, and it is worth being clear about what it is for.\n");
    printf("\n     An explicit scheme does not become inaccurate on a stiff problem, it\n");
    printf("     becomes UNBOUNDED. RK4's stability region reaches |lam h| = 2.78, so on\n\n");
    printf("        xdot1 = -%.0f x1 + u,   xdot2 = x1,   min (1/2) int_0^1 u^2\n\n", STIFF_LAM);
    printf("     with ten segments it needs at least %d steps per segment before the\n",
           (int) ceil(STIFF_LAM/10.0/2.78));
    printf("     propagation is even bounded -- and below that it does not return a poor\n");
    printf("     answer, it returns a meaningless one:\n\n");
    printf("        scheme    steps   J\n");
    stiff_case = 0; STIFF_LAM = 1000.0;
    for (int st : { 8, 16, 36, 72 }) {
        const StiffRow r = solve_stiff(10, st, "RK4");
        if ( r.flag != 0 ) printf("        %-8s %6d   (the solve failed)\n", "RK4", st);
        else               printf("        %-8s %6d   %.6f\n", "RK4", st, r.J);
    }
    for (const char* sch : { "TRBDF2", "ESDIRK3" }) {
        const StiffRow r = solve_stiff(10, 1, sch);
        printf("        %-8s %6d   %.6f\n", sch, 1, r.J);
    }
    printf("\n     One step per segment, for a scheme with no stability limit at all.\n");
    printf("\n     The classical order can only be measured on a NON-STIFF problem, which is\n");
    printf("     a trap rather than an inconvenience: in the stiff limit an L-stable scheme\n");
    printf("     resolves the quasi-steady state almost exactly, and on the problem above\n");
    printf("     TR-BDF2 -- a second-order scheme -- reads as order eight. At lam = 1:\n\n");
    printf("        scheme    steps   |J - J_fine|   ratio\n");
    stiff_case = 1; STIFF_LAM = 1.0;
    for (const char* sch : { "TRBDF2", "ESDIRK3" }) {
        const StiffRow fine = solve_stiff(10, 64, sch);
        double prev = -1.0;
        for (int st : { 2, 4, 8, 16 }) {
            const StiffRow r = solve_stiff(10, st, sch);
            const double e = fabs(r.J - fine.J);
            printf("        %-8s %6d   %.3e", sch, st, e);
            if ( prev > 0.0 ) printf("      %6.2f", prev/e);
            printf("\n");
            prev = e;
        }
    }
    printf("\n     4 = 2^2 and 8 = 2^3, so the two tables are the orders they claim. They were\n");
    printf("     DERIVED rather than transcribed -- gamma from its own defining cubic, then\n");
    printf("     the remaining coefficients from C(2) and B(1..3) -- and checked before use\n");
    printf("     against the row sums, the order conditions, stiff accuracy, the stability\n");
    printf("     function at minus infinity and on the imaginary axis, and the expansion of\n");
    printf("     R(z) against the exponential series.\n");
    printf("\n     algorithm.ms_implicit_iterations is the fixed, unrolled number of\n");
    printf("     modified-Newton iterations each stage system gets; the default of 4 is\n");
    printf("     derived and not tuned. W = I - h gamma J is formed by finite differences\n");
    printf("     once per STEP and reused by every stage -- which is what one repeated\n");
    printf("     diagonal is for -- so nothing in the iteration has to be differentiated and\n");
    printf("     no automatic differentiation is nested. Warm-started from the previous\n");
    printf("     stage, modified Newton contracts by O(h) an iteration, so m of them leave\n");
    printf("     O(h^(m+1)). Measured on a NONLINEAR stiff problem, since modified Newton\n");
    printf("     solves a LINEAR stage system exactly in one iteration and a linear problem\n");
    printf("     would make every count look sufficient:\n\n");
    printf("        iterations   relative change from m = 8\n");
    stiff_case = 2; STIFF_LAM = 1000.0;
    {
        const StiffRow ref = solve_stiff(10, 8, "ESDIRK3", 8);
        for (int m : { 1, 2, 3, 4, 5 }) {
            const StiffRow r = solve_stiff(10, 8, "ESDIRK3", m);
            printf("        %10d   %.3e\n", m, fabs(r.J - ref.J)/fabs(ref.J));
        }
    }
    printf("\n     And the combination that needs points 9 and 11 at once -- STIFF, and a\n");
    printf("     semi-explicit index-1 DAE:\n\n");
    printf("        xdot1 = x2,  xdot2 = u - z - %.0f(x2 - 1),  0 = z^3 + z - x1\n\n", STIFF_LAM);
    printf("        scheme    steps   J                worst |g| at the nodes\n");
    stiff_case = 3;
    for (const char* sch : { "RK4", "ESDIRK3" }) {
        for (int st : { 4, 20, 100 }) {
            const StiffRow r = solve_stiff(10, st, sch);
            if ( r.flag != 0 ) { printf("        %-8s %6d   (the solve failed)\n", sch, st); continue; }
            printf("        %-8s %6d   %.6f     %.3e\n", sch, st, r.J, r.gmax);
        }
    }
    printf("\n     ESDIRK3 at four steps is already within a part in ten thousand, where RK4 is\n");
    printf("     not even stable until a hundred -- and its unstable rows are not reproducible\n");
    printf("     between builds, because a meaningless propagation makes a nonconvex problem\n");
    printf("     land wherever it lands. That is the reading: below its stability limit an\n");
    printf("     explicit scheme neither converges nor respects the algebraic relation it is\n");
    printf("     being solved against, and how badly is not a number worth quoting.\n");
    printf("\n     The implicit column holds |g| at round-off at every step count, which is\n");
    printf("     what STIFF ACCURACY is for: the last stage IS the step, so the constraint\n");
    printf("     holds at the step end, which is where the matching conditions live.\n");
    printf("\n     Two cautions. An implicit step costs several times an explicit one, so on a\n");
    printf("     non-stiff problem these are the wrong choice -- lower order than RK8 and\n");
    printf("     dearer than RK4. And in the stiff limit a Runge-Kutta method converges below\n");
    printf("     its classical order (order reduction), so the reported discretisation error,\n");
    printf("     which scales a Richardson difference by that order, is optimistic there.\n");

    printf("\n--------------------------------------------------------------------------------\n");
    printf("  What is left: higher-index DAEs, where the algebraic relation has to be\n");
    printf("  differentiated more than once and the reduction needs stabilising. Everything\n");
    printf("  else on the accuracy study's improvement list has been built.\n");

    return 0;
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
