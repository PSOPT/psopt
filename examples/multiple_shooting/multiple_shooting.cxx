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
             double err_est; int segments; double t_switch_gap; };

static Row solve_it(int which, const char* transcription, int segments, int steps,
                    const char* upar, int path_samples, bool costates,
                    bool flexible_segments = false, const char* integrator = "RK4",
                    bool automatic = false)
{
    problem_case = which;

    Alg algorithm; Sol solution; Prob problem;
    Row out; out.flag = -1; out.J = 0.0; out.l1_err = 0.0; out.l2_err = 0.0;
    out.u_out_of_bounds = 0.0; out.err_est = 0.0; out.segments = 0;
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

    printf("\n--------------------------------------------------------------------------------\n");
    printf("  What this transcription does not yet have: an IMPLICIT integrator. Semi-explicit\n");
    printf("  index-1 DAEs are reached by point 9 above, so what is left is stiff dynamics,\n");
    printf("  and higher-index systems where the algebraic relation has to be differentiated\n");
    printf("  more than once. A stiffly accurate ESDIRK is the shape that would take.\n");

    return 0;
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
