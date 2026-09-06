//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   Linear tangent steering, and how accurately each    /////
////////          collocation scheme returns its costates             /////
//////// Reference: A. E. Bryson and Y.-C. Ho, "Applied Optimal       /////
////////            Control", Hemisphere, 1975, Section 2.4.          /////
//////////////////////////////////////////////////////////////////////////
////////
//////// Problem:  a point mass thrusting at constant acceleration a in
//////// a uniform field, steered by the thrust angle alpha, is to reach a
//////// given altitude and horizontal speed in minimum time:
////////
////////     xdot  = vx        vxdot = a cos(alpha)
////////     ydot  = vy        vydot = a sin(alpha)
////////
////////     x(0) = y(0) = vx(0) = vy(0) = 0
////////     y(tf) = h,  vx(tf) = vc,  vy(tf) = 0,   minimise tf
////////
//////// with a = 0.02 km/s^2, h = 408 km and vc = 7.66 km/s, which is the
//////// ascent-like parameter set. (examples/lts is the same problem in the
//////// non-dimensional data of Bryson and Ho, without the costate study.)
////////
//////// Maximum principle, with the convention H = L + lambda^T f:
////////
////////     H = l1 vx + l2 vy + l3 a cos(alpha) + l4 a sin(alpha)
////////
////////     l1dot = 0        l3dot = -l1        l1 = 0 here, so l3 is constant
////////     l2dot = 0        l4dot = -l2        l2 constant, so l4 is linear
////////
////////     dH/dalpha = 0  =>  tan(alpha) = l4/l3,
////////
//////// which is the LINEAR TANGENT law that names the problem: the tangent
//////// of the steering angle is linear in time. Imposing the boundary and
//////// transversality conditions leaves two unknowns, the initial steering
//////// angle alpha0 and the final time, fixed by y(tf) = h and vx(tf) = vc:
////////
////////     alpha0 = 0.68558798 rad        tf = 419.58841 s
////////
//////// and then the whole adjoint is known in closed form:
////////
////////     l1*(t) = 0
////////     l2*(t) = -2 sin(alpha0)/(a tf)      = -1.508927e-01
////////     l3*(t) = -cos(alpha0)/a             = -3.870235e+01
////////     l4*(t) = -(sin(alpha0)/a)(1 - 2t/tf), linear, zero at tf/2
////////
//////// WHY THIS EXAMPLE EXISTS. examples/mineng_di checks the covector
//////// mapping on a problem every scheme transcribes exactly, so that any
//////// error left is the mapping's. This example is the complement: a
//////// problem with a genuine discretization error, a free final time, and
//////// an adjoint that is constant in two components and linear in a third,
//////// where the schemes separate by five orders of magnitude. Two of them
//////// are Lobatto schemes, whose covector map is known to leave a
//////// node-to-node oscillation that a smoothing filter only partly
//////// removes; Radau and Gauss have no such difficulty. Running the same
//////// problem through all four is the quickest way to see the difference,
//////// and to see that it is a property of the scheme rather than of the
//////// accuracy of the primal solution, which is excellent in every case.
////////
//////// Usage:   lts_costates [collocation_method] [nodes]
////////
////////     collocation_method   Legendre (default), Chebyshev, Radau,
////////                          Gauss, trapezoidal or Hermite-Simpson
////////     nodes                number of collocation nodes (default 40)
////////
//////// The run writes lts_costates_<method>.dat, whose columns are the node
//////// times, the four computed costates and the four analytical ones, for
//////// plotting.
////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2026                   /////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;
using namespace PSOPT;

//////////////////////////////////////////////////////////////////////////
///////////////////  Problem data and its closed-form solution  //////////
//////////////////////////////////////////////////////////////////////////

static const double A_THRUST = 0.02;            // km/s^2
static const double H_FINAL  = 408.0;           // km
static const double VC_FINAL = 7.66;            // km/s

// The two constants that close the analytical solution, obtained by solving
// y(tf) = h and vx(tf) = vc for alpha0 and tf under the linear tangent law.
static const double ALPHA0   = 0.68558798150;   // rad
static const double TF_STAR  = 419.58840571;    // s

// The analytical costates at time t.
static void analytical_costates(double t, double lam[4])
{
    lam[0] = 0.0;
    lam[1] = -2.0*sin(ALPHA0)/(A_THRUST*TF_STAR);
    lam[2] = -cos(ALPHA0)/A_THRUST;
    lam[3] = -(sin(ALPHA0)/A_THRUST)*(1.0 - 2.0*t/TF_STAR);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return tf;                                  // minimum time
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase,
                       Workspace* workspace)
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
    adouble vx    = states[2];
    adouble vy    = states[3];
    adouble alpha = controls[0];

    derivatives[0] = vx;
    derivatives[1] = vy;
    derivatives[2] = A_THRUST*cos(alpha);
    derivatives[3] = A_THRUST*sin(alpha);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function //////////////////////////
//////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    e[0] = initial_states[0];                   // x(0)   = 0
    e[1] = initial_states[1];                   // y(0)   = 0
    e[2] = initial_states[2];                   // vx(0)  = 0
    e[3] = initial_states[3];                   // vy(0)  = 0
    e[4] = final_states[1];                     // y(tf)  = h
    e[5] = final_states[2];                     // vx(tf) = vc
    e[6] = final_states[3];                     // vy(tf) = 0
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function //////////////////
//////////////////////////////////////////////////////////////////////////

void linkages(adouble* linkages, adouble* xad, Workspace* workspace)
{
    // Single phase: no linkages
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine /////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    string method = (argc > 1) ? argv[1] : "Legendre";
    int    nnodes = (argc > 2) ? atoi(argv[2]) : 40;

    if (nnodes < 4) {
        printf("usage: %s [collocation_method] [nodes >= 4]\n", argv[0]);
        return 1;
    }

    Alg  algorithm;
    Sol  solution;
    Prob problem;

    problem.name        = "Linear tangent steering, costates (" + method + ")";
    problem.outfilename = "lts_costates_" + method + ".txt";

    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 4;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 7;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     = (RowVectorXi(1) << nnodes).finished();
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states << -100.0, -100.0, -20.0, -20.0;
    problem.phases(1).bounds.upper.states << 2000.0, 1000.0,  20.0,  20.0;

    problem.phases(1).bounds.lower.controls << -pi/2.0;
    problem.phases(1).bounds.upper.controls <<  pi/2.0;

    problem.phases(1).bounds.lower.events << 0.0, 0.0, 0.0, 0.0,
                                             H_FINAL, VC_FINAL, 0.0;
    problem.phases(1).bounds.upper.events << 0.0, 0.0, 0.0, 0.0,
                                             H_FINAL, VC_FINAL, 0.0;

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 100.0;
    problem.phases(1).bounds.upper.EndTime   = 1000.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    int ng = 10;
    MatrixXd state_guess(4, ng);
    state_guess << linspace(0.0, 1500.0, ng),
                   linspace(0.0, H_FINAL, ng),
                   linspace(0.0, VC_FINAL, ng),
                   linspace(0.0, 0.0, ng);
    problem.phases(1).guess.controls = 0.5*ones(1, ng);
    problem.phases(1).guess.states   = state_guess;
    problem.phases(1).guess.time     = linspace(0.0, 420.0, ng);

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.collocation_method = method;
    algorithm.nlp_iter_max       = 1000;
    algorithm.nlp_tolerance      = 1.e-8;

    if (psopt(solution, problem, algorithm) != 0 || solution.error_flag) {
        printf("\n lts_costates: the solve did not succeed: %s\n",
               solution.error_msg.c_str());
        return 1;
    }

    //////////////////////////////////////////////////////////////////////
    ///////  What the run found, against what is known  //////////////////
    //////////////////////////////////////////////////////////////////////

    MatrixXd t   = solution.get_time_in_phase(1);
    MatrixXd x   = solution.get_states_in_phase(1);
    MatrixXd u   = solution.get_controls_in_phase(1);
    MatrixXd lam = solution.get_dual_costates_in_phase(1);

    const long   nc = t.cols();
    const double tf = t(0, nc-1);

    double emax[4] = {0.0, 0.0, 0.0, 0.0};
    double tan_err = 0.0;                       // |tan(alpha) - l4/l3|, no closed form needed

    string datname = "lts_costates_" + method + ".dat";
    FILE* fp = fopen(datname.c_str(), "w");
    fprintf(fp, "# linear tangent steering, %s scheme, %ld nodes\n",
            method.c_str(), (long) nnodes);
    fprintf(fp, "# t  lambda1..4 (computed)  lambda1..4 (analytical)\n");

    for (long k = 0; k < nc; k++) {
        const double tk = t(0, k);
        double ana[4];
        analytical_costates(tk, ana);

        fprintf(fp, "%.10e", tk);
        for (int j = 0; j < 4; j++) fprintf(fp, "  %.10e", lam(j, k));
        for (int j = 0; j < 4; j++) fprintf(fp, "  %.10e", ana[j]);
        fprintf(fp, "\n");

        for (int j = 0; j < 4; j++) {
            const double e = fabs(lam(j, k) - ana[j]);
            if (e > emax[j]) emax[j] = e;
        }
        if (k < u.cols() && fabs(lam(2, k)) > 0.0) {
            const double e = fabs(tan(u(0, k)) - lam(3, k)/lam(2, k));
            if (e > tan_err) tan_err = e;
        }
    }
    fclose(fp);

    double ana0[4];
    analytical_costates(0.0, ana0);

    printf("\n=====================================================================\n");
    printf("  Linear tangent steering: %s, %d nodes\n", method.c_str(), nnodes);
    printf("=====================================================================\n");
    printf("  Final time\n");
    printf("    computed   tf = %.8f s\n", tf);
    printf("    analytical tf = %.8f s      error %.3e\n",
           TF_STAR, fabs(tf - TF_STAR));
    printf("    initial steering angle alpha0 = %.8f rad (analytical %.8f)\n",
           u(0,0), ALPHA0);

    printf("\n  Costates, maximum absolute error over the mesh\n");
    printf("    lambda1  %.3e     (analytical 0, constant)\n",          emax[0]);
    printf("    lambda2  %.3e     (analytical %.6e, constant)\n",       emax[1], ana0[1]);
    printf("    lambda3  %.3e     (analytical %.6e, constant)\n",       emax[2], ana0[2]);
    printf("    lambda4  %.3e     (analytical linear, %.6e at t=0)\n",  emax[3], ana0[3]);

    printf("\n  A check that needs no closed form\n");
    printf("    max |tan(alpha) - lambda4/lambda3| = %.3e\n", tan_err);
    printf("    This is the stationarity condition dH/dalpha = 0, formed from\n");
    printf("    PSOPT's own control and its own costates. It is what a user has\n");
    printf("    on a problem whose adjoint is not known.\n");

    printf("\n  Reading these numbers\n");
    printf("    The primal solution is excellent under every scheme: all of them\n");
    printf("    recover tf to a few parts in 1e9. The costates are not, and the\n");
    printf("    schemes separate by orders of magnitude. Legendre and Chebyshev\n");
    printf("    are Lobatto schemes, whose covector map leaves a node-to-node\n");
    printf("    oscillation; PSOPT smooths it, and what survives is largest at\n");
    printf("    the ends of the interval, where the nodes cluster. Radau and\n");
    printf("    Gauss have no such term and return costates as accurate as the\n");
    printf("    states.\n");
    printf("\n");
    printf("    Run all four and compare. Note also that the Lobatto errors are\n");
    printf("    not reproducible in the way an ordinary discretization error is:\n");
    printf("    scale the initial guess by a factor of two, which changes neither\n");
    printf("    the problem nor its solution, and they move by orders of\n");
    printf("    magnitude, because the automatic constraint scaling changes with\n");
    printf("    it. The Radau and Gauss errors move too, and stay small.\n");
    printf("=====================================================================\n");
    printf("\n  wrote %s\n\n", datname.c_str());

    //////////////////////////////////////////////////////////////////////
    ///////  Files and figures  //////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////

    Save(x, "x.dat");
    Save(u, "u.dat");
    Save(t, "t.dat");

    plot(t, x, problem.name + ": states", "time (s)", "states",
         "x y vx vy");
    plot(t, u, problem.name + ": control", "time (s)", "alpha (rad)", "alpha");
    plot(t, lam, problem.name + ": costates", "time (s)", "costates",
         "l1 l2 l3 l4");

    plot(t, x, problem.name + ": states", "time (s)", "states",
         "x y vx vy", "pdf", "lts_costates_states.pdf");
    plot(t, lam, problem.name + ": costates", "time (s)", "costates",
         "l1 l2 l3 l4", "pdf", "lts_costates_adjoint.pdf");

    return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     /////////////////////////////
//////////////////////////////////////////////////////////////////////////
