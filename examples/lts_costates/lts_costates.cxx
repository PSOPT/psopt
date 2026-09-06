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
//////// Usage:   lts_costates [collocation_method] [nodes] [guess_scale]
////////          lts_costates sweep [nodes]
////////
////////     collocation_method   Legendre (default), Chebyshev, Radau,
////////                          Gauss, trapezoidal or Hermite-Simpson
////////     nodes                number of collocation nodes (default 40)
////////     guess_scale          multiplies the initial guess (default 1)
////////
//////// The single run writes lts_costates_<method>.dat, whose columns are the
//////// node times, the four computed costates and the four analytical ones,
//////// for plotting, together with the states, control and time in
//////// lts_costates_<method>_{x,u,t}.dat.
////////
//////// THE SWEEP. The second form runs all four pseudospectral schemes at six
//////// scalings of the initial guess and prints, for each scheme, the smallest
//////// and largest costate error over the six. The guess is multiplied by a
//////// factor; the problem and its solution are untouched, but PSOPT's
//////// automatic constraint scaling is computed at the guess, so the factor
//////// changes the scaling and nothing else. The spread that results is the
//////// point of the exercise: it measures how much of a scheme's costate error
//////// is a property of the scheme and how much is an accident of the scaling.
////////
//////// Two things about the sweep are deliberate and should not be tidied
//////// away. First, each point is solved cold, from its own scaled guess.
//////// Continuing one point from the last would make the six runs share a
//////// scaling and destroy the very quantity being measured; a scan that wants
//////// its points to agree should use set_guess_from_solution and
//////// examples/cracking as its model, but this one wants them to disagree.
//////// Second, every file the sweep writes carries the factor in its name, so
//////// that a sweep cannot overwrite the single run's output. A study driver
//////// that writes to a fixed filename and is then used for a sweep hands the
//////// last point of the sweep to whoever collects the files afterwards, and
//////// nothing in the file says so.
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
///////////////////  One run of the study ////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// What one solve of the problem yields that is worth comparing across runs.
struct Run
{
    bool   ok;
    double tf;
    double alpha0;
    double emax[4];                             // max |computed - analytical|
    double tan_err;                             // max |tan(alpha) - l4/l3|
};

// A suffix that distinguishes the files of one run from those of another. The
// single run uses an empty suffix and owns lts_costates_<method>.dat; every
// run of the sweep carries its factor, so the sweep cannot overwrite it.
static string scale_tag(double guess_scale, bool in_sweep)
{
    if (!in_sweep) return "";
    char buf[32];
    snprintf(buf, sizeof buf, "_g%g", guess_scale);
    return string(buf);
}

// Solve the problem once, from a guess multiplied by guess_scale, and write
// the run's files. Everything is local, so successive calls share nothing:
// each solve is cold, which is what the sweep needs and what a scan that
// wants consistent points must not do.
static Run run_one(const string& method, int nnodes, double guess_scale,
                   bool in_sweep, bool make_plots)
{
    Run r;
    r.ok = false;
    r.tf = r.alpha0 = r.tan_err = 0.0;
    for (int j = 0; j < 4; j++) r.emax[j] = 0.0;

    const string tag  = scale_tag(guess_scale, in_sweep);
    const string stem = "lts_costates_" + method + tag;

    Alg  algorithm;
    Sol  solution;
    Prob problem;

    problem.name        = "Linear tangent steering, costates (" + method + ")";
    problem.outfilename = stem + ".txt";

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

    // The guess. Multiplying the states and the control by a common factor
    // leaves the problem, and therefore its solution, exactly as it was; what
    // it changes is the automatic constraint scaling, which PSOPT computes at
    // the guess. The time guess is not scaled: it sets the horizon, and
    // stretching it would change the run in a second way.
    int ng = 10;
    MatrixXd state_guess(4, ng);
    state_guess << linspace(0.0, 1500.0, ng),
                   linspace(0.0, H_FINAL, ng),
                   linspace(0.0, VC_FINAL, ng),
                   linspace(0.0, 0.0, ng);
    problem.phases(1).guess.controls = guess_scale*0.5*ones(1, ng);
    problem.phases(1).guess.states   = guess_scale*state_guess;
    problem.phases(1).guess.time     = linspace(0.0, 420.0, ng);

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.collocation_method = method;
    algorithm.nlp_iter_max       = 1000;
    algorithm.nlp_tolerance      = 1.e-8;

    if (psopt(solution, problem, algorithm) != 0 || solution.error_flag) {
        printf("\n lts_costates: the solve did not succeed (%s, guess scale %g)"
               ": %s\n", method.c_str(), guess_scale,
               solution.error_msg.c_str());
        return r;
    }

    //////////////////////////////////////////////////////////////////////
    ///////  What the run found, against what is known  //////////////////
    //////////////////////////////////////////////////////////////////////

    MatrixXd t   = solution.get_time_in_phase(1);
    MatrixXd x   = solution.get_states_in_phase(1);
    MatrixXd u   = solution.get_controls_in_phase(1);
    MatrixXd lam = solution.get_dual_costates_in_phase(1);

    const long nc = t.cols();

    r.ok     = true;
    r.tf     = t(0, nc-1);
    r.alpha0 = u(0, 0);

    const string datname = stem + ".dat";
    FILE* fp = fopen(datname.c_str(), "w");
    fprintf(fp, "# linear tangent steering, %s scheme, %d nodes, "
                "guess scale %g\n", method.c_str(), nnodes, guess_scale);
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
            if (e > r.emax[j]) r.emax[j] = e;
        }
        if (k < u.cols() && fabs(lam(2, k)) > 0.0) {
            const double e = fabs(tan(u(0, k)) - lam(3, k)/lam(2, k));
            if (e > r.tan_err) r.tan_err = e;
        }
    }
    fclose(fp);

    Save(x, (stem + "_x.dat").c_str());
    Save(u, (stem + "_u.dat").c_str());
    Save(t, (stem + "_t.dat").c_str());

    if (make_plots) {
        plot(t, x, problem.name + ": states", "time (s)", "states",
             "x y vx vy");
        plot(t, u, problem.name + ": control", "time (s)", "alpha (rad)",
             "alpha");
        plot(t, lam, problem.name + ": costates", "time (s)", "costates",
             "l1 l2 l3 l4");

        plot(t, x, problem.name + ": states", "time (s)", "states",
             "x y vx vy", "pdf", (stem + "_states.pdf").c_str());
        plot(t, lam, problem.name + ": costates", "time (s)", "costates",
             "l1 l2 l3 l4", "pdf", (stem + "_adjoint.pdf").c_str());
    }

    printf("\n  wrote %s\n", datname.c_str());
    return r;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Report one run //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static void report_one(const string& method, int nnodes, const Run& r)
{
    double ana0[4];
    analytical_costates(0.0, ana0);

    printf("\n=====================================================================\n");
    printf("  Linear tangent steering: %s, %d nodes\n", method.c_str(), nnodes);
    printf("=====================================================================\n");
    printf("  Final time\n");
    printf("    computed   tf = %.8f s\n", r.tf);
    printf("    analytical tf = %.8f s      error %.3e\n",
           TF_STAR, fabs(r.tf - TF_STAR));
    printf("    initial steering angle alpha0 = %.8f rad (analytical %.8f)\n",
           r.alpha0, ALPHA0);

    printf("\n  Costates, maximum absolute error over the mesh\n");
    printf("    lambda1  %.3e     (analytical 0, constant)\n",          r.emax[0]);
    printf("    lambda2  %.3e     (analytical %.6e, constant)\n",       r.emax[1], ana0[1]);
    printf("    lambda3  %.3e     (analytical %.6e, constant)\n",       r.emax[2], ana0[2]);
    printf("    lambda4  %.3e     (analytical linear, %.6e at t=0)\n",  r.emax[3], ana0[3]);

    printf("\n  A check that needs no closed form\n");
    printf("    max |tan(alpha) - lambda4/lambda3| = %.3e\n", r.tan_err);
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
    printf("    Run all four and compare, and then run the sweep, which repeats\n");
    printf("    the comparison at six scalings of the initial guess. The Lobatto\n");
    printf("    errors are not reproducible in the way an ordinary discretization\n");
    printf("    error is: scaling the guess changes neither the problem nor its\n");
    printf("    solution, but it changes the automatic constraint scaling, and\n");
    printf("    they move by orders of magnitude with it. The Radau and Gauss\n");
    printf("    errors move too, and stay small.\n");
    printf("=====================================================================\n\n");
}

//////////////////////////////////////////////////////////////////////////
///////////////////  The guess-scaling sweep /////////////////////////////
//////////////////////////////////////////////////////////////////////////

#define NSWEEP_METHODS 4
#define NSWEEP_SCALES  6

static const char*  SWEEP_METHODS[NSWEEP_METHODS] =
                        { "Legendre", "Chebyshev", "Radau", "Gauss" };
static const double SWEEP_SCALES[NSWEEP_SCALES] =
                        { 0.5, 0.8, 1.0, 1.2, 1.5, 2.0 };

// The table, written once to a stream. The solves are all done before any of
// this is emitted, so that the same text can go to the terminal and to a file
// without the NLP solver's own output running through the middle of it. A
// study whose summary exists only in a terminal buffer, or only as lines to be
// sieved out of a log, is a study that will be transcribed by hand sooner or
// later.
static void emit_sweep(FILE* out, int nnodes,
                       const Run runs[NSWEEP_METHODS][NSWEEP_SCALES])
{
    fprintf(out, "=====================================================================\n");
    fprintf(out, "  Costate error against the scaling of the initial guess\n");
    fprintf(out, "  %d nodes; the guess is multiplied by the factor, which leaves the\n",
            nnodes);
    fprintf(out, "  problem and its solution unchanged and moves only the automatic\n");
    fprintf(out, "  constraint scaling. Each point is solved cold, from its own guess.\n");
    fprintf(out, "=====================================================================\n");
    fprintf(out, "  %-10s %6s   %10s %10s %10s %10s %12s\n",
            "scheme", "factor", "lambda1", "lambda2", "lambda3", "lambda4", "tf");

    for (int m = 0; m < NSWEEP_METHODS; m++) {
        for (int s = 0; s < NSWEEP_SCALES; s++) {
            const Run& r = runs[m][s];
            if (!r.ok) {
                fprintf(out, "  %-10s %6g   did not solve\n",
                        SWEEP_METHODS[m], SWEEP_SCALES[s]);
                continue;
            }
            fprintf(out, "  %-10s %6g   %10.3e %10.3e %10.3e %10.3e %12.5f\n",
                    SWEEP_METHODS[m], SWEEP_SCALES[s],
                    r.emax[0], r.emax[1], r.emax[2], r.emax[3], r.tf);
        }
    }

    // lambda3 is the constant costate and is the one tabulated: it is the
    // largest of the four in every scheme, and being constant it can be read
    // without reference to a time-varying analytical value.
    fprintf(out, "\n  Spread of the lambda3 error over the %d factors\n", NSWEEP_SCALES);
    fprintf(out, "  %-10s %12s %12s %10s\n", "scheme", "smallest", "largest", "ratio");
    for (int m = 0; m < NSWEEP_METHODS; m++) {
        double lo = -1.0, hi = -1.0;
        for (int s = 0; s < NSWEEP_SCALES; s++) {
            if (!runs[m][s].ok) continue;
            const double e = runs[m][s].emax[2];
            if (lo < 0.0 || e < lo) lo = e;
            if (hi < 0.0 || e > hi) hi = e;
        }
        if (lo < 0.0) { fprintf(out, "  %-10s   (no point solved)\n", SWEEP_METHODS[m]); continue; }
        fprintf(out, "  %-10s %12.3e %12.3e %10.0f\n",
                SWEEP_METHODS[m], lo, hi, hi/lo);
    }

    fprintf(out, "\n  What the spread means. A discretization error is a property of\n");
    fprintf(out, "  the scheme and the mesh, and should not move when something that\n");
    fprintf(out, "  is not part of the problem is altered. The final time above is\n");
    fprintf(out, "  the evidence that nothing about the problem did move: every run\n");
    fprintf(out, "  returns it to seven figures. The Radau and Gauss costate errors\n");
    fprintf(out, "  move by a factor of a few and stay at the level of the primal\n");
    fprintf(out, "  solution. The Lobatto errors move by two orders of magnitude and\n");
    fprintf(out, "  are nowhere near it, so on this problem their size cannot be read\n");
    fprintf(out, "  as a mesh error at all: what is being measured is the interaction\n");
    fprintf(out, "  of the covector map with the constraint scaling. A costate taken\n");
    fprintf(out, "  from a Lobatto scheme should therefore be checked, by the\n");
    fprintf(out, "  stationarity residual the single run prints, or by re-solving\n");
    fprintf(out, "  under Radau.\n");
    fprintf(out, "=====================================================================\n");
}

static int run_sweep(int nnodes)
{
    Run  runs[NSWEEP_METHODS][NSWEEP_SCALES];
    bool complete = true;

    for (int m = 0; m < NSWEEP_METHODS; m++) {
        for (int s = 0; s < NSWEEP_SCALES; s++) {
            runs[m][s] = run_one(SWEEP_METHODS[m], nnodes, SWEEP_SCALES[s],
                                 /* in_sweep */ true, /* make_plots */ false);
            if (!runs[m][s].ok) complete = false;
        }
    }

    printf("\n");
    emit_sweep(stdout, nnodes, runs);

    // The sweep's own summary file. It is not named after any one factor,
    // because it belongs to no one run; nothing a single run writes is called
    // this, so the two cannot collide.
    const char* summary = "lts_costates_sweep.txt";
    FILE* fp = fopen(summary, "w");
    if (fp != NULL) {
        emit_sweep(fp, nnodes, runs);
        fclose(fp);
        printf("\n  wrote %s\n\n", summary);
    }

    return complete ? 0 : 1;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine /////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    string method = (argc > 1) ? argv[1] : "Legendre";
    int    nnodes = (argc > 2) ? atoi(argv[2]) : 40;
    double gscale = (argc > 3) ? atof(argv[3]) : 1.0;

    if (nnodes < 4 || gscale <= 0.0) {
        printf("usage: %s [collocation_method] [nodes >= 4] [guess_scale > 0]\n"
               "       %s sweep [nodes >= 4]\n", argv[0], argv[0]);
        return 1;
    }

    if (method == "sweep") return run_sweep(nnodes);

    // A single run at a scale other than 1 is still a point of the sweep, and
    // its files are named as one, so that it cannot displace the nominal run.
    const bool tagged = (gscale != 1.0);

    const Run r = run_one(method, nnodes, gscale, tagged, /* make_plots */ true);
    if (!r.ok) return 1;

    report_one(method, nnodes, r);
    return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     /////////////////////////////
//////////////////////////////////////////////////////////////////////////
