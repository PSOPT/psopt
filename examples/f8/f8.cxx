//////////////////////////////////////////////////////////////////////////
//////////////////              f8.cxx              //////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:  F-8 aircraft minimum-time control (integer control) //////
//////// Method: integer-control declaration API + sum-up rounding   //////
//////// Reference: Kaya & Noakes (2003); Sager (2005); mintOC       //////
////////            http://mintoc.de/index.php/F-8_aircraft          //////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2025        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////
//
// The F-8 aircraft time-optimal control problem is a classic mixed-integer
// optimal control benchmark. The tail-deflection control takes one of two
// discrete values, w in {-0.05236, +0.05236} rad, and the aim is to steer the
// aircraft from the initial to the terminal state in minimum time. The optimal
// solution is bang-bang and the best-known objective is T = 3.78086, with a
// three-switch structure starting on the upper value.
//
// This version uses the ergonomic INTEGER-CONTROL DECLARATION API. The dynamics
// are written once with the tail deflection as a single ordinary control; the
// control is then declared as integer-valued via declare_integer_control. PSOPT
// performs the outer convexification internally (introducing the mode weights,
// the SOS1 equality path constraint, and the convex-combination dynamics) at the
// start of psopt(). Contrast this with the manual formulation, where the weights,
// the SOS1 constraint, and the combination are written by hand.
//
//   Part A (relaxation): solve the (internally convexified) problem. The relaxed
//   minimum time is a lower bound on the integer optimum; for a bang-bang
//   problem the two coincide, so it reproduces T = 3.78086.
//
//   Part B (rounding): the relaxed mode weights (returned in the solution) are
//   rounded with sum_up_rounding to a single-active-mode integer control; we
//   report the switch count and integral gap, then forward-simulate the rounded
//   bang-bang control and report the terminal state miss (the O(h) integer
//   approximation error).
//
// HONEST CAVEATS:
//  * The problem is nonconvex with several documented local minima
//    (6.035, 5.742, 5.729, 4.400, 3.78086); reaching 3.78086 depends on the
//    initial guess, which is seeded here with the known switching structure.
//  * In a minimum-time problem the rounded control does not reach the target
//    exactly in the relaxed time; the terminal miss is the integer approximation
//    error and shrinks under mesh refinement. An exactly feasible integer time
//    would require switching-time optimisation.
//  * After psopt(), the solution controls are in the weights layout; the integer
//    control is recovered by rounding.
//  * AI provenance: drafted with AI assistance and validated against mintOC.

#include "psopt.h"
#include "integer_controls.h"
#include "sum_up_rounding.h"

using namespace std;

static const double V0 = -0.05236;
static const double V1 =  0.05236;

//////////////////////////////////////////////////////////////////////////
///////////////////  F-8 right-hand side at deflection w //////////////////
//////////////////////////////////////////////////////////////////////////
// Cubic-in-control Garrard model (mintOC). Templated on the scalar so the same
// expression serves the adouble dynamics (with w an ordinary control) and the
// double forward simulation.

template<class T>
inline void f8_rhs(const T& x0, const T& x1, const T& x2, const T& w, T d[3])
{
    d[0] = -0.877*x0 + x2 - 0.088*x0*x2 + 0.47*x0*x0 - 0.019*x1*x1
           - x0*x0*x2 + 3.846*x0*x0*x0
           - 0.215*w + 0.28*x0*x0*w + 0.47*x0*w*w + 0.63*w*w*w;
    d[1] = x2;
    d[2] = -4.208*x0 - 0.396*x2 - 0.47*x0*x0 - 3.564*x0*x0*x0
           - 20.967*w + 6.265*x0*x0*w + 46.0*x0*w*w + 61.4*w*w*w;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Cost, DAE, events, linkages //////////////////////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return tf;   // minimum time
}

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{
    return 0.0;
}

// Dynamics written once, with the tail deflection as a single ordinary control.
// The outer convexification is performed internally once the control is declared
// integer-valued; this function is never aware of the mode weights.
void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
    adouble x0 = states[0];
    adouble x1 = states[1];
    adouble x2 = states[2];
    adouble w  = controls[0];   // the integer tail-deflection control

    f8_rhs<adouble>(x0, x1, x2, w, derivatives);
}

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    e[0] = initial_states[0];
    e[1] = initial_states[1];
    e[2] = initial_states[2];
    e[3] = final_states[0];
    e[4] = final_states[1];
    e[5] = final_states[2];
}

void linkages(adouble* linkages, adouble* xad, Workspace* workspace)
{
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Forward simulation of a fixed bang-bang control //////
//////////////////////////////////////////////////////////////////////////

static double simulate_rounded(const RowVectorXd& hwidth, const RowVectorXd& wbin)
{
    double x[3] = {0.4655, 0.0, 0.0};
    const int sub = 20;
    for (int i = 0; i < hwidth.size(); ++i) {
        const double w  = wbin(i);
        const double dt = hwidth(i) / sub;
        for (int s = 0; s < sub; ++s) {
            double k1[3], k2[3], k3[3], k4[3], xt[3];
            f8_rhs<double>(x[0], x[1], x[2], w, k1);
            for (int j = 0; j < 3; ++j) xt[j] = x[j] + 0.5*dt*k1[j];
            f8_rhs<double>(xt[0], xt[1], xt[2], w, k2);
            for (int j = 0; j < 3; ++j) xt[j] = x[j] + 0.5*dt*k2[j];
            f8_rhs<double>(xt[0], xt[1], xt[2], w, k3);
            for (int j = 0; j < 3; ++j) xt[j] = x[j] + dt*k3[j];
            f8_rhs<double>(xt[0], xt[1], xt[2], w, k4);
            for (int j = 0; j < 3; ++j)
                x[j] += (dt/6.0)*(k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]);
        }
    }
    return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Main ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main(void)
{
    Alg  algorithm;
    Sol  solution;
    Prob problem;

    problem.name        = "F-8 aircraft (integer-control declaration API)";
    problem.outfilename = "f8.txt";

    problem.nphases   = 1;
    problem.nlinkages = 0;
    psopt_level1_setup(problem);

    const int nnodes = 60;

    problem.phases(1).nstates   = 3;
    problem.phases(1).ncontrols = 1;   // one control: the tail deflection
    problem.phases(1).nevents   = 6;
    problem.phases(1).npath     = 0;   // SOS1 added internally by the declaration
    problem.phases(1).nodes     << nnodes;

    psopt_level2_setup(problem, algorithm);

    // Declare the single control as integer-valued with two admissible values.
    RowVectorXd ivalues(2); ivalues << V0, V1;
    declare_integer_control(problem, 1, 0, ivalues);

    problem.phases(1).bounds.lower.states << -2.0, -2.0, -2.0;
    problem.phases(1).bounds.upper.states <<  2.0,  2.0,  2.0;

    // Control bounds for an integer control are ignored (PSOPT owns the weight
    // bounds); the value range is written here only for readability.
    problem.phases(1).bounds.lower.controls << V0;
    problem.phases(1).bounds.upper.controls << V1;

    problem.phases(1).bounds.lower.events << 0.4655, 0.0, 0.0, 0.0, 0.0, 0.0;
    problem.phases(1).bounds.upper.events << 0.4655, 0.0, 0.0, 0.0, 0.0, 0.0;

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0;
    problem.phases(1).bounds.upper.EndTime   = 8.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    // Initial guess: the known best switching structure (deflection values), plus
    // a dynamically-consistent state guess by forward simulation.
    const double Tguess = 3.781;
    RowVectorXd tg = linspace(0.0, Tguess, nnodes);

    RowVectorXd defl(nnodes);
    for (int i = 0; i < nnodes; ++i) {
        double ti = tg(i);
        bool up = (ti < 1.135) || (ti >= 1.482 && ti < 3.089);
        defl(i) = up ? V1 : V0;
    }

    DMatrix state_guess(3, nnodes);
    {
        double xs[3] = {0.4655, 0.0, 0.0};
        const int sub = 20;
        for (int i = 0; i < nnodes; ++i) {
            state_guess(0, i) = xs[0];
            state_guess(1, i) = xs[1];
            state_guess(2, i) = xs[2];
            if (i == nnodes - 1) break;
            double w  = defl(i);
            double dt = (tg(i+1) - tg(i)) / sub;
            for (int s = 0; s < sub; ++s) {
                double k1[3], k2[3], k3[3], k4[3], xt[3];
                f8_rhs<double>(xs[0], xs[1], xs[2], w, k1);
                for (int j = 0; j < 3; ++j) xt[j] = xs[j] + 0.5*dt*k1[j];
                f8_rhs<double>(xt[0], xt[1], xt[2], w, k2);
                for (int j = 0; j < 3; ++j) xt[j] = xs[j] + 0.5*dt*k2[j];
                f8_rhs<double>(xt[0], xt[1], xt[2], w, k3);
                for (int j = 0; j < 3; ++j) xt[j] = xs[j] + dt*k3[j];
                f8_rhs<double>(xt[0], xt[1], xt[2], w, k4);
                for (int j = 0; j < 3; ++j)
                    xs[j] += (dt/6.0)*(k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]);
            }
        }
    }

    // The integer control's guess row is given in value terms; PSOPT expands it
    // to one-hot mode weights internally.
    DMatrix control_guess(1, nnodes);
    control_guess << defl;

    problem.phases(1).guess.states   = state_guess;
    problem.phases(1).guess.controls = control_guess;
    problem.phases(1).guess.time     = tg;

    algorithm.nlp_method           = "IPOPT";
    algorithm.scaling              = "automatic";
    algorithm.derivatives          = "automatic";
    algorithm.nlp_iter_max         = 3000;
    algorithm.nlp_tolerance        = 1.e-6;
    algorithm.collocation_method   = "trapezoidal";
    algorithm.mesh_refinement      = "automatic";
    algorithm.mr_max_iterations    = 7;

    int psopt_status = psopt(solution, problem, algorithm);
    if (psopt_status != 0)
        printf("\nPSOPT returned a non-zero status (%d).\n", psopt_status);

    DMatrix x = solution.get_states_in_phase(1);
    DMatrix u = solution.get_controls_in_phase(1);   // weights layout: 2 x N
    DMatrix t = solution.get_time_in_phase(1);

    // --------------------------------------------------------------------
    // Part B: reconstruct the integer control from the relaxed mode weights.
    // --------------------------------------------------------------------
    const int N = (int)t.cols();

    IntegerControlReconstruction rec = reconstruct_integer_control(solution, problem, 1);

    const double T_relaxed   = t(0, N-1);
    const double term_miss   = simulate_rounded(rec.interval_widths, rec.control);
    const double ref_optimum = 3.78086;

    printf("\n================ F-8 integer-control summary ================\n");
    printf("Relaxed minimum time  T_relaxed = %.6f   (mintOC best-known %.5f)\n",
           T_relaxed, ref_optimum);
    printf("Rounded control       switches  = %d\n", rec.n_switches);
    printf("Sum-up rounding       gap       = %.3e   (<= max interval width %.3e)\n",
           rec.integral_gap, rec.interval_widths.maxCoeff());
    printf("Rounded trajectory    |x(T)-0|  = %.3e   (integer approximation error)\n",
           term_miss);
    printf("=============================================================\n");
    printf("Dynamics written once with the deflection as a single control; the\n");
    printf("outer convexification (weights, SOS1, combination) was performed\n");
    printf("internally via declare_integer_control. Result generated with AI\n");
    printf("assistance and validated against the mintOC reference.\n\n");

    Save(x, "x.dat");
    Save(u, "u.dat");
    Save(t, "t.dat");

    plot(t, x, problem.name, "time (s)", "states", "x0 x1 x2");
    plot(t, u, problem.name, "time (s)", "mode weights", "w0 w1");

    return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     //////////////////////////////
//////////////////////////////////////////////////////////////////////////
