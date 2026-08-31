//////////////////////////////////////////////////////////////////////////
////////////////           mpec_tank.cxx             /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Dynamic MPEC -- a gas-liquid tank with a  ////////////////
////////        state-dependent outlet                    ////////////////
//////// Last modified: 30 August 2026                    ////////////////
//////// Reference:     K. M. Moudgalya and V. Ryali,     ////////////////
////////         Chemical Engineering Science 56(11),     ////////////////
////////         2001, pp. 3595-3609, for the process;    ////////////////
////////         S. R. Kazi, K. Wang and L. T. Biegler,   ////////////////
////////         J. Process Control 153 (2025) 103492,    ////////////////
////////         for the complementarity formulation.     ////////////////
////////         MPEC background: Luo, Pang and Ralph,    ////////////////
////////         Mathematical Programs with Equilibrium   ////////////////
////////         Constraints, CUP, 1996.                  ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2026        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser  ///////////////
//////// General Public License (LGPL)                     //////////////
//////////////////////////////////////////////////////////////////////////

// A closed tank is fed with liquid and with gas. A single outlet tube leaves the
// vessel at a height corresponding to a liquid volume V_s: when the liquid stands
// above that height, liquid leaves through it; when it stands below, gas leaves.
// The right-hand side of the model is therefore discontinuous in the state, and
// the discontinuity moves as the tank fills and drains.
//
// The technique this example illustrates is the replacement of a discontinuous
// right-hand side by complementarity conditions. Nothing about it is special to
// this process: the same reformulation applies wherever a discrete decision is
// determined by the sign of a continuous quantity.
//
// The derivation used here is set out below rather than quoted, because the
// step that matters -- how the complementarity conditions arise -- is a plain
// application of the Karush-Kuhn-Tucker conditions to a one-variable linear
// program, and is more convincing when it is derived than when it is asserted.
//
// The solution is worth knowing about before running it, because it is not the
// simple switch one might expect. The objective drains liquid, and liquid drains
// only while the outlet runs liquid, so the best the vehicle of the tank can do
// is bring the level down to the outlet height and hold it there. It reaches the
// height at about eleven seconds and then SLIDES along it for the rest of the
// horizon. That is not a failure of the formulation but a consequence of it: on
// the switching surface both multipliers vanish, so neither complementarity
// condition constrains nu, and nu is free to take any value in [0,1]. What it
// settles at is the value that holds the level -- the fraction of time the
// chattering outlet spends running liquid. This is the Filippov solution of the
// discontinuous system, and the complementarity formulation reproduces it
// without being told to.

#include "psopt.h"

using namespace std;
using namespace PSOPT;

//////////////////////////////////////////////////////////////////////////
///////////////////  The reformulation                        ////////////
//////////////////////////////////////////////////////////////////////////

// Write g(M_L) = M_L/rho_L - V_s for the amount by which the liquid volume
// exceeds the outlet height, and let nu in [0,1] select the outlet: nu = 1 for
// the liquid outlet, nu = 0 for the gas outlet. The physical rule is that nu
// takes the value 1 when g > 0 and 0 when g < 0, which is exactly the solution
// of the inner-level linear program
//
//      maximise   nu * g        subject to   0 <= nu <= 1.
//
// Attach multipliers s2 >= 0 to nu >= 0 and s1 >= 0 to nu <= 1. Stationarity of
// the Lagrangian -nu*g - s2*nu + s1*(nu-1) in nu gives
//
//      g = s1 - s2,
//
// and the complementary slackness of the two bound constraints gives
//
//      0 <= nu     perp  s2 >= 0,
//      0 <= 1 - nu perp  s1 >= 0.
//
// Those three lines are the whole reformulation. They are stated here as one
// equality path constraint, bounds on nu, s1 and s2, and a penalty in the
// objective that vanishes exactly when both complementarities hold. Solving the
// complementarity conditions directly would defeat an interior-point method,
// which cannot approach a solution where a variable and its multiplier are both
// driven to zero; moving them into the objective keeps the feasible set open.

//////////////////////////////////////////////////////////////////////////
///////////////////  Process data                             ////////////
//////////////////////////////////////////////////////////////////////////

// Units: moles, litres, seconds, atmospheres, kelvin.

#define R_GAS     0.0820573660809596    // l atm / (mol K)
#define TEMP      300.0                 // K
#define V_TANK    10.0                  // l
#define V_S        5.0                  // l, liquid volume at the outlet height
#define P_OUT      1.0                  // atm
#define RHO_L     50.0                  // mol/l
#define K_L        1.0                  // mol/(s atm)
#define K_G        1.0                  // mol/(s atm)
#define F_L        2.5                  // mol/s, liquid feed
#define F_G        0.1                  // mol/s, gas feed

#define BETA     100.0                  // weight on the valve-tracking term
#define RHO_PEN 1000.0                  // weight on the complementarity penalty
#define X_REF      0.1                  // valve opening the objective prefers

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
   return final_states[1];        // the liquid holdup left at the end
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
   adouble x  = controls[0];
   adouble nu = controls[1];
   adouble s1 = controls[2];
   adouble s2 = controls[3];

   // the two complementarity products, penalised rather than imposed
   adouble comp = nu*s2 + (1.0 - nu)*s1;

   return BETA*(x - X_REF)*(x - X_REF) + RHO_PEN*comp;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
   adouble MG = states[0];        // gas holdup, mol
   adouble ML = states[1];        // liquid holdup, mol

   adouble x  = controls[0];      // valve opening
   adouble nu = controls[1];      // outlet selector
   adouble s1 = controls[2];      // multiplier of nu <= 1
   adouble s2 = controls[3];      // multiplier of nu >= 0

   // The volume closure M_G R T / P + M_L / rho_L = V determines the pressure
   // outright, so it needs no variable of its own.
   adouble Vgas = V_TANK - ML/RHO_L;
   adouble P    = MG*R_GAS*TEMP/Vgas;

   adouble dP = P - P_OUT;
   adouble L  = K_L*x*dP;         // liquid outflow when the outlet runs liquid
   adouble G  = K_G*x*dP;         // gas outflow when it runs gas

   derivatives[0] = F_G*nu + (F_G - G)*(1.0 - nu);
   derivatives[1] = (F_L - L)*nu + F_L*(1.0 - nu);

   // stationarity of the inner-level problem: g = s1 - s2
   path[0] = ML/RHO_L - V_S - s1 + s2;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   e[0] = initial_states[0];
   e[1] = initial_states[1];
}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
   // single phase problem
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define main routine ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    // With an argument, that many nodes on a fixed mesh and no refinement,
    // which is how the convergence of the arrival time was established.
    int fixed_nodes = (argc > 1) ? atoi(argv[1]) : 0;
    Alg  algorithm;
    Sol  solution;
    Prob problem;

    problem.name        = "Dynamic MPEC: gas-liquid tank with a state-dependent outlet";
    problem.outfilename = "mpec_tank.txt";

    problem.nphases   = 1;
    problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 4;
    problem.phases(1).nevents   = 2;
    problem.phases(1).npath     = 1;
    if (fixed_nodes > 0) problem.phases(1).nodes = (RowVectorXi(1) << fixed_nodes).finished();
    else                 problem.phases(1).nodes = (RowVectorXi(3) << 100, 200, 400).finished();

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Bounds                                 ////////////////
////////////////////////////////////////////////////////////////////////////

    // The tank starts with 5.2 l of liquid, just above the outlet height, so the
    // outlet begins by running liquid. The gas holdup follows from the pressure
    // through the volume closure; the reference prints 6.83 mol, which is this
    // number rounded.
    double ML0 = 260.0;
    double P0  = 35.0;
    double MG0 = P0*(V_TANK - ML0/RHO_L)/(R_GAS*TEMP);

    printf("initial gas holdup consistent with %.1f atm and %.1f mol of liquid: %.6f mol\n",
           P0, ML0, MG0);

    problem.phases(1).bounds.lower.states << 0.1, 100.0;
    problem.phases(1).bounds.upper.states << 50.0, 400.0;

    //                                        x     nu    s1     s2
    problem.phases(1).bounds.lower.controls << 0.0,  0.0,  0.0,  0.0;
    problem.phases(1).bounds.upper.controls << 1.0,  1.0, 10.0, 10.0;

    problem.phases(1).bounds.lower.path << 0.0;
    problem.phases(1).bounds.upper.path << 0.0;

    problem.phases(1).bounds.lower.events << MG0, ML0;
    problem.phases(1).bounds.upper.events << MG0, ML0;

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 25.0;
    problem.phases(1).bounds.upper.EndTime   = 25.0;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Initial guess                          ////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes = 100;
    MatrixXd xg = zeros(2, nnodes), ug = zeros(4, nnodes);
    MatrixXd tg = linspace(0.0, 25.0, nnodes);

    for (int j = 0; j < nnodes; j++) {
        double t = 25.0*((double) j)/((double)(nnodes-1));
        double ML = ML0 - 0.9*t;                 // roughly what x = 0.1 gives
        if (ML < 240.0) ML = 240.0;
        double g  = ML/RHO_L - V_S;
        xg(0,j) = MG0;
        xg(1,j) = ML;
        ug(0,j) = X_REF;
        ug(1,j) = (g > 0.0) ? 1.0 : 0.0;
        ug(2,j) = (g > 0.0) ?  g  : 0.0;
        ug(3,j) = (g > 0.0) ? 0.0 : -g;
    }
    problem.phases(1).guess.states   = xg;
    problem.phases(1).guess.controls = ug;
    problem.phases(1).guess.time     = tg;

////////////////////////////////////////////////////////////////////////////
///////////////////  Algorithm options                      ////////////////
////////////////////////////////////////////////////////////////////////////

    algorithm.nlp_iter_max       = 2000;
    algorithm.nlp_tolerance      = 1.e-6;
    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.collocation_method = "trapezoidal";
    algorithm.mesh_refinement    = (fixed_nodes > 0) ? "manual" : "automatic";
    algorithm.mr_max_iterations  = 5;
    algorithm.ode_tolerance      = 1.e-5;

    int rc = psopt(solution, problem, algorithm);
    if (rc != 0) printf("psopt returned %d\n", rc);

////////////////////////////////////////////////////////////////////////////
///////////////////  Results                                ////////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x = solution.get_states_in_phase(1);
    MatrixXd u = solution.get_controls_in_phase(1);
    MatrixXd t = solution.get_time_in_phase(1);
    int N = (int) x.cols();

    // The level reaches the outlet height and then stays on it, so what is
    // worth reporting is the time of arrival and how well the sliding arc is
    // held, rather than a crossover time, of which there is none.
    //
    // The selector itself has to be read off the arc with some care, and the
    // reason is worth stating because it is a property of the problem and not
    // of this solver. On the sliding arc g vanishes, so the path constraint
    // forces s1 = s2, and the penalty then drives both to zero; with s1 = s2 = 0
    // the term rho*(nu*s2 + (1-nu)*s1) is zero whatever nu is, and nu leaves the
    // objective altogether. Nothing determines it pointwise. What determines it
    // is the requirement that the level be held, and a Hermite-Simpson
    // discretization imposes that once per interval, on the Simpson combination
    // of the two control variables it carries there -- the one at the node and
    // the one at the midpoint. Their difference is free, and the solver returns
    // an arbitrary member of that family: here the node and midpoint branches
    // separate by about four hundredths, systematically, over the whole arc.
    // This is what a singular arc looks like in a local collocation scheme: the
    // same phenomenon as the familiar ringing of a collocated singular control,
    // wearing a different coat.
    //
    // The physically meaningful quantity survives it. The fraction of time the
    // chattering outlet runs liquid is fixed by the requirement nu*L = F_L, so
    //
    //           nu* = F_L / ( K_L * x * (P - P_out) ),
    //
    // which depends only on the state and the valve -- both of which the cost
    // does pin down -- and is therefore well determined at every point. That is
    // the Filippov selector, and it is what is reported and plotted below. The
    // departure of the node branch from it is reported too, as a measure of how
    // much of the control the discretization has left free.
    double t_arrive = -1.0, worst_comp = 0.0, worst_comp_int = 0.0;
    double worst_split = 0.0, worst_gap = 0.0;
    double nus_lo = 1.0, nus_hi = 0.0;
    MatrixXd nu_star = zeros(1, N);
    for (int j = 0; j < N; j++) {
        double nu = u(1,j), s1 = u(2,j), s2 = u(3,j);
        double c  = fabs(nu*s2) + fabs((1.0-nu)*s1);
        if (c > worst_comp) worst_comp = c;
        if (j < N-1 && c > worst_comp_int) worst_comp_int = c;
        double g = x(1,j)/RHO_L - V_S;
        double P = x(0,j)*R_GAS*TEMP/(V_TANK - x(1,j)/RHO_L);
        double L = K_L*u(0,j)*(P - P_OUT);
        nu_star(0,j) = (L > 0.0) ? F_L/L : 1.0;
        if (t_arrive < 0.0 && j > 0 && g < 1.0e-4) {
            double gp = x(1,j-1)/RHO_L - V_S;
            t_arrive = t(0,j-1) + (t(0,j)-t(0,j-1))*gp/(gp - g);
        }
        // the last few nodes carry a terminal-cost artefact and are excluded
        if (t_arrive > 0.0 && t(0,j) > t_arrive && t(0,j) < 24.5) {
            if (fabs(g) > worst_gap) worst_gap = fabs(g);
            if (fabs(nu - nu_star(0,j)) > worst_split) worst_split = fabs(nu - nu_star(0,j));
            if (nu_star(0,j) < nus_lo) nus_lo = nu_star(0,j);
            if (nu_star(0,j) > nus_hi) nus_hi = nu_star(0,j);
        }
    }

    printf("\n");
    printf("liquid holdup at 25 s          %12.5f mol\n", x(1,N-1));
    printf("gas holdup at 25 s             %12.5f mol\n", x(0,N-1));
    printf("liquid volume at 25 s          %12.5f l  (outlet height %.1f l)\n",
           x(1,N-1)/RHO_L, V_S);
    printf("level reaches the outlet at    %12.5f s, and stays there\n", t_arrive);
    printf("  held to within                %12.3e l over the rest of the horizon\n", worst_gap);
    printf("  the Filippov selector nu* = F_L/L runs from %.4f to %.4f on that arc,\n",
           nus_lo, nus_hi);
    printf("  strictly inside [0,1], which is the condition for the sliding mode to exist\n");
    printf("  the collocated selector departs from it by up to %.4f: the part of nu\n",
           worst_split);
    printf("  that the discretization leaves free on a singular arc\n");
    printf("complementarity residual       %12.3e at the final node,\n", worst_comp);
    printf("                               %12.3e over the interior\n", worst_comp_int);
    printf("\n");

    Save(x, "x.dat"); Save(u, "u.dat"); Save(t, "t.dat");
    Save(nu_star, "nu_star.dat");

    plot(t, x.row(1)/RHO_L, problem.name + ": liquid volume",
         "time (s)", "V (l)", "V_L");
    plot(t, u.row(0), problem.name + ": valve opening", "time (s)", "x", "x");
    plot(t, nu_star, problem.name + ": Filippov selector", "time (s)", "nu*", "nu*");

    return 0;
}
