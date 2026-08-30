//////////////////////////////////////////////////////////////////////////
////////////////            terrain.cxx              /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Terrain following flight over a ridge line ///////////////
//////// Last modified: 30 August 2026                    ////////////////
//////// Reference:     Original. The vehicle is a         ///////////////
////////         standard linearised longitudinal model    ///////////////
////////         with an elevator actuator; the terrain    ///////////////
////////         and the problem are the author's.         ///////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2026        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser  ///////////////
//////// General Public License (LGPL)                     //////////////
//////////////////////////////////////////////////////////////////////////

// An aircraft flies at constant speed along a fixed ground track and is required
// to stay a set clearance above the terrain while flying as low as it safely can.
// The terrain is a line of four ridges whose widths differ by more than an order
// of magnitude.
//
// The example is a stress test for mesh refinement, and it is built so that the
// difficulty is physical rather than contrived. Two things make it hard.
//
// The dynamics are stiff, and stiff for a reason: the elevator actuator has a
// time constant of 0.05 s, giving a pole at -20 rad/s, while the altitude
// responds over tens of seconds. The short-period mode sits between them at
// -0.95 +/- 1.98i rad/s. The spread of time scales is better than three orders
// of magnitude and none of it was put there to be awkward.
//
// The terrain has features the aircraft cannot follow. Following a ridge of
// height H and width w at speed V needs a flight path angle of about 0.77 H/w.
// The broadest ridge here needs 0.19 rad and can be followed; the narrowest needs
// 1.28 rad and cannot, so the aircraft must climb over it and the constraint
// stops being active at the summit. Where the constraint is active and where it
// is not is therefore decided by the solution, not by the problem statement, and
// the narrowest ridge is crossed in 0.6 s of a 60 s flight. A mesh that does not
// resolve it will cut the corner and violate the clearance.
//
// The ridges are squared hyperbolic secants rather than Gaussians: smooth
// everywhere, no branching for the automatic differentiation to trip over, and
// with a shape whose flanks are straighter, which is what makes the following
// problem interesting.

#include "psopt.h"

using namespace std;
using namespace PSOPT;

//////////////////////////////////////////////////////////////////////////
///////////////////  Vehicle and terrain                      ////////////
//////////////////////////////////////////////////////////////////////////

// SI units: metres, seconds, radians.

#define V_AIR        200.0      // airspeed along the track, m/s
#define ZALPHA_V       0.7      // lift derivative over speed, 1/s
#define M_ALPHA       -4.0      // pitch stiffness, 1/s^2
#define M_Q           -1.2      // pitch damping, 1/s
#define M_DELTA       -8.0      // elevator effectiveness, 1/s^2
#define TAU_ACT        0.05     // elevator actuator time constant, s

#define T_FINAL       60.0      // s, so twelve kilometres of track
#define H_CLEAR      100.0      // required clearance above the terrain, m
#define H_START      250.0      // altitude at the start, m

#define H_SCALE      200.0      // scales the altitude term in the objective
#define W_ELEV       100.0      // weight on elevator activity

// centre, width and height of each ridge
static const double RIDGE_X[4] = { 2000.0, 5000.0, 7500.0, 9800.0 };
static const double RIDGE_W[4] = {  800.0,  400.0,  150.0,   60.0 };
static const double RIDGE_H[4] = {  200.0,  180.0,  140.0,  100.0 };

template <class T> T terrain_height(T x)
{
   T h = 0.0;
   for (int i = 0; i < 4; i++) {
      T z = (x - RIDGE_X[i])/RIDGE_W[i];
      T c = cosh(z);
      h = h + RIDGE_H[i]/(c*c);          // H sech^2(z)
   }
   return h;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
   return 0.0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
   adouble h = states[3];
   adouble u = controls[0];

   adouble x    = V_AIR*time;
   adouble gap  = h - terrain_height(x) - H_CLEAR;

   return (gap/H_SCALE)*(gap/H_SCALE) + W_ELEV*u*u;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
   adouble alpha = states[0];     // angle of attack, rad
   adouble q     = states[1];     // pitch rate, rad/s
   adouble theta = states[2];     // pitch attitude, rad
   adouble h     = states[3];     // altitude, m
   adouble delta = states[4];     // elevator deflection, rad

   adouble u     = controls[0];   // elevator command, rad

   derivatives[0] = q - ZALPHA_V*alpha;
   derivatives[1] = M_ALPHA*alpha + M_Q*q + M_DELTA*delta;
   derivatives[2] = q;
   derivatives[3] = V_AIR*(theta - alpha);
   derivatives[4] = (u - delta)/TAU_ACT;

   // clearance above the terrain
   adouble x = V_AIR*time;
   path[0] = h - terrain_height(x) - H_CLEAR;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   for (int j = 0; j < 5; j++) e[j] = initial_states[j];
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
    // With an argument, that many nodes on a fixed mesh and no refinement.
    int fixed_nodes = (argc > 1) ? atoi(argv[1]) : 0;

    Alg  algorithm;
    Sol  solution;
    Prob problem;

    problem.name        = "Terrain following flight over a ridge line";
    problem.outfilename = "terrain.txt";

    problem.nphases   = 1;
    problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 5;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 5;
    problem.phases(1).npath     = 1;
    if (fixed_nodes > 0) problem.phases(1).nodes = (RowVectorXi(1) << fixed_nodes).finished();
    else                 problem.phases(1).nodes = (RowVectorXi(3) << 80, 160, 320).finished();

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Bounds                                 ////////////////
////////////////////////////////////////////////////////////////////////////

    double amax = 0.30, qmax = 0.6, thmax = 0.45, dmax = 0.35;

    problem.phases(1).bounds.lower.states << -amax, -qmax, -thmax,   0.0, -dmax;
    problem.phases(1).bounds.upper.states <<  amax,  qmax,  thmax, 2000.0, dmax;

    problem.phases(1).bounds.lower.controls << -dmax;
    problem.phases(1).bounds.upper.controls <<  dmax;

    problem.phases(1).bounds.lower.path << 0.0;
    problem.phases(1).bounds.upper.path << 2000.0;

    // level flight at the start
    problem.phases(1).bounds.lower.events << 0.0, 0.0, 0.0, H_START, 0.0;
    problem.phases(1).bounds.upper.events << 0.0, 0.0, 0.0, H_START, 0.0;

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = T_FINAL;
    problem.phases(1).bounds.upper.EndTime   = T_FINAL;

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

    // Level flight at the height of the tallest obstacle: feasible, and far
    // enough from the answer to be an honest starting point.
    int nnodes = 80;
    MatrixXd xg = zeros(5, nnodes), ug = zeros(1, nnodes);
    MatrixXd tg = linspace(0.0, T_FINAL, nnodes);
    for (int j = 0; j < nnodes; j++) {
        xg(3,j) = H_START;
    }
    problem.phases(1).guess.states   = xg;
    problem.phases(1).guess.controls = ug;
    problem.phases(1).guess.time     = tg;

////////////////////////////////////////////////////////////////////////////
///////////////////  Algorithm options                      ////////////////
////////////////////////////////////////////////////////////////////////////

    algorithm.nlp_iter_max       = 3000;
    algorithm.nlp_tolerance      = 1.e-6;
    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.collocation_method = "Hermite-Simpson";
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

    // clearance achieved, and how much of the flight rides the constraint
    double worst = 1.0e30, ride = 0.0;
    MatrixXd clear_(1,N), terr(1,N), track(1,N);
    for (int j = 0; j < N; j++) {
        double xx = V_AIR*t(0,j);
        double te = terrain_height<double>(xx);
        double c  = x(3,j) - te - H_CLEAR;
        track(0,j) = xx/1000.0; terr(0,j) = te; clear_(0,j) = c;
        if (c < worst) worst = c;
        if (c < 5.0 && j > 0) ride += t(0,j) - t(0,j-1);
    }

    printf("\n");
    printf("least clearance above the required floor  %10.4f m\n", worst);
    printf("time spent within 5 m of the floor        %10.4f s of %.0f\n", ride, T_FINAL);
    printf("greatest elevator deflection              %10.4f deg\n",
           u.row(0).cwiseAbs().maxCoeff()*180.0/M_PI);
    printf("greatest flight path angle                %10.4f deg\n",
           (x.row(2)-x.row(0)).cwiseAbs().maxCoeff()*180.0/M_PI);
    printf("altitude range                            %10.1f to %.1f m\n",
           x.row(3).minCoeff(), x.row(3).maxCoeff());
    printf("nodes in the final mesh                   %10d\n", N);
    printf("\n");

    Save(x, "x.dat"); Save(u, "u.dat"); Save(t, "t.dat");
    Save(terr, "terrain.dat"); Save(clear_, "clearance.dat"); Save(track, "track.dat");

    plot(track, x.row(3), problem.name + ": altitude", "track (km)", "h (m)", "h");
    plot(track, clear_,   problem.name + ": clearance", "track (km)", "m", "clearance");
    plot(t, u.row(0)*180.0/M_PI, problem.name + ": elevator",
         "time (s)", "deg", "elevator");

    return 0;
}
