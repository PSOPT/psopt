//////////////////////////////////////////////////////////////////////////
////////////////          vega_ascent.cxx            /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Maximum payload ascent of a four-stage    ////////////////
////////        launch vehicle to a polar orbit           ////////////////
//////// Last modified: 29 August 2026                    ////////////////
//////// Reference:     B. Benedikter, A. Zavoli,         ////////////////
////////         G. Colasurdo, S. Pizzurro and            ////////////////
////////         E. Cavallini, Convex optimization of     ////////////////
////////         launch vehicle ascent trajectory with    ////////////////
////////         heat-flux and splash-down constraints,   ////////////////
////////         J. Spacecraft and Rockets 59(3), 2022,   ////////////////
////////         pp. 900-915, DOI 10.2514/1.A35194.       ////////////////
////////         The vehicle data are theirs; the phase   ////////////////
////////         structure here is simpler than theirs.   ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2026        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser  ///////////////
//////// General Public License (LGPL)                     ///////////////
//////////////////////////////////////////////////////////////////////////

// A four-stage launch vehicle of the Vega class lifts off from the equator on
// the meridian of the Guiana Space Centre and delivers the greatest possible
// payload into a 700 km circular polar orbit. Three solid stages burn in
// succession for fixed durations; the liquid fourth stage burns for a duration
// the optimiser chooses. Dry masses are jettisoned at the stage boundaries and
// the fairing goes with the second stage.
//
// The problem exercises what a multi-phase ascent always exercises: linkages
// between phases, discontinuities in mass across them, a free terminal time and
// terminal constraints on the orbit rather than on the state directly.
//
// Two features of the vehicle model are worth noticing. The vacuum thrust and
// the mass flow of each solid stage vary linearly through the burn, in the ratio
// that keeps the specific impulse constant, which is how a solid motor with a
// regressive grain actually behaves. And the delivered thrust is the vacuum
// thrust less the ambient pressure acting over the nozzle exit area, so the
// first stage gains some 300 kN as it climbs out of the atmosphere.
//
// The phase structure is simpler than the reference's, which uses thirteen
// phases to carry a heat-flux-limited fairing jettison and a splash-down
// constraint on the spent third stage. Neither is modelled here.

#include "psopt.h"
#include "us76.h"

using namespace std;
using namespace PSOPT;

//////////////////////////////////////////////////////////////////////////
///////////////////  Vehicle and environment data  ///////////////////////
//////////////////////////////////////////////////////////////////////////

// SI throughout: metres, seconds, kilograms, newtons.

struct Constants {
  MatrixXd* omega_matrix;
  double mu, Re, g0;
  double rho0, H, p0;
  double cd, sref;
  // per stage: vacuum thrust at start and end, mass flow at start and end,
  // nozzle exit area, burn time, propellant mass, dry mass
  double Tv0[4], Tv1[4], md0[4], md1[4], Ae[4], tb[4], mprop[4], mdry[4];
  double t_start[4];
  double m_fairing;
};

typedef struct Constants Constants_;

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
   if (iphase == 5) return -final_states[6];    // maximise the delivered mass
   else             return 0.0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase, Workspace* workspace)
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
   Constants_& C = *( (Constants_ *) workspace->problem->user_data );
   int j, s = (iphase <= 3) ? iphase-1 : 3;
   bool coasting = (iphase == 4);

   adouble r[3]; r[0]=states[0]; r[1]=states[1]; r[2]=states[2];
   adouble v[3]; v[0]=states[3]; v[1]=states[4]; v[2]=states[5];
   adouble m = states[6];

   adouble rad = sqrt( dot(r, r, 3) );
   adouble altitude = rad - C.Re;

   // Earth-relative velocity
   MatrixXd& W = *C.omega_matrix;
   adouble vrel[3];
   for (j = 0; j < 3; j++)
      vrel[j] = v[j] - W(j,0)*r[0] - W(j,1)*r[1] - W(j,2)*r[2];
   adouble speedrel = sqrt( dot(vrel, vrel, 3) );

   // The atmosphere is the US Standard of 1976, which is what the reference
   // uses. It replaces an exponential fit with a 7200 m scale height that was
   // here before: that fit is within about a quarter of the standard below
   // 70 km, but it is three times too dense at 110 km and four orders of
   // magnitude too thin at 200 km, and the drag it gives through the first
   // stage -- twenty six per cent light at 10 km, where the dynamic pressure
   // peaks -- is large enough to move the payload this problem is maximising.
   adouble rho, pamb, tamb;
   us76<adouble>(altitude, rho, pamb, tamb);

   adouble bc = (rho/(2.0*m))*C.sref*C.cd;
   adouble Drag[3];
   for (j = 0; j < 3; j++) Drag[j] = -vrel[j]*bc*speedrel;

   // thrust and mass flow, linear in the fraction of the burn elapsed
   for (j = 0; j < 3; j++) derivatives[j] = v[j];

   if (coasting) {
      for (j = 0; j < 3; j++)
         derivatives[3+j] = -C.mu*r[j]/(rad*rad*rad) + Drag[j];
      derivatives[6] = 0.0;
      return;
   }

   adouble tau = (s < 3) ? (time - C.t_start[s])/C.tb[s] : 0.0;  // constant for the liquid stage
   adouble Tvac  = C.Tv0[s] + (C.Tv1[s] - C.Tv0[s])*tau;
   adouble mdot  = C.md0[s] + (C.md1[s] - C.md0[s])*tau;
   adouble Thrust = Tvac - pamb*C.Ae[s];       // ambient pressure over the nozzle

   adouble Tovm = Thrust/m;
   adouble* u = controls;

   for (j = 0; j < 3; j++)
      derivatives[3+j] = -C.mu*r[j]/(rad*rad*rad) + Tovm*u[j] + Drag[j];
   derivatives[6] = -mdot;

   path[0] = dot(u, u, 3);          // the thrust direction is a unit vector
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Terminal orbit quantities               /////////////
//////////////////////////////////////////////////////////////////////////

// The target orbit is circular and polar, and both of those are places where
// the classical elements misbehave: the argument of periapsis is undefined at
// zero eccentricity, and computing the inclination through an arccosine is
// needlessly rough. So only the three quantities the terminal constraint needs
// are formed, and each is formed in a way that stays smooth at the target:
// the semi-major axis from vis-viva, the SQUARE of the eccentricity from the
// eccentricity vector, and the cosine of the inclination, which is zero for a
// polar orbit.

static void terminal_orbit(adouble* rv, adouble* vv, double mu,
                           adouble& a, adouble& e2, adouble& cosi)
{
   adouble r  = sqrt( dot(rv, rv, 3) );
   adouble v2 = dot(vv, vv, 3);
   adouble rdotv = dot(rv, vv, 3);

   a = 1.0/( 2.0/r - v2/mu );

   adouble ev[3];
   for (int j = 0; j < 3; j++)
      ev[j] = ( (v2 - mu/r)*rv[j] - rdotv*vv[j] )/mu;
   e2 = dot(ev, ev, 3);

   adouble hv[3];
   cross(rv, vv, hv);
   cosi = hv[2]/sqrt( dot(hv, hv, 3) );
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   Constants_& C = *( (Constants_ *) workspace->problem->user_data );
   int j;

   if (iphase == 1) {
      for (j = 0; j < 7; j++) e[j] = initial_states[j];
   }

   if (iphase == 5) {
      adouble rv[3]; rv[0]=final_states[0]; rv[1]=final_states[1]; rv[2]=final_states[2];
      adouble vv[3]; vv[0]=final_states[3]; vv[1]=final_states[4]; vv[2]=final_states[5];
      adouble a, e2, cosi;
      terminal_orbit( rv, vv, C.mu, a, e2, cosi );
      e[0] = a;          // semi-major axis
      e[1] = e2;         // eccentricity squared
      e[2] = cosi;       // cosine of the inclination: zero for a polar orbit
   }
}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
    Constants_& C = *( (Constants_ *) workspace->problem->user_data );
    int index = 0;

    // auto_link makes the state and the time continuous across the boundary;
    // the mass link is then corrected by whatever is thrown away there.
    auto_link(linkages, &index, xad, 1, 2, workspace);
    linkages[index-2] -= C.mdry[0];
    auto_link(linkages, &index, xad, 2, 3, workspace);
    linkages[index-2] -= C.mdry[1] + C.m_fairing;
    auto_link(linkages, &index, xad, 3, 4, workspace);
    linkages[index-2] -= C.mdry[2];
    auto_link(linkages, &index, xad, 4, 5, workspace);
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define main routine ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{
    Alg  algorithm;
    Sol  solution;
    Prob problem;

    problem.name        = "Maximum payload Vega-class ascent to a polar orbit";
    problem.outfilename = "vega_ascent.txt";

    Constants_ CONSTANTS;
    problem.user_data = (void*) &CONSTANTS;

    problem.nphases   = 5;
    problem.nlinkages = 32;
    psopt_level1_setup(problem);

    for (int i = 1; i <= 5; i++) {
        problem.phases(i).nstates   = 7;
        problem.phases(i).ncontrols = (i == 4) ? 0 : 3;   // phase 4 coasts
        problem.phases(i).npath     = (i == 4) ? 0 : 1;
        problem.phases(i).nevents   = (i == 1) ? 7 : ((i == 5) ? 3 : 0);
        problem.phases(i).nodes     = (RowVectorXi(2) << 15, 20).finished();
    }
    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Vehicle data                          /////////////////
////////////////////////////////////////////////////////////////////////////

    double omega = 7.29211585e-5;
    MatrixXd omega_matrix(3,3);
    omega_matrix << 0.0, -omega, 0.0,
                    omega, 0.0,  0.0,
                    0.0,   0.0,  0.0;

    CONSTANTS.omega_matrix = &omega_matrix;
    CONSTANTS.mu   = 3.986004418e14;
    CONSTANTS.Re   = 6378137.0;
    CONSTANTS.g0   = 9.80665;
    CONSTANTS.rho0 = 1.225;      // retained only for the initial guess
    CONSTANTS.H    = 7200.0;
    CONSTANTS.p0   = 101325.0;
    CONSTANTS.cd   = 0.381;
    CONSTANTS.sref = 9.079;
    CONSTANTS.m_fairing = 535.3;

    //                        stage 1      stage 2      stage 3     stage 4
    double Tv0[4]   = { 2827.37e3,   1075.73e3,   299.81e3,   2.4509e3 };
    double Tv1[4]   = { 1884.91e3,    717.15e3,   221.60e3,   2.4509e3 };
    double md0[4]   = {   1034.09,      382.81,     104.61,     0.7919 };
    double md1[4]   = {    689.40,      255.21,      77.32,     0.7919 };
    double Ae[4]    = {     3.092,       1.697,      1.183,       0.07 };
    double tb[4]    = {     102.0,        75.0,      110.0,      502.1 };
    double mprop[4] = {   87898.0,     23926.0,    10006.0,      397.6 };
    double mdry[4]  = {    8417.7,      2563.8,     1326.5,      813.7 };
    for (int i = 0; i < 4; i++) {
        CONSTANTS.Tv0[i]=Tv0[i]; CONSTANTS.Tv1[i]=Tv1[i];
        CONSTANTS.md0[i]=md0[i]; CONSTANTS.md1[i]=md1[i];
        CONSTANTS.Ae[i]=Ae[i];   CONSTANTS.tb[i]=tb[i];
        CONSTANTS.mprop[i]=mprop[i]; CONSTANTS.mdry[i]=mdry[i];
    }
    CONSTANTS.t_start[0] = 0.0;
    CONSTANTS.t_start[1] = tb[0];
    CONSTANTS.t_start[2] = tb[0] + tb[1];
    CONSTANTS.t_start[3] = tb[0] + tb[1] + tb[2];

    // Lift-off mass, for a payload that the optimiser will improve upon.
    double m_payload_guess = 1400.0;
    double m0 = m_payload_guess + CONSTANTS.m_fairing;
    for (int i = 0; i < 4; i++) m0 += mprop[i] + mdry[i];

    // mass at each phase boundary, for the state bounds and the guess
    double mstart[4], mend[4];
    mstart[0] = m0;
    for (int i = 0; i < 4; i++) {
        mend[i] = mstart[i] - mprop[i];
        if (i < 3) {
            mstart[i+1] = mend[i] - mdry[i];
            if (i == 1) mstart[i+1] -= CONSTANTS.m_fairing;
        }
    }

    printf("lift-off mass %.1f kg; mass at each stage ignition:", m0);
    for (int i = 0; i < 4; i++) printf(" %.1f", mstart[i]);
    printf(" kg\n");

////////////////////////////////////////////////////////////////////////////
///////////////////  Bounds                                /////////////////
////////////////////////////////////////////////////////////////////////////

    double x0 = CONSTANTS.Re, y0 = 0.0, z0 = 0.0;      // equator, Guiana meridian
    MatrixXd r0(3,1); r0 << x0, y0, z0;
    MatrixXd v0 = omega_matrix*r0;                      // at rest on the ground

    double rmax = 2.0*CONSTANTS.Re, rmin = -rmax;
    double vmax = 12000.0, vmin = -vmax;

    // mass range and clock for each phase; phase 4 coasts at constant mass
    double mlo[5] = { mend[0], mend[1], mend[2], mstart[3], mend[3] };
    double mhi[5] = { mstart[0], mstart[1], mstart[2], mstart[3], mstart[3] };

    double t3end = CONSTANTS.t_start[3];        // third stage burns out here
    double coast_max = 2500.0;                  // generous bound on the coast

    for (int i = 1; i <= 5; i++) {
        problem.phases(i).bounds.lower.states << rmin, rmin, rmin, vmin, vmin, vmin, mlo[i-1];
        problem.phases(i).bounds.upper.states << rmax, rmax, rmax, vmax, vmax, vmax, mhi[i-1];
        if (i != 4) {
            problem.phases(i).bounds.lower.controls << -1.0, -1.0, -1.0;
            problem.phases(i).bounds.upper.controls <<  1.0,  1.0,  1.0;
            problem.phases(i).bounds.lower.path << 1.0;
            problem.phases(i).bounds.upper.path << 1.0;
        }
    }

    // the three solid stages burn for their published durations
    for (int i = 1; i <= 3; i++) {
        problem.phases(i).bounds.lower.StartTime = CONSTANTS.t_start[i-1];
        problem.phases(i).bounds.upper.StartTime = CONSTANTS.t_start[i-1];
        problem.phases(i).bounds.lower.EndTime   = CONSTANTS.t_start[i-1] + tb[i-1];
        problem.phases(i).bounds.upper.EndTime   = CONSTANTS.t_start[i-1] + tb[i-1];
    }
    // the coast and the liquid burn are both of a duration the optimiser chooses
    problem.phases(4).bounds.lower.StartTime = t3end;
    problem.phases(4).bounds.upper.StartTime = t3end;
    problem.phases(4).bounds.lower.EndTime   = t3end + 10.0;
    problem.phases(4).bounds.upper.EndTime   = t3end + coast_max;

    problem.phases(5).bounds.lower.StartTime = t3end + 10.0;
    problem.phases(5).bounds.upper.StartTime = t3end + coast_max;
    problem.phases(5).bounds.lower.EndTime   = t3end + 10.0 + 30.0;
    problem.phases(5).bounds.upper.EndTime   = t3end + coast_max + tb[3];

    problem.phases(1).bounds.lower.events << r0(0), r0(1), r0(2), v0(0), v0(1), v0(2), m0;
    problem.phases(1).bounds.upper.events << r0(0), r0(1), r0(2), v0(0), v0(1), v0(2), m0;

    double af    = CONSTANTS.Re + 700000.0;     // 700 km circular
    double e2max = 1.0e-3*1.0e-3;               // eccentricity at most 1e-3
    double cimax = cos(89.99*M_PI/180.0);       // inclination within 0.01 deg of polar
    problem.phases(5).bounds.lower.events << af, 0.0,  -cimax;
    problem.phases(5).bounds.upper.events << af, e2max, cimax;

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

    // A crude guess: the vehicle climbs on a quarter-circle pitch programme
    // through the three solid stages, coasts up to the target radius, and
    // circularises. It is not a trajectory, only something of the right shape
    // and the right magnitudes for the solver to start from.
    {
      double rt = CONSTANTS.Re + 700000.0;
      double vt = sqrt(CONSTANTS.mu/rt);
      double tt[6] = { 0.0, tb[0], tb[0]+tb[1], t3end, t3end+600.0, t3end+600.0+tb[3] };
      double frac_r[6] = { 0.0, 0.02, 0.10, 0.30, 0.95, 1.00 };
      double frac_v[6] = { 0.0, 0.15, 0.45, 0.90, 0.93, 1.00 };

      for (int i = 1; i <= 5; i++) {
        int N = 20;
        MatrixXd xg = zeros(7, N);
        MatrixXd tg = linspace(tt[i-1], tt[i], N);
        MatrixXd ug;
        if (i != 4) ug = zeros(3, N);
        for (int j = 0; j < N; j++) {
            double w  = ((double) j)/((double)(N-1));
            double fr = frac_r[i-1] + (frac_r[i]-frac_r[i-1])*w;
            double fv = frac_v[i-1] + (frac_v[i]-frac_v[i-1])*w;
            double rr = CONSTANTS.Re + (rt - CONSTANTS.Re)*fr;
            double th = 0.5*M_PI*fv;                  // downrange angle, northward
            xg(0,j) = rr*cos(th);
            xg(1,j) = 0.0;
            xg(2,j) = rr*sin(th);
            xg(3,j) = -vt*fv*sin(th);
            xg(4,j) = v0(1)*(1.0 - fv);
            xg(5,j) =  vt*fv*cos(th);
            xg(6,j) = mhi[i-1] + (mlo[i-1]-mhi[i-1])*w;
            if (i != 4) {
                // thrust from vertical at lift-off to horizontal at insertion
                double a = 0.5*M_PI*fv;
                ug(0,j) = cos(a)*cos(th) - sin(a)*sin(th);
                ug(1,j) = 0.0;
                ug(2,j) = cos(a)*sin(th) + sin(a)*cos(th);
                double nn = sqrt(ug(0,j)*ug(0,j)+ug(1,j)*ug(1,j)+ug(2,j)*ug(2,j));
                ug(0,j)/=nn; ug(1,j)/=nn; ug(2,j)/=nn;
            }
        }
        problem.phases(i).guess.states = xg;
        problem.phases(i).guess.time   = tg;
        if (i != 4) problem.phases(i).guess.controls = ug;
      }
    }

////////////////////////////////////////////////////////////////////////////
///////////////////  Algorithm options                      ////////////////
////////////////////////////////////////////////////////////////////////////

    algorithm.nlp_iter_max       = 2000;
    algorithm.nlp_tolerance      = 1.e-6;
    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.defect_scaling     = "jacobian-based";
    algorithm.collocation_method = "Hermite-Simpson";
    algorithm.mesh_refinement    = "automatic";
    algorithm.mr_max_iterations  = 3;
    algorithm.ode_tolerance      = 1.e-5;

    int rc = psopt(solution, problem, algorithm);
    if (rc != 0) printf("psopt returned %d\n", rc);

////////////////////////////////////////////////////////////////////////////
///////////////////  Results                                ////////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x4 = solution.get_states_in_phase(5);
    MatrixXd t4 = solution.get_time_in_phase(5);
    int N4 = (int) x4.cols();

    double mfin = x4(6, N4-1);
    double tfin = t4(0, N4-1);
    printf("\n");
    MatrixXd t5s = solution.get_time_in_phase(5);
    printf("coast before the liquid stage %10.4f s\n", t5s(0,0) - CONSTANTS.t_start[3]);
    printf("burn time of the liquid stage %10.4f s of the %.1f s available\n",
           tfin - t5s(0,0), tb[3]);
    printf("mass at insertion             %10.4f kg\n", mfin);
    printf("payload delivered             %10.4f kg\n", mfin - CONSTANTS.mdry[3]);
    printf("(the reference reports 1400.73 kg for its own thirteen-phase model)\n\n");

    Save(x4, "x4.dat");
    Save(t4, "t4.dat");
    return 0;
}
