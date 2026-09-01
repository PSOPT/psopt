//////////////////////////////////////////////////////////////////////////
////////////////          horus_entry.cxx            /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Minimum heat load re-entry of a winged     ///////////////
////////        vehicle                                    ///////////////
//////// Last modified: 30 August 2026                     ///////////////
//////// Reference:     Vehicle: MBB, Study on Re-entry     //////////////
////////         Guidance and Control, Final Report Vol. 2, //////////////
////////         April 1988, ESA contract 6718/86/NL/AN;    //////////////
////////         aerodynamic tables digitised in E. Mooij,  //////////////
////////         The HORUS-2B Reference Vehicle, TU Delft   //////////////
////////         Memorandum M-692, May 1995, and used with  //////////////
////////         his kind permission. Constraint limits and //////////////
////////         the entry corridor follow M. Sagliano and  //////////////
////////         E. Mooij, Optimal drag-energy entry        //////////////
////////         guidance via pseudospectral convex         //////////////
////////         optimization, Aerospace Science and        //////////////
////////         Technology, 2021.                          //////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2026         //////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser  ///////////////
//////// General Public License (LGPL)                     //////////////
//////////////////////////////////////////////////////////////////////////

// The HORUS-2B was studied as a fully reusable second stage for Ariane 5. It
// enters at 122 km and must reach the terminal area energy management interface
// at 25 km and Mach 2.5, and it is asked here to do so while accumulating as
// little heat as possible, subject to limits on heat flux, load factor and
// dynamic pressure.
//
// Two of its ingredients deserve a word, because neither is quoted from anyone.
//
// THE AERODYNAMICS ARE TRIMMED HERE, not taken as trimmed data. The source is a
// database of UNtrimmed coefficients for the clean configuration together with
// increments for the elevon, the body flap and the rudder. A trimmed coefficient
// is obtained by choosing, at each angle of attack and Mach number, the elevon
// deflection that makes the total pitch coefficient vanish -- falling back to the
// body flap where the elevon saturates, which it does below about Mach 1.5. The
// resulting trimmed lift and drag were then fitted over the region an entry uses,
// 10 to 45 degrees and Mach 2 to 20, by the polynomial in angle of attack and
// inverse Mach number below. The fit is good to 0.027 in either coefficient,
// which is well inside the uncertainty of tables that were themselves measured
// off printed graphs in 1988.
//
// That procedure also settles where the vehicle can fly. Trim is not available
// everywhere: above about 12 degrees at Mach 1.2 the elevon and the body flap
// together cannot balance the vehicle. The trajectory here stays well clear of
// that corner, and so do the published HORUS entries.
//
// THE ATMOSPHERE is us76.h, shared with the Vega ascent example: the US Standard
// Atmosphere of 1976 written as one branch-free expression, the layer model
// below 86 km and the tabulated diffusive branch above it, accurate to 0.29 per
// cent below the join and 0.04 per cent above. A Mach number that means
// something matters here, because the aerodynamic coefficients depend on it and
// the terminal condition is stated in it.

#include "psopt.h"
#include "horus_aero.h"
#include "us76.h"
#include <vector>
#include <array>

using namespace std;
using namespace PSOPT;

//////////////////////////////////////////////////////////////////////////
///////////////////  Vehicle, planet and limits               ////////////
//////////////////////////////////////////////////////////////////////////

// SI units: metres, seconds, kilograms, radians.

#define MU_EARTH   3.986004418e14
#define RE_EARTH   6378137.0
#define OMEGA_E    7.29211585e-5

#define MASS       26029.0      // kg, HORUS-2B landing mass
#define S_REF        110.0      // m^2
#define R_NOSE         0.8      // m

// The three constraint FORMS are those of Sagliano and Mooij, from their
// equations (24) to (26). Two of the three limits are theirs as well. The heat
// flux limit is not, and the reason has to be recorded carefully, because it is
// a statement about our model rather than about theirs.
//
// They limit the heat flux to 5.3e5 W/m^2 and report trajectories that meet it.
// Our model cannot. The peak heat flux on an entry from their tabulated
// interface occurs at the first pull-out, near 68 km and 7140 m/s, before the
// vehicle has slowed appreciably, and the deepest point of that pull-out is
// what sets it. Flying at zero bank -- all the lift used to stay high -- and
// sweeping the angle of attack, the least peak this model can produce is
// 818 kW/m^2, at 45 degrees. It is not the atmosphere: their exponential fit
// and the US Standard 1976 fit used here agree to within three per cent
// through the whole pull-out region. It is not the heat flux correlation:
// theirs and Chapman's agree to within fifteen per cent over the speed and
// density range visited. Reaching 530 kW/m^2 at that pull-out would require
// about 2.8 times the lift our trimmed coefficients give, which is far outside
// any reasonable disagreement between two aerodynamic models.
//
// The likeliest explanation is the aerodynamics. Our coefficients are trimmed
// from the tabulated database by the computation described below; theirs are
// the fitted polynomials of Bergsma and Mooij, which we have not reproduced.
// Rather than print a limit this model demonstrably cannot meet, the heat flux
// limit below is set where this model can work: it binds, and the trajectory
// that meets it is exhibited. The other two limits are theirs unchanged.
//
// Note also that their Table 3 gives a maximum heat flux of 550 kW/m^2 while
// their equation (24) writes the same constraint as 5.3e5 W/m^2. The equation
// and the figure in which the constraint is plotted agree with each other.
#define QDOT_MAX   1.00e6       // W/m^2, ours; see above
#define NLOAD_MAX      2.5      // g,     Sagliano and Mooij eq. (25)
#define QDYN_MAX   1.00e4       // N/m^2, eq. (26)

// THE BANK RATE LIMIT IS ALSO OURS. Sagliano and Mooij bound the bank angle to
// +/- 89 degrees but state no rate limit, and without one the problem is not
// well posed: along the load-factor-constrained arc the constraint fixes the
// dynamic pressure, hence the altitude-speed path, hence the flight path angle
// and its rate, and the bank angle enters only through the second derivative of
// the active constraint. A state constraint of that order leaves the control
// almost undetermined pointwise, and a direct transcription answers by
// chattering the bank between its bounds -- which changes the cost not at all,
// but is an artefact of the discretisation rather than a trajectory a vehicle
// could fly. The remedy is the one the vehicle itself imposes: the bank angle is
// carried as a state and the bank RATE is the control, limited to 5 degrees per
// second, the same order as the Space Shuttle's.
#define SIGDOT_MAX     5.0      // deg/s

#define G0             9.80665

// Convective heat flux, Sagliano and Mooij eq. (24):
//     qdot = (C_Q/sqrt(R_n)) * sqrt(rho) * v^M_Q
// with C_Q in J s^2.15 kg^-0.5 m^-2.15. This is a different correlation from
// Chapman's, but not by much: over the speed and density range this trajectory
// visits, the two agree to within 15 per cent.
#define C_Q        5.28137e-5
#define M_Q             3.15

//////////////////////////////////////////////////////////////////////////
///////////////////  Atmosphere and aerodynamics              ////////////
//////////////////////////////////////////////////////////////////////////

// The atmosphere is us76.h, the same branch-free US Standard of 1976 the Vega
// ascent example uses, rather than the two polynomials this example carried
// before. Those were fitted to the standard over 20 to 125 km and their error
// was quoted here as 1.4 per cent -- but that figure is the root mean square in
// the LOGARITHM of the density, which is 14 per cent in the density itself.
// Measured against the standard: 10 per cent in the root mean square over the
// whole range, 15 per cent at 80 km, 21 per cent at the 122 km entry interface
// and as much as 35 per cent between 110 and 125 km. The speed of sound was
// worse: within 7 m/s below 86 km, but 44 m/s low at 110 km and 141 m/s low at
// 122 km, which put the Mach number at the entry interface at 27 instead of 18.
//
// The speed of sound is taken as sqrt(gamma p / rho), which is exact at every
// altitude and needs no separate model. That is worth noticing, because the
// usual difficulty above 86 km is that the mean molecular weight falls, so a
// speed of sound computed from the temperature and a fixed gas constant is
// wrong -- Bergsma and Mooij fit the molecular weight to 125 km for exactly this
// reason. Taking it from the pressure and the density instead makes the
// molecular weight cancel.
adouble air_density(adouble h)          // h in metres
{
   adouble rho, pres, temp;
   us76<adouble>(h, rho, pres, temp);
   return rho;
}

adouble speed_of_sound(adouble h)
{
   adouble rho, pres, temp;
   us76<adouble>(h, rho, pres, temp);
   return sqrt(1.4*pres/rho);
}

// trimmed coefficients: sum over i of alpha^i times a cubic in 1/M
static void trimmed_coefficients(adouble alpha, adouble mach, adouble& CL, adouble& CD)
{
   adouble w = 1.0/mach;
   adouble ai = 1.0;
   CL = 0.0; CD = 0.0;
   for (int i = 0; i < 5; i++) {
      adouble wj = 1.0;
      for (int j = 0; j < 4; j++) {
         CL = CL + CL_FIT[4*i+j]*ai*wj;
         CD = CD + CD_FIT[4*i+j]*ai*wj;
         wj = wj*w;
      }
      ai = ai*alpha;
   }
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
   adouble h = states[0], v = states[3];
   adouble rho = air_density(h);
   // total convective heat load, scaled so the integral is of order one
   return C_Q/sqrt(R_NOSE)*sqrt(rho)*pow(v, M_Q)/QDOT_MAX;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
   adouble h     = states[0];      // altitude, m
   adouble phi   = states[1];      // longitude, rad
   adouble theta = states[2];      // latitude, rad
   adouble v     = states[3];      // speed, m/s
   adouble gam   = states[4];      // flight path angle, rad
   adouble psi   = states[5];      // heading, rad
   adouble sigma = states[6];      // bank angle, rad

   adouble alpha    = controls[0]; // angle of attack, rad
   adouble sigmadot = controls[1]; // bank rate, rad/s

   adouble r   = RE_EARTH + h;
   adouble g   = MU_EARTH/(r*r);
   adouble rho = air_density(h);
   adouble a   = speed_of_sound(h);
   adouble M   = v/a;

   adouble CL, CD;
   trimmed_coefficients(alpha, M, CL, CD);

   adouble qdyn = 0.5*rho*v*v;
   adouble L    = qdyn*S_REF*CL/MASS;      // specific lift
   adouble D    = qdyn*S_REF*CD/MASS;      // specific drag

   derivatives[0] = v*sin(gam);
   derivatives[1] = v*cos(gam)*sin(psi)/(r*cos(theta));
   derivatives[2] = v*cos(gam)*cos(psi)/r;
   derivatives[3] = -D - g*sin(gam);
   derivatives[4] = L*cos(sigma)/v + cos(gam)*(v/r - g/v);
   derivatives[5] = L*sin(sigma)/(v*cos(gam)) + v*cos(gam)*sin(psi)*tan(theta)/r;
   derivatives[6] = sigmadot;

   // the three constraints that define the entry corridor, in the forms used by
   // Sagliano and Mooij. The load factor is the NORMAL one: the component of the
   // aerodynamic acceleration along the body normal axis, not the magnitude of
   // the whole aerodynamic force.
   path[0] = C_Q/sqrt(R_NOSE)*sqrt(rho)*pow(v, M_Q);       // heat flux, W/m^2
   path[1] = (L*cos(alpha) + D*sin(alpha))/G0;             // normal load factor, g
   path[2] = qdyn;                                         // dynamic pressure, N/m^2
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   for (int j = 0; j < 6; j++) e[j] = initial_states[j];

   adouble hf = final_states[0], vf = final_states[3];
   e[6] = hf;                                   // terminal altitude
   e[7] = vf/speed_of_sound(hf);                // terminal Mach number
   e[8] = final_states[4];                      // terminal flight path angle
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
    // With no argument the mesh is refined automatically; an argument fixes the
    // number of nodes instead, which is how the table of accuracies in the
    // algorithm options below was produced.
    int fixed_nodes = (argc > 1) ? atoi(argv[1]) : 0;

    Alg  algorithm;
    Sol  solution;
    Prob problem;

    problem.name        = "Minimum heat load re-entry of the HORUS-2B";
    problem.outfilename = "horus_entry.txt";

    problem.nphases   = 1;
    problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 7;
    problem.phases(1).ncontrols = 2;
    problem.phases(1).nevents   = 9;
    problem.phases(1).npath     = 3;
    if (fixed_nodes > 0) problem.phases(1).nodes = (RowVectorXi(1) << fixed_nodes).finished();
    else                 problem.phases(1).nodes = (RowVectorXi(3) << 40, 80, 160).finished();

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Bounds                                 ////////////////
////////////////////////////////////////////////////////////////////////////

    double d2r = M_PI/180.0;

    // the entry interface is the one tabulated by Sagliano and Mooij
    double h0 = 122000.0, v0 = 7435.5, gam0 = -1.43*d2r;
    double phi0 = -106.7*d2r, theta0 = -22.3*d2r, psi0 = 70.75*d2r;
    // The terminal interface is a POINT, not a box. Sagliano and Mooij give
    // 26.92 +/- 2 km at Mach 2.5 +/- 0.5, but those tolerances are the
    // dispersions their guidance scheme has to hold against -- their Table 6
    // reports achieved standard deviations of about 1.2 km in altitude, which is
    // what the +/- 2 km is there to cover. A reference trajectory should aim at
    // the nominal point. Treating the tolerances as design freedom instead makes
    // the problem ill-posed: the cost is nearly flat across the box, and the
    // solution then jumps between arriving at the top of it and at the bottom
    // as the mesh is refined, with the heat load differing in the fourth figure.
    double hf_lo = 26920.0, hf_hi = 26920.0;
    double mf_lo = 2.5,     mf_hi = 2.5;
    // The terminal flight path angle is NOT given in the source, and the cost
    // does not settle it: the heat flux at Mach 2 and 25 km is three orders
    // below its peak, so what the vehicle does in the last few seconds is very
    // nearly free. Left to itself the solution simply dives as steeply as the
    // bound on gamma allows -- opening that bound from 20 to 30 to 45 degrees
    // moves the arrival to 20.0, 30.0 and 44.1 degrees nose down while the heat
    // load stays at 163.4 MJ/m^2 to four figures. The angle is therefore
    // imposed rather than optimised, on the physical ground that terminal area
    // guidance takes the vehicle over on a glide and not in a dive.
    double gf_lo = -15.0*d2r, gf_hi = -5.0*d2r;

    // the bound on gamma is a sanity bound only; it is inactive at the solution
    problem.phases(1).bounds.lower.states
        << 20000.0, -M_PI, -70.0*d2r,  300.0, -45.0*d2r, -M_PI, 0.0;
    problem.phases(1).bounds.upper.states
        << 130000.0, M_PI,  70.0*d2r, 8000.0,  45.0*d2r,  M_PI, 89.0*d2r;

    // the angle of attack is held inside the range over which the vehicle can be
    // trimmed and the aerodynamic fit is valid
    // The bank angle is taken non-negative. With no lateral requirement its sign
    // is immaterial -- it enters the vertical channel through cos(sigma) alone --
    // and leaving the sign free makes the problem symmetric and the solution
    // non-unique, which the optimiser finds out the hard way.
    problem.phases(1).bounds.lower.controls << 10.0*d2r, -SIGDOT_MAX*d2r;
    problem.phases(1).bounds.upper.controls << 45.0*d2r,  SIGDOT_MAX*d2r;

    problem.phases(1).bounds.lower.path << 0.0, 0.0, 0.0;
    problem.phases(1).bounds.upper.path << QDOT_MAX, NLOAD_MAX, QDYN_MAX;

    problem.phases(1).bounds.lower.events
        << h0, phi0, theta0, v0, gam0, psi0, hf_lo, mf_lo, gf_lo;
    problem.phases(1).bounds.upper.events
        << h0, phi0, theta0, v0, gam0, psi0, hf_hi, mf_hi, gf_hi;

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   =  400.0;
    problem.phases(1).bounds.upper.EndTime   = 4000.0;

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

    // A crude guess will not do on an entry problem: the constraints are all
    // nearly active and a trajectory that is not dynamically plausible starts
    // the optimiser outside the corridor. The guess is therefore flown, with the
    // angle of attack eased from 40 degrees down to 15 as the vehicle slows --
    // the profile the published HORUS entries use -- and the bank held at the
    // value that brings it to the interface. It is a feasible trajectory, which
    // is what demonstrates that the constraint set is not empty.
    int nnodes = 120;
    MatrixXd xg = zeros(7, nnodes), ug = zeros(2, nnodes);
    MatrixXd tg = zeros(1, nnodes);
    {
        const double dt = 0.5, sigma_g = 58.0*d2r;
        vector<double> th_; vector< array<double,6> > sh_; vector<double> ah_;
        double y[6] = { h0, phi0, theta0, v0, gam0, psi0 }, tt = 0.0;

        while (tt < 4000.0 && y[0] > 25000.0) {
            double k1[6], k2[6], k3[6], k4[6], yt[6];
            double al = 0.0;
            for (int stage = 0; stage < 4; stage++) {
                double* yy = (stage == 0) ? y : yt;
                double* kk = (stage == 0) ? k1 : (stage == 1 ? k2 : (stage == 2 ? k3 : k4));
                double sfrac = (yy[3] - 800.0)/(6000.0 - 800.0);
                if (sfrac < 0.0) sfrac = 0.0; if (sfrac > 1.0) sfrac = 1.0;
                double a_l = (15.0 + 25.0*sfrac)*d2r;
                if (stage == 0) al = a_l;
                adouble hh = yy[0], vv = yy[3];
                double rr = air_density(hh).value(), aa = speed_of_sound(hh).value();
                double M = vv.value()/aa; if (M < 1.5) M = 1.5;
                adouble CLa, CDa; trimmed_coefficients(adouble(a_l), adouble(M), CLa, CDa);
                double r = RE_EARTH + yy[0], g = MU_EARTH/(r*r);
                double q = 0.5*rr*yy[3]*yy[3];
                double L = q*S_REF*CLa.value()/MASS, D = q*S_REF*CDa.value()/MASS;
                kk[0] = yy[3]*sin(yy[4]);
                kk[1] = yy[3]*cos(yy[4])*sin(yy[5])/(r*cos(yy[2]));
                kk[2] = yy[3]*cos(yy[4])*cos(yy[5])/r;
                kk[3] = -D - g*sin(yy[4]);
                kk[4] = L*cos(sigma_g)/yy[3] + cos(yy[4])*(yy[3]/r - g/yy[3]);
                kk[5] = L*sin(sigma_g)/(yy[3]*cos(yy[4]))
                        + yy[3]*cos(yy[4])*sin(yy[5])*tan(yy[2])/r;
                if (stage < 3) {
                    double f = (stage == 2) ? dt : 0.5*dt;
                    for (int i = 0; i < 6; i++) yt[i] = y[i] + f*kk[i];
                }
            }
            array<double,6> ss; for (int i = 0; i < 6; i++) ss[i] = y[i];
            th_.push_back(tt); sh_.push_back(ss); ah_.push_back(al);
            for (int i = 0; i < 6; i++) y[i] += dt/6.0*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
            tt += dt;
        }
        double tf_g = tt;
        printf("shooting guess: %.1f s to 25 km, bank held at %.1f deg\n",
               tf_g, sigma_g/d2r);
        for (int j = 0; j < nnodes; j++) {
            double tj = tf_g*((double) j)/((double)(nnodes-1));
            int idx = (int)(tj/dt); if (idx > (int)th_.size()-2) idx = (int)th_.size()-2;
            double f = (tj - th_[idx])/dt;
            for (int i = 0; i < 6; i++)
                xg(i,j) = sh_[idx][i] + f*(sh_[idx+1][i] - sh_[idx][i]);
            xg(6,j) = sigma_g;
            ug(0,j) = ah_[idx] + f*(ah_[idx+1] - ah_[idx]);
            ug(1,j) = 0.0;
            tg(0,j) = tj;
        }
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
    algorithm.mr_max_iterations  = 4;
    algorithm.ode_tolerance      = 1.e-4;

    // A WORD ON ACCURACY. This example used to report a discretisation error near
    // 0.3, and carried a long note explaining that the figure was dominated by the
    // corners the bang-bang bank rate puts into the bank angle, which a cubic
    // cannot represent, and that the answer was better than the estimate looked.
    // Most of that was wrong, and wrong for a reason worth recording: the
    // estimator was reconstructing the control by splining the node values, which
    // for Hermite-Simpson ignores the control variable at the midpoint of every
    // interval. The residual it measured was largely its own interpolation.
    //
    // The same mistake was in the independent check at the end of main, which
    // interpolated the node controls linearly. Both are fixed, and on fixed
    // meshes the problem now gives
    //
    //     nodes    error estimate    heat load      RK4 check on the heat load
    //       109        1.8e-4        181.152 MJ/m^2        0.013 per cent
    //       200        1.1e-3        180.598               0.002
    //       350        4.3e-4        180.571               0.001
    //
    // where the last column is the block at the end of main: the optimal controls
    // are read off the transcription's own representation -- the interval
    // quadratic through node, midpoint and node -- and the states propagated by
    // Runge-Kutta at a step of 20 ms, two orders below the collocation spacing.
    // The terminal altitude error over that check is -13.0, 0.5 and 0.4 metres,
    // and the terminal speed error -0.70, 0.55 and 2.32 m/s.
    //
    // There was a third instance of the same mistake, on the other side of the
    // comparison. The heat load is this problem's objective, and PSOPT integrates
    // it with the Simpson rule of the transcription over the transcription's own
    // representation; this example accumulated its own trapezoidal sum over the
    // node values of qdot and compared the propagation against that. On the
    // automatic mesh, which is strongly non-uniform, the trapezoidal sum reads
    // 180.904 MJ/m^2 against the objective's 180.679 -- a fifth of a per cent,
    // fifty times the discrepancy being measured. Both figures are printed now.
    //
    // What survives of the old note is more useful than what it said. The check
    // agrees with each collocated solution to a hundredth of a per cent or
    // better, so each mesh's trajectory is a genuine trajectory of the equations
    // of motion. The heat load nevertheless still drifts, from 181.15 to 180.57
    // as the mesh is tripled, because a finer mesh is a richer control
    // parameterisation and buys the optimiser a slightly better optimum. Those
    // are two different questions -- is this a trajectory, and is this the best
    // one -- and only the first is what a propagation check answers.

    int rc = psopt(solution, problem, algorithm);
    if (rc != 0) printf("psopt returned %d\n", rc);

////////////////////////////////////////////////////////////////////////////
///////////////////  Results                                ////////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x = solution.get_states_in_phase(1);
    MatrixXd u = solution.get_controls_in_phase(1);
    MatrixXd t = solution.get_time_in_phase(1);
    int N = (int) x.cols();

    double qmax = 0.0, nmax = 0.0, pmax = 0.0, load = 0.0;
    MatrixXd qd(1,N), nl(1,N), qp(1,N), mach(1,N);
    for (int j = 0; j < N; j++) {
        adouble hh = x(0,j), vv = x(3,j), aa = u(0,j);
        adouble rho = air_density(hh), a = speed_of_sound(hh);
        adouble M = vv/a, CL, CD;
        trimmed_coefficients(aa, M, CL, CD);
        double qdyn = 0.5*rho.value()*x(3,j)*x(3,j);
        double L = qdyn*S_REF*CL.value()/MASS, D = qdyn*S_REF*CD.value()/MASS;
        double q = C_Q/sqrt(R_NOSE)*sqrt(rho.value())*pow(x(3,j), M_Q);
        double n = (L*cos(aa.value()) + D*sin(aa.value()))/G0;
        qd(0,j)=q/1000.0; nl(0,j)=n; qp(0,j)=qdyn/1000.0; mach(0,j)=M.value();
        if (q>qmax) qmax=q; if (n>nmax) nmax=n; if (qdyn>pmax) pmax=qdyn;
        if (j>0) load += 0.5*(qd(0,j)+qd(0,j-1))*(t(0,j)-t(0,j-1));
    }

    printf("\n");
    printf("time of flight             %10.2f s\n", t(0,N-1));
    // The heat load is the problem's own objective, and PSOPT integrates it with the
    // Simpson rule of the transcription, over the transcription's own representation.
    // The trapezoidal sum accumulated in `load` above is second order and is a poor
    // substitute: on this problem it reads 180.90 against 180.68, a fifth of a per cent
    // high, which is the size of the drift the accuracy discussion of the book turns on.
    const double heat_load = solution.get_cost();      // MJ/m^2
    printf("total heat load            %10.3f MJ/m^2  (trapezoidal sum over the nodes: %.3f)\n",
           heat_load, load/1000.0);
    printf("peak heat flux             %10.2f kW/m^2  (limit %.0f)\n", qmax/1000.0, QDOT_MAX/1000.0);
    printf("peak load factor           %10.4f g       (limit %.1f)\n", nmax, NLOAD_MAX);
    printf("peak dynamic pressure      %10.2f kPa     (limit %.1f)\n", pmax/1000.0, QDYN_MAX/1000.0);
    printf("angle of attack range      %10.2f to %.2f deg\n",
           u.row(0).minCoeff()/d2r, u.row(0).maxCoeff()/d2r);
    printf("bank angle range           %10.2f to %.2f deg\n",
           x.row(6).minCoeff()/d2r, x.row(6).maxCoeff()/d2r);
    printf("bank rate range            %10.2f to %.2f deg/s (limit %.1f)\n",
           u.row(1).minCoeff()/d2r, u.row(1).maxCoeff()/d2r, SIGDOT_MAX);
    printf("terminal flight path angle %10.2f deg\n", x(4,N-1)/d2r);
    {   // great-circle range from the entry point to the interface
        double dphi = x(1,N-1) - x(1,0), th0 = x(2,0), thf = x(2,N-1);
        double cd_ = sin(th0)*sin(thf) + cos(th0)*cos(thf)*cos(dphi);
        if (cd_ >  1.0) cd_ =  1.0;
        if (cd_ < -1.0) cd_ = -1.0;
        printf("range flown                %10.2f km\n", acos(cd_)*RE_EARTH/1000.0);
    }
    printf("terminal Mach number       %10.4f\n", mach(0,N-1));
    printf("nodes in the final mesh    %10d\n", N);

    // ---------------------------------------------------------------------
    // An independent check on the transcription. The optimal controls and the
    // bank angle are interpolated and the six physical states are propagated
    // from the entry interface by fourth order Runge-Kutta at a step far below
    // the collocation spacing. If the transcription has done its job the two
    // trajectories agree; the discrepancy is a measure of the discretisation
    // error that owes nothing to the estimate the mesh refinement uses.
    // The control this check propagates has to be the control the transcription
    // used, and under Hermite-Simpson that is not the one get_controls_in_phase
    // returns. That array holds the values at the nodes; the method also carries a
    // control at the midpoint of every interval, and the representation it assumes
    // is the quadratic through the three. Interpolating the node values linearly --
    // which is what this block did until the accessors below existed -- propagates
    // a control the solver never used, and reports the difference as if it were
    // discretisation error. The bank angle is a state, so it is read off the cubic
    // Hermite that the same transcription defines for it, through the node values
    // and the node bank rates.
    MatrixXd uhs = solution.get_hs_controls_in_phase(1);
    MatrixXd ths = solution.get_hs_time_in_phase(1);
    const bool hs = (uhs.cols() == 2*N - 1);
    if (!hs)
        printf("\n(the final mesh is not Hermite-Simpson; the check below "
               "interpolates the node controls linearly)\n");
    {
        const double dts = 0.02;
        double y[6] = { x(0,0), x(1,0), x(2,0), x(3,0), x(4,0), x(5,0) };
        double tt = 0.0, tf = t(0,N-1), qload = 0.0;
        int idx = 0;
        while (tt < tf - 1e-9) {
            double step = (tt + dts > tf) ? (tf - tt) : dts;
            double k[4][6], yt[6];
            for (int stage = 0; stage < 4; stage++) {
                double ts = tt + ((stage == 0) ? 0.0 : (stage == 3 ? step : 0.5*step));
                while (idx < N-2 && t(0,idx+1) < ts) idx++;
                double hk = t(0,idx+1) - t(0,idx);
                double f  = (hk > 0.0) ? (ts - t(0,idx))/hk : 0.0;
                if (f < 0.0) f = 0.0; if (f > 1.0) f = 1.0;
                double a_l, s_l;
                if (hs) {
                    // angle of attack: the interval's quadratic through node,
                    // midpoint and node
                    a_l = (2.0*f - 1.0)*(f - 1.0)*uhs(0, 2*idx)
                        + 4.0*f*(1.0 - f)      *uhs(0, 2*idx+1)
                        + f*(2.0*f - 1.0)      *uhs(0, 2*idx+2);
                    // bank angle: the interval's cubic Hermite, whose end slopes
                    // are the collocated bank rates
                    double f2 = f*f, f3 = f2*f;
                    s_l = ( 2.0*f3 - 3.0*f2 + 1.0)*x(6,idx)
                        + hk*(f3 - 2.0*f2 + f)    *u(1,idx)
                        + (-2.0*f3 + 3.0*f2)      *x(6,idx+1)
                        + hk*(f3 - f2)            *u(1,idx+1);
                }
                else {
                    a_l = u(0,idx) + f*(u(0,idx+1) - u(0,idx));
                    s_l = x(6,idx) + f*(x(6,idx+1) - x(6,idx));
                }
                double* yy = (stage == 0) ? y : yt;
                double* kk = k[stage];
                adouble hh = yy[0];
                double rr = air_density(hh).value(), aa = speed_of_sound(hh).value();
                adouble CLa, CDa;
                trimmed_coefficients(adouble(a_l), adouble(yy[3]/aa), CLa, CDa);
                double r = RE_EARTH + yy[0], g = MU_EARTH/(r*r);
                double qq = 0.5*rr*yy[3]*yy[3];
                double L = qq*S_REF*CLa.value()/MASS, D = qq*S_REF*CDa.value()/MASS;
                kk[0] = yy[3]*sin(yy[4]);
                kk[1] = yy[3]*cos(yy[4])*sin(yy[5])/(r*cos(yy[2]));
                kk[2] = yy[3]*cos(yy[4])*cos(yy[5])/r;
                kk[3] = -D - g*sin(yy[4]);
                kk[4] = L*cos(s_l)/yy[3] + cos(yy[4])*(yy[3]/r - g/yy[3]);
                kk[5] = L*sin(s_l)/(yy[3]*cos(yy[4]))
                        + yy[3]*cos(yy[4])*sin(yy[5])*tan(yy[2])/r;
                if (stage < 3) {
                    double fs = (stage == 2) ? step : 0.5*step;
                    for (int i = 0; i < 6; i++) yt[i] = y[i] + fs*kk[i];
                }
            }
            qload += step*C_Q/sqrt(R_NOSE)
                     *sqrt(air_density(adouble(y[0])).value())*pow(y[3], M_Q);
            for (int i = 0; i < 6; i++)
                y[i] += step/6.0*(k[0][i] + 2*k[1][i] + 2*k[2][i] + k[3][i]);
            tt += step;
        }
        double Mv = y[3]/speed_of_sound(adouble(y[0])).value();
        printf("\n");
        printf("verification by Runge-Kutta propagation of the optimal controls\n");
        printf("   altitude       %10.3f km   (collocation %8.3f, error %6.1f m)\n",
               y[0]/1000.0, x(0,N-1)/1000.0, y[0] - x(0,N-1));
        printf("   speed          %10.3f m/s  (collocation %8.3f, error %6.2f m/s)\n",
               y[3], x(3,N-1), y[3] - x(3,N-1));
        printf("   Mach           %10.4f      (collocation %8.4f)\n", Mv, mach(0,N-1));
        printf("   flight path    %10.3f deg  (collocation %8.3f)\n",
               y[4]/d2r, x(4,N-1)/d2r);
        printf("   heat load      %10.3f MJ/m^2 (collocation %8.3f, %.3f%% apart)\n",
               qload/1.0e6, heat_load, 100.0*fabs(qload/1.0e6 - heat_load)/heat_load);
    }
    printf("\n");

    Save(x,"x.dat"); Save(u,"u.dat"); Save(t,"t.dat");
    Save(qd,"qdot.dat"); Save(nl,"nload.dat"); Save(qp,"qdyn.dat"); Save(mach,"mach.dat");

    MatrixXd km = x.row(0)/1000.0;
    plot(t, km, problem.name+": altitude", "time (s)", "h (km)", "h");
    plot(t, u.row(0)/d2r, problem.name+": angle of attack", "time (s)", "deg", "alpha");
    plot(t, x.row(6)/d2r, problem.name+": bank angle", "time (s)", "deg", "sigma");
    plot(t, qd, problem.name+": heat flux", "time (s)", "kW/m2", "qdot");

    return 0;
}
