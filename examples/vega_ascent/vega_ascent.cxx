//////////////////////////////////////////////////////////////////////////
////////////////          vega_ascent.cxx            /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Maximum payload ascent of a four-stage    ////////////////
////////        launch vehicle to a polar orbit           ////////////////
//////// Last modified: 31 August 2026                    ////////////////
//////// Reference:     B. Benedikter, A. Zavoli,         ////////////////
////////         G. Colasurdo, S. Pizzurro and            ////////////////
////////         E. Cavallini, Convex optimization of     ////////////////
////////         launch vehicle ascent trajectory with    ////////////////
////////         heat-flux and splash-down constraints,   ////////////////
////////         J. Spacecraft and Rockets 59(3), 2022,   ////////////////
////////         pp. 900-915, DOI 10.2514/1.A35194.       ////////////////
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
// succession for their published durations; the liquid fourth stage burns
// twice, on either side of a long coast, for durations the optimiser chooses.
// Dry masses are jettisoned at the stage boundaries and the fairing part way
// through the third-stage burn.
//
// This is a faithful transcription of the reference's phase structure, and it
// is a good deal less free than it looks. Of the 3432 seconds of flight, the
// optimiser controls the attitude during only about 260 of them:
//
//   phase  1   4.10 s   vertical ascent      thrust along the local vertical
//   phase  2   6.60 s   pitch-over           thrust free            <-- control
//   phase  3  91.30 s   gravity turn         thrust along v_rel
//   phase  4   6.60 s   coast, stage 1 spent
//   phase  5  75.00 s   gravity turn         thrust along v_rel
//   phase  6  37.30 s   coast, stage 2 spent
//   phase  7   5.40 s   third stage, fairing on   thrust free       <-- control
//   phase  8 104.60 s   third stage, fairing off  thrust free       <-- control
//   phase  9  15.40 s   coast, stage 3 spent
//   phase 10   free     fourth stage, first firing   thrust free    <-- control
//   phase 11   free     coast to apogee
//   phase 12   free     fourth stage, second firing  thrust free    <-- control
//
// Phases 3 and 5 are zero-lift gravity turns: the vehicle is not steered at
// all, the thrust simply follows the Earth-relative velocity, so the whole of
// the atmospheric flight is set by how far, and in what direction, the vehicle
// pitches over during the 6.6 seconds of phase 2. That single manoeuvre is the
// reason the problem is hard and the reason it is worth solving. Its duration
// is held at the reference's value rather than left free, because a longer
// free-attitude arc would let the optimiser fly at an angle of attack through
// the densest part of the atmosphere -- which is precisely what the gravity
// turn exists to avoid, and what no launch vehicle of this class would do.
//
// PSOPT allows a different number of controls in each phase, and that is used
// here: the seven phases with a prescribed thrust direction carry no controls
// and no path constraint at all, so the NLP is only about a third the size a
// uniform formulation would give.
//
// Two features of the vehicle model are worth noticing. The vacuum thrust and
// the mass flow of each solid stage vary linearly through the burn, in the ratio
// that keeps the specific impulse constant, which is how a solid motor with a
// regressive grain actually behaves. And the delivered thrust is the vacuum
// thrust less the ambient pressure acting over the nozzle exit area, so the
// first stage gains some 300 kN as it climbs out of the atmosphere.
//
// Once the fairing is gone the payload is exposed, and the reference bounds the
// free-molecular heat flux it sees for the rest of the flight, on phases 8 to
// 12. That constraint is imposed here too, and it is the constraint that shapes
// the answer. Left out, the vehicle sheds its fairing far lower and much faster
// than it should, into a flux of tens of kilowatts per square metre; and the
// problem that remains is not merely optimistic but badly posed, since nothing
// then keeps the trajectory out of air dense enough for the drag term to
// dominate -- run without it, the solver does not converge on any mesh.
//
// It is also the constraint that makes the atmosphere model matter. At orbital
// speed 900 W/m^2 is reached at about 140 km, and around 120 km, where the
// jettison happens, the density is a hundredth of a millionth of its sea-level
// value. An exponential extrapolation above 86 km will not do; us76.h therefore
// carries the standard's upper branch as well as its lower one.
//
// The one thing in the reference not modelled here is the splash-down
// constraint on the spent third stage, which needs a thirteenth phase
// propagating the discarded stage to the sea. It is inactive in the case solved
// here, which is the one the reference calls the unconstrained return.

#include "psopt.h"
#include "us76.h"

using namespace std;
using namespace PSOPT;

//////////////////////////////////////////////////////////////////////////
///////////////////  Phase structure                       ///////////////
//////////////////////////////////////////////////////////////////////////

#define NPH  12

// which stage is burning in each phase (0 = none)
static const int PH_STAGE[NPH] = { 1, 1, 1, 0, 2, 0, 3, 3, 0, 4, 0, 4 };

// how the thrust is pointed: 0 coast, 1 free (three controls), 2 along the
// Earth-relative velocity, 3 along the local vertical
#define MODE_COAST   0
#define MODE_FREE    1
#define MODE_ZLGT    2
#define MODE_RADIAL  3
static const int PH_MODE[NPH] = { 3, 1, 2, 0, 2, 0, 1, 1, 0, 1, 0, 1 };

// the heat-flux constraint applies from the fairing jettison onwards
static const int PH_QDOT[NPH] = { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1 };
#define QDOT_MAX  900.0             // W/m^2, the reference's own threshold

// the reference's own phase durations, used for the fixed ones and as the
// starting point for the free ones
static const double PH_DUR[NPH] = { 4.1, 6.6, 91.3, 6.6, 75.0, 37.3,
                                    5.4, 104.6, 15.4, 359.71, 2583.50, 142.39 };

static double PH_T0[NPH+1];        // nominal phase boundary times, filled in main

//////////////////////////////////////////////////////////////////////////
///////////////////  Vehicle and environment data  ///////////////////////
//////////////////////////////////////////////////////////////////////////

// SI throughout: metres, seconds, kilograms, newtons.

struct Constants {
  double mu, Re, g0, omega;
  double cd, sref;
  // per stage: vacuum thrust at start and end, mass flow at start and end,
  // nozzle exit area, burn time, propellant mass, dry mass
  double Tv0[4], Tv1[4], md0[4], md1[4], Ae[4], tb[4], mprop[4], mdry[4];
  double t_ign[4];
  double m_fairing;
};

typedef struct Constants Constants_;

//////////////////////////////////////////////////////////////////////////
///////////////////  The equations of motion               ///////////////
//////////////////////////////////////////////////////////////////////////

// One implementation, used both by the collocation (T = adouble) and by the
// forward integration that builds the initial guess (T = double), so the two
// cannot drift apart.

template <class T>
static void vega_rhs(const Constants_& C, int mode, int stage,
                     const T& t, const T* x, const T* u, T* dx, T* qdot = NULL)
{
   int j;
   T rad = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
   T alt = rad - C.Re;

   // The atmosphere is the US Standard of 1976, which is what the reference
   // uses, and it is needed over the whole of its range: the ambient pressure
   // on the nozzles and the drag through the first stage come from the layer
   // model below 86 km, and the heat-flux constraint after fairing jettison
   // comes from the diffusive branch above it, around 120 km, where the density
   // is eight orders of magnitude smaller.
   T rho, pamb, tamb;
   us76<T>(alt, rho, pamb, tamb);

   // Earth-relative velocity, with the rotation vector along the polar axis
   T vrel[3];
   vrel[0] = x[3] + C.omega*x[1];
   vrel[1] = x[4] - C.omega*x[0];
   vrel[2] = x[5];
   // The relative speed is exactly zero on the launch pad, where its derivative
   // would be unbounded, so it is regularised by a millimetre per second.
   T vr = sqrt( vrel[0]*vrel[0] + vrel[1]*vrel[1] + vrel[2]*vrel[2] + 1.0e-6 );

   T m = x[6];
   T qb = 0.5*rho*vr*C.sref*C.cd/m;         // drag acceleration per unit of v_rel

   T thrust, mdot;
   if (mode == MODE_COAST) {
      thrust = 0.0*rad;                     // keep the type, not the value
      mdot   = 0.0*rad;
   }
   else {
      int s = stage - 1;
      T tau;
      if (s == 3) tau = 0.0*rad;            // the liquid stage is not regressive
      else        tau = (t - C.t_ign[s])/C.tb[s];
      T tvac = C.Tv0[s] + (C.Tv1[s] - C.Tv0[s])*tau;
      mdot   = C.md0[s] + (C.md1[s] - C.md0[s])*tau;
      thrust = tvac - pamb*C.Ae[s];         // ambient pressure over the nozzle
   }

   T d[3];
   if (mode == MODE_FREE)
      { for (j = 0; j < 3; j++) d[j] = u[j]; }
   else if (mode == MODE_ZLGT)
      { for (j = 0; j < 3; j++) d[j] = vrel[j]/vr; }
   else if (mode == MODE_RADIAL)
      { for (j = 0; j < 3; j++) d[j] = x[j]/rad; }
   else
      { for (j = 0; j < 3; j++) d[j] = 0.0*rad; }

   T tovm = thrust/m;
   T g    = C.mu/(rad*rad*rad);

   for (j = 0; j < 3; j++) dx[j] = x[3+j];
   for (j = 0; j < 3; j++)
      dx[3+j] = -g*x[j] - qb*vrel[j] + tovm*d[j];
   dx[6] = -mdot;

   // the free-molecular heat flux the payload sees once the fairing is off
   if (qdot) *qdot = 0.5*rho*vr*vr*vr;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
   // The mass at insertion, less the dry mass of the fourth stage, is the
   // payload -- and it is the payload whether or not any propellant is left in
   // the tanks, because propellant that is not loaded can be replaced kilogram
   // for kilogram by payload without changing the trajectory.
   if (iphase == NPH) return -final_states[6];
   else               return 0.0;
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
   int mode  = PH_MODE[iphase-1];
   int stage = PH_STAGE[iphase-1];

   adouble qdot;
   vega_rhs<adouble>(C, mode, stage, time, states, controls, derivatives, &qdot);

   int ip = 0;
   if (mode == MODE_FREE)
      path[ip++] = controls[0]*controls[0] + controls[1]*controls[1]
                 + controls[2]*controls[2];     // the thrust direction is a unit vector
   if (PH_QDOT[iphase-1])
      path[ip++] = qdot;                        // and the payload must not be cooked
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

template <class T>
static void terminal_orbit(const T* rv, const T* vv, double mu,
                           T& a, T& e2, T& cosi)
{
   T r  = sqrt( rv[0]*rv[0] + rv[1]*rv[1] + rv[2]*rv[2] );
   T v2 = vv[0]*vv[0] + vv[1]*vv[1] + vv[2]*vv[2];
   T rdotv = rv[0]*vv[0] + rv[1]*vv[1] + rv[2]*vv[2];

   a = 1.0/( 2.0/r - v2/mu );

   T ev[3];
   for (int j = 0; j < 3; j++)
      ev[j] = ( (v2 - mu/r)*rv[j] - rdotv*vv[j] )/mu;
   e2 = ev[0]*ev[0] + ev[1]*ev[1] + ev[2]*ev[2];

   T hv[3];
   hv[0] = rv[1]*vv[2] - rv[2]*vv[1];
   hv[1] = rv[2]*vv[0] - rv[0]*vv[2];
   hv[2] = rv[0]*vv[1] - rv[1]*vv[0];
   cosi = hv[2]/sqrt( hv[0]*hv[0] + hv[1]*hv[1] + hv[2]*hv[2] );
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

   // At lift-off the vehicle is at rest on the ground. The mass is NOT fixed:
   // the lift-off mass is what the payload is, and the optimiser chooses it.
   if (iphase == 1)
      for (j = 0; j < 6; j++) e[j] = initial_states[j];

   if (iphase == NPH) {
      adouble rv[3]; rv[0]=final_states[0]; rv[1]=final_states[1]; rv[2]=final_states[2];
      adouble vv[3]; vv[0]=final_states[3]; vv[1]=final_states[4]; vv[2]=final_states[5];
      adouble a, e2, cosi;
      terminal_orbit<adouble>( rv, vv, C.mu, a, e2, cosi );
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

    // auto_link makes the seven states and the time continuous across each of
    // the eleven boundaries; the mass link is then corrected by whatever is
    // thrown away there. auto_link writes the mass link second from last, so
    // linkages[index-2] is the one to correct.
    for (int i = 1; i < NPH; i++) {
        auto_link(linkages, &index, xad, i, i+1, workspace);
        if (i == 3) linkages[index-2] -= C.mdry[0];        // first stage
        if (i == 5) linkages[index-2] -= C.mdry[1];        // second stage
        if (i == 7) linkages[index-2] -= C.m_fairing;      // fairing
        if (i == 8) linkages[index-2] -= C.mdry[2];        // third stage
    }

    // The two firings of the fourth stage cannot together burn more propellant
    // than the stage carries. Written as an inequality rather than as the
    // equality of the reference's eq. (26), so that whether the stage really
    // does burn to depletion comes out as a result and not as an assumption.
    adouble d10 = get_final_time(xad, 10, workspace) - get_initial_time(xad, 10, workspace);
    adouble d12 = get_final_time(xad, 12, workspace) - get_initial_time(xad, 12, workspace);
    linkages[index++] = d10 + d12 - C.tb[3];
}

////////////////////////////////////////////////////////////////////////////
///////////////////  The initial guess                     /////////////////
////////////////////////////////////////////////////////////////////////////

// The guess is a forward integration of the equations of motion themselves.
// The seven uncontrolled phases fly their prescribed steering; of the five
// controlled ones, the pitch-over turns at a constant rate towards a fixed
// direction, the third stage follows the Earth-relative velocity, and the two
// firings of the fourth stage thrust along the inertial velocity.
//
// That leaves five numbers. Three of them -- the two firing durations and the
// coast between them -- are found by shooting, in the order the flight itself
// suggests: the first firing is bisected until the apogee of the resulting
// orbit reaches the target radius, the coast is the Keplerian time from there
// to apogee, and the second firing is chosen to make the orbit as circular as
// it can be made. The remaining two are the size and the azimuth of the
// pitch-over, and they are scanned. The azimuth matters: the pad is moving east
// at 465 m/s, so a vehicle that pitches due north reaches orbit about three and
// a half degrees away from polar, and the pitch-over is the only place in the
// atmospheric flight where that can be taken out.
//
// Starting the solver from a trajectory that satisfies the dynamics -- even if
// it satisfies nothing else -- is what makes a twelve-phase problem with a
// forty-minute coast converge at all.

struct GuessOut {
   MatrixXd xg[NPH], ug[NPH], tg[NPH];
   double   dur[NPH];
   double   qmax;
};

static void guess_direction(const Constants_& C, int iph, double t, double t0,
                            double dur, double delta, double psi,
                            const double* y, double* u)
{
   int j;
   double rad = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
   if (iph == 2) {
      // turn at a constant rate away from the local vertical, towards a fixed
      // compass bearing measured from north
      double rh[3], nn[3], ee[3], sn = 0.0, se = 0.0;
      for (j = 0; j < 3; j++) rh[j] = y[j]/rad;
      ee[0] = -rh[1]; ee[1] = rh[0]; ee[2] = 0.0;            // east
      nn[0] = -rh[2]*rh[0]; nn[1] = -rh[2]*rh[1]; nn[2] = 1.0 - rh[2]*rh[2];
      for (j = 0; j < 3; j++) { sn += nn[j]*nn[j]; se += ee[j]*ee[j]; }
      sn = sqrt(sn) + 1.0e-12; se = sqrt(se) + 1.0e-12;
      double w = (dur > 0.0) ? (t - t0)/dur : 1.0;
      if (w < 0.0) w = 0.0;
      if (w > 1.0) w = 1.0;
      double ca = cos(delta*w), sa = sin(delta*w);
      for (j = 0; j < 3; j++)
         u[j] = ca*rh[j] + sa*( cos(psi)*nn[j]/sn + sin(psi)*ee[j]/se );
   }
   else if (iph == 7 || iph == 8) {
      double vr[3], s = 0.0;
      vr[0] = y[3] + C.omega*y[1];
      vr[1] = y[4] - C.omega*y[0];
      vr[2] = y[5];
      for (j = 0; j < 3; j++) s += vr[j]*vr[j];
      s = sqrt(s) + 1.0e-6;
      for (j = 0; j < 3; j++) u[j] = vr[j]/s;
   }
   else {
      double s = sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5]) + 1.0e-6;
      for (j = 0; j < 3; j++) u[j] = y[3+j]/s;
   }
}

// Integrates one phase for a given duration with a fourth order Runge-Kutta
// rule, optionally storing N equally spaced nodes. Returns false if the vehicle
// falls back into the Earth.
static bool propagate(const Constants_& C, int iph, double t0, double dur,
                      double* y, double delta, double psi,
                      MatrixXd* xg, MatrixXd* ug, MatrixXd* tg, int N,
                      double* qmax = NULL)
{
   int j, k, mode = PH_MODE[iph-1], stage = PH_STAGE[iph-1];
   double u[3], t = t0, qd;

   // half a second through the atmosphere, five seconds in the vacuum
   double hmax = (iph <= 9) ? 0.5 : 5.0;
   int nsub = (int) ceil( dur/((double)(N-1))/hmax );
   if (nsub < 2) nsub = 2;
   double h = dur/((double)(N-1))/((double) nsub);

   if (xg) { *xg = zeros(7,N); *tg = zeros(1,N); if (ug) *ug = zeros(3,N); }

   for (k = 0; k < N; k++) {
      if (xg) {
         for (j = 0; j < 7; j++) (*xg)(j,k) = y[j];
         (*tg)(0,k) = t;
         if (ug) {
            guess_direction(C, iph, t, t0, dur, delta, psi, y, u);
            for (j = 0; j < 3; j++) (*ug)(j,k) = u[j];
         }
      }
      if (k == N-1) break;
      for (int sub = 0; sub < nsub; sub++) {
         double k1[7], k2[7], k3[7], k4[7], yt[7], tt;
         guess_direction(C, iph, t, t0, dur, delta, psi, y, u);
         vega_rhs<double>(C, mode, stage, t, y, u, k1, &qd);
         if (qmax && qd > *qmax) *qmax = qd;
         tt = t + 0.5*h;
         for (j = 0; j < 7; j++) yt[j] = y[j] + 0.5*h*k1[j];
         guess_direction(C, iph, tt, t0, dur, delta, psi, yt, u);
         vega_rhs<double>(C, mode, stage, tt, yt, u, k2);
         for (j = 0; j < 7; j++) yt[j] = y[j] + 0.5*h*k2[j];
         guess_direction(C, iph, tt, t0, dur, delta, psi, yt, u);
         vega_rhs<double>(C, mode, stage, tt, yt, u, k3);
         tt = t + h;
         for (j = 0; j < 7; j++) yt[j] = y[j] + h*k3[j];
         guess_direction(C, iph, tt, t0, dur, delta, psi, yt, u);
         vega_rhs<double>(C, mode, stage, tt, yt, u, k4);
         for (j = 0; j < 7; j++)
            y[j] += (h/6.0)*(k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]);
         t += h;
         double rr = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
         if ( !(rr > C.Re - 1.0) || !(rr < 1.0e9) ) return false;
      }
   }
   return true;
}

// Osculating apoapsis radius, and the Keplerian time from here to apoapsis.
static void kepler_apoapsis(const Constants_& C, const double* y,
                            double& r_apo, double& t_apo)
{
   double a, e2, cosi;
   terminal_orbit<double>(y, y+3, C.mu, a, e2, cosi);
   double e = sqrt(e2 > 0.0 ? e2 : 0.0);
   if ( !(a > 0.0) || e >= 1.0 ) { r_apo = 1.0e12; t_apo = 0.0; return; }
   r_apo = a*(1.0 + e);
   double r  = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
   double rv = y[0]*y[3] + y[1]*y[4] + y[2]*y[5];
   double cE = (1.0 - r/a)/(e + 1.0e-14);
   double sE = rv/( (e + 1.0e-14)*sqrt(C.mu*a) );
   double E  = atan2(sE, cE);
   double M  = E - e*sin(E);
   double n  = sqrt(C.mu/(a*a*a));
   double dM = M_PI - M;
   while (dM < 0.0) dM += 2.0*M_PI;
   t_apo = dM/n;
}

// Flies the whole ascent for a given pitch-over, shooting the three durations
// of the fourth-stage sequence. Returns a score in metres: how far the orbit
// ends from the target, plus a charge for the propellant the fourth stage uses.
static double fly_guess(const Constants_& C, double m0, double delta, double psi,
                        const int* nn, double a_target, GuessOut* out)
{
   int i, j;
   double y[7], t = 0.0;
   double drop[NPH];
   for (i = 0; i < NPH; i++) drop[i] = 0.0;
   drop[2] = C.mdry[0];        // after phase 3
   drop[4] = C.mdry[1];        // after phase 5
   drop[6] = C.m_fairing;      // after phase 7
   drop[7] = C.mdry[2];        // after phase 8

   y[0] = C.Re; y[1] = 0.0; y[2] = 0.0;
   y[3] = 0.0;  y[4] = C.omega*C.Re; y[5] = 0.0;
   y[6] = m0;

   double qmax = 0.0;

   // ---- the nine phases whose durations the vehicle fixes
   for (i = 0; i < 9; i++) {
      bool ok = propagate(C, i+1, t, PH_DUR[i], y, delta, psi,
                          out ? &out->xg[i] : NULL,
                          (out && PH_MODE[i] == MODE_FREE) ? &out->ug[i] : NULL,
                          out ? &out->tg[i] : NULL, nn[i],
                          PH_QDOT[i] ? &qmax : NULL);
      if (!ok) return 1.0e12;
      if (out) out->dur[i] = PH_DUR[i];
      t += PH_DUR[i];
      y[6] -= drop[i];
   }

   double y9[7];
   for (j = 0; j < 7; j++) y9[j] = y[j];
   double t9 = t, tb3 = C.tb[3];

   // ---- first firing: burn until the apoapsis reaches the target radius
   double lo = 0.0, hi = tb3, d10 = tb3, ra, ta, yt[7];
   for (j = 0; j < 7; j++) yt[j] = y9[j];
   if ( propagate(C, 10, t9, hi, yt, delta, psi, NULL, NULL, NULL, 2) ) {
      kepler_apoapsis(C, yt, ra, ta);
      if (ra > a_target) {
         for (int it = 0; it < 50; it++) {
            double md = 0.5*(lo + hi);
            for (j = 0; j < 7; j++) yt[j] = y9[j];
            if ( !propagate(C, 10, t9, md, yt, delta, psi, NULL, NULL, NULL, 2) )
               { lo = md; continue; }
            kepler_apoapsis(C, yt, ra, ta);
            if (ra > a_target) hi = md; else lo = md;
         }
         d10 = 0.5*(lo + hi);
      }
   }
   if (d10 < 5.0) d10 = 5.0;

   // ---- the coast: Keplerian time from the end of the firing to apoapsis
   for (j = 0; j < 7; j++) yt[j] = y9[j];
   if ( !propagate(C, 10, t9, d10, yt, delta, psi, NULL, NULL, NULL, 2) )
      return 1.0e12;
   kepler_apoapsis(C, yt, ra, ta);
   double d11 = ta;
   if ( !(d11 > 10.0) )    d11 = 10.0;
   if (d11 > 7000.0)       d11 = 7000.0;

   double y11[7];
   for (j = 0; j < 7; j++) y11[j] = yt[j];
   if ( !propagate(C, 11, t9+d10, d11, y11, delta, psi, NULL, NULL, NULL, 2) )
      return 1.0e12;

   // ---- second firing: the duration that leaves the orbit most nearly circular
   double a = 0.0, b = tb3 - d10;
   if (b < 1.0) b = 1.0;
   for (int it = 0; it < 40; it++) {
      double m1 = a + 0.382*(b-a), m2 = a + 0.618*(b-a), f1, f2, aa, e2, ci;
      for (j = 0; j < 7; j++) yt[j] = y11[j];
      f1 = propagate(C, 12, t9+d10+d11, m1, yt, delta, psi, NULL, NULL, NULL, 2)
         ? ( terminal_orbit<double>(yt, yt+3, C.mu, aa, e2, ci), sqrt(fabs(e2)) )
         : 1.0e12;
      for (j = 0; j < 7; j++) yt[j] = y11[j];
      f2 = propagate(C, 12, t9+d10+d11, m2, yt, delta, psi, NULL, NULL, NULL, 2)
         ? ( terminal_orbit<double>(yt, yt+3, C.mu, aa, e2, ci), sqrt(fabs(e2)) )
         : 1.0e12;
      if (f1 < f2) b = m2; else a = m1;
   }
   double d12 = 0.5*(a + b);

   // ---- fly the three of them again, this time storing
   double dd[3] = { d10, d11, d12 };
   for (i = 9; i < NPH; i++) {
      bool ok = propagate(C, i+1, t, dd[i-9], y, delta, psi,
                          out ? &out->xg[i] : NULL,
                          (out && PH_MODE[i] == MODE_FREE) ? &out->ug[i] : NULL,
                          out ? &out->tg[i] : NULL, nn[i], &qmax);
      if (!ok) return 1.0e12;
      if (out) out->dur[i] = dd[i-9];
      t += dd[i-9];
   }

   double af, e2, cosi;
   terminal_orbit<double>(y, y+3, C.mu, af, e2, cosi);
   if ( !(af > 0.0) ) return 1.0e12;
   double e = sqrt(e2 > 0.0 ? e2 : 0.0);
   double score = fabs(af - a_target) + a_target*e + a_target*fabs(cosi)
                + 250.0*(d10 + d12);
   if (d10 + d12 > tb3)   score += 1.0e6*(d10 + d12 - tb3);
   if (qmax > QDOT_MAX)   score += 2.0e2*(qmax - QDOT_MAX);
   if (out) out->qmax = qmax;
   return score;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define main routine ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{
    Alg  algorithm;
    Sol  solution;
    Prob problem;
    int i, j;

    problem.name        = "Maximum payload Vega-class ascent to a polar orbit";
    problem.outfilename = "vega_ascent.txt";

    Constants_ CONSTANTS;
    problem.user_data = (void*) &CONSTANTS;

    problem.nphases   = NPH;
    problem.nlinkages = 8*(NPH-1) + 1;         // 88 continuity links, plus the burn budget
    psopt_level1_setup(problem);

    // The nodes are spread by how much happens in each phase, not by its
    // length: the first stage climbs through the whole atmosphere, while the
    // forty-minute coast is a Keplerian arc that a coarse mesh integrates well.
    int nlo[NPH] = {  8,  8, 20,  6, 18,  8,  6, 20,  6, 15, 12, 12 };
    int nhi[NPH] = { 12, 12, 30, 10, 26, 12, 10, 30, 10, 22, 18, 18 };

    for (i = 1; i <= NPH; i++) {
        int mode = PH_MODE[i-1];
        problem.phases(i).nstates   = 7;
        problem.phases(i).ncontrols = (mode == MODE_FREE) ? 3 : 0;
        problem.phases(i).npath     = ((mode == MODE_FREE) ? 1 : 0) + PH_QDOT[i-1];
        problem.phases(i).nevents   = (i == 1) ? 6 : ((i == NPH) ? 3 : 0);
        problem.phases(i).nodes     = (RowVectorXi(2) << nlo[i-1], nhi[i-1]).finished();
    }
    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Vehicle data                          /////////////////
////////////////////////////////////////////////////////////////////////////

    CONSTANTS.mu    = 3.986004418e14;
    CONSTANTS.Re    = 6378137.0;
    CONSTANTS.g0    = 9.80665;
    CONSTANTS.omega = 7.29211585e-5;
    CONSTANTS.cd    = 0.381;
    CONSTANTS.sref  = 9.079;
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
    for (i = 0; i < 4; i++) {
        CONSTANTS.Tv0[i]=Tv0[i]; CONSTANTS.Tv1[i]=Tv1[i];
        CONSTANTS.md0[i]=md0[i]; CONSTANTS.md1[i]=md1[i];
        CONSTANTS.Ae[i]=Ae[i];   CONSTANTS.tb[i]=tb[i];
        CONSTANTS.mprop[i]=mprop[i]; CONSTANTS.mdry[i]=mdry[i];
    }

    // nominal phase boundary times, and from them the stage ignition times
    PH_T0[0] = 0.0;
    for (i = 0; i < NPH; i++) PH_T0[i+1] = PH_T0[i] + PH_DUR[i];
    CONSTANTS.t_ign[0] = PH_T0[0];       // phase 1
    CONSTANTS.t_ign[1] = PH_T0[4];       // phase 5
    CONSTANTS.t_ign[2] = PH_T0[6];       // phase 7
    CONSTANTS.t_ign[3] = PH_T0[9];       // phase 10

////////////////////////////////////////////////////////////////////////////
///////////////////  The mass budget                       /////////////////
////////////////////////////////////////////////////////////////////////////

    // Mass removed from the vehicle -- burnt or jettisoned -- by the start and
    // by the end of each phase. The propellant burnt over part of a solid stage
    // is the integral of a mass flow that is linear in the fraction of the burn
    // elapsed. The bounds built from these are deliberately loose where a phase
    // boundary is free.
    double off_a[NPH], off_b[NPH], removed = 0.0;
    for (i = 0; i < NPH; i++) {
        off_a[i] = removed;
        int s = PH_STAGE[i] - 1;
        if (s >= 0) {
            double ta = (PH_T0[i]   - CONSTANTS.t_ign[s])/tb[s];
            double tbf= (PH_T0[i+1] - CONSTANTS.t_ign[s])/tb[s];
            if (s == 3) removed += md0[3]*PH_DUR[i];            // constant flow
            else removed += tb[s]*( md0[s]*(tbf-ta)
                                  + 0.5*(md1[s]-md0[s])*(tbf*tbf-ta*ta) );
        }
        off_b[i] = removed;
        if (i == 2) removed += mdry[0];
        if (i == 4) removed += mdry[1];
        if (i == 6) removed += CONSTANTS.m_fairing;
        if (i == 7) removed += mdry[2];
    }
    // Where a boundary is free the mass bounds must cover the whole group of
    // phases the boundary moves inside: phases 10 to 12 share the fourth stage.
    off_b[9] = off_b[11];  off_b[10] = off_b[11];
    off_a[10] = off_a[9];  off_a[11] = off_a[9];

    // Everything but the payload: propellant, jettisoned structure, and the dry
    // mass of the fourth stage, which is carried to orbit and never dropped.
    double m_inert = removed + mdry[3];
    double p_lo = 200.0, p_hi = 3000.0;             // payload bracket, kg
    double m0_lo = p_lo + m_inert, m0_hi = p_hi + m_inert;
    double m0_guess = 1400.0 + m_inert;

    printf("inert mass (propellant, structure and fairing) %10.1f kg\n", m_inert);
    printf("lift-off mass for a %.1f kg payload             %10.1f kg\n",
           1400.0, m0_guess);
    printf("mass at fourth-stage ignition                  %10.1f kg\n",
           m0_guess - off_a[9]);

////////////////////////////////////////////////////////////////////////////
///////////////////  Bounds                                /////////////////
////////////////////////////////////////////////////////////////////////////

    double r0[3] = { CONSTANTS.Re, 0.0, 0.0 };      // equator, Guiana meridian
    double v0[3] = { 0.0, CONSTANTS.omega*CONSTANTS.Re, 0.0 };   // at rest on the ground

    double rmax = CONSTANTS.Re + 2.0e6, rmin = -rmax;
    double vmax = 12000.0, vmin = -vmax;

    for (i = 1; i <= NPH; i++) {
        problem.phases(i).bounds.lower.states << rmin, rmin, rmin, vmin, vmin, vmin,
                                                 m0_lo - off_b[i-1];
        problem.phases(i).bounds.upper.states << rmax, rmax, rmax, vmax, vmax, vmax,
                                                 m0_hi - off_a[i-1];
        if (PH_MODE[i-1] == MODE_FREE) {
            problem.phases(i).bounds.lower.controls << -1.0, -1.0, -1.0;
            problem.phases(i).bounds.upper.controls <<  1.0,  1.0,  1.0;
        }
        if (PH_MODE[i-1] == MODE_FREE && PH_QDOT[i-1]) {
            problem.phases(i).bounds.lower.path << 1.0, 0.0;
            problem.phases(i).bounds.upper.path << 1.0, QDOT_MAX;
        }
        else if (PH_MODE[i-1] == MODE_FREE) {
            problem.phases(i).bounds.lower.path << 1.0;
            problem.phases(i).bounds.upper.path << 1.0;
        }
        else if (PH_QDOT[i-1]) {
            problem.phases(i).bounds.lower.path << 0.0;
            problem.phases(i).bounds.upper.path << QDOT_MAX;
        }
    }

    // Times. Every boundary from the end of phase 3 to the start of phase 10 is
    // fixed by the vehicle: the stages burn for their published durations and
    // the separation coasts are part of the flight sequence. The boundary
    // between the pitch-over and the gravity turn is free, and so are the three
    // boundaries in the fourth-stage sequence.
    double lo_t[NPH+1], hi_t[NPH+1];
    for (i = 0; i <= NPH; i++) { lo_t[i] = PH_T0[i]; hi_t[i] = PH_T0[i]; }
    lo_t[10] = PH_T0[9] + 10.0;     hi_t[10] = PH_T0[9] + tb[3]; // end of first firing
    lo_t[11] = PH_T0[9] + 20.0;     hi_t[11] = 8000.0;           // end of the coast
    lo_t[12] = PH_T0[9] + 30.0;     hi_t[12] = 8500.0;           // insertion

    for (i = 1; i <= NPH; i++) {
        problem.phases(i).bounds.lower.StartTime = lo_t[i-1];
        problem.phases(i).bounds.upper.StartTime = hi_t[i-1];
        problem.phases(i).bounds.lower.EndTime   = lo_t[i];
        problem.phases(i).bounds.upper.EndTime   = hi_t[i];
    }

    problem.phases(1).bounds.lower.events << r0[0], r0[1], r0[2], v0[0], v0[1], v0[2];
    problem.phases(1).bounds.upper.events << r0[0], r0[1], r0[2], v0[0], v0[1], v0[2];

    double af    = CONSTANTS.Re + 700000.0;     // 700 km circular
    double e2max = 1.0e-3*1.0e-3;               // eccentricity at most 1e-3
    double cimax = cos(89.99*M_PI/180.0);       // inclination within 0.01 deg of polar
    problem.phases(NPH).bounds.lower.events << af, 0.0,  -cimax;
    problem.phases(NPH).bounds.upper.events << af, e2max, cimax;

    // The eleven continuity linkages are equalities; the burn budget is not.
    problem.bounds.lower.linkage(problem.nlinkages-1) = -tb[3];
    problem.bounds.upper.linkage(problem.nlinkages-1) =  0.0;

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

    {
      GuessOut G;
      double D = 0.0, P = 0.0, best = 1.0e13;

      // a coarse sweep of the size and the bearing of the pitch-over, ...
      for (double d = 1.0; d <= 30.0; d += 1.0)
         for (double p = -10.0; p <= 4.0; p += 1.0) {
            double s = fly_guess(CONSTANTS, m0_guess, d*M_PI/180.0, p*M_PI/180.0,
                                 nhi, af, NULL);
            if (s < best) { best = s; D = d; P = p; }
         }
      // ... then a coordinate descent on the two of them
      for (double step = 0.5; step > 1.0e-4; step *= 0.5)
         for (int pass = 0; pass < 2; pass++)
            for (int k = 0; k < 4; k++) {
               double dd = D + ((k == 0) ? step : (k == 1) ? -step : 0.0);
               double pp = P + ((k == 2) ? step : (k == 3) ? -step : 0.0);
               double s = fly_guess(CONSTANTS, m0_guess, dd*M_PI/180.0, pp*M_PI/180.0,
                                    nhi, af, NULL);
               if (s < best) { best = s; D = dd; P = pp; }
            }

      fly_guess(CONSTANTS, m0_guess, D*M_PI/180.0, P*M_PI/180.0, nhi, af, &G);
      printf("guess: pitch over %6.3f deg on a bearing %+6.3f deg from north;\n"
             "       fourth stage fires for %.1f s, coasts %.1f s, fires %.1f s\n"
             "       peak heat flux after jettison %.1f W/m^2 (limit %.0f)\n",
             D, P, G.dur[9], G.dur[10], G.dur[11], G.qmax, QDOT_MAX);
      {
         double yf[7];
         for (j = 0; j < 7; j++) yf[j] = G.xg[NPH-1](j, nhi[NPH-1]-1);
         double aa, e2g, cig;
         terminal_orbit<double>(yf, yf+3, CONSTANTS.mu, aa, e2g, cig);
         printf("       and reaches a %.1f x %.1f km orbit inclined %.3f deg\n",
                (aa*(1.0-sqrt(fabs(e2g))) - CONSTANTS.Re)/1000.0,
                (aa*(1.0+sqrt(fabs(e2g))) - CONSTANTS.Re)/1000.0,
                acos(cig)*180.0/M_PI);
      }

      for (i = 1; i <= NPH; i++) {
         problem.phases(i).guess.states = G.xg[i-1];
         problem.phases(i).guess.time   = G.tg[i-1];
         if (PH_MODE[i-1] == MODE_FREE) problem.phases(i).guess.controls = G.ug[i-1];
      }
    }

////////////////////////////////////////////////////////////////////////////
///////////////////  Algorithm options                      ////////////////
////////////////////////////////////////////////////////////////////////////

    algorithm.nlp_iter_max       = 3000;
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

    MatrixXd xs[NPH], ts[NPH];
    for (i = 0; i < NPH; i++) {
        xs[i] = solution.get_states_in_phase(i+1);
        ts[i] = solution.get_time_in_phase(i+1);
    }

    int Nl = (int) xs[NPH-1].cols();
    double mfin = xs[NPH-1](6, Nl-1);
    double tfin = ts[NPH-1](0, Nl-1);
    double d10  = ts[9](0, (int)ts[9].cols()-1)  - ts[9](0,0);
    double d11  = ts[10](0,(int)ts[10].cols()-1) - ts[10](0,0);
    double d12  = ts[11](0,(int)ts[11].cols()-1) - ts[11](0,0);
    double d2   = ts[1](0,(int)ts[1].cols()-1)   - ts[1](0,0);
    double d3   = ts[2](0,(int)ts[2].cols()-1)   - ts[2](0,0);

    // where the heat-flux constraint bites, and how hard
    double qpeak = 0.0, tpeak = 0.0, hpeak = 0.0;
    for (i = 7; i < NPH; i++)
       for (int k = 0; k < (int) xs[i].cols(); k++) {
          double yy[7], dd[7], qd, uu[3] = { 0.0, 0.0, 0.0 };
          for (j = 0; j < 7; j++) yy[j] = xs[i](j,k);
          vega_rhs<double>(CONSTANTS, MODE_COAST, 0, ts[i](0,k), yy, uu, dd, &qd);
          if (qd > qpeak) {
             qpeak = qd; tpeak = ts[i](0,k);
             hpeak = sqrt(yy[0]*yy[0]+yy[1]*yy[1]+yy[2]*yy[2]) - CONSTANTS.Re;
          }
       }

    // the state at fairing jettison
    double y8[7];
    for (j = 0; j < 7; j++) y8[j] = xs[7](j, 0);
    double rr = sqrt(y8[0]*y8[0] + y8[1]*y8[1] + y8[2]*y8[2]);
    double vr8[3] = { y8[3] + CONSTANTS.omega*y8[1],
                      y8[4] - CONSTANTS.omega*y8[0], y8[5] };
    double sr = sqrt(vr8[0]*vr8[0] + vr8[1]*vr8[1] + vr8[2]*vr8[2]);
    double rho8, p8, T8;
    us76<double>(rr - CONSTANTS.Re, rho8, p8, T8);

    printf("\n");
    printf("                                     this example        Benedikter et al.\n");
    printf("pitch-over                   %12.2f s %18.2f s\n", d2,   PH_DUR[1]);
    printf("first gravity turn           %12.2f s %18.2f s\n", d3,   PH_DUR[2]);
    printf("fourth stage, first firing   %12.2f s %18.2f s\n", d10,  PH_DUR[9]);
    printf("coast to apogee              %12.2f s %18.2f s\n", d11,  PH_DUR[10]);
    printf("fourth stage, second firing  %12.2f s %18.2f s\n", d12,  PH_DUR[11]);
    printf("time of flight               %12.2f s %18.2f s\n", tfin, PH_T0[NPH]);
    printf("payload                      %12.2f kg%18.2f kg\n",
           mfin - CONSTANTS.mdry[3], 1400.73);
    printf("\n");
    printf("fourth-stage propellant burnt %11.2f kg of the %.1f kg carried\n",
           md0[3]*(d10 + d12), mprop[3]);
    printf("  -- so the reference's eq. (26), imposed there as an equality, comes\n"
           "     out here as an active inequality rather than an assumption\n");
    printf("fairing jettisoned at        %12.2f km, %.1f m/s relative,\n",
           (rr - CONSTANTS.Re)/1000.0, sr);
    printf("  in a heat flux of          %12.1f W/m^2\n", 0.5*rho8*sr*sr*sr);
    printf("peak heat flux over phases 8-12 %9.1f W/m^2 (limit %.0f) at t = %.1f s, %.1f km\n",
           qpeak, QDOT_MAX, tpeak, hpeak/1000.0);
    printf("\n");

    // one file with the whole trajectory, phase by phase
    int ntot = 0;
    for (i = 0; i < NPH; i++) ntot += (int) xs[i].cols();
    MatrixXd traj(9, ntot);
    int c = 0;
    for (i = 0; i < NPH; i++)
       for (int k = 0; k < (int) xs[i].cols(); k++, c++) {
          traj(0,c) = ts[i](0,k);
          for (j = 0; j < 7; j++) traj(1+j,c) = xs[i](j,k);
          traj(8,c) = (double)(i+1);
       }
    Save(traj, "vega_traj.dat");
    return 0;
}
