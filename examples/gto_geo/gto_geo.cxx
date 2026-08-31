//////////////////////////////////////////////////////////////////////////
////////////////            gto_geo.cxx              /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Minimum time GTO to GEO orbit raising      ///////////////
////////        for an all-electric satellite              ///////////////
//////// Last modified: 31 August 2026                     ///////////////
//////// Reference:     Falck and Dankanich (2012),        ///////////////
////////         NASA/TM-2012-217699; Leomanni et al.      ///////////////
////////         (2021), J. Spacecraft and Rockets,        ///////////////
////////         DOI 10.2514/1.A34949. The transfer        ///////////////
////////         orbits and the propulsion model are       ///////////////
////////         theirs; the power level used here is      ///////////////
////////         not (see the note below).                 ///////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2026         ///////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser  ///////////////
//////// General Public License (LGPL)                     ///////////////
//////////////////////////////////////////////////////////////////////////

// The independent variable is the true longitude L rather than the time. That
// choice is what makes the problem tractable. The transfer starts in an orbit
// whose perigee radius is 6564 km, where the spacecraft travels at 10.25 km/s
// and sweeps true longitude at 1.6e-3 rad/s. A mesh uniform in time with thirty
// nodes per revolution would step 113 degrees of true longitude across perigee;
// holding five degrees per step there would need some 677 nodes per revolution,
// and 45000 nodes over the whole transfer. A mesh uniform in L spends its nodes
// evenly around each orbit instead, and thirty per revolution suffice. The
// elapsed time becomes a state, and minimum time means minimising its final
// value.
//
// A spacecraft is released into a geostationary transfer orbit and raises
// itself to a geostationary orbit under continuous low thrust, minimising the
// transfer time. The state is the set of modified equinoctial elements, which
// are free of the singularities the classical elements have at zero
// eccentricity and zero inclination -- both of which this transfer ends at.
//
// The benchmark case in the literature uses a 5 kW system giving 0.31158 N,
// for which the transfer takes about 118 days over some 160 revolutions.
// That problem is large: Leomanni et al. report 130117 nonlinear programming
// variables. This example therefore flies the same mission with a larger power
// system, which changes one number and nothing else, and brings the transfer
// down to a size a reader can run. Setting THRUST_NEWTON to 0.31158 and
// widening the bound on the final time recovers the published case.

#include "psopt.h"
#include <vector>
#include <array>
#include <algorithm>

using namespace std;
using namespace PSOPT;

//////////////////////////////////////////////////////////////////////////
///////////////////  Physical constants and vehicle data  ////////////////
//////////////////////////////////////////////////////////////////////////

// Units throughout: kilometres, seconds, kilograms.

#define MU_EARTH     398600.4418      // km^3/s^2
#define RE_EARTH       6378.137       // km
#define J2_EARTH     1.08262668e-3
#define G0             9.80665        // m/s^2, for the propellant equation

// The thruster is specified by its input power, which is what an electric
// propulsion system is actually sized by. The published benchmark runs at 5 kW,
// which at 55 per cent efficiency and 1800 s of specific impulse is 0.31158 N --
// the value Leomanni and co-workers quote, recovered here to six figures.
//
// The power is settable at run time because it decides how large the problem is.
// Thrust scales with power, the transfer time scales roughly inversely with
// thrust, and the number of revolutions with it: the 5 kW case runs to about 160
// revolutions and several thousand collocation nodes, while 20 kW needs about a
// quarter of that. A reader on a modest machine should start at 20 kW.
#define ISP_SECONDS    1800.0
#define MASS_INITIAL   1200.0         // kg
#define THRUSTER_EFF   0.55

// Steering-law weights, chosen by sweeping them against the published solution.
// At 5 kW these give a guess of 122.3 days over 160.6 revolutions, against the
// published optimum of 118.74 days over about 162 -- so the guess is some three
// per cent slow and has essentially the right revolution count, which is what
// decides the size of the discretised problem. The weight on the semi-major axis
// matters most; the transfer lengthens noticeably if the inclination is given
// too much priority early, and also if it is given too little.
static double LYAP_WA = 6.0, LYAP_WI = 10.0;
static double POWER_KW     = 5.0;
static double THRUST_NEWTON = 2.0*THRUSTER_EFF*5.0e3/(9.80665*1800.0);

// The eclipse model is ON. Every published solution of this transfer includes
// eclipsing, so an eclipse-free result would have nothing to be compared
// against. The model below is the one used by Leomanni and co-workers.
#define USE_ECLIPSE    1

// Solar and shadow constants. The epoch is theirs: JD 2451625.5 is 21 March
// 2000, essentially the vernal equinox, which is the heaviest eclipse season.
#define JD_EPOCH   2451625.5
#define JD_J2000   2451545.0
#define AU_KM      1.495978707e8
#define RSUN_KM    695700.0
#define DEG_RAD    (M_PI/180.0)

// The gain of the logistic that smooths the shadow boundary. This is their
// value, stated as a realistic one for Earth-centred transfers. It sets the
// width of the transition: the shadow function goes from 0.95 to 0.05 over
// about 4/SHADOW_GAIN radians of the geometric margin, which on this transfer
// is a little over one degree of true longitude.
#define SHADOW_GAIN    298.78

// The initial mesh. Coarse intervals carry the slow drift of the elements and
// the once-per-revolution ripple the perturbations put on them; a window of
// fine intervals straddles each eclipse transition. The window is a little wider
// than the transition itself, which spans about one degree of true longitude.
#define MESH_COARSE_DEG    45.0
#define MESH_WINDOW_DEG     2.0
#define MESH_SUBDIV           4
#define MESH_ORDER_COARSE     5
#define MESH_ORDER_FINE       4

// Two safeguards that keep the iterates inside the region where the element set
// is a valid description of the motion.
//
// The equinoctial elements f and g are the components of the eccentricity
// vector, and a box on them individually cannot bound the eccentricity: with
// each free to reach one, the box admits e up to the square root of two, which
// is a hyperbolic orbit. It admits them at every node, and the optimiser will
// take them. Nothing in the equations objects, because they are written with the
// true longitude as the independent variable and every derivative carries a
// factor 1/w with w = 1 + f cos L + g sin L, which merely becomes very large as
// the orbit approaches parabolic; the elapsed time per radian of true longitude
// grows without bound, and the whole first revolution can absorb weeks.
//
// Left unguarded, the twenty kilowatt case did exactly that: the eccentricity
// reached 1.0162 after ten days, the apoapsis reached 1.68 million kilometres,
// and the solver spent an hour of processor time crawling around a transfer of
// 77.66 days that should take about thirty, burning 474 kg of propellant
// against the 180 it needs. The iterate was feasible for the discretised
// problem throughout. It was the continuous problem that had been left open.
//
// So the eccentricity is bounded directly, and the apoapsis with it. Both are
// written as smooth polynomials in the elements -- no square roots, no
// divisions -- and both are set well outside anything the transfer needs: the
// steering-law guess peaks at e = 0.732 and an apoapsis of 67000 km, and the
// starting orbit is itself at e = 0.7306. They are safeguards, not shaping
// constraints, and the solution should be checked against them; the example
// reports how close it came.
#define E_MAX             0.80
#define R_APO_MAX     100000.0        // km

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
   return final_states[6];      // the elapsed time in days; minimum time
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
///////////////////  The shadow function                     /////////////
//////////////////////////////////////////////////////////////////////////

// The spacecraft is shadowed when the angle between the Earth and the Sun, seen
// from the spacecraft, is smaller than the sum of their apparent angular radii.
// This is the condition used by Leomanni and co-workers, and because it adds the
// apparent radius of the Sun rather than subtracting it, it marks the onset of
// penumbra rather than of umbra.
//
// Two departures from the way it is usually written, both deliberate.
//
// The condition is evaluated WITHOUT inverse trigonometric functions. Writing it
// as theta_es <= theta_e + theta_s requires an arc cosine whose derivative is
// unbounded when the Sun is directly behind the Earth -- which is exactly what
// happens in the middle of every eclipse, so it is not a corner case. Comparing
// cosines instead removes the singularity, and dividing by sin(theta_e+theta_s)
// puts the margin back into radians so that their smoothing gain keeps its
// meaning. Against the direct form this shifts the shadow boundary by about a
// tenth of the width of the smoothed transition, and changes the computed
// eclipse fraction over the whole transfer by 0.07 per cent.
//
// The throttle is not carried as a separate control. Their formulation has one,
// bounded above by the shadow function, and reports that at the solution it sits
// at its upper bound everywhere -- which is what a minimum-time transfer must do,
// since thrust with a free direction is never unhelpful. Substituting the shadow
// function directly for the throttle therefore gives the same solution while
// avoiding an extra control and an extra path constraint at every one of several
// thousand nodes.

template <class T> static void sun_vector(T tdays, T& sx, T& sy, T& sz, T& rsun)
{
   // Low-precision solar position, Astronomical Almanac form. Good to about
   // 0.01 degrees, which is far finer than this problem can distinguish.
   T n   = (JD_EPOCH - JD_J2000) + tdays;
   T Lm  = (280.460 + 0.9856474*n)*DEG_RAD;
   T gm  = (357.528 + 0.9856003*n)*DEG_RAD;
   T lam = Lm + (1.915*sin(gm) + 0.020*sin(2.0*gm))*DEG_RAD;
   double eps = 23.439*DEG_RAD;
   rsun = (1.00014 - 0.01671*cos(gm) - 0.00014*cos(2.0*gm))*AU_KM;
   sx = cos(lam);
   sy = cos(eps)*sin(lam);
   sz = sin(eps)*sin(lam);
}

// Inertial position from the equinoctial elements and the true longitude.
template <class T> static void mee_position(T p, T f, T g, T h, T k, T L,
                                            T& rx, T& ry, T& rz)
{
   T s2 = 1.0 + h*h + k*k;
   T w  = 1.0 + f*cos(L) + g*sin(L);
   T r  = p/w;
   T a2 = h*h - k*k;
   rx = r/s2*(cos(L) + a2*cos(L) + 2.0*h*k*sin(L));
   ry = r/s2*(sin(L) - a2*sin(L) + 2.0*h*k*cos(L));
   rz = 2.0*r/s2*(h*sin(L) - k*cos(L));
}

template <class T>
static T sunlit_core(T p, T f, T g, T h, T k, T L, T tdays)
{
   T rx, ry, rz;
   mee_position(p, f, g, h, k, L, rx, ry, rz);

   T sx, sy, sz, rsun;
   sun_vector(tdays, sx, sy, sz, rsun);

   // spacecraft -> Earth, and spacecraft -> Sun
   T ex = -rx,          ey = -ry,          ez = -rz;
   T ux = rsun*sx + ex, uy = rsun*sy + ey, uz = rsun*sz + ez;
   T ne = sqrt(ex*ex + ey*ey + ez*ez);
   T nu = sqrt(ux*ux + uy*uy + uz*uz);

   T cth = (ex*ux + ey*uy + ez*uz)/(ne*nu);      // cos(theta_es)
   T se  = RE_EARTH/ne;                          // sin(theta_e)
   T ss  = RSUN_KM/nu;                           // sin(theta_s)

   // cos(theta_e) is sqrt(1 - se^2), but the optimiser is free to try states
   // whose radius is below the radius of the Earth, where that argument is
   // negative and the square root is not a number. The smooth positive part
   // below agrees with sqrt(q) to better than one part in 1e8 for the values
   // this problem actually visits, and stays real and differentiable if the
   // iterates stray inside the planet.
   T qe = 1.0 - se*se;
   T ce = sqrt(0.5*(qe + sqrt(qe*qe + 1.0e-8)));
   T qs = 1.0 - ss*ss;
   T cs = sqrt(0.5*(qs + sqrt(qs*qs + 1.0e-8)));

   // geometric margin, in radians: positive inside the shadow
   T ell = (cth - (ce*cs - se*ss))/(se*cs + ce*ss);

   // The smoothed step, written with a hyperbolic tangent rather than as
   // 1/(1+exp(c*ell)). The two are identical in exact arithmetic, but the
   // logistic overflows here: deep in shadow c*ell reaches several hundred, and
   // although 1/(1+exp) then merely underflows to zero, its DERIVATIVE carries
   // exp(c*ell) squared, which runs past the largest representable double and
   // returns a not-a-number to the solver. The tangent saturates instead.
   return 0.5*(1.0 - tanh(0.5*SHADOW_GAIN*ell));
}

adouble sunlit_fraction(adouble* states, adouble& L)
{
   // elapsed time is states[6], carried in days
   return sunlit_core<adouble>(states[0], states[1], states[2], states[3],
                               states[4], L, states[6]);
}

// The same function in double precision, for the guess and the mesh builder.
static double sunlit_d(const double* y, double tdays)
{
   return sunlit_core<double>(y[0], y[1], y[2], y[3], y[4], y[5], tdays);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
   // "time" is the independent variable, which here is the true longitude L.
   adouble L = time;

   adouble p = states[0], f = states[1], g = states[2];
   adouble h = states[3], k = states[4], m = states[5];

   // The thrust direction is parameterised by two angles rather than by three
   // components with a unit-norm path constraint. That removes a constraint,
   // removes a control, and makes the direction a unit vector by construction.
   adouble alpha = controls[0];     // in-plane, from the transverse direction
   adouble beta  = controls[1];     // out-of-plane
   adouble ur = cos(beta)*sin(alpha);
   adouble ut = cos(beta)*cos(alpha);
   adouble un = sin(beta);

   adouble sinL = sin(L), cosL = cos(L);
   adouble w    = 1.0 + f*cosL + g*sinL;
   adouble s2   = 1.0 + h*h + k*k;
   adouble r    = p/w;
   adouble q    = h*sinL - k*cosL;
   adouble sq   = sqrt(p/MU_EARTH);

   // Oblateness, expressed directly in the equinoctial variables, in the
   // radial / transverse / normal frame.
   adouble cj2 = -1.5*MU_EARTH*J2_EARTH*RE_EARTH*RE_EARTH/(r*r*r*r);
   adouble aj2_r = cj2*(1.0 - 12.0*q*q/(s2*s2));
   adouble aj2_t = cj2*( 8.0*q*(h*cosL + k*sinL)/(s2*s2) );
   adouble aj2_n = cj2*( 4.0*(1.0 - h*h - k*k)*q/(s2*s2) );

   // Thrust. The factor 1e-3 converts newtons to kg km/s^2.
   adouble duty = 1.0;
#if USE_ECLIPSE
   duty = sunlit_fraction(states, L);
#endif
   adouble athrust = duty*THRUST_NEWTON*1.0e-3/m;

   adouble Dr = aj2_r + athrust*ur;
   adouble Dt = aj2_t + athrust*ut;
   adouble Dn = aj2_n + athrust*un;

   // Gauss variational equations, and the rate of the true longitude
   adouble dp = 2.0*p/w*sq*Dt;
   adouble df = sq*( Dr*sinL + ((w+1.0)*cosL + f)*Dt/w - g*q*Dn/w );
   adouble dg = sq*(-Dr*cosL + ((w+1.0)*sinL + g)*Dt/w + f*q*Dn/w );
   adouble dh = sq*s2*Dn*cosL/(2.0*w);
   adouble dk = sq*s2*Dn*sinL/(2.0*w);
   adouble dm = -duty*THRUST_NEWTON/(G0*ISP_SECONDS);
   adouble dL = sqrt(MU_EARTH*p)*(w/p)*(w/p) + sq*q*Dn/w;

   // Divide through by dL/dt to change the independent variable to L.
   derivatives[0] = dp/dL;
   derivatives[1] = df/dL;
   derivatives[2] = dg/dL;
   derivatives[3] = dh/dL;
   derivatives[4] = dk/dL;
   derivatives[5] = dm/dL;
   derivatives[6] = 1.0/(dL*86400.0);   // elapsed time, carried in DAYS

   // The orbit must stay elliptic, and must not wander out to an apoapsis the
   // transfer has no use for. Both constraints are on the eccentricity: the
   // second is r_apo = p/(1-e) <= R_APO_MAX, squared to remove the square root,
   // which is legitimate because p is bounded well below R_APO_MAX and both
   // sides are therefore positive.
   adouble e2 = f*f + g*g;
   adouble s  = 1.0 - p/R_APO_MAX;
   path[0] = e2;
   path[1] = e2 - s*s;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   for (int i = 0; i < 7; i++) e[i] = initial_states[i];   // p f g h k m t
   e[7]  = final_states[0];      // semilatus rectum
   e[8]  = final_states[1];      // f  )  together these bound the final
   e[9]  = final_states[2];      // g  )  eccentricity
   e[10] = final_states[3];      // h  )  and these the final inclination
   e[11] = final_states[4];      // k  )
}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // single phase problem
}


//////////////////////////////////////////////////////////////////////////
///////////////////  An initial guess from a steering law    /////////////
//////////////////////////////////////////////////////////////////////////

// A crude guess does not converge on a problem of this kind, so the guess is
// generated here by flying the transfer under a Lyapunov steering law: at each
// instant the thrust is pointed in the direction that most rapidly reduces a
// weighted quadratic measure of the distance of the current elements from the
// target ones. The result is a feasible transfer, not an optimal one -- it is
// some thirty per cent slower than the optimum -- but it is a trajectory of the
// right shape, and that is what the optimiser needs to start from.

static void gauss_rows(const double* y, double G[5][3])
{
   double p=y[0], f=y[1], g=y[2], h=y[3], k=y[4], L=y[5];
   double w  = 1.0 + f*cos(L) + g*sin(L);
   double s2 = 1.0 + h*h + k*k;
   double sq = sqrt(p/MU_EARTH);
   double q  = h*sin(L) - k*cos(L);
   G[0][0]=0.0;           G[0][1]=2.0*p/w*sq;                    G[0][2]=0.0;
   G[1][0]=sq*sin(L);     G[1][1]=sq*((w+1.0)*cos(L)+f)/w;       G[1][2]=-sq*g*q/w;
   G[2][0]=-sq*cos(L);    G[2][1]=sq*((w+1.0)*sin(L)+g)/w;       G[2][2]= sq*f*q/w;
   G[3][0]=0.0;           G[3][1]=0.0;                           G[3][2]= sq*s2*cos(L)/(2.0*w);
   G[4][0]=0.0;           G[4][1]=0.0;                           G[4][2]= sq*s2*sin(L)/(2.0*w);
}

// The steering law drives a Lyapunov function of the orbit elements, and which
// elements it uses matters more than one might expect.
//
// The obvious choice is to drive the semi-latus rectum p to its target. That
// converges, but it produces a transfer with far more revolutions than it needs.
// The orbital period depends on the semi-major axis a = p/(1 - e^2), and at the
// eccentricity of a geostationary transfer orbit those two differ by a factor of
// two, so a law that raises p while the orbit is still eccentric leaves the
// period short and spends the transfer grinding through revolutions near the
// starting semi-major axis. Measured here, targeting p gave a guess of 304
// revolutions where the published optimum takes 162.
//
// Targeting a instead raises the period early, which is what the published
// solution does: its apoapsis goes above the geostationary radius before being
// brought back down. Since the size of the discretised problem is proportional
// to the number of revolutions, this is worth more than any amount of tuning
// further downstream.
static void lyapunov_steer(const double* y, double u[3])
{
   const double a_target = 42165.0;
   const double w_a = LYAP_WA, w_e = 1.0, w_i = LYAP_WI;

   double p = y[0], f = y[1], g = y[2], h = y[3], k = y[4];
   double D = 1.0 - f*f - g*g;
   if (D < 1.0e-6) D = 1.0e-6;
   double a = p/D;

   double da_dp = 1.0/D;                  // da/dp
   double da_df = 2.0*p*f/(D*D);          // da/df
   double da_dg = 2.0*p*g/(D*D);          // da/dg

   double G[5][3]; gauss_rows(y, G);
   double ea = (a - a_target)/a_target;

   double grad[3] = {0.0, 0.0, 0.0};
   for (int j = 0; j < 3; j++) {
      double da = da_dp*G[0][j] + da_df*G[1][j] + da_dg*G[2][j];
      grad[j] = w_a*ea*da/a_target
              + w_e*( f*G[1][j] + g*G[2][j] )
              + w_i*( h*G[3][j] + k*G[4][j] );
   }
   double n = sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);
   if (n > 1.0e-14) { for (int j=0;j<3;j++) u[j] = -grad[j]/n; }
   else             { u[0]=0.0; u[1]=1.0; u[2]=0.0; }
}

static void guess_rates(const double* y, double m, const double* u, double* dy,
                        double duty)
{
   double p=y[0], f=y[1], g=y[2], h=y[3], k=y[4], L=y[5];
   double sinL=sin(L), cosL=cos(L);
   double w=1.0+f*cosL+g*sinL, s2=1.0+h*h+k*k, r=p/w;
   double q=h*sinL-k*cosL, sq=sqrt(p/MU_EARTH);
   double cj2=-1.5*MU_EARTH*J2_EARTH*RE_EARTH*RE_EARTH/(r*r*r*r);
   double a = duty*THRUST_NEWTON*1.0e-3/m;
   double Dr=cj2*(1.0-12.0*q*q/(s2*s2)) + a*u[0];
   double Dt=cj2*(8.0*q*(h*cosL+k*sinL)/(s2*s2)) + a*u[1];
   double Dn=cj2*(4.0*(1.0-h*h-k*k)*q/(s2*s2)) + a*u[2];
   dy[0]=2.0*p/w*sq*Dt;
   dy[1]=sq*( Dr*sinL + ((w+1.0)*cosL+f)*Dt/w - g*q*Dn/w );
   dy[2]=sq*(-Dr*cosL + ((w+1.0)*sinL+g)*Dt/w + f*q*Dn/w );
   dy[3]=sq*s2*Dn*cosL/(2.0*w);
   dy[4]=sq*s2*Dn*sinL/(2.0*w);
   dy[5]=sqrt(MU_EARTH*p)*(w/p)*(w/p) + sq*q*Dn/w;
}

// Builds a multi-interval mesh on (0,1) that is coarse through the smooth arcs
// and clustered on the eclipse transitions. PSOPT accepts such a mesh directly
// through hp_breakpoints and hp_orders, which is what makes this problem
// tractable at all: a mesh uniform in true longitude would need of the order of
// a quarter of a million nodes to see transitions this narrow.
static void build_hp_mesh(const vector<double>& Ltrans, double Lf,
                          double coarse_deg, double window_deg, int nsub,
                          int order_coarse, int order_fine,
                          RowVectorXd& breaks, RowVectorXi& orders)
{
   const double W = window_deg*M_PI/180.0;      // half-width of a fine window
   vector<double> b;                            // interior breakpoints, in L

   double step = coarse_deg*M_PI/180.0;
   for (double L = step; L < Lf - 1.0e-9; L += step) b.push_back(L);

   for (size_t i = 0; i < Ltrans.size(); i++) {
      double Lc = Ltrans[i];
      for (int j = 0; j <= nsub; j++) {
         double Lb = Lc - W + 2.0*W*((double) j)/((double) nsub);
         if (Lb > 1.0e-9 && Lb < Lf - 1.0e-9) b.push_back(Lb);
      }
   }

   sort(b.begin(), b.end());
   // Drop near-duplicates, which would make a degenerate interval. The tolerance
   // has to be small against the fine window, not against the whole transfer:
   // scaled to Lf it would be some ten degrees here and would collapse every
   // cluster back to a single point, quietly undoing the clustering.
   vector<double> c;
   double tol = 0.05*W;
   for (size_t i = 0; i < b.size(); i++)
      if (c.empty() || b[i] - c.back() > tol) c.push_back(b[i]);

   int K = (int) c.size() + 1;                  // number of intervals
   breaks.resize(K - 1);
   orders.resize(K);
   for (int i = 0; i < K - 1; i++) breaks(i) = c[i]/Lf;      // normalise to (0,1)

   // an interval is "fine" if its midpoint lies inside one of the windows
   for (int i = 0; i < K; i++) {
      double a = (i == 0)     ? 0.0 : c[i-1];
      double bb= (i == K - 1) ? Lf  : c[i];
      double mid = 0.5*(a + bb);
      bool fine = false;
      for (size_t j = 0; j < Ltrans.size(); j++)
         if (fabs(mid - Ltrans[j]) < W) { fine = true; break; }
      orders(i) = fine ? order_fine : order_coarse;
   }
}

// Flies the steering law and samples it onto nnodes points.
// Returns the transfer time it achieved.
static double build_guess(int nnodes, MatrixXd& x_guess, MatrixXd& u_guess, MatrixXd& t_guess,
                          const double* y0, double m0, double Lf_target,
                          vector<double>* Ltrans = NULL, double* ecl_frac = NULL)
{
   const double dt = 60.0, tmax = 400.0*86400.0;
   const double mdot = THRUST_NEWTON/(G0*ISP_SECONDS);
   int cap = (int) (tmax/dt) + 2;
   vector<double> th; th.reserve(cap);
   vector< array<double,7> > sh; sh.reserve(cap);
   vector< array<double,3> > uh; uh.reserve(cap);

   double y[6], m = m0, t = 0.0;
   for (int i=0;i<6;i++) y[i]=y0[i];

   while (t < tmax) {
      double u[3]; lyapunov_steer(y, u);
      array<double,7> ss; for(int i=0;i<6;i++) ss[i]=y[i]; ss[6]=m;
      array<double,3> uu = {u[0],u[1],u[2]};
      th.push_back(t); sh.push_back(ss); uh.push_back(uu);

      // The guess must fly the same dynamics the optimiser will see, eclipses
      // included: a guess that thrusts through the shadow is inconsistent with
      // the model by the eclipse fraction, which is enough to start the solver
      // badly on a transfer of this length.
      double duty = 1.0;
#if USE_ECLIPSE
      duty = sunlit_d(y, t/86400.0);
#endif
      double k1[6],k2[6],k3[6],k4[6],yt[6],ut_[3];
      lyapunov_steer(y,ut_);           guess_rates(y,m,ut_,k1,duty);
      for(int i=0;i<6;i++) yt[i]=y[i]+0.5*dt*k1[i];
      lyapunov_steer(yt,ut_);          guess_rates(yt,m,ut_,k2,duty);
      for(int i=0;i<6;i++) yt[i]=y[i]+0.5*dt*k2[i];
      lyapunov_steer(yt,ut_);          guess_rates(yt,m,ut_,k3,duty);
      for(int i=0;i<6;i++) yt[i]=y[i]+dt*k3[i];
      lyapunov_steer(yt,ut_);          guess_rates(yt,m,ut_,k4,duty);
      for(int i=0;i<6;i++) y[i]+= dt/6.0*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
      m -= duty*mdot*dt; t += dt;

      double ecc = sqrt(y[1]*y[1]+y[2]*y[2]);
      double sma = y[0]/(1.0-ecc*ecc);
      double inc = 2.0*atan(sqrt(y[3]*y[3]+y[4]*y[4]))*180.0/M_PI;
      if (fabs(sma-42164.0) < 200.0 && ecc < 2.0e-3 && inc < 0.05) break;
   }
   double tf = t, Lf = y[5];
   // If a revolution count has been imposed, squeeze the same trajectory into it.
   double squeeze = 1.0;
   if (Lf_target > 0.0) { squeeze = Lf_target/Lf; Lf = Lf_target; tf *= squeeze; }

   // Walk the generated trajectory and record where the shadow function crosses
   // one half. These are the places the mesh has to resolve: the transition is
   // barely a degree of true longitude wide, so a mesh that does not know about
   // them will step straight over the eclipses without ever seeing one.
   if (Ltrans || ecl_frac) {
      double shaded = 0.0, span = 0.0, prev = -1.0;
      for (size_t i = 0; i < sh.size(); i++) {
         double yi[6] = { sh[i][0], sh[i][1], sh[i][2], sh[i][3], sh[i][4], sh[i][5] };
         double psi = sunlit_d(yi, squeeze*th[i]/86400.0);
         shaded += (1.0 - psi)*dt;   span += dt;
         if (Ltrans && i > 0 && (prev - 0.5)*(psi - 0.5) < 0.0) {
            double fr = (0.5 - prev)/(psi - prev);       // in [0,1]
            Ltrans->push_back(squeeze*(sh[i-1][5] + fr*(sh[i][5] - sh[i-1][5])));
         }
         prev = psi;
      }
      if (ecl_frac) *ecl_frac = (span > 0.0) ? shaded/span : 0.0;
   }

   x_guess = zeros(7, nnodes);
   u_guess = zeros(2, nnodes);
   t_guess = linspace(0.0, Lf, nnodes);      // the independent variable is L

   // sh[] holds (p,f,g,h,k,L,m) at each step; resample it uniformly in L
   double prev_alpha = 0.0;
   size_t idx = 0;
   for (int j = 0; j < nnodes; j++) {
      double Lj = Lf*((double) j)/((double)(nnodes-1));
      double Lsrc = Lj/squeeze;      // the same point on the generated trajectory
      while (idx + 2 < sh.size() && sh[idx+1][5] < Lsrc) idx++;
      double L0_ = sh[idx][5], L1_ = sh[idx+1][5];
      double frac = (L1_ > L0_) ? (Lsrc - L0_)/(L1_ - L0_) : 0.0;
      if (frac < 0.0) frac = 0.0; if (frac > 1.0) frac = 1.0;
      for (int i = 0; i < 5; i++)                     // p f g h k
         x_guess(i,j) = sh[idx][i] + frac*(sh[idx+1][i] - sh[idx][i]);
      x_guess(5,j) = sh[idx][6] + frac*(sh[idx+1][6] - sh[idx][6]);      // mass
      x_guess(6,j) = squeeze*(th[idx] + frac*(th[idx+1] - th[idx]))/86400.0; // elapsed time, days

      double ur = uh[idx][0] + frac*(uh[idx+1][0]-uh[idx][0]);
      double ut = uh[idx][1] + frac*(uh[idx+1][1]-uh[idx][1]);
      double un = uh[idx][2] + frac*(uh[idx+1][2]-uh[idx][2]);
      double nn = sqrt(ur*ur+ut*ut+un*un); ur/=nn; ut/=nn; un/=nn;
      double beta  = asin(un);
      double alpha = atan2(ur, ut);
      while (alpha - prev_alpha >  M_PI) alpha -= 2.0*M_PI;
      while (alpha - prev_alpha < -M_PI) alpha += 2.0*M_PI;
      prev_alpha = alpha;
      u_guess(0,j) = alpha;
      u_guess(1,j) = beta;
   }
   return tf;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define main routine ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    // With an argument, the number of revolutions is fixed at that value and the
    // problem becomes minimum time for that many revolutions. Sweeping it traces
    // the minimum time against revolution count, whose least value is the answer
    // being sought; with no argument the revolution count is free and the solver
    // finds only the local minimum nearest its initial guess.
    double nrev_fixed = (argc > 1) ? atof(argv[1]) : 0.0;
    bool guess_only = (argc > 3 && string(argv[3]) == "guess");
    if (argc > 2) POWER_KW = atof(argv[2]);
    if (argc > 4) LYAP_WA = atof(argv[4]);
    if (argc > 5) LYAP_WI = atof(argv[5]);
    THRUST_NEWTON = 2.0*THRUSTER_EFF*POWER_KW*1.0e3/(G0*ISP_SECONDS);
    printf("\nthruster: %.1f kW at %.0f per cent efficiency, Isp %.0f s"
           "  ->  %.6f N\n", POWER_KW, 100.0*THRUSTER_EFF, ISP_SECONDS, THRUST_NEWTON);
    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        = "Minimum time GTO to GEO orbit raising";
    problem.outfilename = "gto_geo.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   = 1;
    problem.nlinkages = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates     = 7;
    problem.phases(1).ncontrols   = 2;
    problem.phases(1).nevents     = 12;
    problem.phases(1).npath       = 2;
    problem.phases(1).nodes       << 1000;

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    // the geostationary transfer orbit: 185.5 km by 35786.2 km at 28.5 degrees
    double p0 = 11359.07;      // km
    double f0 = 0.7306;
    double g0_ = 0.0;
    double h0 = 0.2539676;     // tan(i/2) with i = 28.5 deg, right ascension zero
    double k0 = 0.0;
    double L0 = 0.0;           // released at perigee

    double pf = 42165.0;       // km, geostationary
    double e_tol = 1.0e-3;                 // final eccentricity
    double i_tol = tan(0.5*0.01*M_PI/180.0);   // final inclination, 0.01 degrees

    problem.phases(1).bounds.lower.states << 10000.0, -1.0, -1.0, -1.0, -1.0, 600.0,          0.0;
    problem.phases(1).bounds.upper.states << 46000.0,  1.0,  1.0,  1.0,  1.0, MASS_INITIAL, 120.0;

    problem.phases(1).bounds.lower.controls << -4.0*M_PI, -M_PI/2.0;
    problem.phases(1).bounds.upper.controls <<  4.0*M_PI,  M_PI/2.0;

    problem.phases(1).bounds.lower.path << 0.0, -1.0;
    problem.phases(1).bounds.upper.path << E_MAX*E_MAX, 0.0;

    problem.phases(1).bounds.lower.events
        << p0, f0, g0_, h0, k0, MASS_INITIAL, 0.0,
           pf, -e_tol/sqrt(2.0), -e_tol/sqrt(2.0), -i_tol/sqrt(2.0), -i_tol/sqrt(2.0);
    problem.phases(1).bounds.upper.events
        << p0, f0, g0_, h0, k0, MASS_INITIAL, 0.0,
           pf,  e_tol/sqrt(2.0),  e_tol/sqrt(2.0),  i_tol/sqrt(2.0),  i_tol/sqrt(2.0);

    // the independent variable is the true longitude, in radians
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;

    // set below, once the guess has been flown and its own revolution count is
    // known -- a fixed floor of twenty revolutions is wrong at every power but
    // the smallest
    problem.phases(1).bounds.lower.EndTime = 0.0;
    problem.phases(1).bounds.upper.EndTime = 200.0*2.0*M_PI;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes = 1000;

    MatrixXd x_guess, u_guess, t_guess;
    vector<double> Ltrans;
    double ecl_frac = 0.0;
    double y0[6] = {p0, f0, g0_, h0, k0, L0};
    double tf_guess = build_guess(nnodes, x_guess, u_guess, t_guess, y0, MASS_INITIAL,
                                  nrev_fixed > 0.0 ? nrev_fixed*2.0*M_PI : 0.0,
                                  &Ltrans, &ecl_frac);

    double Lf_guess = t_guess(0, nnodes-1);
    printf("\n");
    printf("guess: %.2f days over %.1f revolutions\n",
           tf_guess/86400.0, Lf_guess/(2.0*M_PI));
#if USE_ECLIPSE
    printf("eclipse model on:\n");
    printf("   shadow transitions found along the guess   %d\n", (int) Ltrans.size());
    printf("   eclipse fraction of the guess              %.4f\n", ecl_frac);

    // An eclipse-aware mesh, clustered on the transitions.
    build_hp_mesh(Ltrans, Lf_guess, MESH_COARSE_DEG, MESH_WINDOW_DEG,
                  MESH_SUBDIV, MESH_ORDER_COARSE, MESH_ORDER_FINE,
                  problem.phases(1).hp_breakpoints, problem.phases(1).hp_orders);
    {
        int K = (int) problem.phases(1).hp_orders.size();
        long nodes_total = 0;
        int nfine = 0;
        for (int i = 0; i < K; i++) {
            nodes_total += problem.phases(1).hp_orders(i);
            if (problem.phases(1).hp_orders(i) == MESH_ORDER_FINE) nfine++;
        }
        printf("   mesh intervals                             %d  (%d fine)\n", K, nfine);
        printf("   collocation nodes                          %ld\n", nodes_total);
    }
#endif
    // The final true longitude. With the eclipse model on it is FIXED, and that
    // is not a detail: the clustered mesh is anchored to the true longitudes at
    // which the guess crosses the shadow boundary, and PSOPT's breakpoints are
    // fractions of the phase. Were the final longitude free, those fractions
    // would stretch with it and the fine windows would slide off the eclipses --
    // not by a little, but by whole revolutions, since ten per cent of the final
    // longitude is four revolutions at twenty kilowatts. The problem solved is
    // therefore minimum time for a given number of revolutions; sweeping that
    // number traces the envelope whose least value is the answer sought.
    //
    // What does still move is the elapsed time, and with it the Sun: the guess
    // is a few per cent slow, so the solution arrives at each revolution at a
    // slightly different date and its shadow crossings sit a little away from
    // the guess's. Whether that little is smaller than the half-width of the
    // fine windows is not something to assume, so the example measures it after
    // the solve and reports it.
#if USE_ECLIPSE
    double Lf_fixed = (nrev_fixed > 0.0) ? nrev_fixed*2.0*M_PI : Lf_guess;
    problem.phases(1).bounds.lower.EndTime = Lf_fixed;
    problem.phases(1).bounds.upper.EndTime = Lf_fixed;
#else
    if (nrev_fixed > 0.0) {
        problem.phases(1).bounds.lower.EndTime = nrev_fixed*2.0*M_PI;
        problem.phases(1).bounds.upper.EndTime = nrev_fixed*2.0*M_PI;
    } else {
        problem.phases(1).bounds.lower.EndTime = 0.60*Lf_guess;
        problem.phases(1).bounds.upper.EndTime = 1.60*Lf_guess;
    }
#endif

    printf("\n");
    printf("steering-law guess: transfer time %.3f days, %.1f revolutions\n",
           tf_guess/86400.0, t_guess(0,nnodes-1)/(2.0*M_PI));
    {
       double emax = 0.0, ramax = 0.0;
       for (int j = 0; j < nnodes; j++) {
          double fj = x_guess(1,j), gj = x_guess(2,j);
          double ej = sqrt(fj*fj + gj*gj);
          double ra = x_guess(0,j)/(1.0 - ej);
          if (ej > emax) emax = ej;
          if (ra > ramax) ramax = ra;
       }
       printf("   peak eccentricity %.4f, peak apoapsis %.0f km\n", emax, ramax);
    }
    if (guess_only) return 0;

    MatrixXd param_guess;
    auto_phase_guess(problem, u_guess, x_guess, param_guess, t_guess);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////

    algorithm.nlp_iter_max        = 5000;
    algorithm.nlp_tolerance       = 1.e-6;
    algorithm.nlp_method          = "IPOPT";
    algorithm.scaling             = "automatic";
    algorithm.derivatives         = "automatic";
    algorithm.defect_scaling      = "jacobian-based";
#if USE_ECLIPSE
    // An explicit multi-interval mesh requires a pseudospectral method. Radau is
    // the natural choice: it collocates at one endpoint, which suits a problem
    // whose intervals are joined end to end, and its costate map is the accurate
    // one of the four (see the comparison in Chapter 5 of the book).
    algorithm.collocation_method  = "Radau";
    algorithm.mesh_refinement     = "manual";
    algorithm.mr_max_iterations   = 1;
#else
    algorithm.collocation_method  = "Hermite-Simpson";
    algorithm.mesh_refinement     = "automatic";
    algorithm.mr_max_iterations   = 3;
#endif
    // Three refinements take the mesh from 1000 nodes to 1413 and the maximum
    // discretisation error from 4.2e-3 to 1.4e-3, still short of the tolerance
    // asked for here. The transfer time moves from 27.887 to 27.843 days across
    // them, so quote it to three significant figures and no more. Relaxing the
    // tolerance instead makes the refinement lazy: at 1e-3 it adds only a
    // hundred nodes and leaves the error at 4.5e-3, above the value requested.
    algorithm.ode_tolerance       = 1.e-4;

    // PSOPT gives IPOPT an hour of processor time by default. That is ample for
    // every other example in this distribution and nowhere near enough for this
    // one: the twenty kilowatt case has some 28000 variables and takes rather
    // longer than that. A run that ends with IPOPT return code -4 has hit this
    // limit and not converged, which is easy to mistake for a failure of the
    // formulation.
    algorithm.ipopt_max_cpu_time  = 6.0*3600.0;

////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    int rc = psopt(solution, problem, algorithm);
    if (rc != 0) printf("psopt returned %d\n", rc);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x = solution.get_states_in_phase(1);
    MatrixXd u = solution.get_controls_in_phase(1);
    MatrixXd t = solution.get_time_in_phase(1);

    int N = (int) x.cols();
    double Lfin = t(0, N-1);                 // independent variable at the end
    double pfin = x(0, N-1);
    double ffin = x(1, N-1), gfin = x(2, N-1);
    double hfin = x(3, N-1), kfin = x(4, N-1);
    double mfin = x(5, N-1), tfin = x(6, N-1);
    double efin = sqrt(ffin*ffin + gfin*gfin);
    double afin = pfin/(1.0 - efin*efin);
    double ifin = 2.0*atan(sqrt(hfin*hfin + kfin*kfin))*180.0/M_PI;

    printf("\n");
    printf("transfer time      %12.4f days\n", tfin);
    printf("revolutions        %12.2f\n",      Lfin/(2.0*M_PI));
    printf("propellant used    %12.4f kg\n",   MASS_INITIAL - mfin);
    printf("final semi-major   %12.4f km\n",   afin);
    printf("final eccentricity %12.3e\n",      efin);
    printf("final inclination  %12.3e deg\n",  ifin);

    // How near the solution came to the two safeguards. If either is close to
    // its bound the bound is shaping the answer and should be re-examined; if
    // both are comfortably inside, as they should be, they did nothing but keep
    // the iterates out of the region where the element set stops meaning
    // anything, which is what they are for.
    {
       double emax = 0.0, ramax = 0.0;
       for (int j = 0; j < N; j++) {
          double ej = sqrt(x(1,j)*x(1,j) + x(2,j)*x(2,j));
          double ra = x(0,j)/(1.0 - ej);
          if (ej > emax) emax = ej;
          if (ra > ramax) ramax = ra;
       }
       printf("peak eccentricity  %12.4f   (bound %.2f)\n", emax, E_MAX);
       printf("peak apoapsis      %12.0f km (bound %.0f km)\n", ramax, R_APO_MAX);
    }

#if USE_ECLIPSE
    // Did the fine windows stay over the eclipses? The mesh was clustered on the
    // true longitudes at which the GUESS crossed the shadow boundary, and the
    // solution takes a different time to get to each revolution, so the Sun is
    // somewhere slightly different when it arrives and the crossing moves. If
    // that movement exceeds the half-width of the windows they are no longer
    // doing their job, and the transitions are being stepped over.
    {
       double worst = 0.0;
       int nfound = 0, noutside = 0;
       double prev = -1.0;
       for (int j = 0; j < N; j++) {
          double y[6] = { x(0,j), x(1,j), x(2,j), x(3,j), x(4,j), t(0,j) };
          double psi = sunlit_d(y, x(6,j));
          if (j > 0 && (prev - 0.5)*(psi - 0.5) < 0.0) {
             double fr = (0.5 - prev)/(psi - prev);
             double Lc = t(0,j-1) + fr*(t(0,j) - t(0,j-1));
             double best = 1.0e30;
             for (size_t q = 0; q < Ltrans.size(); q++)
                if (fabs(Lc - Ltrans[q]) < best) best = fabs(Lc - Ltrans[q]);
             nfound++;
             if (best > worst) worst = best;
             if (best > MESH_WINDOW_DEG*M_PI/180.0) noutside++;
          }
          prev = psi;
       }
       printf("shadow crossings found along the solution  %d (guess had %d)\n",
              nfound, (int) Ltrans.size());
       printf("   furthest one has moved   %8.3f deg of true longitude"
              "  (window half-width %.1f deg)\n", worst*180.0/M_PI, MESH_WINDOW_DEG);
       printf("   crossings now outside a fine window      %d\n", noutside);
    }
#endif
    printf("\n");

    Save(x, "x.dat");
    Save(u, "u.dat");
    Save(t, "t.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results                                  ////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd days = x.row(6);

    plot(days, x.row(0), problem.name + ": semilatus rectum", "time (days)", "p (km)", "p");
    plot(days, x.row(5), problem.name + ": mass", "time (days)", "m (kg)", "m");
    plot(days, u, problem.name + ": steering angles", "time (days)", "rad", "alpha beta");

    return 0;
}
