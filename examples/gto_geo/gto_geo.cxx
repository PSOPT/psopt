//////////////////////////////////////////////////////////////////////////
////////////////            gto_geo.cxx              /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Minimum time GTO to GEO orbit raising      ///////////////
////////        for an all-electric satellite              ///////////////
//////// Last modified: 29 August 2026                     ///////////////
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

// 20 kW at 55 per cent efficiency and 1800 s of specific impulse.
// The published benchmark is 5 kW, giving 0.31158 N; everything else is the same.
#define THRUST_NEWTON  1.24632
#define ISP_SECONDS    1800.0
#define MASS_INITIAL   1200.0         // kg

// The eclipse model is off by default. The book's results are the eclipse-free
// ones; the model is here so that the published benchmark, which does include
// eclipsing, can be reproduced as a check on the dynamics.
#define USE_ECLIPSE    0

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

// A cylindrical shadow, smoothed so that the transition is differentiable.
// Returns approximately 1 in sunlight and 0 in shadow. The Sun direction is
// taken in the equatorial plane, rotating at the mean rate of the Earth about
// the Sun, which is accurate enough over a transfer of a few weeks.

adouble sunlit_fraction(adouble* states, adouble& L)
{
   adouble p = states[0], f = states[1], g = states[2];
   adouble h = states[3], k = states[4];
   adouble time = states[6]*86400.0;         // elapsed time is a state, in days

   adouble s2 = 1.0 + h*h + k*k;
   adouble w  = 1.0 + f*cos(L) + g*sin(L);
   adouble r  = p/w;
   adouble a2 = h*h - k*k;

   // inertial position
   adouble rx = r/s2*(cos(L) + a2*cos(L) + 2.0*h*k*sin(L));
   adouble ry = r/s2*(sin(L) - a2*sin(L) + 2.0*h*k*cos(L));
   adouble rz = 2.0*r/s2*(h*sin(L) - k*cos(L));

   double  wsun = 2.0*M_PI/(365.25*86400.0);       // rad/s
   adouble sx = cos(wsun*time), sy = sin(wsun*time), sz = 0.0;

   adouble proj = rx*sx + ry*sy + rz*sz;            // along the Sun direction
   adouble d2   = rx*rx + ry*ry + rz*rz - proj*proj;   // perpendicular distance^2

   // in shadow when proj < 0 and d < RE
   double  sharp = 40.0;
   adouble behind = 0.5*(1.0 + tanh(-sharp*proj/RE_EARTH));
   adouble inside = 0.5*(1.0 + tanh(sharp*(RE_EARTH*RE_EARTH - d2)/(RE_EARTH*RE_EARTH)));
   return 1.0 - behind*inside;
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

static void lyapunov_steer(const double* y, double u[3])
{
   static const double target[5] = {42165.0, 0.0, 0.0, 0.0, 0.0};
   static const double weight[5] = {1.0, 1.0, 1.0, 30.0, 30.0};
   static const double scale [5] = {1.0/42165.0, 1.0, 1.0, 1.0, 1.0};
   double G[5][3]; gauss_rows(y, G);
   double grad[3] = {0.0, 0.0, 0.0};
   for (int i = 0; i < 5; i++) {
      double d = (y[i] - target[i])*scale[i];
      for (int j = 0; j < 3; j++) grad[j] += weight[i]*d*G[i][j]*scale[i];
   }
   double n = sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);
   if (n > 1.0e-14) { for (int j=0;j<3;j++) u[j] = -grad[j]/n; }
   else             { u[0]=0.0; u[1]=1.0; u[2]=0.0; }
}

static void guess_rates(const double* y, double m, const double* u, double* dy)
{
   double p=y[0], f=y[1], g=y[2], h=y[3], k=y[4], L=y[5];
   double sinL=sin(L), cosL=cos(L);
   double w=1.0+f*cosL+g*sinL, s2=1.0+h*h+k*k, r=p/w;
   double q=h*sinL-k*cosL, sq=sqrt(p/MU_EARTH);
   double cj2=-1.5*MU_EARTH*J2_EARTH*RE_EARTH*RE_EARTH/(r*r*r*r);
   double a = THRUST_NEWTON*1.0e-3/m;
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

// Flies the steering law and samples it onto nnodes points.
// Returns the transfer time it achieved.
static double build_guess(int nnodes, MatrixXd& x_guess, MatrixXd& u_guess, MatrixXd& t_guess,
                          const double* y0, double m0, double Lf_target)
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

      double k1[6],k2[6],k3[6],k4[6],yt[6],ut_[3];
      lyapunov_steer(y,ut_);           guess_rates(y,m,ut_,k1);
      for(int i=0;i<6;i++) yt[i]=y[i]+0.5*dt*k1[i];
      lyapunov_steer(yt,ut_);          guess_rates(yt,m,ut_,k2);
      for(int i=0;i<6;i++) yt[i]=y[i]+0.5*dt*k2[i];
      lyapunov_steer(yt,ut_);          guess_rates(yt,m,ut_,k3);
      for(int i=0;i<6;i++) yt[i]=y[i]+dt*k3[i];
      lyapunov_steer(yt,ut_);          guess_rates(yt,m,ut_,k4);
      for(int i=0;i<6;i++) y[i]+= dt/6.0*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
      m -= mdot*dt; t += dt;

      double ecc = sqrt(y[1]*y[1]+y[2]*y[2]);
      double sma = y[0]/(1.0-ecc*ecc);
      double inc = 2.0*atan(sqrt(y[3]*y[3]+y[4]*y[4]))*180.0/M_PI;
      if (fabs(sma-42164.0) < 200.0 && ecc < 2.0e-3 && inc < 0.05) break;
   }
   double tf = t, Lf = y[5];
   // If a revolution count has been imposed, squeeze the same trajectory into it.
   double squeeze = 1.0;
   if (Lf_target > 0.0) { squeeze = Lf_target/Lf; Lf = Lf_target; tf *= squeeze; }

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
    problem.phases(1).npath       = 0;
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

    problem.phases(1).bounds.lower.events
        << p0, f0, g0_, h0, k0, MASS_INITIAL, 0.0,
           pf, -e_tol/sqrt(2.0), -e_tol/sqrt(2.0), -i_tol/sqrt(2.0), -i_tol/sqrt(2.0);
    problem.phases(1).bounds.upper.events
        << p0, f0, g0_, h0, k0, MASS_INITIAL, 0.0,
           pf,  e_tol/sqrt(2.0),  e_tol/sqrt(2.0),  i_tol/sqrt(2.0),  i_tol/sqrt(2.0);

    // the independent variable is the true longitude, in radians
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;

    if (nrev_fixed > 0.0) {
        problem.phases(1).bounds.lower.EndTime = nrev_fixed*2.0*M_PI;
        problem.phases(1).bounds.upper.EndTime = nrev_fixed*2.0*M_PI;
    } else {
        problem.phases(1).bounds.lower.EndTime =  20.0*2.0*M_PI;
        problem.phases(1).bounds.upper.EndTime = 150.0*2.0*M_PI;
    }

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
    double y0[6] = {p0, f0, g0_, h0, k0, L0};
    double tf_guess = build_guess(nnodes, x_guess, u_guess, t_guess, y0, MASS_INITIAL,
                                  nrev_fixed > 0.0 ? nrev_fixed*2.0*M_PI : 0.0);
    printf("steering-law guess: transfer time %.3f days, %.1f revolutions\n",
           tf_guess/86400.0, t_guess(0,nnodes-1)/(2.0*M_PI));

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
    algorithm.collocation_method  = "Hermite-Simpson";
    algorithm.mesh_refinement     = "automatic";
    algorithm.mr_max_iterations   = 3;
    // Three refinements take the mesh from 1000 nodes to 1413 and the maximum
    // discretisation error from 4.2e-3 to 1.4e-3, still short of the tolerance
    // asked for here. The transfer time moves from 27.887 to 27.843 days across
    // them, so quote it to three significant figures and no more. Relaxing the
    // tolerance instead makes the refinement lazy: at 1e-3 it adds only a
    // hundred nodes and leaves the error at 4.5e-3, above the value requested.
    algorithm.ode_tolerance       = 1.e-4;

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
