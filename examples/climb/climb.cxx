//////////////////////////////////////////////////////////////////////////
//////////////////            climb.cxx              /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Minimum time to climb of a supersonic aircraft ///////////
//////// Last modified: 12 January 2009                   ////////////////
//////// Reference:     GPOPS Manual			  ////////////////
//////// (See PSOPT handbook for full reference)           ///////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

using namespace PSOPT;




//////////////////////////////////////////////////////////////////////////
/////////  Declare an auxiliary structure to hold local constants ////////
//////////////////////////////////////////////////////////////////////////

struct Constants {
  double g0;
  double S;
  double Re;
  double Isp;
  double mu;
  MatrixXd* CLa_table;
  MatrixXd* CD0_table;
  MatrixXd* eta_table;
  MatrixXd* T_table;
  MatrixXd* M1;
  MatrixXd* M2;
  MatrixXd* h1;
  MatrixXd* htab;
  MatrixXd* ttab;
  MatrixXd* ptab;
  MatrixXd* gtab;
};

typedef struct Constants Constants_;

void atmosphere(adouble* alt,adouble* sigma,adouble* delta,adouble* theta, Constants_& CONSTANTS);

void atmosphere_model(adouble* rho, adouble* M, adouble v, adouble h, Constants_& CONSTANTS);


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
	return tf;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)

{
    return  0.0;
}



//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

  Constants_& CONSTANTS = *( (Constants_ *) workspace->problem->user_data );

  adouble alpha = controls[ 0]; // Angle of attack (rad)


  adouble h     =   states[ 0 ]; // Altitude (ft)
  adouble v     =   states[ 1 ]; // Velocity (ft/s)
  adouble gamma =   states[ 2 ]; // Flight path angle (rad)
  adouble w     =   states[ 3 ]; // weight (lb)


  double g0   = CONSTANTS.g0;
  double S    = CONSTANTS.S;
  double Re   = CONSTANTS.Re;
  double Isp  = CONSTANTS.Isp;
  double mu   = CONSTANTS.mu;

  MatrixXd& M1         = *CONSTANTS.M1;
  MatrixXd& M2         = *CONSTANTS.M2;
  MatrixXd& h1         = *CONSTANTS.h1;
  MatrixXd& CLa_table  = *CONSTANTS.CLa_table;
  MatrixXd& CD0_table  = *CONSTANTS.CD0_table;
  MatrixXd& eta_table  = *CONSTANTS.eta_table;
  MatrixXd& T_table    = *CONSTANTS.T_table;

  int lM1 = length(M1);

  adouble rho;
  adouble   m = w/g0;
  adouble   M;

  atmosphere_model( &rho, &M, v, h, CONSTANTS);

  adouble CL_a, CD0, eta, T;


  spline_interpolation( &CL_a, M, M1, CLa_table,  lM1);
  spline_interpolation( &CD0,  M, M1, CD0_table, lM1);
  spline_interpolation( &eta,  M, M1, eta_table, lM1);
  spline_2d_interpolation(&T, M, h, M2, h1, T_table, workspace);

//  smooth_linear_interpolation( &CL_a, M, M1, CLa_table,  lM1);
//  smooth_linear_interpolation( &CD0,  M, M1, CD0_table, lM1);
//  smooth_linear_interpolation( &eta,  M, M1, eta_table, lM1);
//  smooth_bilinear_interpolation(&T, M, h, M2, h1, T_table);


//  linear_interpolation( &CL_a, M, M1, CLa_table,  lM1);
//  linear_interpolation( &CD0,  M, M1, CD0_table, lM1);
//  linear_interpolation( &eta,  M, M1, eta_table, lM1);
//  bilinear_interpolation(&T, M, h, M2, h1, T_table);


   adouble CL = CL_a*alpha;
   adouble CD = CD0 + eta*CL_a*alpha*alpha;

   adouble D = 0.5*CD*S*rho*v*v;
   adouble L = 0.5*CL*S*rho*v*v;


   adouble hdot     = v*sin(gamma);
   adouble vdot     = 1.0/m*(T*cos(alpha)-D) - mu/pow(Re+h,2.0)*sin(gamma);
   adouble gammadot = (1.0/(m*v))*(T*sin(alpha)+L) + cos(gamma)*(v/(Re+h)-mu/(v*pow(Re+h,2.0)));
   adouble wdot     = -T/Isp;



   derivatives[ 0 ] = hdot;
   derivatives[ 1 ] = vdot;
   derivatives[ 2 ] = gammadot;
   derivatives[ 3 ] = wdot;


}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{



      adouble h0     = initial_states[0];
      adouble v0     = initial_states[1];
      adouble gamma0 = initial_states[2];
      adouble w0     = initial_states[3];

      adouble hf     = final_states[0];
      adouble vf     = final_states[1];
      adouble gammaf = final_states[2];

      e[ 0 ] = h0;
      e[ 1 ] = v0;
      e[ 2 ] = gamma0;
      e[ 3 ] = w0;
      e[ 4 ] = hf;
      e[ 5 ] = vf;
      e[ 6 ] = gammaf;

}


///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{

  // Single phase problem

}



////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{
////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name = "Minimum time to climb for a supersonic aircraft";
    problem.outfilename = "climb.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   = 4;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 7;  
    problem.phases(1).npath     = 0;

    problem.phases(1).nodes     = (RowVectorXi(1) << 30).finished();


    psopt_level2_setup(problem, algorithm);
    
////////////////////////////////////////////////////////////////////////////
///////////////////  Declare an instance of Constants structure /////////////
////////////////////////////////////////////////////////////////////////////

    
    Constants_ CONSTANTS;    
    
    problem.user_data = (void*) &CONSTANTS;

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare MatrixXd objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x, u, t, H;



////////////////////////////////////////////////////////////////////////////
///////////////////  Initialize CONSTANTS and //////////////////////////////
///////////////////  declare local variables  //////////////////////////////
////////////////////////////////////////////////////////////////////////////

   CONSTANTS.g0   = 32.174; // ft/s^2
   CONSTANTS.S    = 530.0; // ft^2
   CONSTANTS.Re   = 20902900.0; // ft
   CONSTANTS.Isp  = 1600.00; //s
   CONSTANTS.mu   = 0.14076539E17; // ft^3/s^2


   MatrixXd M1(1,9);
   M1 << 0.E0, .4E0, .8E0, .9E0, 1.E0, 1.2E0, 1.4E0, 1.6E0, 1.8E0;

   MatrixXd M2(1,10);
   M2 << 0.E0, .2E0, .4E0, .6E0, .8E0, 1.E0, 1.2E0, 1.4E0, 1.6E0, 1.8E0;

   MatrixXd h1(1,10);
   h1 << 0.E0, 5E3, 10.E3, 15.E3, 20.E3, 25.E3, 30.E3, 40.E3, 50.E3, 70.E3;

   MatrixXd CLa_table(1,9);
   CLa_table << 3.44E0, 3.44E0, 3.44E0, 3.58E0, 4.44E0, 3.44E0, 3.01E0, 2.86E0, 2.44E0;

   MatrixXd CD0_table(1,9); 
   CD0_table << .013E0, .013E0, .013E0, .014E0, .031E0, 0.041E0, .039E0, .036E0, .035E0;

   MatrixXd eta_table(1,9); 
   eta_table << .54E0, .54E0, .54E0, .75E0, .79E0, .78E0, .89E0, .93E0, .93E0;

   MatrixXd T_table(10,10);
   T_table << 24200., 24000.,  20300.,  17300.,14500.,12200.,10200.,5700.,3400.,100.,
   28000., 24600.,  21100.,  18100.,15200.,12800.,10700.,6500.,3900.,200.,
   28300., 25200.,  21900.,  18700.,15900.,13400.,11200.,7300.,4400.,400.,
   30800., 27200.,  23800.,  20500.,17300.,14700.,12300.,8100.,4900.,800.,
   34500., 30300.,  26600.,  23200.,19800.,16800.,14100.,9400.,5600.,1100.,
   37900., 34300.,  30400.,  26800.,23300.,19800.,16800.,11200.,6800.,1400.,
   36100., 38000.,  34900.,  31300.,27300.,23600.,20100.,13400.,8300.,1700.,
   34300., 36600.,  38500.,  36100.,31600.,28100.,24200.,16200.,10000.,2200.,
   32500., 35200.,  42100.,  38700.,35700.,32000.,28100.,19300.,11900.,2900.,
   30700., 33800.,  45700.,  41300.,39800.,34600.,31100.,21700.,13300.,3100. ;

   MatrixXd htab(1,8);
   htab <<      0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852;
   MatrixXd ttab(1,8);
   ttab <<      288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946;
   MatrixXd ptab(1,8);
   ptab <<      1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 1.0945601E-3,
                                     6.6063531E-4, 3.9046834E-5, 3.68501E-6;
   MatrixXd gtab(1,8);
   gtab << -6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0;

//   M1.Print("M1");

//   M2.Print("M2");

//   h1.Print("h1");

//   CLa_table.Print("CLa_table");

//   CD0_table.Print("CD0_table");

//   eta_table.Print("eta_table");

//   T_table.Print("T_table");


   CONSTANTS.M1         = &M1;
   CONSTANTS.M2         = &M2;
   CONSTANTS.h1         = &h1;
   CONSTANTS.CLa_table  = &CLa_table;
   CONSTANTS.CD0_table  = &CD0_table;
   CONSTANTS.eta_table  = &eta_table;
   CONSTANTS.T_table    = &T_table;
   CONSTANTS.htab       = &htab;
   CONSTANTS.ttab       = &ttab;
   CONSTANTS.ptab       = &ptab;
   CONSTANTS.gtab       = &gtab;

   double h0     = 0.0;
   double hf     = 65600.0;
   double v0     = 424.26;
   double vf     = 968.148;
   double gamma0 = 0.0;
   double gammaf = 0.0;
   double w0     = 42000.0;

   double hmin =  0;
   double hmax =  69000.0;
   double vmin =  1.0;
   double vmax =  2000.0;
   double gammamin = -89.0*pi/180.0; // -89.0*pi/180.0;
   double gammamax =  89.0*pi/180.0; //  89.0*pi/180.0;
   double wmin     = 0.0;
   double wmax     = 45000.0;
   double alphamin = -20.0*pi/180.0;
   double alphamax =  20.0*pi/180.0;

   double t0min = 0.0;
   double t0max = 0.0;
   double tfmin = 200.0;
   double tfmax = 500.0;


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    int iphase =  1;

    problem.phases(iphase).bounds.lower.StartTime    = t0min;
    problem.phases(iphase).bounds.upper.StartTime    = t0max;

    problem.phases(iphase).bounds.lower.EndTime      = tfmin;
    problem.phases(iphase).bounds.upper.EndTime      = tfmax;


    problem.phases(iphase).bounds.lower.states(0) = hmin;
    problem.phases(iphase).bounds.upper.states(0) = hmax;
    problem.phases(iphase).bounds.lower.states(1) = vmin;
    problem.phases(iphase).bounds.upper.states(1) = vmax;
    problem.phases(iphase).bounds.lower.states(2) = gammamin;
    problem.phases(iphase).bounds.upper.states(2) = gammamax;
    problem.phases(iphase).bounds.lower.states(3) = wmin;
    problem.phases(iphase).bounds.upper.states(3) = wmax;


    problem.phases(iphase).bounds.lower.controls(0) = alphamin;
    problem.phases(iphase).bounds.upper.controls(0) = alphamax;

    // The following bounds fix the initial and final state conditions

    problem.phases(iphase).bounds.lower.events(0) = h0;
    problem.phases(iphase).bounds.upper.events(0) = h0;
    problem.phases(iphase).bounds.lower.events(1) = v0;
    problem.phases(iphase).bounds.upper.events(1) = v0;
    problem.phases(iphase).bounds.lower.events(2) = gamma0;
    problem.phases(iphase).bounds.upper.events(2) = gamma0;
    problem.phases(iphase).bounds.lower.events(3) = w0;
    problem.phases(iphase).bounds.upper.events(3) = w0;
    problem.phases(iphase).bounds.lower.events(4) = hf;
    problem.phases(iphase).bounds.upper.events(4) = hf;
    problem.phases(iphase).bounds.lower.events(5) = vf;
    problem.phases(iphase).bounds.upper.events(5) = vf;
    problem.phases(iphase).bounds.lower.events(6) = gammaf;
    problem.phases(iphase).bounds.upper.events(6) = gammaf;


////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes = problem.phases(iphase).nodes(0);


    MatrixXd stateGuess(4,nnodes);

    stateGuess.row(0) = linspace(h0,hf,nnodes);
    stateGuess.row(1) = linspace(v0,vf,nnodes);
    stateGuess.row(2) = linspace(gamma0,gammaf,nnodes);
    stateGuess.row(3) = linspace(w0,0.8*w0,nnodes);



    problem.phases(1).guess.controls       = zeros(1,nnodes);
    problem.phases(1).guess.states         = stateGuess;
    problem.phases(1).guess.time           = linspace(t0min, tfmax, nnodes);


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae 		= &dae;
    problem.events 		= &events;
    problem.linkages		= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_method                  	= "IPOPT";
    algorithm.scaling                     	= "automatic";
    algorithm.derivatives                 	= "numerical";
    algorithm.collocation_method            = "trapezoidal";
    algorithm.nlp_iter_max                	= 1000;
    algorithm.nlp_tolerance               	= 1.e-6;
    algorithm.mesh_refinement               = "automatic";
    algorithm.mr_max_iterations             = 4;
    algorithm.defect_scaling                = "jacobian-based";




////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////


    psopt(solution, problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    x 		= solution.get_states_in_phase(1);
    u 		= solution.get_controls_in_phase(1);
    t 		= solution.get_time_in_phase(1);
    H       = solution.get_dual_hamiltonian_in_phase(1);

    MatrixXd h     = x.row(0); 
    MatrixXd v     = x.row(1); 
    MatrixXd gamma = x.row(2); 
    MatrixXd w     = x.row(3); 

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x,"x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////


    plot(t,h/1000.0,problem.name + ": altitude", "time (s)", "altitude (x1,000 ft)", "h");
    plot(t,v/100.0,problem.name + ": velocity", "time (s)", "velocity (x100 ft/s)", "v");
    plot(t,gamma*180/pi,problem.name + ": flight path angle", "time (s)", "gamma (deg)", "gamma");
    plot(t,w/10000.0,problem.name + ": weight", "time (s)", "w (x10,000 lb)", "w");
    plot(t,u*180/pi,problem.name + ": angle of attack", "time (s)", "alpha (deg)", "alpha");

    plot(t,h/1000.0,problem.name + ": altitude", "time (s)", "altitude (x1,000 ft", "h",
                                   "pdf","climb_altitude.pdf");
    plot(t,v/100.0,problem.name + ": velocity", "time (s)", "velocity (x100 ft/s)", "v",
                            "pdf","climb_velocity.pdf");
    plot(t,gamma*180/pi,problem.name + ": flight path angle", "time (s)", "gamma (deg)", "gamma",
                                     "pdf","climb_fpa.pdf");
    plot(t,w/10000.0,problem.name + ": weight", "time (s)", "w (x10,000 lb)", "w", "pdf",
                                      "weight.pdf");
    plot(t,u*180/pi,problem.name + ": angle of attack", "time (s)", "alpha (deg)", "alpha",
                                      "pdf", "alpha.pdf");


}


void atmosphere(adouble* alt,adouble* sigma,adouble* delta,adouble* theta, Constants_& CONSTANTS)
// US Standard Atmosphere Model 1976
// Adopted from original Fortran 90 code by Ralph Carmichael
// Fortran code located at: http://www.pdas.com/programs/atmos.f90
{
/*!   -------------------------------------------------------------------------
! PURPOSE - Compute the properties of the 1976 standard atmosphere to 86 km.
! AUTHOR - Ralph Carmichael, Public Domain Aeronautical Software
! NOTE - If alt > 86, the values returned will not be correct, but they will
!   not be too far removed from the correct values for density.
!   The reference document does not use the terms pressure and temperature
!   above 86 km.
  IMPLICIT NONE
!============================================================================
!     A R G U M E N T S                                                     |
!============================================================================
  alt        ! geometric altitude, km.
  sigma      ! density/sea-level standard density
  delta      ! pressure/sea-level standard pressure
  theta      ! temperature/sea-level standard temperature
*/

/*!============================================================================
!     L O C A L   C O N S T A N T S                                         |
!============================================================================
*/
  double REARTH = 6369.0;                 // radius of the Earth (km)
  double GMR = 34.163195;                 // hydrostatic constant
  int NTAB=8;       // number of entries in the defining tables
/*!============================================================================
!     L O C A L   V A R I A B L E S                                         |
!============================================================================
*/
  int i,j,k;                                                  // counters
  adouble h;                                      // geopotential altitude (km)
  adouble tgrad, tbase;      // temperature gradient and base temp of this layer
  adouble tlocal;                                           // local temperature
  adouble deltah;                             // height above base of this layer
/*!============================================================================
!     L O C A L   A R R A Y S   ( 1 9 7 6   S T D.  A T M O S P H E R E )   |
!============================================================================
*/

  MatrixXd& htab = *CONSTANTS.htab;
  MatrixXd& ttab = *CONSTANTS.ttab;
  MatrixXd& ptab = *CONSTANTS.ptab;
  MatrixXd& gtab = *CONSTANTS.gtab;


//!----------------------------------------------------------------------------
  h=(*alt)*REARTH/((*alt)+REARTH);      //convert geometric to geopotential altitude

  i=1;
  j=NTAB;                                       // setting up for binary search
  while (j<=i+1) {
    k=(i+j)/2;                                              // integer division
    if (h < htab(k-1)) {
      j=k;
    } else {
       i=k;
    }
  }

  tgrad=gtab(i-1);                                     // i will be in 1...NTAB-1
  tbase=ttab(i-1);
  deltah=h-htab(i-1);
  tlocal=tbase+tgrad*deltah;
  *theta=tlocal/ttab(0);                                    // temperature ratio

  if (tgrad == 0.0) {                                  //  pressure ratio
    *delta=ptab(i-1)*exp(-GMR*deltah/tbase);
  } else {
    *delta=ptab(i-1)*pow(tbase/tlocal, GMR/tgrad);
  }

  *sigma=(*delta)/(*theta);                                           // density ratio
  return;

}


void atmosphere_model(adouble* rho, adouble* M, adouble v, adouble h, Constants_& CONSTANTS)
{
   double feet2meter = 0.3048;
   double kgperm3_to_slug_per_feet3 = 0.062427960841/32.174049;
   adouble alt, sigma, delta, theta;
   alt = h.value()*feet2meter/1000.0;

   // Call the standard atmosphere model 1976

   atmosphere(&alt, &sigma, &delta, &theta, CONSTANTS);

   adouble rho1 = 1.22521 * sigma; // Multiply by standard density at zero altitude and 15 deg C.

   rho1 = rho1*kgperm3_to_slug_per_feet3;

   *rho = rho1;

   adouble T;
   adouble mach;

   double TempStandardSeaLevel = 288.15; // in K, or 15 deg C.

   T = theta*TempStandardSeaLevel;

   adouble a = 20.0468 * sqrt(T); // Speed of sound in m/s.

   a = a/feet2meter;  // Speed of sound in ft/s

   mach = v/a;

   *M = mach;


  return;
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
